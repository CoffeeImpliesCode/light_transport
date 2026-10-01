//! Headless render benchmark.
//!
//! Times the render without a window, a compositor or a GPU:
//!
//! ```text
//! cargo run --release -- --bench
//! cargo run --release -- --bench --serial --out before.ppm
//! ```
//!
//! Two modes, because they answer different questions.
//!
//! The default drives the real worker pool, so it reports the throughput a
//! user actually gets. It is not reproducible: each worker draws from its own
//! thread RNG, so two runs differ in the last bits.
//!
//! `--serial` renders on one seeded thread. That makes it reproducible, which
//! is what makes it the oracle for "did this change alter the image?" — diff
//! the two `.ppm` files. It is also what to point a profiler at, because a
//! single thread gives a single stack instead of twelve.
//!
//! Options, all optional:
//!
//! | flag | default | meaning |
//! |------|---------|---------|
//! | `--size WxH` | `512x512` | image size in pixels |
//! | `--samples N` | `8` | samples per pixel |
//! | `--bounces N` | `3` | path depth |
//! | `--iters N` | `5` | measured passes, after the warm-up |
//! | `--warmup N` | `1` | unmeasured passes |
//! | `--reps N` | `1` | whole measurements to repeat, to get a spread |
//! | `--spheres N` | `5` | spheres in the scene; past five, hashed from the index |
//! | `--seed N` | `1` | seed for `--serial` |
//! | `--serial` | off | single threaded, reproducible |
//! | `--compare PATH` | none | paired A/B against another `light_transport` |
//! | `--out PATH` | none | write the last render as a binary PPM |
//! | `--quiet` | off | print only the summary line |
//!
//! ## Where the time goes, and what was already tried
//!
//! `perf record -e cycles:u` over a `--serial` pass at 512x512, 8 samples and
//! 3 bounces, on a Core i7-9750H:
//!
//! | share | symbol |
//! |-------|--------|
//! | 51.6% | `Scene::intersect` |
//! | 41.1% | `Renderer::cast` |
//! |  2.4% | `chacha20`, the RNG |
//! |  2.2% | `Scene::get_info` |
//! |  2.0% | `RenderWorkload::handle` |
//! |  0.4% | `powf`, the gamma curve |
//!
//! Four changes aimed at that table were built, measured against this binary
//! in interleaved pairs, and **reverted, because none of them beat it**:
//!
//! - Hoisting `|direction|^2` onto `Ray` so the sphere test stops recomputing
//!   it five times per ray, and computing the far root only when the near one
//!   is not already positive. Both genuinely remove work.
//! - Replacing the `Vec` index loops with iterator zips, to drop the 3.8% the
//!   profile attributes to the slice bounds-check panic path.
//! - Carrying the hit `Material` by reference instead of cloning it per hit.
//! - A backface reject in the sphere test: `c > 0 && b > 0` means the ray can
//!   never enter, so the discriminant, the sqrt and both divisions can be
//!   skipped. **This one measured as an 8.5% regression.** In this scene
//!   nearly every sphere is in front of the ray, so the test almost never
//!   fires and the extra compares are pure cost on the hottest function in the
//!   program. It is correct, and it is not worth it here.
//!
//! The part worth keeping: the profile says 51.6% is in `intersect`, but that
//! time is the arithmetic itself — three dots, a sqrt and two divisions per
//! sphere — not bookkeeping. Clearing the bookkeeping bought nothing. The wins
//! left in `intersect` are changes to the traversal rather than the tidying:
//! a bounding-volume hierarchy, or SIMD across the five spheres. Both are far
//! bigger than anything above, and both have to keep the render bit-exact,
//! which the `--serial --out` oracle is there to check.
//!
//! To reproduce the ranking:
//!
//! ```text
//! cargo build --profile profiling
//! perf record -F 999 -e cycles:u -g -- \
//!     target/profiling/light_transport --bench --serial
//! perf report -i perf.data --stdio --no-children -g none --sort srcline
//! ```
//!
//! ## Deciding whether a change is faster
//!
//! The four reverts above were not a wrong profile. They were a missing error
//! bar: the same unchanged binary, measured twice, did not give the same
//! number to the precision anyone was arguing over, on a box that is usually
//! busy. A single mean reports no spread, so a difference of a few percent
//! could be read as a result when it was the machine. Three flags exist to
//! close that hole.
//!
//! `--reps N` repeats the whole warm-up-plus-iters measurement N times and
//! reports the spread across the rep means, as a coefficient of variation,
//! `sd / median`. That percentage is the bar: a later A/B difference smaller
//! than it is noise, and `--compare` says "inconclusive" rather than naming a
//! winner. A single rep has no spread to report, so the line is only printed
//! when `--reps` is above 1 and the default output is unchanged.
//!
//! `--compare PATH` is the paired form. It alternates this binary and the
//! other one, one rep each, rather than running all of one and then all of the
//! other: load drifts over seconds, and a blocked layout hands that drift to
//! whichever side ran second. The child is spawned with the flags it would
//! have received, minus `--compare` and its own `--reps`, and its
//! `Mrays/s mean` is parsed back out. A child that exits non-zero or does not
//! print that field is a hard error carrying its stderr — a zero there would
//! silently turn into "the other build is infinitely fast".
//!
//! `--spheres N` scales the intersection loop, which is where the profile
//! above puts most of the time, so a traversal change can be measured against
//! a scene that is actually sphere-bound. The demo's own five stay at indices
//! 0..5, so `--spheres 5` is the scene the earlier rounds were measured
//! against, and everything past the fifth is derived from its index with an
//! integer hash — never an RNG — so the render stays reproducible and
//! `--serial --out` is still the oracle.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Duration, Instant};

use rand::SeedableRng;

use crate::geometry::Sphere;
use crate::image::{Color, Image};
use crate::material::Material;
use crate::math::{Vec3, F};
use crate::renderer::{
    blit_chunk, chunk_ranges, demo_scene, RenderChunk, RenderWorkload, Renderer, Scene,
};

/// Tile edge used by `Renderer::render`. Kept in step with it so both modes
/// hand out the same work.
const TILE: usize = 8;

/// A pass that never finishes is a bug in the renderer, not a slow scene.
const PASS_TIMEOUT: Duration = Duration::from_secs(300);

/// `Scene::get_info` reads an id at or above this as a plane index, so a
/// sphere numbered that high would look up the wrong geometry. `--spheres`
/// is capped below it; the demo's five planes sit at 100..105.
const PLANE_ID_BASE: usize = 100;

const USAGE: &str = "\
usage: light_transport --bench [--serial] [--size WxH] [--samples N] [--bounces N]
                       [--iters N] [--warmup N] [--reps N] [--spheres N] [--seed N]
                       [--compare PATH] [--out PATH] [--quiet]";

pub struct Options {
    serial: bool,
    quiet: bool,
    out: Option<PathBuf>,
    compare: Option<PathBuf>,
    size: [usize; 2],
    samples: usize,
    bounces: usize,
    iters: usize,
    warmup: usize,
    reps: usize,
    /// `None` leaves `demo_scene` alone. `Some(0)` is a mistake and is
    /// rejected by name, so the unset case cannot share a value with it.
    spheres: Option<usize>,
    seed: u64,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            serial: false,
            quiet: false,
            out: None,
            compare: None,
            size: [512, 512],
            samples: 8,
            bounces: 3,
            iters: 5,
            warmup: 1,
            reps: 1,
            spheres: None,
            seed: 1,
        }
    }
}

impl Options {
    /// `Err` names the argument that was rejected, so a typo cannot silently
    /// benchmark the wrong workload and report a plausible number.
    fn parse(args: &[String]) -> Result<Self, String> {
        fn value(args: &[String], i: &mut usize, flag: &str) -> Result<String, String> {
            *i += 1;
            args.get(*i)
                .cloned()
                .ok_or_else(|| format!("{flag} needs a value"))
        }
        fn num(args: &[String], i: &mut usize, flag: &str) -> Result<usize, String> {
            let raw = value(args, i, flag)?;
            raw.parse()
                .map_err(|_| format!("{flag} needs a number, got `{raw}`"))
        }
        fn seed(args: &[String], i: &mut usize) -> Result<u64, String> {
            let raw = value(args, i, "--seed")?;
            raw.parse()
                .map_err(|_| format!("--seed needs a number, got `{raw}`"))
        }

        let mut opts = Self::default();
        let mut i = 0;

        while i < args.len() {
            match args[i].as_str() {
                "--serial" => opts.serial = true,
                "--quiet" => opts.quiet = true,
                "--out" => opts.out = Some(PathBuf::from(value(args, &mut i, "--out")?)),
                "--compare" => {
                    opts.compare = Some(PathBuf::from(value(args, &mut i, "--compare")?))
                }
                "--seed" => opts.seed = seed(args, &mut i)?,
                "--samples" => opts.samples = num(args, &mut i, "--samples")?,
                "--bounces" => opts.bounces = num(args, &mut i, "--bounces")?,
                "--iters" => opts.iters = num(args, &mut i, "--iters")?,
                "--warmup" => opts.warmup = num(args, &mut i, "--warmup")?,
                "--reps" => opts.reps = num(args, &mut i, "--reps")?,
                "--spheres" => opts.spheres = Some(num(args, &mut i, "--spheres")?),
                "--size" => {
                    let raw = value(args, &mut i, "--size")?;
                    let (w, h) = raw
                        .split_once(['x', 'X'])
                        .ok_or_else(|| format!("--size wants WxH, got `{raw}`"))?;
                    opts.size = [
                        w.parse()
                            .map_err(|_| format!("--size width needs a number, got `{w}`"))?,
                        h.parse()
                            .map_err(|_| format!("--size height needs a number, got `{h}`"))?,
                    ];
                }
                other => return Err(format!("unknown argument `{other}`")),
            }
            i += 1;
        }

        if opts.size[0] == 0 || opts.size[1] == 0 {
            return Err("--size needs a non-zero width and height".into());
        }
        if opts.samples == 0 {
            return Err("--samples must be at least 1".into());
        }
        if opts.iters == 0 {
            return Err("--iters must be at least 1".into());
        }
        if opts.reps == 0 {
            return Err("--reps must be at least 1".into());
        }
        if let Some(count) = opts.spheres {
            if count == 0 {
                return Err("--spheres must be at least 1".into());
            }
            if count > PLANE_ID_BASE {
                return Err(format!(
                    "--spheres must be at most {PLANE_ID_BASE}, ids from there up read as planes"
                ));
            }
        }
        Ok(opts)
    }
}

/// Some if the command line asked for the benchmark, so `main` hands over.
pub fn take_command_line() -> Option<Options> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if !args.iter().any(|a| a == "--bench") {
        return None;
    }
    let rest: Vec<String> = args.into_iter().filter(|a| a != "--bench").collect();

    match Options::parse(&rest) {
        Ok(opts) => Some(opts),
        Err(e) => {
            eprintln!("light_transport: {e}\n{USAGE}");
            std::process::exit(2);
        }
    }
}

/// Render every tile on this thread, walked in a fixed order from a fixed
/// seed, so the image does not depend on scheduling.
fn render_serial(scene: &Scene, size: [usize; 2], seed: u64) -> Image {
    let mut image = Image::new(size);

    for &(j0, j1) in &chunk_ranges(size[1], TILE) {
        for &(i0, i1) in &chunk_ranges(size[0], TILE) {
            let mut workload = RenderWorkload::new((i0, j0), (i1, j1), size);
            let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
            workload.handle(scene, &mut rng);
            blit_chunk(
                &mut image,
                &RenderChunk {
                    start: workload.start,
                    end: workload.end,
                    pixels: workload.pixels,
                },
            );
        }
    }

    image
}

/// Drive the real worker pool the way the Render button does.
///
/// `render` pushes every tile into the work channel synchronously, so the tile
/// count is known before the workers start. Counting the tiles that come back
/// is therefore an exact completion signal, which pixel inspection is not: a
/// ray that lands on a non-emissive surface legitimately resolves to black.
fn render_parallel(scene: &Scene, size: [usize; 2]) -> Image {
    let mut renderer = Renderer::new(scene);
    *renderer.size.lock() = size;

    let expected = chunk_ranges(size[0], TILE).len() * chunk_ranges(size[1], TILE).len();
    // A clone rather than `&renderer.rx`: a `while let` scrutinee keeps its
    // temporary borrow alive for the whole body, which would collide with the
    // mutable borrow of `renderer.image`. Nothing else drains the receiver in
    // this mode, so the clone is the only consumer.
    let rx = renderer.rx.clone();
    let start = Instant::now();
    renderer.render(scene);

    let mut received = 0usize;
    while received < expected {
        while let Ok(chunk) = rx.try_recv() {
            blit_chunk(&mut renderer.image, &chunk);
            received += 1;
        }
        assert!(
            start.elapsed() < PASS_TIMEOUT,
            "only {received}/{expected} tiles came back"
        );
        std::thread::sleep(Duration::from_millis(1));
    }

    // `Renderer` has a `Drop` impl, so its fields cannot be moved out.
    // Swapping in an empty image leaves the join to run against a renderer
    // every worker already finished with.
    std::mem::replace(&mut renderer.image, Image::new([0, 0]))
}

/// One full pass, timed. The image comes back for the `--out` oracle; keeping
/// it is the caller's choice.
fn pass(opts: &Options, scene: &Scene) -> (Image, Duration) {
    let start = Instant::now();
    let image = if opts.serial {
        render_serial(scene, opts.size, opts.seed)
    } else {
        render_parallel(scene, opts.size)
    };
    (image, start.elapsed())
}

/// Write the image as a binary PPM (P6). Trivial to produce and trivial to
/// compare, so a before/after diff needs no image library.
fn write_ppm(path: &Path, image: &Image) -> std::io::Result<()> {
    let mut out = std::fs::File::create(path)?;
    write!(out, "P6\n{} {}\n255\n", image.size[0], image.size[1])?;

    let bytes = image.bytes();
    let mut row = Vec::with_capacity(image.size[0] * 3);
    for y in 0..image.size[1] {
        row.clear();
        for x in 0..image.size[0] {
            let p = (y * image.size[0] + x) * 4;
            row.extend_from_slice(&bytes[p..p + 3]);
        }
        out.write_all(&row)?;
    }
    out.flush()
}

/// The `k`-th hashed value of `index`, in `[0, 1)`.
///
/// splitmix64's finalizer over `index * 8 + k`. Integer only, so a given
/// index always yields the same sphere, and eight values per index stop the
/// axes, radius and colour from repeating in step with each other.
fn unit(index: usize, k: usize) -> F {
    let mut x = (index as u64)
        .wrapping_mul(8)
        .wrapping_add(k as u64)
        .wrapping_add(0x9e37_79b9_7f4a_7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^= x >> 31;
    (x >> 40) as F / (1u64 << 24) as F
}

/// The `index`-th added sphere, derived from `index` and nothing else.
///
/// `Scene::get_info` looks a sphere up by its `id`, so the id has to be the
/// position in `Scene::spheres` — which is why `--spheres` is capped below
/// the plane id base. Radii and colours vary because a lattice of identical
/// spheres is a uniform workload, and this exists to make the loop
/// non-uniform the way a real scene would.
fn extra_sphere(index: usize) -> Sphere {
    let radius = 0.05 + 0.03 * unit(index, 3);
    let reflecting = unit(index, 4);
    Sphere {
        id: index,
        // The demo camera sits at (-5, 0, 0) and sees the whole one-unit room,
        // so anything inside the box the demo spheres occupy is on frame and
        // inside the walls. Outside it, spheres cost time and change nothing.
        center: Vec3::new([
            2.0 * unit(index, 0) - 1.0,
            2.0 * unit(index, 1) - 1.0,
            2.0 * unit(index, 2) - 1.0,
        ]) * 0.35,
        radius,
        material: Material {
            color: Color::new(
                0.25 + 0.75 * unit(index, 5),
                0.25 + 0.75 * unit(index, 6),
                0.25 + 0.75 * unit(index, 7),
                1.0,
            ),
            // The demo's five keep these summing to 1, so raising `--spheres`
            // does not quietly change how much light the scene returns.
            emmission: 0.0,
            reflecting,
            diffuse: 1.0 - reflecting,
        },
    }
}

/// The demo scene with `--spheres` applied.
///
/// Indices 0..5 stay the demo's own, so `--spheres 5` is exactly the scene
/// the earlier rounds were measured against and those numbers stay
/// comparable. Below five the list is truncated, which keeps the remaining
/// ids equal to their positions in the shortened list.
fn build_scene(opts: &Options) -> Scene {
    let mut scene = demo_scene(opts.samples, opts.bounces);
    let Some(want) = opts.spheres else {
        return scene;
    };

    let demo = scene.spheres.len();
    scene.spheres.truncate(want.min(demo));
    for index in demo..want {
        scene.spheres.push(extra_sphere(index));
    }
    scene
}

/// One rep: the warm-up passes, then the measured ones.
///
/// A rep is the unit that `--reps` repeats and that `--compare` alternates,
/// so it reports its own mean and hands back its own image. Pooling is the
/// caller's business, because the two callers pool differently.
struct Rep {
    image: Image,
    mean: f64,
    per_pass: Vec<f64>,
}

fn measure(opts: &Options, scene: &Scene) -> Rep {
    // Warm-up faults the tiles in and lets the allocator settle, so the
    // measured passes are not paying first-touch costs.
    for _ in 0..opts.warmup {
        let _ = pass(opts, scene);
    }

    let mut per_pass = Vec::with_capacity(opts.iters);
    let mut last: Option<Image> = None;
    for _ in 0..opts.iters {
        let (image, elapsed) = pass(opts, scene);
        per_pass.push(elapsed.as_secs_f64());
        last = Some(image);
    }

    let mean = per_pass.iter().sum::<f64>() / per_pass.len() as f64;
    Rep {
        image: last.expect("--iters is at least 1, so a pass ran"),
        mean,
        per_pass,
    }
}

/// Middle value, averaging the two middle ones when the count is even.
fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("pass times are finite"));
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    }
}

/// Sample standard deviation. One value cannot show spread, so it reports
/// none rather than dividing by zero.
fn std_dev(values: &[f64]) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let sum: f64 = values.iter().map(|v| (v - mean) * (v - mean)).sum();
    (sum / (values.len() - 1) as f64).sqrt()
}

/// How much of a median a spread is, in percent.
///
/// This is the bar a difference has to clear. A few percent of swing between
/// two runs of the same binary is what the coefficient measures, so a claim
/// smaller than it is a claim about the machine, not about the code.
fn spread_pct(sd: f64, med: f64) -> f64 {
    100.0 * sd / med
}

/// The flags the other binary needs to measure the same workload.
///
/// A flag and its value are two separate entries in `argv`, which is how a
/// shell hands them over. Packed into one entry, `--size 512x512` reaches
/// the child as a single token matching no arm of its parser: it exits 2
/// having measured nothing, and the parent's error names that exit.
///
/// `--compare` is left out so the child does not try to compare itself, and
/// `--reps` is left out because the parent owns the alternation: a child
/// running its own reps would do N times the work behind one table row.
/// `--out` is left out for a different reason — the sides alternate writes to
/// one path, so the file would end up holding whichever ran last, which is
/// not what the caller asked for. A writes its own image from its own last
/// pass instead.
fn child_args(opts: &Options) -> Vec<String> {
    let mut args = Vec::with_capacity(16);
    args.push("--bench".to_string());
    args.push("--size".to_string());
    args.push(format!("{}x{}", opts.size[0], opts.size[1]));
    args.push("--samples".to_string());
    args.push(opts.samples.to_string());
    args.push("--bounces".to_string());
    args.push(opts.bounces.to_string());
    args.push("--iters".to_string());
    args.push(opts.iters.to_string());
    args.push("--warmup".to_string());
    args.push(opts.warmup.to_string());
    args.push("--seed".to_string());
    args.push(opts.seed.to_string());
    if opts.serial {
        args.push("--serial".to_string());
    }
    if let Some(count) = opts.spheres {
        args.push("--spheres".to_string());
        args.push(count.to_string());
    }
    if opts.quiet {
        args.push("--quiet".to_string());
    }
    args
}

/// The child's stderr, or a note that it said nothing, so a failure names
/// something the reader can act on.
fn because(stderr: &str) -> String {
    if stderr.is_empty() {
        return " with no stderr".to_string();
    }
    format!(": {stderr}")
}

/// Seconds per pass, as the other binary reported it.
///
/// The summary line ends `... 12.34 Mrays/s mean  13.56 Mrays/s best`, so the
/// mean is the token before the *first* `Mrays/s`. Anything else — a non-zero
/// exit, a missing field, a value that is not a finite time — is an error
/// carrying the child's own stderr. Substituting a zero would turn a broken
/// child into "the other build is infinitely fast" and quietly hand the
/// verdict to A.
fn child_seconds(path: &Path, opts: &Options) -> Result<f64, String> {
    let output = Command::new(path)
        .args(child_args(opts))
        .output()
        .map_err(|e| format!("cannot run {}: {e}", path.display()))?;
    let stderr = String::from_utf8_lossy(&output.stderr);
    let stderr = stderr.trim();

    if !output.status.success() {
        return Err(format!(
            "{} exited with {}{}",
            path.display(),
            output.status,
            because(stderr)
        ));
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let mray_per_s = stdout
        .split_whitespace()
        .collect::<Vec<_>>()
        .windows(2)
        .find(|pair| pair[1] == "Mrays/s")
        .and_then(|pair| pair[0].parse::<f64>().ok())
        .ok_or_else(|| {
            format!(
                "{} printed no `Mrays/s mean` field{}",
                path.display(),
                because(stderr)
            )
        })?;

    if !mray_per_s.is_finite() || mray_per_s <= 0.0 {
        return Err(format!(
            "{} reported {mray_per_s} Mrays/s, which is not a time{}",
            path.display(),
            because(stderr)
        ));
    }

    let rays = (opts.size[0] * opts.size[1] * opts.samples) as f64;
    Ok(rays / (mray_per_s * 1e6))
}

/// Paired A/B against another build.
///
/// The two sides alternate one rep at a time instead of running in blocks.
/// Load on a busy box drifts over seconds, and a blocked layout gives all of
/// that drift to whichever side ran second — which is how a change that does
/// nothing gets declared a winner. A verdict is only printed when the two
/// observed ranges are disjoint; otherwise the honest answer is that the
/// medians differ by less than the spread of a single side, and saying so is
/// worth more than picking a side.
fn compare(opts: &Options, scene: &Scene, other: &Path) {
    let here = std::env::current_exe()
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| "<this binary>".to_string());
    println!("bench compare A={here} B={}", other.display());

    let mut a_means = Vec::with_capacity(opts.reps);
    let mut b_means = Vec::with_capacity(opts.reps);
    let mut last: Option<Image> = None;

    for rep in 1..=opts.reps {
        let a = measure(opts, scene);
        let b = child_seconds(other, opts).unwrap_or_else(|e| {
            eprintln!("bench compare: {e}");
            std::process::exit(1);
        });
        if !opts.quiet {
            println!(
                "bench compare rep {rep}   A {:.3}   B {b:.3}   A/B {:.3}",
                a.mean,
                a.mean / b
            );
        }
        a_means.push(a.mean);
        b_means.push(b);
        last = Some(a.image);
    }

    let a_median = median(&a_means);
    let b_median = median(&b_means);
    let a_lo = a_means.iter().copied().fold(f64::INFINITY, f64::min);
    let a_hi = a_means.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let b_lo = b_means.iter().copied().fold(f64::INFINITY, f64::min);
    let b_hi = b_means.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    let diff = 100.0 * (a_median - b_median) / b_median;
    let verdict = if a_lo <= b_hi && b_lo <= a_hi {
        format!("inconclusive: the spreads overlap, the medians differ by {diff:.1}%")
    } else if diff >= 0.0 {
        format!("A is {diff:.1}% slower")
    } else {
        format!("A is {:.1}% faster", -diff)
    };

    println!(
        "bench compare A {a_median:.3}s median  B {b_median:.3}s median  {verdict}  \
         (A sd {:.1}%, B sd {:.1}%)",
        spread_pct(std_dev(&a_means), a_median),
        spread_pct(std_dev(&b_means), b_median),
    );

    if let (Some(path), Some(image)) = (&opts.out, &last) {
        write_ppm(path, image).unwrap_or_else(|e| {
            eprintln!("bench: cannot write {}: {e}", path.display());
            std::process::exit(1);
        });
        let bytes = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
        println!("bench wrote {} ({bytes} bytes)", path.display());
    }
}

pub fn run(opts: &Options) {
    let scene = build_scene(opts);

    if let Some(other) = &opts.compare {
        compare(opts, &scene, other);
        return;
    }

    let label = if opts.serial { "serial" } else { "parallel" };

    if !opts.quiet {
        println!(
            "bench {label} {}x{} samples={} bounces={} iters={} warmup={} seed={} \
             spheres={} threads={}",
            opts.size[0],
            opts.size[1],
            opts.samples,
            opts.bounces,
            opts.iters,
            opts.warmup,
            opts.seed,
            scene.spheres.len(),
            num_cpus::get(),
        );
    }

    let mut per_pass: Vec<f64> = Vec::with_capacity(opts.iters * opts.reps);
    let mut rep_means: Vec<f64> = Vec::with_capacity(opts.reps);
    let mut last: Option<Image> = None;
    for _ in 0..opts.reps {
        let rep = measure(opts, &scene);
        per_pass.extend_from_slice(&rep.per_pass);
        rep_means.push(rep.mean);
        last = Some(rep.image);
    }

    // The error bar the pooled line cannot have. One rep has no spread to
    // report, so the line is only printed when there is more than one: the
    // default output has to stay what the earlier rounds compared against.
    if opts.reps > 1 {
        let med = median(&rep_means);
        let lo = rep_means.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = rep_means.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let sd = std_dev(&rep_means);
        println!(
            "bench {label} reps={} median {med:.3}s  min {lo:.3}s  max {hi:.3}s  \
             sd {sd:.3}s  ({:.1}% of median)",
            opts.reps,
            spread_pct(sd, med),
        );
    }

    let passes = per_pass.len();
    let rays = (opts.size[0] * opts.size[1] * opts.samples * passes) as f64;
    let total: f64 = per_pass.iter().sum();
    let best = per_pass.iter().copied().fold(f64::INFINITY, f64::min);
    let mean = total / passes as f64;
    let median_pass = median(&per_pass);

    println!(
        "bench {label} {best:.4}s best  {mean:.4}s mean  {median_pass:.4}s median  \
         {:.2} Mrays/s mean  {:.2} Mrays/s best",
        rays / total / 1e6,
        rays / (best * passes as f64) / 1e6,
    );
    if !opts.quiet {
        println!("bench per-pass seconds: {per_pass:?}");
    }

    if let (Some(path), Some(image)) = (&opts.out, &last) {
        write_ppm(path, image).unwrap_or_else(|e| {
            eprintln!("bench: cannot write {}: {e}", path.display());
            std::process::exit(1);
        });
        let bytes = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
        println!("bench wrote {} ({bytes} bytes)", path.display());
    }
}
