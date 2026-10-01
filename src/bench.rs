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
//! | `--seed N` | `1` | seed for `--serial` |
//! | `--serial` | off | single threaded, reproducible |
//! | `--out PATH` | none | write the last render as a binary PPM |
//! | `--quiet` | off | print only the summary line |

use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use rand::SeedableRng;

use crate::image::Image;
use crate::renderer::{
    blit_chunk, chunk_ranges, demo_scene, RenderChunk, RenderWorkload, Renderer, Scene,
};

/// Tile edge used by `Renderer::render`. Kept in step with it so both modes
/// hand out the same work.
const TILE: usize = 8;

/// A pass that never finishes is a bug in the renderer, not a slow scene.
const PASS_TIMEOUT: Duration = Duration::from_secs(300);

const USAGE: &str = "\
usage: light_transport --bench [--serial] [--size WxH] [--samples N] [--bounces N]
                        [--iters N] [--warmup N] [--seed N] [--out PATH] [--quiet]";

pub struct Options {
    serial: bool,
    quiet: bool,
    out: Option<PathBuf>,
    size: [usize; 2],
    samples: usize,
    bounces: usize,
    iters: usize,
    warmup: usize,
    seed: u64,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            serial: false,
            quiet: false,
            out: None,
            size: [512, 512],
            samples: 8,
            bounces: 3,
            iters: 5,
            warmup: 1,
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
                "--seed" => opts.seed = seed(args, &mut i)?,
                "--samples" => opts.samples = num(args, &mut i, "--samples")?,
                "--bounces" => opts.bounces = num(args, &mut i, "--bounces")?,
                "--iters" => opts.iters = num(args, &mut i, "--iters")?,
                "--warmup" => opts.warmup = num(args, &mut i, "--warmup")?,
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

pub fn run(opts: &Options) {
    let scene = demo_scene(opts.samples, opts.bounces);
    let label = if opts.serial { "serial" } else { "parallel" };

    if !opts.quiet {
        println!(
            "bench {label} {}x{} samples={} bounces={} iters={} warmup={} seed={} threads={}",
            opts.size[0],
            opts.size[1],
            opts.samples,
            opts.bounces,
            opts.iters,
            opts.warmup,
            opts.seed,
            num_cpus::get(),
        );
    }

    // Warm-up faults the tiles in and lets the allocator settle, so the
    // measured passes are not paying first-touch costs.
    for _ in 0..opts.warmup {
        let _ = pass(opts, &scene);
    }

    let mut per_pass = Vec::with_capacity(opts.iters);
    let mut last: Option<Image> = None;
    for _ in 0..opts.iters {
        let (image, elapsed) = pass(opts, &scene);
        per_pass.push(elapsed.as_secs_f64());
        last = Some(image);
    }

    let rays = (opts.size[0] * opts.size[1] * opts.samples * opts.iters) as f64;
    let total: f64 = per_pass.iter().sum();
    let best = per_pass.iter().copied().fold(f64::INFINITY, f64::min);
    let mean = total / per_pass.len() as f64;
    let mut sorted = per_pass.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("pass times are finite"));
    let median = sorted[sorted.len() / 2];

    println!(
        "bench {label} {best:.4}s best  {mean:.4}s mean  {median:.4}s median  \
         {:.2} Mrays/s mean  {:.2} Mrays/s best",
        rays / total / 1e6,
        rays / (best * per_pass.len() as f64) / 1e6,
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
