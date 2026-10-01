use num_cpus;
use std::{
    sync::{Arc, RwLock},
    thread::JoinHandle,
};

use eframe::{egui, epaint::ColorImage};
use egui::mutex::Mutex;

use crate::image::{Color, RGBA};
use crate::Image;

use crate::material::{schlick_fresnel, Material};
use crate::math::{Constants, Vec3, F};

use crate::geometry::*;

pub struct Renderer {
    pub scene: Arc<RwLock<Scene>>,
    pub workers: Vec<JoinHandle<()>>,
    pub send: Option<crossbeam::channel::Sender<RenderWorkload>>,
    /// Finished chunks coming back from the workers. The UI thread owns
    /// `image` and is its only reader, so a chunk has to travel as data.
    pub rx: crossbeam::channel::Receiver<RenderChunk>,
    chunks: Option<crossbeam::channel::Sender<RenderChunk>>,
    // pub image: Arc<Mutex<Option<ColorImage>>>,
    pub image: Image,
    pub size: Arc<Mutex<[usize; 2]>>,
    // Kept current by the worker threads; no call site reads it yet.
    #[allow(dead_code)]
    pub avg_rps: Arc<Mutex<(f64, usize)>>,
}

/// A finished tile of one render pass, owned by whoever holds it.
#[derive(Debug)]
pub struct RenderChunk {
    /// Inclusive top left corner, in pixels.
    pub start: (usize, usize),
    /// Inclusive bottom right corner, in pixels.
    pub end: (usize, usize),
    /// Row major, `width * height` pixels of the tile.
    pub pixels: Vec<RGBA>,
}

/// Split `0..len` into inclusive `(start, end)` ranges of at most `size`
/// elements. The ranges tile the whole range exactly: no gaps, no overlap
/// and nothing past `len`, whatever the remainder is.
pub fn chunk_ranges(len: usize, size: usize) -> Vec<(usize, usize)> {
    assert!(size > 0, "chunk size must be non-zero");
    if len == 0 {
        return Vec::new();
    }

    (0..len)
        .step_by(size)
        .map(|start| (start, start.saturating_add(size).min(len) - 1))
        .collect()
}

#[derive(Debug, Clone)]
pub struct Camera {
    pub origin: Vec3,
    pub right: Vec3,
    pub up: Vec3,
    pub width: F,
    pub height: F,
}

// trait Element<T: ?Sized> {
//     fn insert(&mut self, t: T);
//     fn get(&self, id: Id) -> Option<&T>;
//     fn get_mut(&mut self, id: Id) -> Option<&mut T>;
// }

// impl Element<Sphere> for Shapes {
//     fn insert(&mut self, t: Sphere) {
//         todo!()
//     }

//     fn get(&self, id: Id) -> Option<&Sphere> {
//         None
//     }

//     fn get_mut(&mut self, id: Id) -> Option<&mut Sphere> {
//         None
//     }
// }

// impl Element<Plane> for Shapes {
//     fn insert(&mut self, t: Plane) {
//         todo!()
//     }

//     fn get(&self, id: Id) -> Option<&Plane> {
//         None
//     }

//     fn get_mut(&mut self, id: Id) -> Option<&mut Plane> {
//         None
//     }
// }

// impl Element<Box<dyn Intersect>> for Shapes {
//     fn insert(&mut self, t: Box<dyn Intersect>) {
//         panic!()
//     }

//     fn get(&self, id: Id) -> Option<&Box<dyn Intersect>> {
//         None
//     }

//     fn get_mut(&mut self, id: Id) -> Option<&mut Box<dyn Intersect>> {
//         None
//     }
// }

#[derive(Debug, Clone)]
pub struct Scene {
    pub camera: Camera,
    pub spheres: Vec<Sphere>,
    pub planes: Vec<Plane>,
    // Directional light. The BRDF does not read it yet; the emitters carry
    // their own `emmission` term instead.
    #[allow(dead_code)]
    pub light: Vec3,
    pub ambient: Material,
    pub num_samples: usize,
    pub num_bounces: usize,
}

impl Scene {
    fn get_info(&self, ray: &Ray, inters: &Intersection) -> HitInfo {
        if inters.id >= 100 {
            let plane = &self.planes[inters.id - 100];

            HitInfo {
                normal: plane.normal,
                material: plane.material.clone(),
            }
        } else {
            let sphere = &self.spheres[inters.id];

            HitInfo {
                normal: ((ray.origin + ray.direction * inters.distance) - sphere.center)
                    .normalized(),
                material: sphere.material.clone(),
            }
        }
    }
}

impl Intersect for Scene {
    fn intersect(&self, ray: &Ray) -> Option<Intersection> {
        let mut closest: Intersection = Intersection {
            distance: F::INFINITY,
            id: 100,
        };
        for s in &self.spheres {
            if let Some(inter @ Intersection { distance, .. }) = s.intersect(ray) {
                if distance < closest.distance {
                    closest = inter
                }
            }
        }

        for p in &self.planes {
            if let Some(inter @ Intersection { distance, .. }) = p.intersect(ray) {
                if distance < closest.distance {
                    closest = inter
                }
            }
        }

        if closest.distance == F::INFINITY {
            None
        } else {
            Some(closest)
        }
    }
}
/// One tile of a render pass. The worker owns `pixels` outright and hands it
/// back to the UI thread, so no two threads ever touch the same buffer.
#[derive(Debug)]
pub struct RenderWorkload {
    pub start: (usize, usize),
    pub end: (usize, usize),
    /// Full image size, needed to map pixels to camera coordinates.
    pub size: [usize; 2],
    pub pixels: Vec<RGBA>,
}

impl RenderWorkload {
    pub fn new(start: (usize, usize), end: (usize, usize), size: [usize; 2]) -> Self {
        let pixels = vec![RGBA::new(0, 0, 0, 255); (end.0 - start.0 + 1) * (end.1 - start.1 + 1)];
        Self {
            start,
            end,
            size,
            pixels,
        }
    }

    /// Render this tile into `self.pixels`. Returns the number of rays cast.
    pub fn handle<R: rand::Rng + ?Sized>(&mut self, scene: &Scene, rng: &mut R) -> usize {
        let (direction, dx, dy) = {
            let direction = scene.camera.up.cross(scene.camera.right).normalized();
            let dx = scene.camera.right * scene.camera.width * 0.5;
            let dy = scene.camera.up * scene.camera.height * 0.5;
            (direction, dx, dy)
        };

        let pre_render: RGBA = Color::new(0.1, 0.1, 0.1, 1.0).into();

        let width = self.end.0 - self.start.0 + 1;
        debug_assert_eq!(self.pixels.len(), width * (self.end.1 - self.start.1 + 1));

        // Tile local, row major offset of a pixel.
        let at = |i: usize, j: usize| (j - self.start.1) * width + (i - self.start.0);

        let camera_dx = 2.0 / self.size[0] as F;
        let camera_dy = 2.0 / self.size[1] as F;

        for j in self.start.1..=self.end.1 {
            for i in self.start.0..=self.end.0 {
                self.pixels[at(i, j)] = pre_render;
            }
        }

        for i in self.start.0..=self.end.0 {
            self.pixels[at(i, self.start.1)] = Color::BLACK.into();
            self.pixels[at(i, self.end.1)] = Color::BLACK.into();
        }

        for j in self.start.1..=self.end.1 {
            self.pixels[at(self.start.0, j)] = Color::BLACK.into();
            self.pixels[at(self.end.0, j)] = Color::BLACK.into();
        }

        let mut cast_rays = 0;

        for j in self.start.1..=self.end.1 {
            let camera_y = j as F * camera_dy - 1.0;
            for i in self.start.0..=self.end.0 {
                let camera_x = i as F * camera_dx - 1.0;

                let ray_direction = direction - dy * camera_y + dx * camera_x;

                let ray = Ray {
                    origin: scene.camera.origin,
                    direction: ray_direction.normalized(),
                };

                let mut color = Color::BLACK;

                for _ in 0..scene.num_samples {
                    let c = Renderer::cast(&scene, &ray, scene.num_bounces, rng);
                    color += c;
                    cast_rays += 1;
                }
                color *= 1.0 / (scene.num_samples as F);
                self.pixels[at(i, j)] = color.gamma().into();
                // stage.store(j * image.dimension[1] + i, Ordering::Release);
            }
        }
        return cast_rays;
    }

    /// Hand the rendered tile over to the UI thread.
    pub fn into_chunk(self) -> RenderChunk {
        RenderChunk {
            start: self.start,
            end: self.end,
            pixels: self.pixels,
        }
    }
}

impl Renderer {
    pub fn new(scene: &Scene) -> Self {
        let (s, r) = crossbeam::channel::unbounded::<RenderWorkload>();
        let (tx, rx) = crossbeam::channel::unbounded::<RenderChunk>();

        let scene = Arc::new(RwLock::new(scene.clone()));
        let image = Image::new([crate::DEFAULT_IMAGE_WIDTH, crate::DEFAULT_IMAGE_HEIGHT]);
        let size = Arc::new(Mutex::new([
            crate::DEFAULT_IMAGE_WIDTH,
            crate::DEFAULT_IMAGE_HEIGHT,
        ]));

        let avg_rps = Arc::new(Mutex::new((0.0 as f64, 0 as usize)));

        let workers: Vec<_> = (0..num_cpus::get())
            .into_iter()
            .map(|_worker| {
                let r = r.clone();
                // let image = image.clone();
                let scene = scene.clone();
                let tx = tx.clone();
                let avg_rps = avg_rps.clone();
                // let stage = stage.clone();

                std::thread::spawn(move || {
                    // One generator per worker. Fetching `rng()` per sample
                    // made rand's reseeding-counter atomic the hottest single
                    // instruction in the binary under perf.
                    let mut rng = rand::rng();
                    for mut workload in r.iter() {
                        let start = std::time::Instant::now();
                        // running.store(true, Ordering::Relaxed);

                        let num_cast_rays = {
                            // A poisoned lock means some other render pass
                            // panicked. Stop this worker rather than add a
                            // panic of its own.
                            let Ok(scene) = scene.read() else {
                                return;
                            };

                            // workload.store(image.dimension[0] * image.dimension[1], Ordering::Relaxed);
                            // stage.store(0, Ordering::Relaxed);

                            /*let (direction, dx, dy) = {
                                                       let direction = scene.camera.up.cross(scene.camera.right).normalized();
                                                       let dx = scene.camera.right * scene.camera.width * 0.5;
                                                       let dy = scene.camera.up * scene.camera.height * 0.5;
                                                       (direction, dx, dy)
                                                   };

                                                   let camera_dx = 2.0 / img.size[0] as f64;
                                                   let camera_dy = 2.0 / img.size[1] as f64;
                            */
                            workload.handle(&scene, &mut rng)
                            /*
                            for j in 0..img.size[1] {
                                let camera_y = j as f64 * camera_dy - 1.0;
                                for i in 0..img.size[0] {
                                    let camera_x = i as f64 * camera_dx - 1.0;

                                    let ray_direction = direction - dy * camera_y + dx * camera_x;

                                    let ray = Ray {
                                        origin: scene.camera.origin,
                                        direction: ray_direction.normalized(),
                                    };

                                    let mut color = Color::BLACK;
                                    for _ in 0..scene.num_samples {
                                        color += Renderer::cast(&scene, &ray, scene.num_bounces);
                                    }
                                    color *= 1.0 / (scene.num_samples as f32);
                                    img[(i, j)] = color.into();
                                    // stage.store(j * image.dimension[1] + i, Ordering::Release);
                                }
                            }
                            */
                        };

                        // *image.lock() = Some(img);

                        // The finished tile travels back as owned data: no
                        // worker ever writes into the UI thread's image.
                        let chunk = workload.into_chunk();

                        let end = std::time::Instant::now();
                        // stage.store(2, Ordering::Release);
                        // stage.store(3, Ordering::Release);
                        // running.store(false, Ordering::Release);

                        let millis = ((end - start).as_micros() as f64) / 1000.0;
                        let secs = millis / 1000.0;

                        // A clock that coarse reads 0 for a fast chunk.
                        // Skipping beats poisoning the average with infinity.
                        if secs > 0.0 {
                            let rps = (num_cast_rays as f64) / secs;

                            {
                                let mut lock = avg_rps.lock();
                                let new_count = lock.1 + 1;
                                let new_avg_rps = lock.0 + (rps - lock.0) / (new_count as f64);
                                lock.0 = new_avg_rps;
                                lock.1 = new_count;
                            }
                        }

                        /*println!(
                            "Worker {:>2}: cast {:>8} rays in {:>8.2}ms ({:>14.2} rps)",
                            worker, num_cast_rays, millis, rps
                        );*/
                        /*println!(
                            "Worker {}: Chunk ({} {})..({} {}) rendered in {:.3}ms",
                            worker,
                            chunk.start.0,
                            chunk.start.1,
                            chunk.end.0,
                            chunk.end.1,
                            ((end - start).as_micros() as f32) / 1000.0
                        );*/

                        // Nobody left to hand the tile to, so stop working.
                        if tx.send(chunk).is_err() {
                            return;
                        }
                    }
                })
            })
            .collect();

        Renderer {
            image,
            size,
            scene,
            workers,
            avg_rps,
            send: Some(s),
            rx,
            chunks: Some(tx),
        }
    }

    pub fn render(&mut self, scene: &Scene) {
        // println!("Starting Render");
        /*{
            let lock = self.avg_rps.lock();
            println!("Avg rps per thread: {} ({} samples)", lock.0, lock.1);
        }*/

        // The UI thread owns `image`. Resetting it here, before any tile is
        // queued, is what keeps a second Render click from pulling the
        // rug out from under the workers of the first pass.
        let size = *self.size.lock();
        self.image = Image::new(size);

        match self.scene.write() {
            Ok(mut lock) => {
                *lock = scene.clone();
            }
            Err(e) => println!("Error: {:?}", e),
        }

        // let image = Arc::new(UnsafeCell::new(Image::new([0, 0])));

        const SIZE_X: usize = 8;
        const SIZE_Y: usize = 8;

        let xs = chunk_ranges(size[0], SIZE_X);
        let ys = chunk_ranges(size[1], SIZE_Y);

        // A tile runs from (x0, y0) to (x1, y1), both corners inclusive.
        let mut workloads = Vec::with_capacity(xs.len() * ys.len());
        for &(j0, j1) in &ys {
            for &(i0, i1) in &xs {
                workloads.push(RenderWorkload::new((i0, j0), (i1, j1), size));
            }
        }

        // println!("Workloads: {:#?}", workloads);
        // panic!();

        // self.send.clone().unwrap().try_send(()).unwrap();
        let s = self.send.as_ref().unwrap();

        workloads.into_iter().for_each(|w| s.try_send(w).unwrap());
    }

    /// Copy every finished tile the workers have handed back into the image.
    /// UI thread only, so the blit never races a worker.
    pub fn drain(&mut self) {
        // `try_recv` hands over an owned chunk, so the borrow of the
        // receiver ends before the image is touched.
        while let Ok(chunk) = self.rx.try_recv() {
            self.blit(&chunk);
        }
    }

    fn blit(&mut self, chunk: &RenderChunk) {
        blit_chunk(&mut self.image, chunk);
    }

    pub fn take_image(&mut self) -> ColorImage {
        ColorImage::from_rgba_unmultiplied(self.image.size, self.image.bytes())
    }

    #[inline(always)]
    pub fn cast<R: rand::Rng + ?Sized>(scene: &Scene, ray: &Ray, n: usize, rng: &mut R) -> Color {
        if let Some(closest) = scene.intersect(ray) {
            let pos = ray.origin + ray.direction * closest.distance * 0.9999;
            let hit = scene.get_info(ray, &closest);

            // assert_eq!(hit.normal, Vec3::new(0.0, 0.0, 1.0));
            let emmission = hit.material.color * hit.material.emmission;
            let reflecting_direction = (-ray.direction).reflect(hit.normal);
            if n == 0 {
                // terminate recursion. The path still ends on a lit surface,
                // so the emission collected here is what the eye should see.
                emmission
                /*scene.ambient.color
                * hit.material.color
                * Renderer::brdf(
                    reflecting_direction,
                    -ray.direction,
                    hit.normal,
                    hit.material.reflecting as f64,
                    hit.material.diffuse as f64,
                )
                + emmission*/
            } else {
                // const NUM_CASTS: usize = 10;
                let mut average_color = Color::BLACK;

                if hit.material.reflecting > 0.0 {
                    let r = Ray {
                        origin: pos,
                        direction: reflecting_direction,
                    };
                    let color_incoming = Renderer::cast(scene, &r, n - 1, rng);

                    // The BRDF already carries the reflectance, so the
                    // weight is just the cosine. The diffuse estimator next
                    // door keeps its own albedo and solid angle terms.
                    average_color += color_incoming
                        * Self::reflecting_brdf(
                            reflecting_direction,
                            -ray.direction,
                            hit.normal,
                            hit.material.reflecting as F,
                            hit.material.diffuse as F,
                        )
                        * (reflecting_direction * hit.normal);
                }

                if hit.material.diffuse > 0.0 {
                    let r = Ray {
                        origin: pos,
                        direction: Vec3::random_on_hemisphere(hit.normal, rng).normalized(),
                    };
                    let color_incoming = Renderer::cast(scene, &r, n - 1, rng);

                    average_color += color_incoming
                        * Self::diffuse_brdf(
                            r.direction,
                            -ray.direction,
                            hit.normal,
                            hit.material.reflecting as F,
                            hit.material.diffuse as F,
                        )
                        * (r.direction * hit.normal)
                        * hit.material.diffuse as F
                        * F::TAU;
                }

                /*let casts = [
                    Vec3::random_on_hemisphere(hit.normal)
                        .lerp(reflecting_direction, hit.material.reflecting as f64)
                        .normalized(),
                    // Vec3::random_on_hemisphere(hit.normal),
                    // Vec3::random_on_hemisphere(hit.normal),
                    // Vec3::random_on_hemisphere(hit.normal),
                    // Vec3::random_on_hemisphere(hit.normal),
                    // Vec3::random_on_hemisphere(hit.normal),
                    // Vec3::random_on_hemisphere(hit.normal),
                    // Vec3::random_on_hemisphere(hit.normal),
                    // Vec3::random_on_hemisphere(hit.normal),
                ];

                let dw = 2.0 * std::f64::consts::PI / (casts.len() as f64);

                for cast in casts {
                    /*let diffuse_direction = Vec3::random_on_hemisphere(hit.normal);
                    let r = Ray {
                        origin: pos,
                        direction: (diffuse_direction * (hit.material.diffuse as f64)
                            + reflecting_direction * (hit.material.reflecting as f64))
                            .normalized(),
                    };*/

                    // let attenuation = (ray.direction * ray.direction.reflect(hit.normal)).max(0.0);
                    // let hit_color = Renderer::cast(scene, &r, n - 1);
                    let r = Ray {
                        origin: pos,
                        direction: cast,
                    };
                    let color_incoming = Renderer::cast(scene, &r, n - 1);
                    average_color += color_incoming
                        * Renderer::brdf(
                            cast,
                            -ray.direction,
                            hit.normal,
                            hit.material.reflecting as f64,
                            hit.material.diffuse as f64,
                        )
                        * (cast * hit.normal)
                        * dw;
                }*/

                // let attenuation = primary_dir * r.direction;

                average_color * hit.material.color + emmission
            }
        } else {
            // direct ambient hit
            scene.ambient.color * scene.ambient.emmission
        }
    }

    /// Schlick Fresnel reflectance. `incoming` points away from the surface
    /// and `outgoing` points into it, so the angle of incidence is the
    /// absolute cosine between `outgoing` and the normal.
    pub fn reflecting_brdf(
        _incoming: Vec3,
        outgoing: Vec3,
        normal: Vec3,
        reflecting: F,
        _diffuse: F,
    ) -> F {
        let f0 = reflecting.clamp(0.0, 1.0);
        let cos_theta = (outgoing * normal).abs().clamp(0.0, 1.0);
        schlick_fresnel(f0, cos_theta)
    }

    pub fn diffuse_brdf(
        _incoming: Vec3,
        _outgoing: Vec3,
        _normal: Vec3,
        _reflecting: F,
        _diffuse: F,
    ) -> F {
        F::FRAC_1_PI
    }

    /*pub fn brdf(
        incoming: Vec3,
        outgoing: Vec3,
        normal: Vec3,
        reflecting: f64,
        diffuse: f64,
    ) -> f64 {
        // let direct = incoming.reflect(normal) * outgoing >= 0.8;
        /*let reflecting_component = if direct {
            // let fresnel = reflecting + (1.0 - reflecting)
            // println!("REFLECTING!");
            let val = 1.0 / ((incoming * outgoing).abs());
            assert!(!val.is_infinite());
            assert!(!val.is_nan());
            assert!(val > 0.0);
            val
        } else {
            0.0
        };*/
        let reflecting_component: f64 = 1.0 / ((incoming * outgoing).abs());
        let diffuse_component = /*if direct {
            0.0
        } else {*/
            1.0 / std::f64::consts::PI
        ;
        return reflecting_component * reflecting + diffuse_component * diffuse;
    }*/
}

/// Copy one finished tile into the image.
///
/// A tile that does not fit belongs to an earlier pass rendered at a larger
/// size. Dropping it keeps the stale pixels out of the image.
pub fn blit_chunk(image: &mut Image, chunk: &RenderChunk) {
    let (x0, y0) = chunk.start;
    let (x1, y1) = chunk.end;
    let (width, height) = (image.size[0], image.size[1]);

    if x0 >= width || y0 >= height || x1 >= width || y1 >= height {
        return;
    }

    let chunk_width = x1 - x0 + 1;
    for (row, y) in (y0..=y1).enumerate() {
        let src = row * chunk_width;
        let dst = y * width + x0;
        image.pixels[dst..dst + chunk_width].copy_from_slice(&chunk.pixels[src..src + chunk_width]);
    }
}

/// The scene the app builds and the benchmark renders: five spheres, five
/// planes, an emitter at the origin and a 512x512 camera. Kept in one place so
/// the GUI and the headless benchmark cannot drift apart.
pub fn demo_scene(num_samples: usize, num_bounces: usize) -> Scene {
    fn mat(r: F, g: F, b: F, e: F, refl: F, diff: F) -> Material {
        Material {
            color: Color::new(r, g, b, 1.0),
            emmission: e,
            reflecting: refl,
            diffuse: diff,
        }
    }

    Scene {
        camera: Camera {
            origin: Vec3::new([-5.0, 0.0, 0.0]),
            right: Vec3::new([0.0, -1.0, 0.0]),
            up: Vec3::new([0.0, 0.0, 1.0]),
            width: 0.5,
            height: 0.5,
        },
        spheres: vec![
            Sphere {
                id: 0,
                center: Vec3::new([-0.2, 0.1, 0.3]),
                radius: 0.25,
                material: mat(0.8, 0.2, 0.2, 0.0, 0.3, 0.7),
            },
            Sphere {
                id: 1,
                center: Vec3::new([0.3, -0.2, 0.1]),
                radius: 0.2,
                material: mat(0.2, 0.8, 0.3, 0.0, 0.5, 0.5),
            },
            Sphere {
                id: 2,
                center: Vec3::new([0.1, 0.35, -0.15]),
                radius: 0.18,
                material: mat(0.2, 0.3, 0.9, 0.0, 0.1, 0.9),
            },
            Sphere {
                id: 3,
                center: Vec3::new([-0.35, -0.3, -0.25]),
                radius: 0.22,
                material: mat(0.9, 0.9, 0.2, 0.0, 0.8, 0.2),
            },
            Sphere {
                id: 4,
                center: Vec3::new([0.0, 0.0, 0.0]),
                radius: 0.2,
                material: mat(1.0, 1.0, 1.0, 1.0, 0.5, 0.5),
            },
        ],
        planes: vec![
            Plane {
                id: 100,
                support: Vec3::new([0.0, 0.0, -0.5]),
                normal: Vec3::new([0.0, 0.0, 1.0]),
                material: mat(1.0, 1.0, 1.0, 0.0, 0.0, 1.0),
            },
            Plane {
                id: 101,
                support: Vec3::new([0.0, -1.0, 0.0]),
                normal: Vec3::new([0.0, 1.0, 0.0]),
                material: mat(1.0, 1.0, 1.0, 0.0, 0.0, 1.0),
            },
            Plane {
                id: 102,
                support: Vec3::new([1.0, 0.0, 0.0]),
                normal: Vec3::new([-1.0, 0.0, 0.0]),
                material: mat(1.0, 1.0, 1.0, 0.0, 0.0, 1.0),
            },
            Plane {
                id: 103,
                support: Vec3::new([-10.0, 0.0, 0.0]),
                normal: Vec3::new([1.0, 0.0, 0.0]),
                material: mat(1.0, 1.0, 1.0, 0.0, 0.0, 1.0),
            },
            Plane {
                id: 104,
                support: Vec3::new([0.0, 1.0, 0.0]),
                normal: Vec3::new([0.0, -1.0, 0.0]),
                material: mat(1.0, 1.0, 1.0, 0.0, 0.0, 1.0),
            },
        ],
        light: Vec3::new([0.0, 0.0, 0.0]),
        ambient: mat(0.1, 0.1, 0.1, 1.0, 0.0, 0.0),
        num_bounces,
        num_samples,
    }
}

impl Drop for Renderer {
    fn drop(&mut self) {
        // Closing the channels ends every worker loop. A join can still fail
        // if a worker panicked, and a panic inside a destructor aborts the
        // process, so the result is deliberately ignored.
        drop(self.send.take());
        drop(self.chunks.take());
        for worker in self.workers.drain(..) {
            let _ = worker.join();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    const LENGTHS: [usize; 17] = [
        1, 7, 8, 9, 63, 64, 100, 127, 128, 255, 511, 512, 513, 519, 1000, 1023, 1024,
    ];

    const CHUNK_SIZES: [usize; 7] = [1, 2, 3, 7, 8, 16, 64];

    #[test]
    fn chunk_ranges_tile_exactly() {
        for &len in &LENGTHS {
            for &size in &CHUNK_SIZES {
                let ranges = chunk_ranges(len, size);

                assert!(
                    !ranges.is_empty(),
                    "len {len} size {size} produced no ranges"
                );
                assert_eq!(ranges[0].0, 0, "len {len} size {size}");

                let mut next = 0;
                for &(start, end) in &ranges {
                    assert!(start <= end, "len {len} size {size}: empty range");
                    assert!(end < len, "len {len} size {size}: end {end} past len");
                    assert_eq!(
                        start, next,
                        "len {len} size {size}: gap or overlap at {start}"
                    );
                    next = end + 1;
                }

                assert_eq!(next, len, "len {len} size {size}: incomplete tiling");
                assert_eq!(
                    ranges.len(),
                    len.div_ceil(size),
                    "len {len} size {size}: wrong range count"
                );
            }
        }
    }

    #[test]
    fn chunk_ranges_of_empty_is_empty() {
        assert!(chunk_ranges(0, 8).is_empty());
    }

    #[test]
    fn chunk_ranges_keeps_the_tail_of_an_odd_length() {
        // 100 is 12 full chunks of 8 plus a 4 pixel tail, so the old
        // inclusive range dropped the last 4 rows and columns.
        let ranges = chunk_ranges(100, 8);
        assert_eq!(ranges.len(), 13);
        assert_eq!(*ranges.last().unwrap(), (96, 99));
    }

    /// Writing every tile must stay inside the image, for the sizes that
    /// used to overrun the buffer (511) and to leave black strips (100).
    #[test]
    fn tiles_write_in_bounds() {
        for &[width, height] in &[[511, 511], [100, 100]] {
            let mut image = Image::new([width, height]);
            let mut written = 0;

            for &(i0, i1) in &chunk_ranges(width, 8) {
                for &(j0, j1) in &chunk_ranges(height, 8) {
                    for j in j0..=j1 {
                        for i in i0..=i1 {
                            image[(i, j)] = RGBA::new(1, 2, 3, 255);
                            written += 1;
                        }
                    }
                }
            }

            assert_eq!(written, width * height, "{width}x{height}: coverage");
            assert_eq!(image.pixels.len(), width * height);
            assert!(
                image.pixels.iter().all(|p| *p == RGBA::new(1, 2, 3, 255)),
                "{width}x{height}: untouched pixel left behind"
            );
        }
    }

    /// A tile must own exactly the pixels it covers. 511 leaves a 7 pixel
    /// tail, so the old fixed 8 by 8 tile wrote 15 pixels past the buffer.
    #[test]
    fn a_tile_owns_exactly_its_own_pixels() {
        let tail = RenderWorkload::new((504, 504), (510, 510), [511, 511]);
        assert_eq!(tail.pixels.len(), 49);

        let last = RenderWorkload::new((96, 96), (99, 99), [100, 100]);
        assert_eq!(last.pixels.len(), 16);

        let full = RenderWorkload::new((0, 0), (7, 7), [100, 100]);
        assert_eq!(full.pixels.len(), 64);
    }

    /// An empty scene, so every camera ray misses and returns the ambient
    /// color. That makes one expected pixel value stand for "fully covered".
    fn test_scene(ambient: Color) -> Scene {
        Scene {
            camera: Camera {
                origin: Vec3::new([0.0, 0.0, 0.0]),
                right: Vec3::new([1.0, 0.0, 0.0]),
                up: Vec3::new([0.0, 1.0, 0.0]),
                width: 1.0,
                height: 1.0,
            },
            spheres: Vec::new(),
            planes: Vec::new(),
            light: Vec3::new([0.0, 0.0, 1.0]),
            ambient: Material {
                color: ambient,
                emmission: 1.0,
                reflecting: 0.0,
                diffuse: 0.0,
            },
            num_samples: 1,
            num_bounces: 1,
        }
    }

    /// Drain until the image is fully covered, or give up loudly.
    fn drain_until(renderer: &mut Renderer, expected: RGBA) {
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(30);
        loop {
            renderer.drain();
            if renderer.image.pixels.iter().all(|p| *p == expected) {
                return;
            }
            assert!(
                std::time::Instant::now() < deadline,
                "tiles never covered every pixel"
            );
            std::thread::sleep(std::time::Duration::from_millis(2));
        }
    }

    /// Every pixel of every tile must land in the image, at the sizes that
    /// used to overrun the buffer (511) or to leave black strips (100).
    #[test]
    fn render_covers_every_pixel() {
        let ambient = Color::rgb(0.5, 0.25, 0.125);
        let expected: RGBA = ambient.gamma().into();

        for &[width, height] in &[[100usize, 100usize], [511, 511], [1023, 1024]] {
            let mut renderer = Renderer::new(&test_scene(ambient));
            *renderer.size.lock() = [width, height];

            renderer.render(&test_scene(ambient));

            drain_until(&mut renderer, expected);
            assert_eq!(renderer.image.size, [width, height]);
            assert_eq!(renderer.image.pixels.len(), width * height);
        }
    }

    /// The old shared image was freed by the second click while the workers
    /// of the first pass still held references into it.
    #[test]
    fn a_second_render_does_not_trip_over_the_first() {
        let first = Color::rgb(0.5, 0.25, 0.125);
        let second = Color::rgb(0.125, 0.25, 0.5);
        let mut renderer = Renderer::new(&test_scene(first));
        *renderer.size.lock() = [100, 100];

        renderer.render(&test_scene(first));
        renderer.render(&test_scene(second));

        drain_until(&mut renderer, second.gamma().into());
    }

    #[test]
    fn drain_drops_a_tile_from_a_larger_pass() {
        let ambient = Color::rgb(0.5, 0.25, 0.125);
        let mut renderer = Renderer::new(&test_scene(ambient));
        *renderer.size.lock() = [8, 8];
        renderer.image = Image::new([8, 8]);

        // Cloned, so the borrow does not outlive the drain calls below.
        let tx = renderer.chunks.clone().unwrap();

        // Left over from a pass rendered at a larger size.
        tx.send(RenderChunk {
            start: (8, 8),
            end: (15, 15),
            pixels: vec![RGBA::new(1, 2, 3, 255); 64],
        })
        .unwrap();
        renderer.drain();
        assert!(renderer
            .image
            .pixels
            .iter()
            .all(|p| *p == RGBA::new(0, 0, 0, 255)));

        tx.send(RenderChunk {
            start: (0, 0),
            end: (1, 1),
            pixels: vec![RGBA::new(9, 9, 9, 255); 4],
        })
        .unwrap();
        renderer.drain();

        assert_eq!(renderer.image[(0, 0)], RGBA::new(9, 9, 9, 255));
        assert_eq!(renderer.image[(1, 1)], RGBA::new(9, 9, 9, 255));
        assert_eq!(renderer.image[(2, 2)], RGBA::new(0, 0, 0, 255));
    }

    /// The same workload as the headless `--bench` mode, driven through
    /// `cargo test` instead of the binary. Ignored by default because it is a
    /// benchmark, not a correctness test.
    ///
    /// Prefer the binary for measurement — it reports best/mean/median and
    /// can write the image out for a before/after diff:
    ///
    ///     cargo run --release -- --bench
    ///     cargo run --release -- --bench --serial --out before.ppm
    ///
    /// This form stays for `cargo test`-driven profiling:
    ///
    ///     cargo test --release -- --ignored --nocapture bench_render_action
    ///
    /// `BENCH_ITERS` and `BENCH_SAMPLES` tune the workload.
    #[test]
    #[ignore = "profiling harness; run explicitly with --ignored"]
    fn bench_render_action() {
        let samples: usize = std::env::var("BENCH_SAMPLES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(8);
        let iters: usize = std::env::var("BENCH_ITERS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(5);

        // The headless `--bench` mode renders the same scene, so the two
        // cannot drift apart.
        let scene = demo_scene(samples, 3);

        let [w, h] = [512usize, 512];
        let mut renderer = Renderer::new(&scene);
        *renderer.size.lock() = [w, h];

        // `render()` pushes every tile into the work channel synchronously, so
        // the total is known up front and counting returned tiles is an exact
        // completion signal. Pixel inspection is not: a ray that bottoms out on
        // a non-emissive surface legitimately resolves to pure black.
        //
        // This receives and blits in one step rather than calling `drain()`:
        // a cloned crossbeam receiver is a *second consumer*, so counting from
        // one while draining the other races and loses tiles.
        let expected_tiles = chunk_ranges(w, 8).len() * chunk_ranges(h, 8).len();
        let drain_full = |r: &mut Renderer| {
            let rx = r.rx.clone();
            let begin = std::time::Instant::now();
            let mut received = 0usize;
            let mut last_report = std::time::Instant::now();
            loop {
                while let Ok(chunk) = rx.try_recv() {
                    r.blit(&chunk);
                    received += 1;
                }
                if received >= expected_tiles {
                    eprintln!(
                        "BENCH   {received}/{expected_tiles} tiles in {:.3}s",
                        begin.elapsed().as_secs_f64()
                    );
                    return;
                }
                if last_report.elapsed() > std::time::Duration::from_secs(2) {
                    eprintln!(
                        "BENCH   progress {received}/{expected_tiles} after {:.1}s",
                        begin.elapsed().as_secs_f64()
                    );
                    last_report = std::time::Instant::now();
                }
                assert!(
                    begin.elapsed() < std::time::Duration::from_secs(300),
                    "only {received}/{expected_tiles} tiles came back"
                );
                std::thread::sleep(std::time::Duration::from_millis(1));
            }
        };

        // Warm up so the profile is not dominated by first-touch page faults.
        renderer.render(&scene);
        drain_full(&mut renderer);

        let start = std::time::Instant::now();
        for _ in 0..iters {
            renderer.render(&scene);
            drain_full(&mut renderer);
        }
        let elapsed = start.elapsed();

        let per = elapsed.as_secs_f64() / iters as f64;
        let rays = (w * h * samples * iters) as f64;
        eprintln!(
            "BENCH render_action {w}x{h} samples={samples} bounces={} iters={iters} \
             total={elapsed:?} per_render={per:.4}s throughput={:.2} Mrays/s",
            scene.num_bounces,
            rays / elapsed.as_secs_f64() / 1e6,
        );
        // Keep the result observable so the work cannot be optimized away.
        assert!(renderer
            .image
            .pixels
            .iter()
            .any(|p| *p != RGBA::new(0, 0, 0, 255)));
    }

    /// SHA-256 (FIPS 180-4), written out because the crate has no hashing
    /// dependency and no dependency may be added. `std::hash::Hasher` is not
    /// an option either: it is not stable across rustc versions, so the test
    /// would fail on a toolchain bump for no real reason.
    fn sha256(data: &[u8]) -> [u8; 32] {
        // Packed seven per line. rustfmt defaults to one per line here, and 64
        // lines of round constants bury the algorithm they belong to.
        #[rustfmt::skip]
        const K: [u32; 64] = [
            0x428a_2f98, 0x7137_4491, 0xb5c0_fbcf, 0xe9b5_dba5, 0x3956_c25b, 0x59f1_11f1,
            0x923f_82a4, 0xab1c_5ed5, 0xd807_aa98, 0x1283_5b01, 0x2431_85be, 0x550c_7dc3,
            0x72be_5d74, 0x80de_b1fe, 0x9bdc_06a7, 0xc19b_f174, 0xe49b_69c1, 0xefbe_4786,
            0x0fc1_9dc6, 0x240c_a1cc, 0x2de9_2c6f, 0x4a74_84aa, 0x5cb0_a9dc, 0x76f9_88da,
            0x983e_5152, 0xa831_c66d, 0xb003_27c8, 0xbf59_7fc7, 0xc6e0_0bf3, 0xd5a7_9147,
            0x06ca_6351, 0x1429_2967, 0x27b7_0a85, 0x2e1b_2138, 0x4d2c_6dfc, 0x5338_0d13,
            0x650a_7354, 0x766a_0abb, 0x81c2_c92e, 0x9272_2c85, 0xa2bf_e8a1, 0xa81a_664b,
            0xc24b_8b70, 0xc76c_51a3, 0xd192_e819, 0xd699_0624, 0xf40e_3585, 0x106a_a070,
            0x19a4_c116, 0x1e37_6c08, 0x2748_774c, 0x34b0_bcb5, 0x391c_0cb3, 0x4ed8_aa4a,
            0x5b9c_ca4f, 0x682e_6ff3, 0x748f_82ee, 0x78a5_636f, 0x84c8_7814, 0x8cc7_0208,
            0x90be_fffa, 0xa450_6ceb, 0xbef9_a3f7, 0xc671_78f2,
        ];

        let mut state: [u32; 8] = [
            0x6a09_e667,
            0xbb67_ae85,
            0x3c6e_f372,
            0xa54f_f53a,
            0x510e_527f,
            0x9b05_688c,
            0x1f83_d9ab,
            0x5be0_cd19,
        ];

        // Append 0x80, then zeroes, then the message length in bits as a big
        // endian u64, so the padded message is a whole number of 64 byte
        // blocks.
        let mut padded = Vec::with_capacity(data.len() + 72);
        padded.extend_from_slice(data);
        padded.push(0x80);
        while padded.len() % 64 != 56 {
            padded.push(0);
        }
        padded.extend_from_slice(&((data.len() as u64) * 8).to_be_bytes());

        for block in padded.chunks_exact(64) {
            let mut w = [0u32; 64];
            for (word, quad) in w.iter_mut().zip(block.chunks_exact(4)) {
                *word = u32::from_be_bytes([quad[0], quad[1], quad[2], quad[3]]);
            }
            for i in 16..64 {
                let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
                let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
                w[i] = w[i - 16]
                    .wrapping_add(s0)
                    .wrapping_add(w[i - 7])
                    .wrapping_add(s1);
            }

            let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = state;
            for i in 0..64 {
                let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
                let choose = (e & f) ^ (!e & g);
                let temp1 = h
                    .wrapping_add(s1)
                    .wrapping_add(choose)
                    .wrapping_add(K[i])
                    .wrapping_add(w[i]);
                let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
                let majority = (a & b) ^ (a & c) ^ (b & c);
                let temp2 = s0.wrapping_add(majority);

                h = g;
                g = f;
                f = e;
                e = d.wrapping_add(temp1);
                d = c;
                c = b;
                b = a;
                a = temp1.wrapping_add(temp2);
            }

            for (slot, word) in state.iter_mut().zip([a, b, c, d, e, f, g, h]) {
                *slot = slot.wrapping_add(word);
            }
        }

        let mut digest = [0u8; 32];
        for (chunk, word) in digest.chunks_exact_mut(4).zip(state) {
            chunk.copy_from_slice(&word.to_be_bytes());
        }
        digest
    }

    fn digest_hex(data: &[u8]) -> String {
        sha256(data)
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect()
    }

    /// The render, byte for byte, at the size and seed the `--bench --out`
    /// run in the failure message produces.
    const GOLDEN_SHA256: &str = "8187bc4541908aecaf0f28515cf17c192cfe944d3c2176001760a6642b672338";

    /// Every optimization round in this project is held to one bar: the
    /// render must not change by a single byte. That is checkable by hand,
    ///
    ///     light_transport --bench --serial --warmup 0 --iters 1 --size 48x48
    ///       --samples 2 --bounces 1 --seed 1 --out golden.ppm
    ///
    /// piped into `cmp`, and it is exactly the check that gets skipped when
    /// someone is in a hurry. This pins it to `cargo test` so it cannot be.
    ///
    /// The render and the serialisation below are the benchmark's serial path
    /// and its writer, copied rather than called: `bench::render_serial` and
    /// `bench::write_ppm` are private to that module, and a test that shared
    /// them would stop noticing a change to either one.
    ///
    /// The test has a sensitivity floor, and it sits under the resolution of
    /// an 8-bit byte. At 2 samples and 1 bounce a perturbation too small to
    /// move a single byte leaves the image identical and the test correctly
    /// stays green. Probing it with a change under that floor proves nothing:
    /// nudge the surface epsilon `0.9999 -> 0.99999` and the digest does not
    /// move. To show this test still bites, perturb coarsely, `0.9999 -> 0.95`.
    #[test]
    fn render_is_bit_identical_to_the_golden() {
        let scene = demo_scene(2, 1);
        let size = [48usize, 48];

        // The walk order and the per tile seeding are load bearing: they are
        // what makes the image independent of thread scheduling.
        let mut image = Image::new(size);
        for &(j0, j1) in &chunk_ranges(size[1], 8) {
            for &(i0, i1) in &chunk_ranges(size[0], 8) {
                let mut workload = RenderWorkload::new((i0, j0), (i1, j1), size);
                let mut rng = rand::rngs::StdRng::seed_from_u64(1);
                workload.handle(&scene, &mut rng);
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

        // Byte for byte what `write_ppm` puts in the file: a 13 byte header
        // and 48 * 48 * 3 bytes of colour, so 6925 bytes.
        let mut ppm = Vec::with_capacity(6925);
        ppm.extend_from_slice(format!("P6\n{} {}\n255\n", size[0], size[1]).as_bytes());
        let bytes = image.bytes();
        for y in 0..size[1] {
            for x in 0..size[0] {
                let p = (y * size[0] + x) * 4;
                ppm.extend_from_slice(&bytes[p..p + 3]);
            }
        }
        assert_eq!(
            ppm.len(),
            6925,
            "the serialised PPM is the wrong length, so it is not what \
             `write_ppm` produces and the digest below means nothing"
        );

        assert_eq!(
            digest_hex(&ppm),
            GOLDEN_SHA256,
            "the render changed. This is only correct if the change was intended to alter \
             the image; otherwise an optimization has changed a pixel. If it was intended, \
             re-derive the digest with:\n  \
             light_transport --bench --serial --warmup 0 --iters 1 --size 48x48 \
             --samples 2 --bounces 1 --seed 1 --out /tmp/golden.ppm\n  \
             and sha256sum /tmp/golden.ppm"
        );
    }
}
