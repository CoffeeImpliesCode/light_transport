use num_cpus;
use std::{
    cell::UnsafeCell,
    sync::{
        atomic::{AtomicBool, AtomicUsize, Ordering},
        Arc, RwLock,
    },
    thread::JoinHandle,
};

use std::cmp::PartialOrd;

use eframe::{
    egui,
    epaint::{Color32, ColorImage},
};
use egui::mutex::Mutex;

use crate::Image;
use crate::{
    image::{Color, GlobalMutImage},
    DEFAULT_IMAGE_HEIGHT, DEFAULT_IMAGE_WIDTH,
};

use crate::material::Material;
use crate::math::{Constants, Vec3, F};

use crate::geometry::*;
pub struct Renderer {
    pub present: Arc<AtomicBool>,
    pub running: Arc<AtomicBool>,
    pub stage: Arc<AtomicUsize>,
    pub workload: Arc<AtomicUsize>,
    pub scene: Arc<RwLock<Scene>>,
    pub workers: Vec<JoinHandle<()>>,
    pub send: Option<crossbeam::channel::Sender<RenderWorkload>>,
    // pub image: Arc<Mutex<Option<ColorImage>>>,
    pub image: GlobalMutImage,
    pub size: Arc<Mutex<[usize; 2]>>,
    pub avg_rps: Arc<Mutex<(f64, usize)>>,
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
#[derive(Debug)]
pub struct RenderWorkload {
    pub image: crate::image::GlobalMutImage,
    pub start: (usize, usize),
    pub end: (usize, usize),
}

unsafe impl Send for RenderWorkload {}

impl RenderWorkload {
    pub fn handle<'a>(&self, scene: &'a Scene) -> usize {
        let (direction, dx, dy) = {
            let direction = scene.camera.up.cross(scene.camera.right).normalized();
            let dx = scene.camera.right * scene.camera.width * 0.5;
            let dy = scene.camera.up * scene.camera.height * 0.5;
            (direction, dx, dy)
        };

        let pre_render: crate::image::RGBA = Color::new(0.1, 0.1, 0.1, 1.0).into();

        let image = unsafe { &mut *self.image.get() as &'a mut Image };

        let camera_dx = 2.0 / image.size[0] as F;
        let camera_dy = 2.0 / image.size[1] as F;

        for j in self.start.1..=self.end.1 {
            for i in self.start.0..=self.end.0 {
                image[(i, j)] = pre_render;
            }
        }

        for i in self.start.0..=self.end.0 {
            image[(i, self.start.1)] = Color::BLACK.into();
            image[(i, self.end.1)] = Color::BLACK.into();
        }

        for j in self.start.1..=self.end.1 {
            image[(self.start.0, j)] = Color::BLACK.into();
            image[(self.end.0, j)] = Color::BLACK.into();
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
                    let c = Renderer::cast(&scene, &ray, scene.num_bounces);
                    color += c;
                    cast_rays += 1;
                }
                color *= 1.0 / (scene.num_samples as F);
                image[(i, j)] = color.gamma().into();
                // stage.store(j * image.dimension[1] + i, Ordering::Release);
            }
        }
        return cast_rays;
    }
}

impl Renderer {
    pub fn new(scene: &Scene) -> Self {
        let (s, r) = crossbeam::channel::unbounded::<RenderWorkload>();

        let present = Arc::new(AtomicBool::new(true));
        let running = Arc::new(AtomicBool::new(false));
        let workload = Arc::new(AtomicUsize::new(0));
        let stage = Arc::new(AtomicUsize::new(0));
        let scene = Arc::new(RwLock::new(scene.clone()));
        let image = Arc::new(UnsafeCell::new(Image::new([
            crate::DEFAULT_IMAGE_WIDTH,
            crate::DEFAULT_IMAGE_HEIGHT,
        ])));
        let size = Arc::new(Mutex::new([
            crate::DEFAULT_IMAGE_WIDTH,
            crate::DEFAULT_IMAGE_HEIGHT,
        ]));

        let avg_rps = Arc::new(Mutex::new((0.0 as f64, 0 as usize)));

        let workers: Vec<_> = (0..num_cpus::get())
            .into_iter()
            .map(|worker| {
                let r = r.clone();
                let present = present.clone();
                // let image = image.clone();
                let scene = scene.clone();
                let running = running.clone();
                let avg_rps = avg_rps.clone();
                // let stage = stage.clone();

                std::thread::spawn(move || {
                    for workload in r.iter() {
                        let start = std::time::Instant::now();
                        present.store(false, Ordering::SeqCst);
                        // running.store(true, Ordering::Relaxed);

                        let scene = scene.read().unwrap();

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
                        let num_cast_rays = workload.handle(&scene);
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

                        // *image.lock() = Some(img);

                        let end = std::time::Instant::now();
                        // stage.store(2, Ordering::Release);
                        // stage.store(3, Ordering::Release);
                        present.store(true, Ordering::SeqCst);
                        // running.store(false, Ordering::Release);

                        let millis = ((end - start).as_micros() as f64) / 1000.0;
                        let secs = millis / 1000.0;
                        let rps = (num_cast_rays as f64) / secs;

                        {
                            let mut lock = avg_rps.lock();
                            let new_count = lock.1 + 1;
                            let new_avg_rps = lock.0 + (rps - lock.0) / (new_count as f64);
                            lock.0 = new_avg_rps;
                            lock.1 = new_count;
                        }

                        /*println!(
                            "Worker {:>2}: cast {:>8} rays in {:>8.2}ms ({:>14.2} rps)",
                            worker, num_cast_rays, millis, rps
                        );*/
                        /*println!(
                            "Worker {}: Chunk ({} {})..({} {}) rendered in {:.3}ms",
                            worker,
                            workload.start.0,
                            workload.start.1,
                            workload.end.0,
                            workload.end.1,
                            ((end - start).as_micros() as f32) / 1000.0
                        );*/
                    }
                })
            })
            .collect();

        Renderer {
            image,
            size,
            present,
            running,
            workload,
            stage,
            scene,
            workers,
            avg_rps,
            send: Some(s),
        }
    }

    pub fn render(&mut self, scene: &Scene) {
        // println!("Starting Render");
        /*{
            let lock = self.avg_rps.lock();
            println!("Avg rps per thread: {} ({} samples)", lock.0, lock.1);
        }*/

        unsafe {
            *(self.image.get()) = Image::new(*self.size.lock());
        }

        match self.scene.write() {
            Ok(mut lock) => {
                *lock = scene.clone();
            }
            Err(e) => println!("Error: {:?}", e),
        }

        // let image = Arc::new(UnsafeCell::new(Image::new([0, 0])));

        let mut workloads = Vec::with_capacity(10);
        let (width, height) = {
            let lock = self.size.lock();
            (lock[0], lock[1])
        };

        const SIZE_X: usize = 8;
        const SIZE_Y: usize = 8;
        for j in (0..=height).array_chunks::<SIZE_Y>() {
            for i in (0..=width).array_chunks::<SIZE_X>() {
                workloads.push(RenderWorkload {
                    image: self.image.clone(),
                    start: (*i.first().unwrap(), *j.first().unwrap()),
                    end: (*i.last().unwrap(), *j.last().unwrap()),
                });
            }
        }

        // println!("Workloads: {:#?}", workloads);
        // panic!();

        // self.send.clone().unwrap().try_send(()).unwrap();
        let s = self.send.as_ref().unwrap();

        workloads.into_iter().for_each(|w| s.try_send(w).unwrap());
    }

    pub fn take_image(&mut self) -> ColorImage {
        unsafe {
            let image = self.image.as_ref().get();
            ColorImage::from_rgba_unmultiplied((*image).size, (*image).bytes())
        }
    }

    #[inline(always)]
    pub fn cast(scene: &Scene, ray: &Ray, n: usize) -> Color {
        if let Some(closest) = scene.intersect(ray) {
            let pos = ray.origin + ray.direction * closest.distance * 0.9999;
            let hit = scene.get_info(ray, &closest);

            // assert_eq!(hit.normal, Vec3::new(0.0, 0.0, 1.0));
            let emmission = hit.material.color * hit.material.emmission;
            let reflecting_direction = (-ray.direction).reflect(hit.normal);
            if n == 0 {
                // terminate recursion
                // hit.material.color * scene.ambient.color * scene.ambient.emmission + emmission
                Color::BLACK
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
                    let color_incoming = Renderer::cast(scene, &r, n - 1);

                    average_color += color_incoming
                        * Self::reflecting_brdf(
                            reflecting_direction,
                            -ray.direction,
                            hit.normal,
                            hit.material.reflecting as F,
                            hit.material.diffuse as F,
                        )
                        * (reflecting_direction * hit.normal)
                        * hit.material.reflecting as F
                        * F::TAU;
                }

                if hit.material.diffuse > 0.0 {
                    let r = Ray {
                        origin: pos,
                        direction: Vec3::random_on_hemisphere(hit.normal).normalized(),
                    };
                    let color_incoming = Renderer::cast(scene, &r, n - 1);

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

    pub fn reflecting_brdf(
        incoming: Vec3,
        outgoing: Vec3,
        normal: Vec3,
        reflecting: F,
        diffuse: F,
    ) -> F {
        if incoming.reflect(normal) * outgoing >= 0.99 {
            1.0
        } else {
            0.0
        }
    }

    pub fn diffuse_brdf(
        incoming: Vec3,
        outgoing: Vec3,
        normal: Vec3,
        reflecting: F,
        diffuse: F,
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

impl Drop for Renderer {
    fn drop(&mut self) {
        drop(self.send.take());
        println!("Dropping renderer.");
        self.workers.drain(..).for_each(|w| {
            w.join().unwrap();
        });
        println!("Done.");
    }
}
