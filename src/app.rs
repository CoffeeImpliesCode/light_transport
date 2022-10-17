use std::sync::{atomic::Ordering, Arc};
use std::time::Instant;

use eframe::egui::{self, Sense};
use eframe::egui::{Separator, Slider, Visuals};
use eframe::epaint::color::Hsva;
use eframe::epaint::{Color32, ColorImage, Vec2};
use egui::mutex::Mutex;
use egui::Key;
use rand::Rng;

use crate::image::{Color, Image};
use crate::math::Vec3;
use crate::renderer::{Camera, Material, Plane, Renderer, Scene, Sphere};

pub struct LightTransport {
    pub renderer: Renderer,
    pub scene: Arc<Mutex<Scene>>,

    pub texture_handle: egui::TextureHandle,
    // pub rotating_triangle: Arc<Mutex<RotatingTriangle>>,
    // pub angle: f32,
    //
    pub begin: Instant,
    pub last: Instant,
}

impl LightTransport {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let gl = cc.gl.as_ref(); // .expect("Run with glow backend!");

        let image = ColorImage::new(
            [crate::IMAGE_WIDTH, crate::IMAGE_HEIGHT],
            Color32::from_rgb(0, 0, 0),
        );

        // let image = Image::random([1000, 1000]);
        let texture_handle = cc.egui_ctx.load_texture(
            "render",
            image, // egui::ColorImage::from_rgba_unmultiplied(image.dimension, image.bytes()),
            egui::TextureFilter::Nearest,
        );

        let mut rng = rand::thread_rng();

        let mut spheres: Vec<Sphere> = (0..4)
            .into_iter()
            .map(|i| Sphere {
                center: Vec3::new(
                    rng.gen_range(-0.5..0.5),
                    rng.gen_range(-0.5..0.5),
                    rng.gen_range(-0.5..0.5),
                ),
                radius: rng.gen_range(0.1..0.3),
                id: i,
                material: Material {
                    color: Color::new(
                        rng.gen_range(0.5..1.0),
                        rng.gen_range(0.5..1.0),
                        rng.gen_range(0.5..1.0),
                        1.0,
                    ),
                    emmission: 0.0,
                    reflecting: 0.0,
                    diffuse: 1.0,
                },
            })
            .collect();

        spheres.push(Sphere {
            center: Vec3::zero(),
            radius: 0.2,
            id: spheres.len(),
            material: Material {
                color: Color::RED + Color::GREEN,
                emmission: 1.0,
                diffuse: 0.5,
                reflecting: 0.5,
            },
        });

        let scene = Arc::new(Mutex::new(Scene {
            camera: Camera {
                origin: Vec3::new(-5.0, 0.0, 0.0),
                right: Vec3::new(0.0, -1.0, 0.0),
                up: Vec3::new(0.0, 0.0, 1.0),
                width: 0.5, //0.5 * (800.0 / 600.0),
                height: 0.5,
            },
            spheres,
            planes: vec![Plane {
                support: Vec3::new(0.0, 0.0, -0.5),
                normal: Vec3::new(0.0, 0.0, 1.0),
                id: 100,
                material: Material {
                    color: Color::new(1.0, 1.0, 1.0, 1.0),
                    emmission: 0.0,
                    diffuse: 0.0,
                    reflecting: 1.0,
                },
            }],
            light: Vec3::new(0.0, 0.0, 0.0),
            ambient: Material {
                color: Color::new(0.1, 0.1, 0.1, 1.0),
                diffuse: 0.0,
                reflecting: 0.0,
                emmission: 1.0,
            },
            num_bounces: 1,
            num_samples: 1,
        }));

        let renderer = Renderer::new(&*scene.lock());
        Self {
            renderer,
            texture_handle,
            scene,
            begin: Instant::now(),
            last: Instant::now(),
            // rotating_triangle: Arc::new(Mutex::new(RotatingTriangle::new(gl))),
            // angle: 0.0,
        }
    }

    /*fn custom_painting(&mut self, ui: &mut egui::Ui) {
        let (rect, response) =
            ui.allocate_exact_size(egui::Vec2::splat(1000.0), egui::Sense::drag());

        self.angle += response.drag_delta().x * 0.01;

        let angle = self.angle;
        let rotating_triangle = self.rotating_triangle.clone();

        /*let callback = egui::PaintCallback {
            rect,
            callback: Arc::new(egui_glow::CallbackFn::new(move |_info, painter| {
                rotating_triangle.lock().paint(painter.gl(), angle);
            })),
        };
        ui.painter().add(callback);

        */
    }*/
}

impl eframe::App for LightTransport {
    fn update(&mut self, ctx: &egui::Context, frame: &mut eframe::Frame) {
        let mut vis = Visuals::default();
        vis.dark_mode = true;
        ctx.set_visuals(vis);

        let secs_since_start = Instant::now().duration_since(self.begin).as_secs_f64();
        let secs_since_last = Instant::now().duration_since(self.last).as_secs_f64();

        /*if secs_since_start > 3.0 {
            frame.quit();
        }*/

        // puffin_egui::profiler_window(ctx);
        {
            let image = self.renderer.take_image();
            self.texture_handle.set(image, egui::TextureFilter::Linear);
            /*if self.renderer.send.as_ref().unwrap().is_empty() {
                self.renderer.render(&self.scene.lock());
            }*/
        }
        /*
        if self.renderer.present.swap(false, Ordering::SeqCst) {
            println!("PRESENT");
            if let Some(image) = self.renderer.take_image() {
                //egui::ColorImage::from_rgba_unmultiplied(image.size, image.pixels);
                self.texture_handle
                    .set(image /*, egui::TextureFilter::Nearest*/);
            } else {
            }
            self.renderer.render(&self.scene.lock());
        }*/

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.horizontal(|ui| {
                let drag_delta = ui
                    .image(self.texture_handle.id(), Vec2::new(1000.0, 1000.0))
                    .interact(Sense::drag())
                    .drag_delta();
                if ui.button("Render").clicked() {
                    self.renderer.render(&self.scene.lock());
                }
                egui::ScrollArea::vertical().show(ui, |ui| {
                    let mut scene = self.scene.lock();
                    let dx = drag_delta.x as f64 * 0.0005;
                    let dy = drag_delta.y as f64 * 0.0005;

                    let u = scene.camera.up;
                    let r = scene.camera.right;
                    // let d = u.cross(r).normalized();

                    let r = r * dx.cos() + u.cross(r) * dx.sin() + u * (u * r) * (1.0 - dx.cos());

                    scene.camera.right = r;

                    let u = u * dy.cos() + r.cross(u) * dy.sin() + r * (r * u) * (1.0 - dy.cos());

                    scene.camera.up = u;

                    let right = scene.camera.right;
                    let up = scene.camera.up;
                    let dir = up.cross(right).normalized();

                    {
                        let input = ctx.input();

                        if input.key_down(Key::W) {
                            scene.camera.origin += dir * 0.02;
                        }
                        if input.key_down(Key::S) {
                            scene.camera.origin -= dir * 0.02;
                        }
                        if input.key_down(Key::A) {
                            scene.camera.origin -= right * 0.02;
                        }
                        if input.key_down(Key::D) {
                            scene.camera.origin += right * 0.02;
                        }
                        if input.key_down(Key::Q) {
                            scene.camera.origin += up * 0.02;
                        }
                        if input.key_down(Key::E) {
                            scene.camera.origin -= up * 0.02;
                        }
                        if input.key_down(Key::Escape) {
                            frame.quit();
                        }
                        if input.key_down(Key::R) {
                            let mut lock = self.renderer.avg_rps.lock();
                            lock.0 = 0.0;
                            lock.1 = 0;
                        }
                    }

                    /*scene.light = Vec3::new(
                        -1.0 + secs_since_start.sin() * 0.3,
                        1.0 + secs_since_start.cos() * 0.3,
                        1.0,
                    );*/
                    ui.vertical(|ui| {
                        ui.horizontal(|ui| {
                            ui.label("Samples");
                            ui.add(Slider::new(&mut scene.num_samples, 1..=1000));
                        });
                        ui.horizontal(|ui| {
                            ui.label("bounces");
                            ui.add(Slider::new(&mut scene.num_bounces, 1..=10));
                        });
                        ui.horizontal(|ui| {
                            ui.label("camera_x:");
                            ui.add(Slider::new(&mut scene.camera.origin[0], -1.0..=1.0));
                        });
                        ui.horizontal(|ui| {
                            ui.label("camera_y:");
                            ui.add(Slider::new(&mut scene.camera.origin[1], -1.0..=1.0));
                        });
                        ui.horizontal(|ui| {
                            ui.label("camera_z:");
                            ui.add(Slider::new(&mut scene.camera.origin[2], -1.0..=1.0));
                        });

                        ui.label("scene color");
                        ui.add(Slider::new(&mut scene.ambient.color[0], 0.0..=1.0));
                        ui.add(Slider::new(&mut scene.ambient.color[1], 0.0..=1.0));
                        ui.add(Slider::new(&mut scene.ambient.color[2], 0.0..=1.0));
                        ui.label("scene emmission");
                        ui.add(Slider::new(&mut scene.ambient.emmission, 0.0..=1.0));

                        for (idx, sphere) in scene.spheres.iter_mut().enumerate() {
                            ui.add(Separator::default());
                            ui.label(format!("Sphere {}", idx));
                            ui.label("Position");
                            ui.add(Slider::new(&mut sphere.center[0], -1.0..=1.0));
                            ui.add(Slider::new(&mut sphere.center[1], -1.0..=1.0));
                            ui.add(Slider::new(&mut sphere.center[2], -1.0..=1.0));
                            ui.label("Radius");
                            ui.add(Slider::new(&mut sphere.radius, 0.0..=1.0));

                            ui.horizontal(|ui| {
                                ui.label("Color");
                                let mut c = sphere.material.color.as_rgb_slice_mut();
                                let mut color = Hsva::from_rgb(*c);
                                if egui::widgets::color_picker::color_picker_hsva_2d(
                                    ui,
                                    &mut color,
                                    egui::color_picker::Alpha::Opaque,
                                ) {
                                    *c = color.to_rgb();
                                }
                            });
                            ui.label("Diffuse | Reflecting | Emmission");
                            ui.add(Slider::new(&mut sphere.material.diffuse, 0.0..=1.0));
                            ui.add(Slider::new(&mut sphere.material.reflecting, 0.0..=1.0));
                            ui.add(Slider::new(&mut sphere.material.emmission, 0.0..=1.0));
                        }

                        for (idx, plane) in scene.planes.iter_mut().enumerate() {
                            ui.add(Separator::default());
                            ui.label(format!("Plane {}", idx));
                            ui.label("Diffuse | Reflecting | Emmission");
                            ui.add(Slider::new(&mut plane.material.diffuse, 0.0..=1.0));
                            ui.add(Slider::new(&mut plane.material.reflecting, 0.0..=1.0));
                            ui.add(Slider::new(&mut plane.material.emmission, 0.0..=1.0));
                        }
                    });
                });

                /*if self.renderer.send.clone().unwrap().is_empty() {
                    self.renderer.render(&self.scene.lock());
                }*/

                /*if self.renderer.running.load(Ordering::Relaxed) {
                    ui.add_enabled(false, egui::Button::new("Render"));
                    ui.label(format!(
                        "Rendering... ({}/{})",
                        self.renderer.stage.load(Ordering::Relaxed),
                        self.renderer.workload.load(Ordering::Relaxed)
                    ));
                } else {
                    self.renderer.render(&self.scene.lock());
                }*/ /*
                    else if ui.add(egui::Button::new("Render")).clicked() {
                        self.renderer.render();
                    }*/
            });

            self.last = Instant::now();
        });
        ctx.request_repaint();
    }

    fn on_exit(&mut self, gl: Option<&glow::Context>) {}
}
