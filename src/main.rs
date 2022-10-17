#![feature(iter_array_chunks)]
#![allow(dead_code)]

use eframe::egui;

mod app;
mod image;
mod math;
mod renderer;

use app::LightTransport;
use image::Image;
use math::Vec3;
use renderer::Renderer;

const IMAGE_WIDTH: usize = 256;
const IMAGE_HEIGHT: usize = 256;

fn main() {
    let up = Vec3::new(0.0, 0.0, 1.0);

    let outer = 10;
    let inner = 100000;

    let samples = (0..outer)
        .into_iter()
        .map(|_| {
            let outgoing = Vec3::random_on_hemisphere(up);
            let dw = (2.0 * std::f64::consts::PI) / (inner as f64);
            let brdf_samples = (0..inner)
                .into_iter()
                .map(|_| {
                    let incoming = Vec3::random_on_hemisphere(up);
                    Renderer::brdf(incoming, outgoing, up, 1.0, 0.0, dw) * (incoming * up) * dw
                })
                .collect::<Vec<_>>();
            brdf_samples.iter().sum::<f64>() / (outer as f64)
        })
        .collect::<Vec<_>>();

    for sample in &samples {
        // assert!(*sample <= 1.0);
    }

    println!("{:#?}", samples);

    let options = eframe::NativeOptions {
        initial_window_size: Some(egui::vec2(1300.0, 1020.0)),
        multisampling: 0,
        vsync: false,

        // renderer: eframe::Renderer::Glow,
        ..Default::default()
    };

    puffin::set_scopes_on(false);

    eframe::run_native(
        "Light Transport",
        options,
        Box::new(|cc| Box::new(LightTransport::new(cc))),
    );
}
