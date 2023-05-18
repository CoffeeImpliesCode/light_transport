#![feature(iter_array_chunks)]
#![allow(dead_code)]

use std::io::Write;

use eframe::egui;

mod app;
mod geometry;
mod image;
mod material;
mod math;
mod renderer;

use app::LightTransport;
use image::Image;
use math::Vec3;
use renderer::Renderer;

const DEFAULT_IMAGE_WIDTH: usize = 512;
const DEFAULT_IMAGE_HEIGHT: usize = 512;

fn main() {
    let mut random_on_sphere = std::fs::File::create("tests/random_on_hemisphere.csv").unwrap();
    let norm = Vec3::new([0.0, 1.0, 1.0]);
    writeln!(&mut random_on_sphere, "X,Y,Z").unwrap();
    for _ in 0..1000 {
        let v = Vec3::random_on_hemisphere(norm);
        writeln!(&mut random_on_sphere, "{},{},{}", v[0], v[1], v[2]).unwrap();
    }

    let up = Vec3::new([0.0, 0.0, 1.0]);

    let outer = 10;
    let inner = 100000;

    /*let samples = (0..outer)
            .into_iter()
            .map(|_| {
                let outgoing = Vec3::random_on_hemisphere(up);
                let dw = (2.0 * std::f64::consts::PI) / (inner as f64);
                let brdf_samples = (0..inner)
                    .into_iter()
                    .map(|_| {
                        let incoming = Vec3::random_on_hemisphere(up);
                        Renderer::brdf(incoming, outgoing, up, 1.0, 0.0) * (incoming * up) * dw
                    })
                    .collect::<Vec<_>>();
                brdf_samples.iter().sum::<f64>() / (outer as f64)
            })
            .collect::<Vec<_>>();

        for sample in &samples {
            // assert!(*sample <= 1.0);
        }

        println!("{:#?}", samples);
    */
    let options = eframe::NativeOptions {
        initial_window_size: Some(egui::vec2(1300.0, 1020.0)),
        multisampling: 0,
        vsync: false,

        // renderer: eframe::Renderer::Glow,
        ..Default::default()
    };

    eframe::run_native(
        "Light Transport",
        options,
        Box::new(|cc| Box::new(LightTransport::new(cc))),
    );
}
