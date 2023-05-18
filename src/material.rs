use crate::geometry::Local;
use crate::image::Color;
use crate::math::{Constants, Vec3, F};

fn schlick_fresnel(f0: F, cos_theta: F) -> F {
    f0 + (1.0 - f0) * (1.0 - cos_theta).powi(5)
}

fn phong_pdf(theta: F, phi: F, s: F) -> F {
    (s + 1.0) / (F::TAU) * theta.cos().powf(s)
}

fn phong_sample(theta: F, phi: F) -> F {
    todo!()
}

#[derive(Debug, Clone)]
pub struct Material {
    pub color: Color,
    pub emmission: f32,
    pub reflecting: f32,
    pub diffuse: f32,
}

pub trait BSDF {
    fn eval(&self, incoming: Vec3, outgoing: Vec3, local: Local) -> F;
    fn sample(&self, incoming: Vec3, local: Local) -> Vec3;
    // fn pdf(&self) -> ?;
}

pub struct Lambertian {
    pub albedo: F,
}

impl BSDF for Lambertian {
    fn eval(&self, incoming: Vec3, outgoing: Vec3, local: Local) -> F {
        self.albedo / F::PI * (incoming * local.normal)
    }

    fn sample(&self, incoming: Vec3, local: Local) -> Vec3 {
        Vec3::random_on_hemisphere(local.normal)
    }
}

pub struct Reflecting {}

impl BSDF for Reflecting {
    fn eval(&self, incoming: Vec3, outgoing: Vec3, local: Local) -> F {
        1.0
    }

    fn sample(&self, incoming: Vec3, local: Local) -> Vec3 {
        incoming.reflect(local.normal)
    }
}
