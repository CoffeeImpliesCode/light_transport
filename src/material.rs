use crate::image::Color;
use crate::math::F;

pub(crate) fn schlick_fresnel(f0: F, cos_theta: F) -> F {
    f0 + (1.0 - f0) * (1.0 - cos_theta).powi(5)
}

#[derive(Debug, Clone)]
pub struct Material {
    pub color: Color,
    pub emmission: f32,
    pub reflecting: f32,
    pub diffuse: f32,
}
