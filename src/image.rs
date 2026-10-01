use core::slice;
use std::{
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign},
    sync::Arc,
};

use eframe::epaint::Color32;
use rand::RngExt;

pub type GlobalImage = Arc<Image>;

// The render path hands the RGBA buffer to egui as raw bytes, so that view is
// only sound while RGBA is four tightly packed, byte-aligned bytes. Pin the
// layout at compile time: a change becomes a build error instead of undefined
// behavior. `Color32` is no longer reinterpreted — it is `align(4)` since
// ecolor 0.36, so the two layouts are no longer compatible.
const _: () = assert!(std::mem::size_of::<RGBA>() == 4);
const _: () = assert!(std::mem::align_of::<RGBA>() == 1);
const _: () = assert!(std::mem::size_of::<Color32>() == 4);

#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(transparent)]
pub struct Color([f32; 4]);
impl Color {
    pub const RED: Color = Color([1.0, 0.0, 0.0, 1.0]);
    pub const GREEN: Color = Color([0.0, 1.0, 0.0, 1.0]);
    pub const BLUE: Color = Color([0.0, 0.0, 1.0, 1.0]);
    pub const BLACK: Color = Color([0.0, 0.0, 0.0, 1.0]);
    pub const WHITE: Color = Color([1.0, 1.0, 1.0, 1.0]);

    #[inline(always)]
    pub fn new(r: f32, g: f32, b: f32, a: f32) -> Self {
        Color([r, g, b, a])
    }

    #[inline(always)]
    pub fn rgb(r: f32, g: f32, b: f32) -> Self {
        Color([r, g, b, 1.0])
    }

    /// Opaque color from a hue in [0, 1), a saturation in [0, 1] and a
    /// lightness in [0, 1]. The hue wraps, so 1.0 is the same as 0.0.
    #[inline(always)]
    pub fn hsl(h: f32, g: f32, b: f32) -> Self {
        Self::hsla(h, g, b, 1.0)
    }

    /// Same as [`Color::hsl`], with `a` as the opacity in [0, 1].
    #[inline(always)]
    pub fn hsla(h: f32, g: f32, b: f32, a: f32) -> Self {
        let s = g.clamp(0.0, 1.0);
        let l = b.clamp(0.0, 1.0);
        let hue = h.rem_euclid(1.0) * 6.0;

        let c = (1.0 - (2.0 * l - 1.0).abs()) * s;
        let x = c * (1.0 - (hue % 2.0 - 1.0).abs());
        let m = l - c / 2.0;

        // `hue` is in [0, 6), so the cast truncates to a sector in 0..=5.
        let (r, gg, bb) = match hue as u32 {
            0 => (c, x, 0.0),
            1 => (x, c, 0.0),
            2 => (0.0, c, x),
            3 => (0.0, x, c),
            4 => (x, 0.0, c),
            _ => (c, 0.0, x),
        };

        // Rounding can push a channel a hair outside [0, 1] at the corners
        // of the cube, so clamp rather than hand out a negative color.
        Color::new(
            (r + m).clamp(0.0, 1.0),
            (gg + m).clamp(0.0, 1.0),
            (bb + m).clamp(0.0, 1.0),
            a,
        )
    }

    /// Apply the gamma curve to the color channels. Alpha is opacity, not
    /// radiance, so it passes through untouched.
    #[inline(always)]
    pub fn gamma(&self) -> Self {
        Color([
            self[0].powf(2.2),
            self[1].powf(2.2),
            self[2].powf(2.2),
            self[3],
        ])
    }

    /// Inverse of [`Color::gamma`]. Alpha passes through untouched.
    #[inline(always)]
    pub fn ungamma(&self) -> Self {
        Color([
            self[0].powf(1.0 / 2.2),
            self[1].powf(1.0 / 2.2),
            self[2].powf(1.0 / 2.2),
            self[3],
        ])
    }

    pub fn as_rgb_slice(&self) -> &[f32; 3] {
        self.0[..3].try_into().unwrap()
    }

    pub fn as_rgb_slice_mut(&mut self) -> &mut [f32; 3] {
        (&mut self.0[..3]).try_into().unwrap()
    }
}

impl Add<Color> for Color {
    type Output = Color;

    #[inline(always)]
    fn add(self, other: Color) -> Self::Output {
        return Color::new(
            self[0] + other[0],
            self[1] + other[1],
            self[2] + other[2],
            self[3] + other[3],
        );
    }
}

impl AddAssign<Color> for Color {
    #[inline(always)]
    fn add_assign(&mut self, other: Color) {
        self[0] += other[0];
        self[1] += other[1];
        self[2] += other[2];
        self[3] += other[3];
    }
}

impl Sub<Color> for Color {
    type Output = Color;

    #[inline(always)]
    fn sub(self, other: Color) -> Self::Output {
        return Color::new(
            self[0] - other[0],
            self[1] - other[1],
            self[2] - other[2],
            self[3] - other[3],
        );
    }
}

impl SubAssign<Color> for Color {
    #[inline(always)]
    fn sub_assign(&mut self, other: Color) {
        self[0] -= other[0];
        self[1] -= other[1];
        self[2] -= other[2];
        self[3] -= other[3];
    }
}

impl Mul<f32> for Color {
    type Output = Color;

    #[inline(always)]
    fn mul(self, scale: f32) -> Self::Output {
        return Color::new(self[0] * scale, self[1] * scale, self[2] * scale, self[3]);
    }
}

impl MulAssign<f32> for Color {
    #[inline(always)]
    fn mul_assign(&mut self, scale: f32) {
        self[0] *= scale;
        self[1] *= scale;
        self[2] *= scale;
    }
}

impl Mul<f64> for Color {
    type Output = Color;

    #[inline(always)]
    fn mul(self, scale: f64) -> Self::Output {
        let scale = scale as f32;
        return Color::new(self[0] * scale, self[1] * scale, self[2] * scale, self[3]);
    }
}

impl MulAssign<f64> for Color {
    #[inline(always)]
    fn mul_assign(&mut self, scale: f64) {
        let scale = scale as f32;
        self[0] *= scale;
        self[1] *= scale;
        self[2] *= scale;
    }
}

impl Mul<Color> for Color {
    type Output = Color;

    #[inline(always)]
    fn mul(self, other: Color) -> Self::Output {
        Color::new(
            self[0] * other[0],
            self[1] * other[1],
            self[2] * other[2],
            self[3],
        )
    }
}

impl Index<usize> for Color {
    type Output = f32;

    #[inline(always)]
    fn index(&self, idx: usize) -> &f32 {
        &self.0[idx]
    }
}

impl IndexMut<usize> for Color {
    #[inline(always)]
    fn index_mut(&mut self, idx: usize) -> &mut f32 {
        &mut self.0[idx]
    }
}

impl Into<Color32> for Color {
    #[inline(always)]
    fn into(self) -> Color32 {
        Color32::from_rgb(
            (self[0] * 255.0) as u8,
            (self[1] * 255.0) as u8,
            (self[2] * 255.0) as u8,
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(transparent)]
pub struct RGBA([u8; 4]);

impl RGBA {
    const BLACK: RGBA = RGBA([0, 0, 0, 255]);
    const WHITE: RGBA = RGBA([255, 255, 255, 255]);
    const TRANSPARENT: RGBA = RGBA([0, 0, 0, 0]);

    pub fn new(r: u8, g: u8, b: u8, a: u8) -> Self {
        RGBA([r, g, b, a])
    }

    pub fn rgb(r: u8, g: u8, b: u8) -> Self {
        RGBA([r, g, b, 255])
    }
}

impl From<[u8; 4]> for RGBA {
    fn from(f: [u8; 4]) -> Self {
        RGBA(f)
    }
}

impl From<&[u8]> for RGBA {
    fn from(s: &[u8]) -> Self {
        RGBA(s.try_into().unwrap())
    }
}

impl From<Color> for RGBA {
    #[inline(always)]
    fn from(c: Color) -> Self {
        RGBA::new(
            (c[0].min(1.0) * 255.0) as u8,
            (c[1].min(1.0) * 255.0) as u8,
            (c[2].min(1.0) * 255.0) as u8,
            255, // (c[3] * 255.0) as u8,
        )
    }
}

pub struct Image {
    pub size: [usize; 2],
    pub pixels: Vec<RGBA>,
}

impl Image {
    pub fn new(dimension: [usize; 2]) -> Self {
        let mut data: Vec<RGBA> = Vec::with_capacity(dimension[0] * dimension[1]);
        for _ in 0..dimension[0] * dimension[1] {
            data.push(RGBA::BLACK);
        }

        Self {
            size: dimension,
            pixels: data,
        }
    }

    pub fn random(dimension: [usize; 2]) -> Self {
        let mut rng = rand::rng();
        let mut data: Vec<RGBA> = Vec::with_capacity(dimension[0] * dimension[1]);

        for _ in 0..dimension[0] * dimension[1] {
            data.push(RGBA::rgb(rng.random(), rng.random(), rng.random()));
        }

        Self {
            size: dimension,
            pixels: data,
        }
    }

    pub fn bytes(&self) -> &[u8] {
        unsafe {
            std::slice::from_raw_parts(self.pixels.as_ptr() as *const u8, self.pixels.len() * 4)
        }
    }

    pub fn bytes_mut(&mut self) -> &mut [u8] {
        unsafe {
            std::slice::from_raw_parts_mut(
                self.pixels.as_mut_ptr() as *mut u8,
                self.pixels.len() * 4,
            )
        }
    }
}

impl Index<(usize, usize)> for Image {
    type Output = RGBA;
    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let offset = index.1 * self.size[0] + index.0;
        &self.pixels[offset]
    }
}

impl IndexMut<(usize, usize)> for Image {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let offset = index.1 * self.size[0] + index.0;
        &mut self.pixels[offset]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: Color, b: Color) -> bool {
        (0..4).all(|i| (a[i] - b[i]).abs() < 1e-5)
    }

    #[test]
    fn gamma_leaves_alpha_alone() {
        let c = Color::new(0.5, 0.5, 0.5, 0.25);
        assert!(close(
            c.gamma(),
            Color::new(0.5f32.powf(2.2), 0.5f32.powf(2.2), 0.5f32.powf(2.2), 0.25)
        ));
        assert!(close(c.gamma().ungamma(), Color::new(0.5, 0.5, 0.5, 0.25)));
    }

    #[test]
    fn hsl_builds_primaries_and_wraps() {
        assert!(close(Color::hsl(0.0, 1.0, 0.5), Color::RED));
        assert!(close(Color::hsl(1.0 / 3.0, 1.0, 0.5), Color::GREEN));
        assert!(close(Color::hsl(2.0 / 3.0, 1.0, 0.5), Color::BLUE));
        assert!(close(Color::hsl(1.0, 1.0, 0.5), Color::RED));
        assert!(close(Color::hsl(0.0, 0.0, 0.0), Color::BLACK));
        assert!(close(Color::hsl(0.0, 0.0, 1.0), Color::WHITE));
    }

    #[test]
    fn hsla_keeps_opacity() {
        let c = Color::hsla(0.0, 1.0, 0.5, 0.5);
        assert!(close(c, Color::new(1.0, 0.0, 0.0, 0.5)));
    }

    #[test]
    fn hsl_stays_in_the_unit_cube() {
        for i in 0..=60 {
            for j in 0..=20 {
                for k in 0..=20 {
                    let c = Color::hsl(i as f32 / 60.0, j as f32 / 20.0, k as f32 / 20.0);
                    for ch in 0..4 {
                        assert!(
                            (0.0..=1.0).contains(&c[ch]),
                            "channel {ch} out of range: {} at {c:?}",
                            c[ch]
                        );
                    }
                }
            }
        }
    }
}
