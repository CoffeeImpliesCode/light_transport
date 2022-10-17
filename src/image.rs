use core::slice;
use std::{
    cell::UnsafeCell,
    ops::{Add, AddAssign, Deref, DerefMut, Index, IndexMut, Mul, MulAssign, Sub, SubAssign},
    sync::Arc,
};

use eframe::epaint::Color32;
use rand::Rng;

pub type GlobalImage = Arc<Image>;
pub type GlobalMutImage = Arc<UnsafeCell<Image>>;

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

    #[inline(always)]
    pub fn hsl(h: f32, g: f32, b: f32) -> Self {
        unimplemented!()
    }

    #[inline(always)]
    pub fn hsla(h: f32, g: f32, b: f32, a: f32) -> Self {
        unimplemented!()
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
        let mut rng = rand::thread_rng();
        let mut data: Vec<RGBA> = Vec::with_capacity(dimension[0] * dimension[1]);

        for _ in 0..dimension[0] * dimension[1] {
            data.push(RGBA::rgb(rng.gen(), rng.gen(), rng.gen()));
        }

        Self {
            size: dimension,
            pixels: data,
        }
    }

    pub fn color_32(&self) -> &[Color32] {
        unsafe {
            std::slice::from_raw_parts(self.pixels.as_ptr() as *const Color32, self.pixels.len())
        }
    }

    pub fn color_32_mut(&mut self) -> &mut [Color32] {
        unsafe {
            std::slice::from_raw_parts_mut(
                self.pixels.as_mut_ptr() as *mut Color32,
                self.pixels.len(),
            )
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
