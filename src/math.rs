use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
};

use rand::Rng;

#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(transparent)]
pub struct Vec3([f64; 3]);

impl Vec3 {
    #[inline(always)]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3([x, y, z])
    }

    #[inline(always)]
    pub fn zero() -> Self {
        Vec3([0.0, 0.0, 0.0])
    }

    #[inline(always)]
    pub fn x(&self) -> f64 {
        self.0[0]
    }

    #[inline(always)]
    pub fn y(&self) -> f64 {
        self.0[1]
    }

    #[inline(always)]
    pub fn z(&self) -> f64 {
        self.0[2]
    }

    #[inline(always)]
    pub fn dot(&self, other: Vec3) -> f64 {
        self[0] * other[0] + self[1] * other[1] + self[2] * other[2]
    }

    #[inline(always)]
    pub fn len(&self) -> f64 {
        self.dot(*self).sqrt()
    }

    #[inline(always)]
    pub fn len_sq(&self) -> f64 {
        self.dot(*self)
    }

    #[inline(always)]
    pub fn normalize(&mut self) {
        *self *= 1.0 / self.len();
    }

    #[inline(always)]
    pub fn normalized(&self) -> Vec3 {
        *self / self.len()
    }

    #[inline(always)]
    pub fn cross(&self, other: Vec3) -> Vec3 {
        Vec3::new(
            self[1] * other[2] - self[2] * other[1],
            self[2] * other[0] - self[0] * other[2],
            self[0] * other[1] - self[1] * other[0],
        )
    }

    #[inline(always)]
    pub fn random_on_sphere() -> Vec3 {
        let mut rng = rand::thread_rng();
        let theta = rng.gen_range(0.0..2.0 * std::f64::consts::PI);
        let phi = rng.gen_range(-1.0f64..1.0f64).acos();
        let sin_theta = theta.sin();
        let sin_phi = phi.sin();
        let x = sin_phi * theta.cos();
        let y = sin_phi * sin_theta;
        let z = phi.cos();

        Vec3::new(x, y, z)
    }

    #[inline(always)]
    pub fn from_spherical(r: f64, theta: f64, phi: f64) -> Vec3 {
        unimplemented!()
    }

    #[inline(always)]
    pub fn from_spherical_unit(theta: f64, phi: f64) -> Vec3 {
        let sin_theta = theta.sin();
        let cos_theta = theta.cos();
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();

        Vec3::new(sin_phi * cos_theta, sin_phi * sin_theta, cos_phi)
    }

    #[inline(always)]
    pub fn random_on_hemisphere(norm: Vec3) -> Vec3 {
        let r = Vec3::random_on_sphere();
        if r * norm < 0.0 {
            -r
        } else {
            r
        }
    }

    #[inline(always)]
    /// rho, theta, phi
    pub fn spherical(&self) -> (f64, f64, f64) {
        let rho = self.len();
        (rho, self[0].atan2(self[1]), (self[2] / rho).acos())
    }

    #[inline(always)]
    /// theta, phi
    pub fn norm_spherical(&self) -> (f64, f64) {
        (self[1].atan2(self[0]), self[2].acos())
    }

    #[inline(always)]
    pub fn reflect(&self, norm: Vec3) -> Vec3 {
        norm * (*self * norm) * 2.0 - *self
    }
}

impl Index<usize> for Vec3 {
    type Output = f64;
    fn index(&self, idx: usize) -> &Self::Output {
        &self.0[idx]
    }
}

impl IndexMut<usize> for Vec3 {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        &mut self.0[idx]
    }
}

impl Add<Vec3> for Vec3 {
    type Output = Vec3;
    fn add(self, other: Vec3) -> Vec3 {
        Vec3::new(self[0] + other[0], self[1] + other[1], self[2] + other[2])
    }
}

impl AddAssign<Vec3> for Vec3 {
    fn add_assign(&mut self, other: Vec3) {
        self[0] += other[0];
        self[1] += other[1];
        self[2] += other[2];
    }
}

impl Add<f64> for Vec3 {
    type Output = Vec3;
    fn add(self, offset: f64) -> Vec3 {
        Vec3::new(self[0] + offset, self[1] + offset, self[2] + offset)
    }
}

impl AddAssign<f64> for Vec3 {
    fn add_assign(&mut self, offset: f64) {
        self[0] += offset;
        self[1] += offset;
        self[2] += offset;
    }
}

impl Sub<Vec3> for Vec3 {
    type Output = Vec3;
    fn sub(self, other: Vec3) -> Vec3 {
        Vec3::new(self[0] - other[0], self[1] - other[1], self[2] - other[2])
    }
}

impl Neg for Vec3 {
    type Output = Vec3;
    fn neg(self) -> Vec3 {
        Vec3::new(-self[0], -self[1], -self[2])
    }
}

impl SubAssign<Vec3> for Vec3 {
    fn sub_assign(&mut self, other: Vec3) {
        self[0] -= other[0];
        self[1] -= other[1];
        self[2] -= other[2];
    }
}

impl Sub<f64> for Vec3 {
    type Output = Vec3;
    fn sub(self, offset: f64) -> Vec3 {
        Vec3::new(self[0] - offset, self[1] - offset, self[2] - offset)
    }
}

impl SubAssign<f64> for Vec3 {
    fn sub_assign(&mut self, offset: f64) {
        self[0] -= offset;
        self[1] -= offset;
        self[2] -= offset;
    }
}

impl Mul<Vec3> for Vec3 {
    type Output = f64;
    fn mul(self, other: Vec3) -> f64 {
        self.dot(other)
    }
}

impl Mul<f64> for Vec3 {
    type Output = Vec3;
    fn mul(self, scale: f64) -> Vec3 {
        Vec3::new(self[0] * scale, self[1] * scale, self[2] * scale)
    }
}

impl MulAssign<f64> for Vec3 {
    fn mul_assign(&mut self, scale: f64) {
        self[0] *= scale;
        self[1] *= scale;
        self[2] *= scale;
    }
}

impl Div<f64> for Vec3 {
    type Output = Vec3;
    fn div(self, scale: f64) -> Vec3 {
        let over = 1.0 / scale;
        Vec3::new(self[0] * over, self[1] * over, self[2] * over)
    }
}

impl DivAssign<f64> for Vec3 {
    fn div_assign(&mut self, scale: f64) {
        let over = 1.0 / scale;
        self[0] *= over;
        self[1] *= over;
        self[2] *= over;
    }
}
