use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
};

use num::{Float, Num};

use rand::Rng;

pub trait Constants {
    const E: Self;
    const FRAC_1_PI: Self;
    const FRAC_1_SQRT_2: Self;
    const FRAC_2_PI: Self;
    const FRAC_2_SQRT_PI: Self;
    const FRAC_PI_2: Self;
    const FRAC_PI_3: Self;
    const FRAC_PI_4: Self;
    const FRAC_PI_6: Self;
    const FRAC_PI_8: Self;
    const LN_2: Self;
    const LN_10: Self;
    const LOG2_10: Self;
    const LOG2_E: Self;
    const LOG10_2: Self;
    const LOG10_E: Self;
    const PI: Self;
    const SQRT_2: Self;
    const TAU: Self;
}

impl Constants for f32 {
    const E: Self = std::f32::consts::E;
    const FRAC_1_PI: Self = std::f32::consts::FRAC_1_PI;
    const FRAC_1_SQRT_2: Self = std::f32::consts::FRAC_1_SQRT_2;
    const FRAC_2_PI: Self = std::f32::consts::FRAC_2_PI;
    const FRAC_2_SQRT_PI: Self = std::f32::consts::FRAC_2_SQRT_PI;
    const FRAC_PI_2: Self = std::f32::consts::FRAC_PI_2;
    const FRAC_PI_3: Self = std::f32::consts::FRAC_PI_3;
    const FRAC_PI_4: Self = std::f32::consts::FRAC_PI_4;
    const FRAC_PI_6: Self = std::f32::consts::FRAC_PI_6;
    const FRAC_PI_8: Self = std::f32::consts::FRAC_PI_8;
    const LN_2: Self = std::f32::consts::LN_2;
    const LN_10: Self = std::f32::consts::LN_10;
    const LOG2_10: Self = std::f32::consts::LOG2_10;
    const LOG2_E: Self = std::f32::consts::LOG2_E;
    const LOG10_2: Self = std::f32::consts::LOG10_2;
    const LOG10_E: Self = std::f32::consts::LOG10_E;
    const PI: Self = std::f32::consts::PI;
    const SQRT_2: Self = std::f32::consts::SQRT_2;
    const TAU: Self = std::f32::consts::TAU;
}

impl Constants for f64 {
    const E: Self = std::f64::consts::E;
    const FRAC_1_PI: Self = std::f64::consts::FRAC_1_PI;
    const FRAC_1_SQRT_2: Self = std::f64::consts::FRAC_1_SQRT_2;
    const FRAC_2_PI: Self = std::f64::consts::FRAC_2_PI;
    const FRAC_2_SQRT_PI: Self = std::f64::consts::FRAC_2_SQRT_PI;
    const FRAC_PI_2: Self = std::f64::consts::FRAC_PI_2;
    const FRAC_PI_3: Self = std::f64::consts::FRAC_PI_3;
    const FRAC_PI_4: Self = std::f64::consts::FRAC_PI_4;
    const FRAC_PI_6: Self = std::f64::consts::FRAC_PI_6;
    const FRAC_PI_8: Self = std::f64::consts::FRAC_PI_8;
    const LN_2: Self = std::f64::consts::LN_2;
    const LN_10: Self = std::f64::consts::LN_10;
    const LOG2_10: Self = std::f64::consts::LOG2_10;
    const LOG2_E: Self = std::f64::consts::LOG2_E;
    const LOG10_2: Self = std::f64::consts::LOG10_2;
    const LOG10_E: Self = std::f64::consts::LOG10_E;
    const PI: Self = std::f64::consts::PI;
    const SQRT_2: Self = std::f64::consts::SQRT_2;
    const TAU: Self = std::f64::consts::TAU;
}

#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(transparent)]
pub struct BiVec3<T: Copy>([T; 3]);

impl<T: Copy> BiVec3<T> {
    pub fn new(components: [T; 3]) -> Self {
        BiVec3(components)
    }
}

impl<T: Copy> Into<BiVec3<T>> for Vec<T, 3> {
    fn into(self) -> BiVec3<T> {
        BiVec3(self.0)
    }
}

impl<T: Copy> Into<Vec<T, 3>> for BiVec3<T> {
    fn into(self) -> Vec<T, 3> {
        Vec(self.0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rotor3<T: Copy>(T, BiVec3<T>);

impl<T: Copy + Float> Rotor3<T> {
    fn rot_ab(a: &Vec<T, 3>, b: &Vec<T, 3>) -> Self {
        let s = T::one() + b.dot(*a);
        let b = b.outer(*a);
        let res = Rotor3(s, b);
        res.normalized()
    }

    fn rot_angle_plane(plane: BiVec3<T>, angle: T) -> Self {
        let sin_s = (angle / (T::one() + T::one())).sin();
        let s = (angle / (T::one() + T::one())).cos();
        let b = BiVec3([
            -sin_s * plane.0[0],
            -sin_s * plane.0[1],
            -sin_s * plane.0[2],
        ]);
        Self(s, b)
    }

    fn normalize(&mut self) {}

    fn normalized(&self) -> Self {
        todo!()
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(transparent)]
pub struct Vec<T: Copy, const Dim: usize>([T; Dim]);

pub type F = f32;
pub type Vec3 = Vec<F, 3>;

pub type Vec2f = Vec<f32, 2>;
pub type Vec2d = Vec<f32, 3>;
pub type Vec3f = Vec<f32, 3>;
pub type Vec3d = Vec<f64, 3>;
pub type Vec4f = Vec<f32, 4>;
pub type Vec4d = Vec<f64, 4>;

impl<T: Copy, const Dim: usize> Vec<T, Dim> {
    #[inline(always)]
    pub const fn new(values: [T; Dim]) -> Self {
        Vec(values)
    }
}

impl<T: Copy + Num, const Dim: usize> Vec<T, Dim> {
    #[inline(always)]
    pub fn dot(&self, other: Vec<T, Dim>) -> T {
        let mut res = T::zero();
        for i in 0..Dim {
            res = res + self[i] * other[i];
        }
        res
    }
}

impl<T: Copy + Num, const Dim: usize> Vec<T, Dim> {
    #[inline(always)]
    pub fn zero() -> Self {
        Vec([T::zero(); Dim])
    }

    pub fn one() -> Self {
        Vec([T::one(); Dim])
    }

    #[inline(always)]
    pub fn len_sq(&self) -> T {
        self.dot(*self)
    }
}

impl<T: Float> Vec<T, 3> {
    pub fn cross(&self, other: Vec<T, 3>) -> Vec<T, 3> {
        Vec::new([
            self[1] * other[2] - self[2] * other[1],
            self[2] * other[0] - self[0] * other[2],
            self[0] * other[1] - self[1] * other[0],
        ])
    }

    pub fn outer(&self, other: Vec<T, 3>) -> BiVec3<T> {
        BiVec3::new([
            self[0] * other[1] - other[0] * self[1],
            self[2] * other[0] - other[2] * self[0],
            self[1] * other[2] - other[1] * self[2],
        ])
    }

    pub fn prod(&self, other: Vec<T, 3>) -> Rotor3<T> {
        Rotor3(self.dot(other), self.outer(other))
    }

    pub fn rot_ab(&self, a: Vec<T, 3>, b: Vec<T, 3>) -> Vec<T, 3> {
        let ab = a.prod(b);
        let ba = b.prod(a);
        todo!()
    }
}

impl<T: Float + Constants + From<f32>> Vec<T, 3> {
    /// Uniform direction on the unit sphere via Archimedes' cylinder
    /// projection: uniform `(x, y)` in the unit disc map to a unit vector with
    /// `z` distributed exactly as the sphere's polar density demands.
    ///
    /// The previous formulation sampled an angle pair and spent five
    /// transcendental calls per direction (`acos` plus two `sin`/`cos` pairs),
    /// which perf attributed to ~30% of total render time. This needs one
    /// `sqrt` and a single `gen` per component. Acceptance is pi/4, so the
    /// loop averages about 1.27 iterations.
    ///
    /// The generator is passed in rather than fetched per call. `rand`'s
    /// `ThreadRng` checks a reseeding counter on every draw, and that atomic
    /// load was the hottest single instruction in the whole binary under
    /// perf. One generator per worker pays it once rather than per sample.
    pub fn random_on_sphere<R: rand::Rng + ?Sized>(rng: &mut R) -> Vec<T, 3> {
        loop {
            let x: f32 = rng.gen::<f32>() * 2.0 - 1.0;
            let y: f32 = rng.gen::<f32>() * 2.0 - 1.0;
            let s = x * x + y * y;
            if s < 1.0 && s > 0.0 {
                // (x*f)^2 + (y*f)^2 + (1 - 2s)^2 == 4s(1-s) + 1 - 4s + 4s^2 == 1.
                let f = 2.0 * (1.0 - s).sqrt();
                return Vec::new([
                    <f32 as Into<T>>::into(x * f),
                    <f32 as Into<T>>::into(y * f),
                    <f32 as Into<T>>::into(1.0 - 2.0 * s),
                ]);
            }
        }
    }

    #[inline(always)]
    pub fn from_spherical(r: T, theta: T, phi: T) -> Vec<T, 3> {
        unimplemented!()
    }

    #[inline(always)]
    pub fn from_spherical_unit(theta: T, phi: T) -> Vec<T, 3> {
        let sin_theta = theta.sin();
        let cos_theta = theta.cos();
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();

        Vec::new([sin_phi * cos_theta, sin_phi * sin_theta, cos_phi])
    }

    #[inline(always)]
    pub fn random_on_hemisphere<R: rand::Rng + ?Sized>(norm: Vec<T, 3>, rng: &mut R) -> Vec<T, 3> {
        let r = Vec::random_on_sphere(rng);
        if r * norm < T::zero() {
            -r
        } else {
            r
        }
    }

    #[inline(always)]
    /// rho, theta, phi
    pub fn spherical(&self) -> (T, T, T) {
        let rho = self.len();
        (rho, self[0].atan2(self[1]), (self[2] / rho).acos())
    }

    #[inline(always)]
    /// theta, phi
    pub fn norm_spherical(&self) -> (T, T) {
        (self[1].atan2(self[0]), self[2].acos())
    }
}

impl<T: Float + From<f32>, const Dim: usize> Vec<T, Dim> {
    #[inline(always)]
    pub fn new_normalized(values: [T; Dim]) -> Self {
        Self::new(values).normalized()
    }

    #[inline(always)]
    pub fn len(&self) -> T {
        self.dot(*self).sqrt()
    }

    #[inline(always)]
    pub fn normalized(&self) -> Self {
        let len = self.len();
        // A zero or non-finite length reciprocates to inf, so the product would be NaN.
        if len.is_finite() && len > T::zero() {
            *self * len.recip()
        } else {
            Self::zero()
        }
    }

    #[inline(always)]
    pub fn reflect(&self, norm: Vec<T, Dim>) -> Vec<T, Dim> {
        norm * (*self * norm) * <f32 as Into<T>>::into(2.0) - *self
    }
}

impl<T: Copy, const Dim: usize> Index<usize> for Vec<T, Dim> {
    type Output = T;
    fn index(&self, idx: usize) -> &Self::Output {
        &self.0[idx]
    }
}

impl<T: Copy, const Dim: usize> IndexMut<usize> for Vec<T, Dim> {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        &mut self.0[idx]
    }
}

impl<T: Copy + Num, const Dim: usize> Add<Vec<T, Dim>> for Vec<T, Dim> {
    type Output = Vec<T, Dim>;
    fn add(self, other: Vec<T, Dim>) -> Vec<T, Dim> {
        let mut res = self.clone();
        for i in 0..res.0.len() {
            res[i] = res[i] + other[i];
        }
        res
    }
}

impl<T: Copy + Num, const Dim: usize> AddAssign<Vec<T, Dim>> for Vec<T, Dim> {
    fn add_assign(&mut self, other: Vec<T, Dim>) {
        for i in 0..self.0.len() {
            self[i] = self[i] + other[i];
        }
    }
}

impl<T: Copy + Num, const Dim: usize> Add<T> for Vec<T, Dim> {
    type Output = Vec<T, Dim>;
    fn add(self, offset: T) -> Vec<T, Dim> {
        let mut res = self.clone();
        for i in 0..self.0.len() {
            res[i] = res[i] + offset
        }
        res
    }
}

impl<T: Copy + Num, const Dim: usize> AddAssign<T> for Vec<T, Dim> {
    fn add_assign(&mut self, offset: T) {
        for i in 0..self.0.len() {
            self[i] = self[i] + offset;
        }
    }
}

impl<T: Copy + Num, const Dim: usize> Sub<Vec<T, Dim>> for Vec<T, Dim> {
    type Output = Vec<T, Dim>;
    fn sub(self, other: Vec<T, Dim>) -> Vec<T, Dim> {
        let mut res = self.clone();
        for i in 0..self.0.len() {
            res[i] = res[i] - other[i];
        }
        res
    }
}

impl<T: Copy + Num, const Dim: usize> Neg for Vec<T, Dim> {
    type Output = Vec<T, Dim>;
    fn neg(self) -> Vec<T, Dim> {
        let mut res = self.clone();
        for i in 0..self.0.len() {
            res[i] = T::zero() - res[i];
        }
        res
    }
}

impl<T: Copy + Num, const Dim: usize> SubAssign<Vec<T, Dim>> for Vec<T, Dim> {
    fn sub_assign(&mut self, other: Vec<T, Dim>) {
        for i in 0..self.0.len() {
            self[i] = self[i] - other[i];
        }
    }
}

impl<T: Copy + Num, const Dim: usize> Sub<T> for Vec<T, Dim> {
    type Output = Vec<T, Dim>;
    fn sub(self, offset: T) -> Vec<T, Dim> {
        let mut res = self.clone();
        for i in 0..self.0.len() {
            res[i] = res[i] - offset;
        }
        res
    }
}

impl<T: Copy + Num, const Dim: usize> SubAssign<T> for Vec<T, Dim> {
    fn sub_assign(&mut self, offset: T) {
        for i in 0..self.0.len() {
            self[i] = self[i] - offset;
        }
    }
}

impl<T: Copy + Num, const Dim: usize> Mul<Vec<T, Dim>> for Vec<T, Dim> {
    type Output = T;
    fn mul(self, other: Vec<T, Dim>) -> T {
        self.dot(other)
    }
}

impl<T: Copy + Num, const Dim: usize> Mul<T> for Vec<T, Dim> {
    type Output = Vec<T, Dim>;
    fn mul(self, scale: T) -> Vec<T, Dim> {
        let mut res = self.clone();
        for i in 0..self.0.len() {
            res[i] = res[i] * scale;
        }
        res
    }
}

impl<T: Copy + Num, const Dim: usize> MulAssign<T> for Vec<T, Dim> {
    fn mul_assign(&mut self, scale: T) {
        for i in 0..self.0.len() {
            self[i] = self[i] * scale;
        }
    }
}

impl<T: Copy + Float, const Dim: usize> Div<T> for Vec<T, Dim> {
    type Output = Vec<T, Dim>;
    fn div(self, scale: T) -> Vec<T, Dim> {
        let over = scale.recip();
        let mut ret = self.clone();
        for i in 0..self.0.len() {
            ret[i] = ret[i] * over;
        }
        ret
    }
}

impl<T: Copy + Float, const Dim: usize> DivAssign<T> for Vec<T, Dim> {
    fn div_assign(&mut self, scale: T) {
        let over = scale.recip();
        for i in 0..self.0.len() {
            self[i] = self[i] * over;
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    /// The old angle-pair sampler was replaced because it burned five
    /// transcendental calls per direction. Correctness of the replacement
    /// matters more than the speed: these pin down unit length and the
    /// spherical distribution, which a trig-free method could easily get wrong.
    #[test]
    fn random_on_sphere_is_unit_length() {
        let mut rng = rand::thread_rng();
        for _ in 0..20_000 {
            let v = Vec3::random_on_sphere(&mut rng);
            let len = v.len();
            assert!(
                (len - 1.0).abs() < 1e-4,
                "direction not unit length: {len} for {v:?}"
            );
        }
    }

    /// Uniform on the sphere means z is uniform on [-1, 1], and each of the six
    /// signed axis directions is hit one eighth of the time.
    #[test]
    fn random_on_sphere_is_uniform() {
        const N: usize = 120_000;
        let mut z_sum = 0.0f32;
        let mut octants = [0u32; 8];
        let mut rng = rand::thread_rng();
        for _ in 0..N {
            let v = Vec3::random_on_sphere(&mut rng);
            z_sum += v[2];
            let oct =
                (v[0] > 0.0) as usize | ((v[1] > 0.0) as usize) << 1 | ((v[2] > 0.0) as usize) << 2;
            octants[oct] += 1;
        }

        // Mean of z must be ~0. The standard error of the mean here is
        // about 1/sqrt(3N), so 0.01 is a loose but non-vacuous bound.
        let mean_z = z_sum / N as f32;
        assert!(mean_z.abs() < 0.01, "z mean drifted: {mean_z}");

        // Every octant within 2% of the expected eighth is far outside the
        // sampling noise at this N (expected 15000, sd ~37).
        for (i, count) in octants.iter().enumerate() {
            let expected = N as f32 / 8.0;
            assert!(
                (*count as f32 - expected).abs() < expected * 0.02,
                "octant {i} skewed: {count} vs {expected}"
            );
        }
    }

    /// The hemisphere helper must never return a direction behind the plane.
    #[test]
    fn random_on_hemisphere_stays_on_the_near_side() {
        let normals = [
            Vec3::new([0.0, 0.0, 1.0]),
            Vec3::new([1.0, 0.0, 0.0]),
            Vec3::new([0.0, -1.0, 0.0]),
            Vec3::new([0.577, 0.577, 0.577]),
        ];
        let mut rng = rand::thread_rng();
        for n in normals {
            let n = n.normalized();
            for _ in 0..10_000 {
                let d = Vec3::random_on_hemisphere(n, &mut rng);
                assert!(d * n >= 0.0, "direction {d:?} behind normal {n:?}");
            }
        }
    }
}
