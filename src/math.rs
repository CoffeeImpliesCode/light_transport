// WIP linear-algebra module. The rotor, the `Vec2/3/4` aliases, the spherical
// conversions and the `Constants` impls are written ahead of their call sites.
// Scoped here rather than at the crate root so that dead code in the
// application modules still warns.
#![allow(dead_code)]

use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
};

use num::{Float, Num};

use rand::RngExt;

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
pub struct Vec<T: Copy, const DIM: usize>([T; DIM]);

pub type F = f32;
pub type Vec3 = Vec<F, 3>;

pub type Vec2f = Vec<f32, 2>;
pub type Vec2d = Vec<f32, 3>;
pub type Vec3f = Vec<f32, 3>;
pub type Vec3d = Vec<f64, 3>;
pub type Vec4f = Vec<f32, 4>;
pub type Vec4d = Vec<f64, 4>;

impl<T: Copy, const DIM: usize> Vec<T, DIM> {
    #[inline(always)]
    pub const fn new(values: [T; DIM]) -> Self {
        Vec(values)
    }
}

impl<T: Copy + Num, const DIM: usize> Vec<T, DIM> {
    /// Sum of the componentwise products, accumulated left to right.
    ///
    /// Zipping the two component arrays keeps the index inside the slice
    /// iterators, so the bounds check is gone by construction instead of by
    /// relying on the optimiser to fold it out of a `DIM`-trip loop. Trip
    /// count and summation order are unchanged, so the result is bit for bit
    /// what the index loop produced. Called from every primitive intersection.
    #[inline(always)]
    pub fn dot(&self, other: Vec<T, DIM>) -> T {
        let mut res = T::zero();
        for (s, o) in self.0.iter().zip(other.0.iter()) {
            res = res + *s * *o;
        }
        res
    }
}

impl<T: Copy + Num, const DIM: usize> Vec<T, DIM> {
    #[inline(always)]
    pub fn zero() -> Self {
        Vec([T::zero(); DIM])
    }

    pub fn one() -> Self {
        Vec([T::one(); DIM])
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
        let _ab = a.prod(b);
        let _ba = b.prod(a);
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
    /// `sqrt` and a single `random` per component. Acceptance is pi/4, so the
    /// loop averages about 1.27 iterations.
    ///
    /// The generator is passed in rather than fetched per call. `rand`'s
    /// `ThreadRng` checks a reseeding counter on every draw, and that atomic
    /// load was the hottest single instruction in the whole binary under
    /// perf. One generator per worker pays it once rather than per sample.
    ///
    /// Force-inline: `random_on_hemisphere` is already `#[inline(always)]` and
    /// calls this on every scatter, so leaving the sampler out of line plants a
    /// call and a register spill in the middle of the bounce loop. Pure
    /// codegen change; the draw sequence is untouched.
    #[inline(always)]
    pub fn random_on_sphere<R: rand::Rng + ?Sized>(rng: &mut R) -> Vec<T, 3> {
        loop {
            let x: f32 = rng.random::<f32>() * 2.0 - 1.0;
            let y: f32 = rng.random::<f32>() * 2.0 - 1.0;
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
    pub fn from_spherical(_r: T, _theta: T, _phi: T) -> Vec<T, 3> {
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

impl<T: Float + From<f32>, const DIM: usize> Vec<T, DIM> {
    #[inline(always)]
    pub fn new_normalized(values: [T; DIM]) -> Self {
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
        //
        // Positivity is tested first on purpose. Both operands are pure
        // predicates, so `&&` short-circuits them in either order and the
        // result is identical, but `len > 0` rejects zero, every negative and
        // NaN in one ordinary compare. `is_finite` is an abs-and-compare
        // against the bit pattern and only then has to run, for the `+inf` a
        // huge vector can still produce. A direction reaching here is unit
        // length to within rounding, so `len` is near 1 and the common path
        // pays one compare instead of two.
        if len > T::zero() && len.is_finite() {
            *self * len.recip()
        } else {
            Self::zero()
        }
    }

    #[inline(always)]
    pub fn reflect(&self, norm: Vec<T, DIM>) -> Vec<T, DIM> {
        norm * (*self * norm) * <f32 as Into<T>>::into(2.0) - *self
    }
}

impl<T: Copy, const DIM: usize> Index<usize> for Vec<T, DIM> {
    type Output = T;
    fn index(&self, idx: usize) -> &Self::Output {
        &self.0[idx]
    }
}

impl<T: Copy, const DIM: usize> IndexMut<usize> for Vec<T, DIM> {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        &mut self.0[idx]
    }
}

// Every operator below walks the component arrays with `zip` or with a plain
// `iter_mut` rather than `for i in 0..DIM { self[i] = ... }`. The index form
// makes each read and write go through `Index`, which is bounds-checked and
// which this module keeps that way for its callers; here the index lives
// inside a slice iterator, so the check has nowhere to be generated. The
// operation, its operand order and its evaluation order per component are
// exactly what they were, so results are bit for bit unchanged.
impl<T: Copy + Num, const DIM: usize> Add<Vec<T, DIM>> for Vec<T, DIM> {
    type Output = Vec<T, DIM>;
    #[inline(always)]
    fn add(self, other: Vec<T, DIM>) -> Vec<T, DIM> {
        let mut res = self.0;
        for (r, o) in res.iter_mut().zip(other.0.iter()) {
            *r = *r + *o;
        }
        Vec(res)
    }
}

impl<T: Copy + Num, const DIM: usize> AddAssign<Vec<T, DIM>> for Vec<T, DIM> {
    #[inline(always)]
    fn add_assign(&mut self, other: Vec<T, DIM>) {
        for (r, o) in self.0.iter_mut().zip(other.0.iter()) {
            *r = *r + *o;
        }
    }
}

impl<T: Copy + Num, const DIM: usize> Add<T> for Vec<T, DIM> {
    type Output = Vec<T, DIM>;
    #[inline(always)]
    fn add(self, offset: T) -> Vec<T, DIM> {
        let mut res = self.0;
        for r in res.iter_mut() {
            *r = *r + offset;
        }
        Vec(res)
    }
}

impl<T: Copy + Num, const DIM: usize> AddAssign<T> for Vec<T, DIM> {
    #[inline(always)]
    fn add_assign(&mut self, offset: T) {
        for r in self.0.iter_mut() {
            *r = *r + offset;
        }
    }
}

impl<T: Copy + Num, const DIM: usize> Sub<Vec<T, DIM>> for Vec<T, DIM> {
    type Output = Vec<T, DIM>;
    #[inline(always)]
    fn sub(self, other: Vec<T, DIM>) -> Vec<T, DIM> {
        let mut res = self.0;
        for (r, o) in res.iter_mut().zip(other.0.iter()) {
            *r = *r - *o;
        }
        Vec(res)
    }
}

impl<T: Copy + Num, const DIM: usize> Neg for Vec<T, DIM> {
    type Output = Vec<T, DIM>;
    #[inline(always)]
    fn neg(self) -> Vec<T, DIM> {
        // `T::zero() - x` and `-x` are not the same function: on floats they
        // disagree in the sign of zero (`0 - 0` is `+0`, `-0` is `-0`) and a
        // negated zero can reach the renderer through
        // `random_on_hemisphere`. The subtract is the instruction that was
        // already being emitted, so it costs nothing to keep and it keeps the
        // operator bit-exact.
        let zero = T::zero();
        let mut res = self.0;
        for r in res.iter_mut() {
            *r = zero - *r;
        }
        Vec(res)
    }
}

impl<T: Copy + Num, const DIM: usize> SubAssign<Vec<T, DIM>> for Vec<T, DIM> {
    #[inline(always)]
    fn sub_assign(&mut self, other: Vec<T, DIM>) {
        for (r, o) in self.0.iter_mut().zip(other.0.iter()) {
            *r = *r - *o;
        }
    }
}

impl<T: Copy + Num, const DIM: usize> Sub<T> for Vec<T, DIM> {
    type Output = Vec<T, DIM>;
    #[inline(always)]
    fn sub(self, offset: T) -> Vec<T, DIM> {
        let mut res = self.0;
        for r in res.iter_mut() {
            *r = *r - offset;
        }
        Vec(res)
    }
}

impl<T: Copy + Num, const DIM: usize> SubAssign<T> for Vec<T, DIM> {
    #[inline(always)]
    fn sub_assign(&mut self, offset: T) {
        for r in self.0.iter_mut() {
            *r = *r - offset;
        }
    }
}

impl<T: Copy + Num, const DIM: usize> Mul<Vec<T, DIM>> for Vec<T, DIM> {
    type Output = T;
    #[inline(always)]
    fn mul(self, other: Vec<T, DIM>) -> T {
        self.dot(other)
    }
}

impl<T: Copy + Num, const DIM: usize> Mul<T> for Vec<T, DIM> {
    type Output = Vec<T, DIM>;
    #[inline(always)]
    fn mul(self, scale: T) -> Vec<T, DIM> {
        let mut res = self.0;
        for r in res.iter_mut() {
            *r = *r * scale;
        }
        Vec(res)
    }
}

impl<T: Copy + Num, const DIM: usize> MulAssign<T> for Vec<T, DIM> {
    #[inline(always)]
    fn mul_assign(&mut self, scale: T) {
        for r in self.0.iter_mut() {
            *r = *r * scale;
        }
    }
}

impl<T: Copy + Float, const DIM: usize> Div<T> for Vec<T, DIM> {
    type Output = Vec<T, DIM>;
    #[inline(always)]
    fn div(self, scale: T) -> Vec<T, DIM> {
        // Reciprocal then multiply, kept on purpose. This is what the operator
        // has always computed and the render is pinned bit-for-bit against it;
        // a real `r / scale` rounds differently in the last bit and that
        // difference is visible in the golden image. Not a wart to fix here.
        let over = scale.recip();
        let mut ret = self.0;
        for r in ret.iter_mut() {
            *r = *r * over;
        }
        Vec(ret)
    }
}

impl<T: Copy + Float, const DIM: usize> DivAssign<T> for Vec<T, DIM> {
    #[inline(always)]
    fn div_assign(&mut self, scale: T) {
        let over = scale.recip();
        for r in self.0.iter_mut() {
            *r = *r * over;
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    /// A fixed-seed generator, so a statistical assertion fails the same way
    /// on every run instead of once in twenty. See `random_on_sphere_is_uniform`.
    fn seeded_rng(seed: u64) -> rand::rngs::StdRng {
        use rand::SeedableRng;
        rand::rngs::StdRng::seed_from_u64(seed)
    }

    /// The old angle-pair sampler was replaced because it burned five
    /// transcendental calls per direction. Correctness of the replacement
    /// matters more than the speed: these pin down unit length and the
    /// spherical distribution, which a trig-free method could easily get wrong.
    #[test]
    fn random_on_sphere_is_unit_length() {
        let mut rng = seeded_rng(0x5EED_0000_0000_0001);
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
        let mut rng = seeded_rng(0x5EED_0000_0000_0002);
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

        // Each octant count is Binomial(N, 1/8), so its standard deviation is
        // sqrt(N * p * (1 - p)) — about 115 counts at this N, not the 37 an
        // earlier comment claimed. A 2% band was 2.6 sigma, which failed once
        // in twenty runs on a correct sampler. Bound at 5 sigma instead: the
        // chance of a false failure across all eight octants drops below 1e-5.
        let expected = N as f32 / 8.0;
        let sigma = (N as f32 * 0.125 * 0.875).sqrt();
        for (i, count) in octants.iter().enumerate() {
            assert!(
                (*count as f32 - expected).abs() < 5.0 * sigma,
                "octant {i} skewed: {count} vs {expected} (5 sigma = {})",
                5.0 * sigma
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
        let mut rng = seeded_rng(0x5EED_0000_0000_0003);
        for n in normals {
            let n = n.normalized();
            for _ in 0..10_000 {
                let d = Vec3::random_on_hemisphere(n, &mut rng);
                assert!(d * n >= 0.0, "direction {d:?} behind normal {n:?}");
            }
        }
    }
}
