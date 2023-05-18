use crate::math::{Vec3, F};
use crate::material::Material;

type Id = usize;

#[derive(Debug)]
pub struct HitInfo {
    pub normal: Vec3,
    pub material: Material,
}

pub trait Hit {
    fn hit(&self, closest: Vec3) -> HitInfo;
}

#[derive(Debug, Clone)]
pub struct Ray {
    pub origin: Vec3,
    pub direction: Vec3,
}

impl Ray {
    fn at(&self, t: F) -> Vec3 {
        return self.origin + self.direction * t;
    }
}

pub struct Local {
    pub normal: Vec3,
}

#[derive(Debug, Clone)]
pub struct Sphere {
    pub id: Id,
    pub center: Vec3,
    pub radius: F,
    pub material: Material,
}

impl Intersect for Sphere {
    #[inline(always)]
    fn intersect(&self, ray: &Ray) -> Option<Intersection> {
        let m = ray.origin - self.center;
        let b = m * ray.direction;
        let c = m * m - self.radius * self.radius;

        if c > 0.0 && b > 0.0 {
            return None;
        }
        let discriminant = b * b - c;

        if discriminant < 0.0 {
            return None;
        }

        return Some(Intersection {
            distance: -b - discriminant.sqrt(),
            id: self.id,
        });
    }
}

impl Hit for Sphere {
    #[inline(always)]
    fn hit(&self, closest: Vec3) -> HitInfo {
        todo!()
    }
}

#[derive(Debug, Clone)]
pub struct Triangle {
    pub id: Id,
    pub verts: [Vec3; 3],
    pub material: Material,
}

impl Intersect for Triangle {
    #[inline(always)]
    fn intersect(&self, ray: &Ray) -> Option<Intersection> {
        let edge1 = self[1] - self[0];
        let edge2 = self[2] - self[0];

        let h = ray.direction.cross(edge2);
        let a = ray.direction.dot(h);
        if a.abs() < 0.0001 {
            return None;
        }

        let f = 1.0 / a;
        let s = ray.origin - self[0];
        let u = f * (s * h);
        if u < 0.0 || u > 1.0 {
            return None;
        }

        let q = s.cross(edge1);
        let v = f * (ray.direction * q);
        if v < 0.0 || u + v > 1.0 {
            return None;
        }

        let t = f * (self[2] * q);
        if t < 0.0001 {
            return None;
        } else {
            return Some(Intersection {
                distance: t,
                id: self.id,
            });
        }
    }
}

impl std::ops::Index<usize> for Triangle {
    type Output = Vec3;
    fn index(&self, idx: usize) -> &Self::Output {
        &self.verts[idx]
    }
}

#[derive(Debug, Clone)]
pub struct Plane {
    pub id: Id,
    pub support: Vec3,
    pub normal: Vec3,
    pub material: Material,
}

impl Intersect for Plane {
    #[inline(always)]
    fn intersect(&self, ray: &Ray) -> Option<Intersection> {
        let denom = ray.direction * self.normal;

        if denom.abs() < 0.0001 {
            return None;
        }
        // let d = -self.support * self.normal;
        let t = (self.support - ray.origin) * self.normal / denom;

        if t < 0.0001 {
            return None;
        } else {
            return Some(Intersection {
                distance: t,
                id: self.id,
            });
        }
    }
}

#[derive(Debug, Clone)]
pub struct Cuboid {}

impl Intersect for Cuboid {
    fn intersect(&self, ray: &Ray) -> Option<Intersection> {
        todo!()
    }
}

#[derive(Debug, Clone)]
pub struct Shapes {
    spheres: Vec<Sphere>,
    planes: Vec<Plane>,
    triangles: Vec<Triangle>,
    cuboids: Vec<Cuboid>,
}

impl Intersect for Shapes {
    fn intersect(&self, ray: &Ray) -> Option<Intersection> {
        let mut closest = Intersection::invalid();

        for sphere in &self.spheres {
            if let Some(inters) = sphere.intersect(ray) {
                if inters < closest {
                    closest = inters;
                }
            }
        }

        for plane in &self.planes {
            if let Some(inters) = plane.intersect(ray) {
                if inters < closest {
                    closest = inters;
                }
            }
        }

        for triangle in &self.triangles {
            if let Some(inters) = triangle.intersect(ray) {
                if inters < closest {
                    closest = inters;
                }
            }
        }

        for cuboid in &self.cuboids {
            if let Some(inters) = cuboid.intersect(ray) {
                if inters < closest {
                    closest = inters;
                }
            }
        }

        closest.wrap()
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Intersection {
    pub distance: F,
    pub id: Id,
}

impl PartialOrd for Intersection {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.distance.total_cmp(&other.distance))
    }
}

impl Intersection {
    pub fn new(distance: F, id: Id) -> Self {
        Intersection { distance, id }
    }

    pub fn invalid() -> Self {
        Intersection {
            distance: std::f32::INFINITY as F,
            id: 0,
        }
    }

    pub fn is_invalid(&self) -> bool {
        self.distance == F::INFINITY
    }

    pub fn wrap(self) -> Option<Intersection> {
        if self.is_invalid() {
            None
        } else {
            Some(self)
        }
    }
}

pub trait Intersect {
    fn intersect(&self, ray: &Ray) -> Option<Intersection>;
}
