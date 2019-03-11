use crate::geometry::distance;
use crate::geometry::Vector3;
use crate::math;
use crate::ray::Ray;
use crate::types::Float;
use crate::types::Int;
use crate::types::Number;
use crate::Point2;
use crate::Point2f;
use crate::Point3;
use crate::Point3f;
use crate::Vector3f;
use std::ops;

#[derive(Clone, Default)]
pub struct BoundingSphere {
    pub center: Point3f,
    pub radius: Float,
}

#[derive(Clone, Default)]
pub struct Bounds2 {
    pub min: Point2f,
    pub max: Point2f,
}

#[derive(Clone, Default)]
pub struct Bounds3 {
    pub min: Point3f,
    pub max: Point3f,
}

impl Bounds3 {
    pub fn bounding_sphere(&self) -> BoundingSphere {
        let center = self.min + self.max / 2.;
        BoundingSphere {
            center,
            radius: if self.inside(&center) {
                distance(center, self.max)
            } else {
                0.
            },
        }
    }

    pub fn corner(&self, corner: i32) -> Point3f {
        Point3 {
            x: self[corner & 1].x,
            y: self[(corner & 2) / 2].y,
            z: self[(corner & 4) / 4].z,
        }
    }

    pub fn diagonal(&self) -> Vector3f {
        &self.max - &self.min
    }

    pub fn surface_area(&self) -> Float {
        let d = self.diagonal();
        (d.x * d.y + d.x * d.z + d.y * d.z) * 2.
    }

    pub fn maximum_extent(&self) -> Int {
        let d = self.diagonal();
        if d.x > d.y && d.x > d.z {
            0
        } else if d.y > d.z {
            1
        } else {
            2
        }
    }

    pub fn intersect_p(
        &self,
        r: &Ray,
        inv_dir: &Vector3f,
        direction_is_negative: [i32; 3],
    ) -> bool {
        // check for ray intersection against x and y slabs
        let mut t_min = (self[direction_is_negative[0]].x - r.origin.x) * inv_dir.x;
        let mut t_max = (self[1 - direction_is_negative[0]].x - r.origin.x) * inv_dir.x;
        let ty_min = (self[direction_is_negative[1]].y - r.origin.y) * inv_dir.y;
        let mut ty_max = (self[1 - direction_is_negative[1]].y - r.origin.y) * inv_dir.y;

        // update tMax and tyMax to ensure robust bounds intersection
        t_max *= 1. + 2. * math::gamma(3.);
        ty_max *= 1. + 2. * math::gamma(3.);
        if t_min > ty_max || ty_min > t_max {
            return false;
        }
        if ty_min > t_min {
            t_min = ty_min;
        }
        if ty_max < t_max {
            t_max = ty_max;
        }

        // check for ray intersection against z slab
        let tz_min = (self[direction_is_negative[2]].z - r.origin.z) * inv_dir.z;
        let mut tz_max = (self[1 - direction_is_negative[2]].z - r.origin.z) * inv_dir.z;

        // update tzMax to ensure robust bounds intersection
        tz_max *= 1. + 2. * math::gamma(3.);
        if t_min > tz_max || tz_min > t_max {
            return false;
        }
        if tz_min > t_min {
            t_min = tz_min;
        }
        if tz_max < t_max {
            t_max = tz_max;
        }
        return t_min < r.time_max && t_max > 0.;
    }

    pub fn lerp(&self, other: &Point3f) -> Point3f {
        Point3f {
            x: math::lerp(other.x, self.min.x, self.max.x),
            y: math::lerp(other.y, self.min.y, self.max.y),
            z: math::lerp(other.z, self.min.z, self.max.z),
        }
    }

    pub fn offset(&self, other: &Point3f) -> Vector3f {
        let mut o = other - &self.min;
        if self.max.x > self.min.x {
            o.x /= self.max.x - self.min.x;
        }
        if self.max.y > self.min.y {
            o.y /= self.max.y - self.min.y;
        }
        if self.max.z > self.min.z {
            o.z /= self.max.z - self.min.z;
        }
        return o;
    }

    pub fn union_point(&mut self, p: &Point3f) {
        self.min = Point3::min(&self.min, p);
        self.max = Point3::max(&self.max, p);
    }

    pub fn union(&mut self, b: &Bounds3) {
        self.min = Point3::min(&self.min, &b.min);
        self.max = Point3::max(&self.max, &b.max);
    }

    pub fn inside(&self, p: &Point3f) -> bool {
        p.x >= self.min.x
            && p.x <= self.max.x
            && p.y >= self.min.y
            && p.y <= self.max.y
            && p.z >= self.min.z
            && p.z <= self.max.z
    }
}

impl ops::Index<i32> for Bounds3 {
    type Output = Point3f;
    fn index(&self, index: i32) -> &Self::Output {
        if index == 0 {
            &self.min
        } else {
            &self.max
        }
    }
}
