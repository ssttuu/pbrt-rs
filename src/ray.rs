use crate::geometry::vector::dot;
use crate::math::next_float_down;
use crate::math::next_float_up;
use crate::medium::Medium;
use crate::types::Float;
use crate::Normal3f;
use crate::Point3f;
use crate::Vector3f;
use num::abs;
use num::Signed;

#[derive(Clone)]
pub struct Ray {
    pub origin: Point3f,
    pub direction: Vector3f,
    pub time: Float,
    pub time_max: Float,

    pub origin_error: Option<Vector3f>,
    pub direction_error: Option<Vector3f>,

    //    pub medium: Option<&'a Medium>,
    pub has_differentials: bool,
    pub rx_origin: Option<Point3f>,
    pub ry_origin: Option<Point3f>,
    pub rx_direction: Option<Vector3f>,
    pub ry_direction: Option<Vector3f>,
}

impl Ray {
    pub fn new(
        origin: Point3f,
        direction: Vector3f,
        time: Float,
        time_max: Float,
        //        medium: Option<&'a Medium>,
    ) -> Ray {
        Ray {
            origin,
            direction,
            time,
            time_max,

            origin_error: None,
            direction_error: None,

            //            medium,
            has_differentials: false,
            rx_origin: None,
            ry_origin: None,
            rx_direction: None,
            ry_direction: None,
        }
    }
}

pub fn offset_ray_origin(p: &Point3f, p_error: &Vector3f, n: &Normal3f, w: &Vector3f) -> Point3f {
    let d = dot(&n.abs(), &p_error);
    let mut offset = *n * d;
    if dot(w, n) < 0.0 {
        offset = -offset;
    }

    let mut po = *p + offset;
    for i in 0..3 {
        if offset[i] > 0.0 {
            po[i] = next_float_up(&po[i]);
        } else if offset[i] < 0.0 {
            po[i] = next_float_down(&po[i]);
        }
    }

    po
}
