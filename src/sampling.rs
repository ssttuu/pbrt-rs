use crate::types::Float;
use crate::math::consts;
use crate::Point2f;
use crate::Point3f;
use crate::math;

pub fn uniform_cone_pdf(cos_theta_max: Float) -> Float {
    1. / (2. * consts::PI * (1. - cos_theta_max))
}

pub fn uniform_sample_sphere(u: &Point2f) -> Point3f {
    let z = 1. - 2. * u.x;
    let r = math::sqrt(math::max(0., 1. - z * z));
    let phi = 2. * consts::PI * u.y;
    Point3f {
        x: r * math::cos(phi),
        y: r * math::sin(phi),
        z,
    }
}