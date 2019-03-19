pub mod vector;
pub mod vector2;
pub mod vector3;

pub use self::vector::*;
pub use self::vector2::*;
pub use self::vector3::*;
use crate::math;
use crate::types::Float;
use crate::Vector3f;

pub fn coordinate_system(v1: Vector3f) -> (Vector3f, Vector3f) {
    let v2 = if v1.x.abs() > v1.y.abs() {
        let v = v1.x * v1.x + v1.z * v1.z;
        Vector3f::new(-v1.z, 0., v1.x) / v
    } else {
        let v = v1.y * v1.y + v1.z * v1.z;
        Vector3f::new(0., v1.z, -v1.y) / v
    };

    (v2, cross(&v1, &v2))
}

pub fn spherical_direction_xyz(
    sin_theta: Float,
    cos_theta: Float,
    phi: Float,
    x: &Vector3f,
    y: &Vector3f,
    z: &Vector3f,
) -> Vector3f {
    x * sin_theta * math::cos(phi) + y * sin_theta * math::sin(phi) + z * cos_theta
}
