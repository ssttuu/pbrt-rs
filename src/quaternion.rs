use crate::geometry::vector::Dot;
use crate::math;
use crate::transform::Matrix4x4;
use crate::transform::Transform;
use crate::types::Float;
use crate::Vector3f;
use std::ops;

#[derive(Clone, Copy, Default)]
pub struct Quaternion {
    pub v: Vector3f,
    pub w: Float,
}

impl Quaternion {
    pub fn normalized(self) -> Quaternion {
        &self / self.dot(&self).sqrt()
    }

    pub fn to_transform(&self) -> Transform {
        let (xx, yy, zz) = (
            self.v.x * self.v.x,
            self.v.y * self.v.y,
            self.v.z * self.v.z,
        );
        let (xy, xz, yz) = (
            self.v.x * self.v.y,
            self.v.x * self.v.z,
            self.v.y * self.v.z,
        );
        let (wx, wy, wz) = (self.v.x * self.w, self.v.y * self.w, self.v.z * self.w);

        let mut m = Matrix4x4::new();
        m[0][0] = 1. - 2. * (yy + zz);
        m[0][1] = 2. * (xy + wz);
        m[0][2] = 2. * (xz - wy);
        m[1][0] = 2. * (xy - wz);
        m[1][1] = 1. - 2. * (xx + zz);
        m[1][2] = 2. * (yz + wx);
        m[2][0] = 2. * (xz + wy);
        m[2][1] = 2. * (yz - wx);
        m[2][2] = 1. - 2. * (xx + yy);

        // Transpose since we are left-handed
        return Transform {
            matrix: m.transpose(),
            inverse: m,
        };
    }

    pub fn slerp(t: Float, q1: &Quaternion, q2: &Quaternion) -> Quaternion {
        let cos_theta = q1.dot(q2);
        if cos_theta > 0.9995 {
            return (q1 * (1.0 - t) + q2 * t).normalized();
        }

        let theta = math::acos(math::clamp(cos_theta, -1., 1.));
        let thetap = theta * t;
        let qperp = q2 - &(q1 * cos_theta);

        q1 * math::cos(thetap) + &qperp * math::sin(thetap)
    }
}

impl Dot<Float> for Quaternion {
    fn dot(&self, other: &Self) -> Float {
        self.v.dot(&other.v) + self.w * other.w
    }
    fn abs_dot(&self, other: &Self) -> Float {
        self.dot(other).abs()
    }
}

impl ops::Add for &Quaternion {
    type Output = Quaternion;
    fn add(self, rhs: Self) -> Self::Output {
        Self::Output {
            v: self.v + rhs.v,
            w: self.w + rhs.w,
        }
    }
}

impl ops::Add for Quaternion {
    type Output = Quaternion;
    fn add(self, rhs: Self) -> Self::Output {
        Self::Output {
            v: self.v + rhs.v,
            w: self.w + rhs.w,
        }
    }
}

impl ops::Div<Float> for &Quaternion {
    type Output = Quaternion;
    fn div(self, rhs: Float) -> Self::Output {
        Self::Output {
            v: self.v / rhs,
            w: self.w / rhs,
        }
    }
}

impl ops::Mul<Float> for &Quaternion {
    type Output = Quaternion;
    fn mul(self, rhs: Float) -> Self::Output {
        Self::Output {
            v: self.v * rhs,
            w: self.w * rhs,
        }
    }
}

impl ops::MulAssign<Float> for &mut Quaternion {
    fn mul_assign(&mut self, rhs: Float) {
        self.v *= rhs;
        self.w *= rhs;
    }
}

impl ops::Neg for Quaternion {
    type Output = Quaternion;
    fn neg(self) -> Self::Output {
        Self::Output {
            v: -self.v,
            w: -self.w,
        }
    }
}

impl ops::Sub for &Quaternion {
    type Output = Quaternion;
    fn sub(self, rhs: Self) -> Self::Output {
        Self::Output {
            v: self.v - rhs.v,
            w: self.w - rhs.w,
        }
    }
}

impl ops::Sub for Quaternion {
    type Output = Quaternion;
    fn sub(self, rhs: Self) -> Self::Output {
        Self::Output {
            v: self.v - rhs.v,
            w: self.w - rhs.w,
        }
    }
}
