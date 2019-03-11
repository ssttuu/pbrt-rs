use crate::bounds::Bounds3;
use crate::geometry::vector::face_forward;
use crate::geometry::vector::Cross;
use crate::geometry::vector::Dot;
use crate::geometry::vector::Length;
use crate::geometry::vector3::Vector3;
use crate::interaction::Shading;
use crate::interaction::SurfaceInteraction;
use crate::math;
use crate::quaternion::Quaternion;
use crate::ray::Ray;
use crate::types::Float;
use crate::types::Int;
use crate::types::Number;
use crate::Point3f;
use crate::Vector3f;
use num::abs;
use std::fmt;
use std::mem;
use std::num::ParseIntError;
use std::ops;
use std::ptr;
use std::rc::Rc;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Matrix4x4 {
    data: [[Float; 4]; 4],
}

#[derive(Clone, Debug)]
pub enum MatrixError {
    SingularMatrix,
    SameDirection,
}

impl Matrix4x4 {
    pub fn new() -> Matrix4x4 {
        Matrix4x4 {
            data: [
                [1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.],
            ],
        }
    }

    pub fn inverse(&self) -> Result<Matrix4x4, MatrixError> {
        //        var indxc, indxr, ipiv [4]int
        let mut indxc = [0, 0, 0, 0];
        let mut indxr = [0, 0, 0, 0];
        let mut ipiv = [0, 0, 0, 0];

        let mut minv = self.clone();

        for i in 0..4 {
            let (mut irow, mut icol) = (0, 0);
            let mut big = 0.;

            // Choose pivot
            for j in 0..4 {
                if ipiv[j] != 1 {
                    for k in 0..4 {
                        if ipiv[k] == 0 {
                            if abs(minv[j][k]) >= big {
                                big = abs(minv[j][k]);
                                irow = j;
                                icol = k;
                            }
                        } else if ipiv[k] > 1 {
                            return Err(MatrixError::SingularMatrix);
                        }
                    }
                }
            }

            ipiv[icol] += 1;

            // Swap rows _irow_ and _icol_ for pivot
            if irow != icol {
                for k in 0..4 {
                    minv.swap((irow, k), (icol, k));
                    //                    mem::swap(&mut minv[irow][k], &mut minv[icol][k])
                }
            }

            indxr[i] = irow;
            indxc[i] = icol;
            if minv[icol][icol] == 0. {
                return Err(MatrixError::SingularMatrix);
            }

            // Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
            let pivinv = 1. / minv[icol][icol];
            minv[icol][icol] = 1.;
            for j in 0..4 {
                minv[icol][j] *= pivinv;
            }

            // Subtract this row from others to zero out their columns
            for j in 0..4 {
                if j != icol {
                    let save = minv[j][icol];
                    minv[j][icol] = 0.;
                    for k in 0..4 {
                        minv[j][k] -= minv[icol][k] * save;
                    }
                }
            }
        }

        // Swap columns to reflect permutation
        for j in (0..4).rev() {
            if indxr[j] != indxc[j] {
                for k in 0..4 {
                    minv.swap((k, indxr[j]), (k, indxc[j]));
                    //                    mem::swap(&mut minv[k][indxr[j]], &mut minv[k][indxc[j]]);
                }
            }
        }

        return Ok(minv);
    }

    // swap (row a, column a) with (row b, column b)
    fn swap(&mut self, a: (usize, usize), b: (usize, usize)) {
        unsafe { ptr::swap(&mut self[a.0][a.1], &mut self[b.0][b.1]) }
    }

    pub fn transpose(&self) -> Self {
        Self {
            data: [
                [self[0][0], self[1][0], self[2][0], self[3][0]],
                [self[0][1], self[1][1], self[2][1], self[3][1]],
                [self[0][2], self[1][2], self[2][2], self[3][2]],
                [self[0][3], self[1][3], self[2][3], self[3][3]],
            ],
        }
    }
}

impl ops::Index<usize> for Matrix4x4 {
    type Output = [Float; 4];

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl ops::IndexMut<usize> for Matrix4x4 {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

impl ops::Mul for Matrix4x4 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        let mut out = Self::new();
        for i in 0..4 {
            for j in 0..4 {
                out[i][j] = self[i][0] * rhs[0][j]
                    + self[i][1] * rhs[1][j]
                    + self[i][2] * rhs[2][j]
                    + self[i][3] * rhs[3][j]
            }
        }
        out
    }
}

#[derive(Clone, Default, PartialEq)]
pub struct Transform {
    pub matrix: Matrix4x4,
    pub inverse: Matrix4x4,
}

impl Transform {
    pub fn new(m: Matrix4x4) -> Transform {
        Transform {
            inverse: m.inverse().unwrap(),
            matrix: m,
        }
    }

    pub fn is_identity(&self) -> bool {
        self.matrix[0][0] == 1.
            && self.matrix[0][1] == 0.
            && self.matrix[0][2] == 0.
            && self.matrix[0][3] == 0.
            && self.matrix[1][0] == 0.
            && self.matrix[1][1] == 1.
            && self.matrix[1][2] == 0.
            && self.matrix[1][3] == 0.
            && self.matrix[2][0] == 0.
            && self.matrix[2][1] == 0.
            && self.matrix[2][2] == 1.
            && self.matrix[2][3] == 0.
            && self.matrix[3][0] == 0.
            && self.matrix[3][1] == 0.
            && self.matrix[3][2] == 0.
            && self.matrix[3][3] == 1.
    }

    pub fn invert(&self) -> Self {
        Self {
            matrix: self.inverse.clone(),
            inverse: self.matrix.clone(),
        }
    }

    pub fn transform_point(&self, p: &Point3f, p_error: &Vector3f) -> (Point3f, Vector3f) {
        // compute transformed points from p
        let xp = self.matrix[0][0] * p.x
            + self.matrix[0][1] * p.y
            + self.matrix[0][2] * p.z
            + self.matrix[0][3];
        let yp = self.matrix[1][0] * p.x
            + self.matrix[1][1] * p.y
            + self.matrix[1][2] * p.z
            + self.matrix[1][3];
        let zp = self.matrix[2][0] * p.x
            + self.matrix[2][1] * p.y
            + self.matrix[2][2] * p.z
            + self.matrix[2][3];
        let wp = self.matrix[3][0] * p.x
            + self.matrix[3][1] * p.y
            + self.matrix[3][2] * p.z
            + self.matrix[3][3];

        let abs_error = Vector3f {
            x: (math::gamma(3.0) + 1.)
                * (abs(self.matrix[0][0]) * p_error.x
                    + abs(self.matrix[0][1]) * p_error.y
                    + abs(self.matrix[0][2]) * p_error.z)
                + (math::gamma(3.0)
                    * (abs(self.matrix[0][0] * p.x)
                        + abs(self.matrix[0][1]) * p.y
                        + abs(self.matrix[0][2] * p.z + abs(self.matrix[0][3])))),
            y: (math::gamma(3.0) + 1.)
                * (abs(self.matrix[1][0]) * p_error.x
                    + abs(self.matrix[1][1]) * p_error.y
                    + abs(self.matrix[1][2]) * p_error.z)
                + (math::gamma(3.0)
                    * (abs(self.matrix[1][0] * p.x)
                        + abs(self.matrix[1][1]) * p.y
                        + abs(self.matrix[1][2] * p.z + abs(self.matrix[1][3])))),
            z: (math::gamma(3.0) + 1.)
                * (abs(self.matrix[2][0]) * p_error.x
                    + abs(self.matrix[2][1]) * p_error.y
                    + abs(self.matrix[2][2]) * p_error.z)
                + (math::gamma(3.0)
                    * (abs(self.matrix[2][0] * p.x)
                        + abs(self.matrix[2][1]) * p.y
                        + abs(self.matrix[2][2] * p.z + abs(self.matrix[2][3])))),
        };

        let new_point = Point3f::new(xp, yp, zp);

        if wp == 1. {
            return (new_point, abs_error);
        }

        (new_point / wp, abs_error)
    }

    pub fn transform_vector(&self, v: &Vector3f) -> Vector3f {
        Vector3f {
            x: self.matrix[0][0] * v.x + self.matrix[0][1] * v.y + self.matrix[0][2] * v.z,
            y: self.matrix[1][0] * v.x + self.matrix[1][1] * v.y + self.matrix[1][2] * v.z,
            z: self.matrix[2][0] * v.x + self.matrix[2][1] * v.y + self.matrix[2][2] * v.z,
        }
    }

    pub fn transform_vector_with_abs_error(&self, v: &Vector3f) -> (Vector3f, Vector3f) {
        let new_vector = Vector3f {
            x: self.matrix[0][0] * v.x + self.matrix[0][1] * v.y + self.matrix[0][2] * v.z,
            y: self.matrix[1][0] * v.x + self.matrix[1][1] * v.y + self.matrix[1][2] * v.z,
            z: self.matrix[2][0] * v.x + self.matrix[2][1] * v.y + self.matrix[2][2] * v.z,
        };
        let abs_error = Vector3f {
            x: math::gamma(3.0)
                * (abs(self.matrix[0][0] * v.x)
                    + abs(self.matrix[0][1] * v.y)
                    + abs(self.matrix[0][2] * v.z)),
            y: math::gamma(3.0)
                * (abs(self.matrix[1][0] * v.x)
                    + abs(self.matrix[1][1] * v.y)
                    + abs(self.matrix[1][2] * v.z)),
            z: math::gamma(3.0)
                * (abs(self.matrix[2][0] * v.x)
                    + abs(self.matrix[2][1] * v.y)
                    + abs(self.matrix[2][2] * v.z)),
        };
        (new_vector, abs_error)
    }

    pub fn transform_normal(&self, n: &Vector3f) -> Vector3f {
        Vector3f {
            x: self.inverse[0][0] * n.x + self.inverse[1][0] * n.y + self.inverse[2][0] * n.z,
            y: self.inverse[0][1] * n.x + self.inverse[1][1] * n.y + self.inverse[2][1] * n.z,
            z: self.inverse[0][2] * n.x + self.inverse[1][2] * n.y + self.inverse[2][2] * n.z,
        }
    }

    pub fn transform_ray(&self, ray: &Ray) -> Ray {
        let (mut origin, origin_error) = self.transform_point(&ray.origin, &Vector3f::default());
        let (direction, direction_error) = self.transform_vector_with_abs_error(&ray.direction);
        let length_squared = direction.length_squared();
        if length_squared > 0. {
            let dt = abs(direction).dot(&origin_error) / length_squared;
            origin += direction * dt;
        }

        Ray {
            origin,
            direction,
            origin_error: Some(origin_error),
            direction_error: Some(direction_error),
            ..ray.clone()
        }
    }

    pub fn transform_surface_interaction<'prim>(
        &self,
        si: &SurfaceInteraction<'prim>,
    ) -> SurfaceInteraction<'prim> {
        let (point, point_error) = self.transform_point(&si.point, &si.point_error);

        SurfaceInteraction {
            point,
            point_error,
            normal: self.transform_normal(&si.normal),
            wo: self.transform_vector(&si.wo).normalized(),
            dpdu: self.transform_vector(&si.dpdu),
            dpdv: self.transform_vector(&si.dpdv),
            dndu: self.transform_normal(&si.dndu),
            dndv: self.transform_normal(&si.dndv),
            shading: Shading {
                normal: face_forward(si.shading.normal, si.normal),
                ..si.shading.clone()
            },
            ..si.clone()
        }
    }

    pub fn transform_bounds(&self, b: &Bounds3) -> Bounds3 {
        let (corner, _) = self.transform_point(&b.min, &Vector3f::default());
        let mut transformed_bounds = Bounds3 {
            min: corner,
            max: corner,
        };
        for i in 0..8 {
            let (corner, _) = self.transform_point(&b.corner(i), &Vector3f::default());
            transformed_bounds.union_point(&corner);
        }
        transformed_bounds
    }
}

impl ops::Mul for Transform {
    type Output = Transform;
    fn mul(self, rhs: Transform) -> Self {
        Self {
            matrix: self.matrix * rhs.matrix,
            inverse: self.inverse * rhs.inverse,
        }
    }
}

pub fn translate(delta: &Vector3<Float>) -> Transform {
    Transform {
        matrix: Matrix4x4 {
            data: [
                [1., 0., 0., delta.x],
                [0., 1., 0., delta.y],
                [0., 0., 1., delta.z],
                [0., 0., 0., 1.],
            ],
        },
        inverse: Matrix4x4 {
            data: [
                [1., 0., 0., -delta.x],
                [0., 1., 0., -delta.y],
                [0., 0., 1., -delta.z],
                [0., 0., 0., 1.],
            ],
        },
    }
}

pub fn scale(x: Float, y: Float, z: Float) -> Transform {
    Transform {
        matrix: Matrix4x4 {
            data: [
                [x, 0., 0., 0.],
                [0., y, 0., 0.],
                [0., 0., z, 0.],
                [0., 0., 0., 1.],
            ],
        },
        inverse: Matrix4x4 {
            data: [
                [1. / x, 0., 0., 0.],
                [0., 1. / y, 0., 0.],
                [0., 0., 1. / z, 0.],
                [0., 0., 0., 1.],
            ],
        },
    }
}

pub fn rotate_x(degrees: Float) -> Transform {
    let sin_theta = math::sin(math::radians(degrees));
    let cos_theta = math::cos(math::radians(degrees));
    let m = Matrix4x4 {
        data: [
            [1., 0., 0., 0.],
            [0., cos_theta, -sin_theta, 0.],
            [0., sin_theta, cos_theta, 0.],
            [0., 0., 0., 1.],
        ],
    };
    return Transform {
        inverse: m.transpose(),
        matrix: m,
    };
}

pub fn rotate_y(degrees: Float) -> Transform {
    let sin_theta = math::sin(math::radians(degrees));
    let cos_theta = math::cos(math::radians(degrees));
    let m = Matrix4x4 {
        data: [
            [cos_theta, 0., sin_theta, 0.],
            [0., 1., 0., 0.],
            [-sin_theta, 0., cos_theta, 0.],
            [0., 0., 0., 1.],
        ],
    };
    return Transform {
        inverse: m.transpose(),
        matrix: m,
    };
}

pub fn rotate_z(degrees: Float) -> Transform {
    let sin_theta = math::sin(math::radians(degrees));
    let cos_theta = math::cos(math::radians(degrees));
    let m = Matrix4x4 {
        data: [
            [cos_theta, -sin_theta, 0., 0.],
            [sin_theta, cos_theta, 0., 0.],
            [0., 0., 1., 0.],
            [0., 0., 0., 1.],
        ],
    };
    return Transform {
        inverse: m.transpose(),
        matrix: m,
    };
}

pub fn rotate(theta: Float, axis: &Vector3f) -> Transform {
    let a = axis.normalized();
    let sin_theta = math::sin(math::radians(theta));
    let cos_theta = math::cos(math::radians(theta));

    let mut m = Matrix4x4::new();

    // Compute rotation of first basis vector
    m[0][0] = a.x * a.x + (1. - a.x * a.x) * cos_theta;
    m[0][1] = a.x * a.y * (1. - cos_theta) - a.z * sin_theta;
    m[0][2] = a.x * a.z * (1. - cos_theta) + a.y * sin_theta;
    m[0][3] = 0.;

    // Compute rotations of second and third basis vectors
    m[1][0] = a.x * a.y * (1. - cos_theta) + a.z * sin_theta;
    m[1][1] = a.y * a.y + (1. - a.y * a.y) * cos_theta;
    m[1][2] = a.y * a.z * (1. - cos_theta) - a.x * sin_theta;
    m[1][3] = 0.;

    m[2][0] = a.x * a.z * (1. - cos_theta) - a.y * sin_theta;
    m[2][1] = a.y * a.z * (1. - cos_theta) + a.x * sin_theta;
    m[2][2] = a.z * a.z + (1. - a.z * a.z) * cos_theta;
    m[2][3] = 0.;

    return Transform {
        inverse: m.transpose(),
        matrix: m,
    };
}

pub fn look_at(pos: &Point3f, look: &Point3f, up: &Vector3f) -> Result<Transform, String> {
    let mut camera_to_world = Matrix4x4::new();
    // Initialize fourth column of viewing Matrix
    camera_to_world[0][3] = pos.x;
    camera_to_world[1][3] = pos.y;
    camera_to_world[2][3] = pos.z;
    camera_to_world[3][3] = 1.;

    // Initialize first three columns of viewing Matrix
    let dir = (look - pos).normalized();
    if up.normalized().cross(&dir).length() == 0. {
        return Err(String::from("same direction"));
    }
    let right = up.normalized().cross(&dir).normalized();
    let new_up = dir.cross(&right);
    camera_to_world[0][0] = right.x;
    camera_to_world[1][0] = right.y;
    camera_to_world[2][0] = right.z;
    camera_to_world[3][0] = 0.;
    camera_to_world[0][1] = new_up.x;
    camera_to_world[1][1] = new_up.y;
    camera_to_world[2][1] = new_up.z;
    camera_to_world[3][1] = 0.;
    camera_to_world[0][2] = dir.x;
    camera_to_world[1][2] = dir.y;
    camera_to_world[2][2] = dir.z;
    camera_to_world[3][2] = 0.;

    match camera_to_world.inverse() {
        Ok(inv) => Ok(Transform {
            matrix: camera_to_world,
            inverse: inv,
        }),
        Err(e) => Err(String::from("inverting camera")),
    }
}

pub fn orthographic(z_near: Float, z_far: Float) -> Transform {
    scale(1., 1., 1. / (z_far - z_near)) * translate(&Vector3f::new(0., 0., -z_near))
}

pub fn perspective(fov: Float, n: Float, f: Float) -> Transform {
    let persp = Matrix4x4 {
        data: [
            [1., 0., 0., 0.],
            [0., 1., 0., 0.],
            [0., 0., 1., 0.],
            [0., 0., 0., 1.],
        ],
    };
    let inv_tan_ang = 1. / math::tan(math::radians(fov) / 2.0);

    scale(inv_tan_ang, inv_tan_ang, 1.) * Transform::new(persp)
}

#[derive(Clone, Copy, Default)]
pub struct DerivativeTerm {
    pub kc: Float,
    pub kx: Float,
    pub ky: Float,
    pub kz: Float,
}

impl DerivativeTerm {
    pub fn eval(&self, p: &Point3f) -> Float {
        self.kc + self.kx * p.x + self.ky * p.y + self.kz * p.z
    }
}

#[derive(Clone, Default)]
pub struct AnimatedTransform {
    pub start: Transform,
    pub end: Transform,
    pub start_time: Float,
    pub end_time: Float,
    pub actually_animated: bool,
    pub start_t: Vector3f,
    pub end_t: Vector3f,
    pub start_r: Quaternion,
    pub end_r: Quaternion,
    pub start_s: Matrix4x4,
    pub end_s: Matrix4x4,
    pub has_rotation: bool,
    pub c1: [DerivativeTerm; 3],
    pub c2: [DerivativeTerm; 3],
    pub c3: [DerivativeTerm; 3],
    pub c4: [DerivativeTerm; 3],
    pub c5: [DerivativeTerm; 3],
}

impl AnimatedTransform {
    pub fn new(start: Transform, end: Transform, start_time: Float, end_time: Float) -> Self {
        let actually_animated = &start == &end;
        let mut at = Self {
            start,
            end,
            start_time,
            end_time,
            actually_animated,
            ..Self::default()
        };

        if !actually_animated {
            return at;
        }

        // TODO
        //at.startT, at.startR, at.startS = Decompose(start.Matrix)
        //at.endT, at.endR, at.endS = Decompose(end.Matrix)

        // Flip r if needed to select shortest path
        if at.start_r.dot(&at.end_r) < 0. {
            at.end_r = &at.end_r * -1.
        }

        at.has_rotation = at.start_r.dot(&at.end_r) < 0.9995;

        // compute terms of motion derivative function
        if at.has_rotation {
            // TODO
        }

        return at;
    }

    pub fn interpolate(&self, time: Float) -> Transform {
        if !self.actually_animated || time <= self.start_time {
            return self.start.clone();
        }
        if time >= self.end_time {
            return self.end.clone();
        }

        let dt = (time - self.start_time) / (self.end_time - self.start_time);
        let t = self.start_t * (1. - dt) + self.end_t * dt;
        let r = Quaternion::slerp(dt, &self.start_r, &self.end_r);

        let mut s = Matrix4x4::new();
        for i in 0..3 {
            for j in 0..3 {
                s[i][j] = math::lerp(dt, self.start_s[i][j], self.end_s[i][j]);
            }
        }

        translate(&t) * r.to_transform() * Transform::new(s)
    }

    pub fn motion_bounds(&self, b: &Bounds3) -> Bounds3 {
        if !self.actually_animated {
            return self.start.transform_bounds(b);
        }

        return Bounds3::default();
    }

    pub fn transform_ray(&self, ray: &Ray) -> Ray {
        if !self.actually_animated || ray.time <= self.start_time {
            return self.start.transform_ray(ray);
        } else if ray.time >= self.end_time {
            return self.end.transform_ray(ray);
        }
        self.interpolate(ray.time).transform_ray(ray)
    }

    pub fn transform_point_at_time(&self, p: &Point3f, time: Float) -> Point3f {
        let mut point_error = Vector3f::default();
        if !self.actually_animated || time <= self.start_time {
            let (point, _) = self.start.transform_point(p, &point_error);
            return point;
        } else if time >= self.end_time {
            let (point, _) = self.end.transform_point(p, &point_error);
        }

        let (point, _) = self.interpolate(time).transform_point(p, &point_error);
        return point;
    }

    pub fn transform_vector_at_time(&self, v: &Vector3f, time: Float) -> Vector3f {
        let mut point_error = Vector3f::default();
        if !self.actually_animated || time <= self.start_time {
            let (point, _) = self.start.transform_point(v, &point_error);
            return point;
        } else if time >= self.end_time {
            let (point, _) = self.end.transform_point(v, &point_error);
        }

        let (point, _) = self.interpolate(time).transform_point(v, &point_error);
        return point;
    }
}
