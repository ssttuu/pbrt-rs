use crate::geometry::vector3::Vector3;
use crate::types::Float;
use crate::types::Int;
use crate::types::Number;
use crate::Vector3f;
use num::abs;
use std::fmt;
use std::mem;
use std::num::ParseIntError;
use std::ops;
use std::ptr;

#[derive(Clone, Debug)]
pub struct Matrix4x4 {
    data: [[Float; 4]; 4],
}

#[derive(Clone, Debug)]
pub enum MatrixError {
    SingularMatrix,
}

impl Matrix4x4 {
    pub fn new() -> Matrix4x4 {
        Matrix4x4 {
            data: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
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
            let mut big = 0.0;

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
            if minv[icol][icol] == 0.0 {
                return Err(MatrixError::SingularMatrix);
            }

            // Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
            let pivinv = 1.0 / minv[icol][icol];
            minv[icol][icol] = 1.0;
            for j in 0..4 {
                minv[icol][j] *= pivinv;
            }

            // Subtract this row from others to zero out their columns
            for j in 0..4 {
                if j != icol {
                    let save = minv[j][icol];
                    minv[j][icol] = 0.0;
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
        unsafe {
            let pa: *mut Float = &mut self[a.0][a.1];
            let pb: *mut Float = &mut self[b.0][b.1];
            ptr::swap(pa, pb)
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

pub struct Transform {
    matrix: Matrix4x4,
    inverse: Matrix4x4,
}

impl Transform {
    pub fn new(m: Matrix4x4) -> Transform {
        Transform {
            inverse: m.inverse().unwrap(),
            matrix: m,
        }
    }
}

pub fn translate<T: Number>(delta: Vector3<T>) -> Transform {
    Transform {
        matrix: Matrix4x4 {
            data: [
                [1.0, 0.0, 0.0, delta.x.to_float()],
                [0.0, 1.0, 0.0, delta.y.to_float()],
                [0.0, 0.0, 1.0, delta.z.to_float()],
                [0.0, 0.0, 0.0, 1.0],
            ],
        },
        inverse: Matrix4x4 {
            data: [
                [1.0, 0.0, 0.0, -delta.x.to_float()],
                [0.0, 1.0, 0.0, -delta.y.to_float()],
                [0.0, 0.0, 1.0, -delta.z.to_float()],
                [0.0, 0.0, 0.0, 1.0],
            ],
        },
    }
}
