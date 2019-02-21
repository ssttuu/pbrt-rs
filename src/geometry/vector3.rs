use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::AddAssign;
use std::ops::Sub;
use std::ops::SubAssign;
use std::ops::Div;
use std::ops::DivAssign;
use std::ops::Add;
use std::ops::Index;
use std::ops::IndexMut;
use num::Num;
use num;

use crate::types::{Int, Float};
use crate::math;

use num::Signed;
use num::ToPrimitive;
use std::ops::Neg;
use num::NumCast;
use num::Integer;
use num::Zero;
use std::ops::Rem;
use num::One;

use crate::geometry::vector::Length;
use crate::geometry::vector::Cross;
use crate::geometry::vector::Dot;
use crate::types::Number;
use crate::types::to_float;

pub struct ParseVectorError;

//fn to_float<T: Copy + ToPrimitive>(x: T) -> Float {
//    x.to_f64().unwrap()
//}

#[derive(Debug, Eq, PartialEq, Clone, Copy)]
pub struct Vector3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}


impl<T> Vector3<T> where T: Number {
    pub fn new(x: T, y: T, z: T) -> Self {
        Vector3 {
            x,
            y,
            z,
        }
    }

    pub fn is_negative(&self) -> Vector3<bool> {
        Vector3 {
            x: self.x < T::zero(),
            y: self.y < T::zero(),
            z: self.z < T::zero(),
        }
    }

    pub fn normalize(v: Self) -> Vector3<Float> {
        v.normalized()
    }

    pub fn normalized(&self) -> Vector3<Float> {
        let mul = to_float(1.0) / self.length();

        Vector3 {
            x: self.x.to_float() * mul,
            y: self.y.to_float() * mul,
            z: self.z.to_float() * mul,
        }
    }
}

impl<T> Length<T> for Vector3<T> where T: Number {
    fn length_squared(&self) -> T {
        self.x * self.x + self.y * self.y
    }

    fn length(&self) -> Float {
        to_float(self.length_squared()).sqrt()
    }
}

impl<T: Add<Output=T>> Add for Vector3<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<T: Add<Output=T> + Copy> Add<T> for Vector3<T> {
    type Output = Self;

    fn add(self, other: T) -> Self {
        Self {
            x: self.x + other,
            y: self.y + other,
            z: self.z + other,
        }
    }
}

impl<T: AddAssign + Copy> AddAssign<Self> for Vector3<T> {
    fn add_assign(&mut self, other: Self) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl<T: AddAssign + Copy> AddAssign<T> for Vector3<T> {
    fn add_assign(&mut self, other: T) {
        self.x += other;
        self.y += other;
        self.z += other;
    }
}

impl<T> Div<T> for Vector3<T> where T: Num + ToPrimitive + Copy + Div<Output=T> + Mul<Output=T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self {
        let inv = (T::one() / rhs) as T;
        Self {
            x: self.x * inv,
            y: self.y * inv,
            z: self.z * inv,
        }
    }
}

impl<T> Div for Vector3<T> where T: Div<Output=T> + Copy {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        Self {
            x: self.x / rhs.x,
            y: self.y / rhs.y,
            z: self.z / rhs.z,
        }
    }
}

impl<T: DivAssign + Copy> DivAssign<T> for Vector3<T> {
    fn div_assign(&mut self, other: T) {
        self.x /= other;
        self.y /= other;
        self.z /= other;
    }
}

impl<T> Dot<T> for Vector3<T> where T: Num + Copy + Signed {
    fn dot(&self, other: Self) -> T {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn abs_dot(&self, other: Self) -> T {
        self.dot(other).abs()
    }
}

impl<T> Index<i32> for Vector3<T> {
    type Output = T;

    fn index(&self, index: i32) -> &T {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("index out of range"),
        }
    }
}

impl<T> IndexMut<i32> for Vector3<T> {
    fn index_mut(&mut self, index: i32) -> &mut T {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("index out of range"),
        }
    }
}

impl<T: Mul<Output=T> + Copy> Mul for Vector3<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
            z: self.z * rhs.z,
        }
    }
}

impl<T: Mul<Output=T> + Copy> Mul<T> for Vector3<T> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl<T: Neg<Output=T> + Copy> Neg for Vector3<T> {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: self.x.neg(),
            y: self.y.neg(),
            z: self.z.neg(),
        }
    }
}

impl<T: Num + Copy> Num for Vector3<T> {
    type FromStrRadixErr = ParseVectorError;

    fn from_str_radix(str: &str, radix: u32) -> Result<Self, ParseVectorError> {
        Err(ParseVectorError)
    }
}

impl<T: One + Copy> One for Vector3<T> {
    fn one() -> Self {
        Self {
            x: T::one(),
            y: T::one(),
            z: T::one(),
        }
    }
}

impl<T: MulAssign + Copy> MulAssign<T> for Vector3<T> {
    fn mul_assign(&mut self, rhs: T) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl<T: Rem<Output=T> + Copy> Rem for Vector3<T> {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x % rhs.x,
            y: self.y % rhs.y,
            z: self.z % rhs.z,
        }
    }
}

impl<T: Signed + Copy> Signed for Vector3<T> {
    fn abs(&self) -> Self {
        Self {
            x: self.x.abs(),
            y: self.y.abs(),
            z: self.z.abs(),
        }
    }

    fn abs_sub(&self, other: &Self) -> Self {
        Self {
            x: self.x.abs_sub(&other.x),
            y: self.y.abs_sub(&other.y),
            z: self.z.abs_sub(&other.z),
        }
    }

    fn signum(&self) -> Self {
        Self {
            x: self.x.signum(),
            y: self.y.signum(),
            z: self.z.signum(),
        }
    }

    fn is_positive(&self) -> bool {
        self.x.is_positive() && self.y.is_positive() && self.z.is_positive()
    }

    fn is_negative(&self) -> bool {
        self.x.is_negative() && self.y.is_negative() && self.z.is_negative()
    }
}

impl<T: Sub<Output=T> + Copy> Sub<Self> for Vector3<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl<T: Sub<Output=T> + Copy> Sub<T> for Vector3<T> {
    type Output = Self;

    fn sub(self, rhs: T) -> Self {
        Self {
            x: self.x - rhs,
            y: self.y - rhs,
            z: self.z - rhs,
        }
    }
}

impl<T: SubAssign + Copy> SubAssign<T> for Vector3<T> {
    fn sub_assign(&mut self, rhs: T) {
        self.x -= rhs;
        self.y -= rhs;
        self.z -= rhs;
    }
}

impl<T: Zero + Copy> Zero for Vector3<T> {
    fn zero() -> Self {
        Self {
            x: T::zero(),
            y: T::zero(),
            z: T::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y.is_zero() && self.z.is_zero()
    }
}

