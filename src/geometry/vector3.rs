use num;
use num::Num;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::DivAssign;
use std::ops::Index;
use std::ops::IndexMut;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Sub;
use std::ops::SubAssign;

use crate::math;
use crate::types::{Float, Int};

use num::Integer;
use num::NumCast;
use num::One;
use num::Signed;
use num::ToPrimitive;
use num::Zero;
use std::ops::Neg;
use std::ops::Rem;

use crate::geometry::vector::Cross;
use crate::geometry::vector::Distance;
use crate::geometry::vector::Dot;
use crate::geometry::vector::Length;
use crate::types::to_float;
use crate::types::Number;
use std::cmp;
use std::ops::Deref;

pub struct ParseVectorError;

//fn to_float<T: Copy + ToPrimitive>(x: T) -> Float {
//    x.to_f64().unwrap()
//}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Vector3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Vector3<T>
where
    T: Number,
{
    pub fn fill(v: T) -> Self {
        Self { x: v, y: v, z: v }
    }

    pub fn new(x: T, y: T, z: T) -> Self {
        Vector3 { x, y, z }
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

    pub fn max(v1: &Self, v2: &Self) -> Self {
        Self {
            x: math::max(v1.x, v2.x),
            y: math::max(v1.y, v2.y),
            z: math::max(v1.z, v2.z),
        }
    }

    pub fn min(v1: &Self, v2: &Self) -> Self {
        Self {
            x: math::min(v1.x, v2.x),
            y: math::min(v1.y, v2.y),
            z: math::min(v1.z, v2.z),
        }
    }
}

impl<T> Cross for Vector3<T>
where
    T: Number,
{
    fn cross(&self, other: &Self) -> Self {
        Self {
            x: (self.y * other.z) - (self.z * other.y),
            y: (self.z * other.x) - (self.x * other.z),
            z: (self.x * other.y) - (self.y * other.x),
        }
    }
}

impl<T> Length<T> for Vector3<T>
where
    T: Number,
{
    fn length_squared(&self) -> T {
        self.x * self.x + self.y * self.y
    }

    fn length(&self) -> Float {
        to_float(self.length_squared()).sqrt()
    }
}

impl<T> Add for Vector3<T>
where
    T: Number,
{
    type Output = Vector3<T>;
    fn add(self, other: Self) -> Self::Output {
        Self::Output {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<T> Add for &Vector3<T>
where
    T: Number,
{
    type Output = Vector3<T>;
    fn add(self, other: Self) -> Self::Output {
        Self::Output {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<T> Add<T> for Vector3<T>
where
    T: Number,
{
    type Output = Self;

    fn add(self, other: T) -> Self {
        Self {
            x: self.x + other,
            y: self.y + other,
            z: self.z + other,
        }
    }
}

impl<T> AddAssign for Vector3<T>
where
    T: Number,
{
    fn add_assign(&mut self, other: Self) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl<T> AddAssign<T> for Vector3<T>
where
    T: Number,
{
    fn add_assign(&mut self, other: T) {
        self.x += other;
        self.y += other;
        self.z += other;
    }
}

impl<T> Distance for Vector3<T>
where
    T: Number,
{
    fn distance(&self, other: &Self) -> Float {
        (self - other).length()
    }
}

impl<T> Div<T> for Vector3<T>
where
    T: Number + From<Float>,
{
    type Output = Vector3<T>;

    fn div(self, rhs: T) -> Self::Output {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl<T> Div for Vector3<T>
where
    T: Number,
{
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        Self {
            x: self.x / rhs.x,
            y: self.y / rhs.y,
            z: self.z / rhs.z,
        }
    }
}

impl<T> DivAssign<T> for Vector3<T>
where
    T: Number,
{
    fn div_assign(&mut self, other: T) {
        self.x /= other;
        self.y /= other;
        self.z /= other;
    }
}

impl<T> Dot<T> for Vector3<T>
where
    T: Number,
{
    fn dot(&self, other: &Self) -> T {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn abs_dot(&self, other: &Self) -> T {
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

impl<T> Mul for Vector3<T>
where
    T: Number,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
            z: self.z * rhs.z,
        }
    }
}

//impl<T> Mul<T> for Vector3<T>
//where
//    T: Number,
//{
//    type Output = Self;
//
//    fn mul(self, rhs: T) -> Self {
//        Self {
//            x: self.x * rhs,
//            y: self.y * rhs,
//            z: self.z * rhs,
//        }
//    }
//}

impl<T> Mul<T> for Vector3<T>
where
    T: Number,
{
    type Output = Self;

    fn mul(self, rhs: T) -> Self {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl<T> Neg for Vector3<T>
where
    T: Number,
{
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: self.x.neg(),
            y: self.y.neg(),
            z: self.z.neg(),
        }
    }
}

impl<T> Num for Vector3<T>
where
    T: Number,
{
    type FromStrRadixErr = ParseVectorError;

    fn from_str_radix(str: &str, radix: u32) -> Result<Self, ParseVectorError> {
        Err(ParseVectorError)
    }
}

impl<T> One for Vector3<T>
where
    T: Number,
{
    fn one() -> Self {
        Self {
            x: T::one(),
            y: T::one(),
            z: T::one(),
        }
    }
}

impl<T> MulAssign<T> for Vector3<T>
where
    T: Number,
{
    fn mul_assign(&mut self, rhs: T) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl<T> Rem for Vector3<T>
where
    T: Number,
{
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x % rhs.x,
            y: self.y % rhs.y,
            z: self.z % rhs.z,
        }
    }
}

impl<T> Signed for Vector3<T>
where
    T: Number,
{
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

impl<T> Sub for &Vector3<T>
where
    T: Number,
{
    type Output = Vector3<T>;
    fn sub(self, rhs: Self) -> Self::Output {
        Self::Output {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl<T> Sub for Vector3<T>
where
    T: Number,
{
    type Output = Vector3<T>;
    fn sub(self, rhs: Self) -> Self::Output {
        Self::Output {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl<T> Sub<T> for Vector3<T>
where
    T: Number,
{
    type Output = Self;

    fn sub(self, rhs: T) -> Self {
        Self {
            x: self.x - rhs,
            y: self.y - rhs,
            z: self.z - rhs,
        }
    }
}

impl<T> SubAssign<T> for Vector3<T>
where
    T: Number,
{
    fn sub_assign(&mut self, rhs: T) {
        self.x -= rhs;
        self.y -= rhs;
        self.z -= rhs;
    }
}

impl<T> Zero for Vector3<T>
where
    T: Number,
{
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
