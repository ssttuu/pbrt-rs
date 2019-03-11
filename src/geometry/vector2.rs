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
use crate::geometry::vector::Dot;
use crate::geometry::vector::Length;
use crate::types::to_float;
use crate::types::Number;

pub struct ParseVectorError;

//fn to_float<T: Copy + ToPrimitive>(x: T) -> Float {
//    x.to_f64().unwrap()
//}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Vector2<T> {
    pub x: T,
    pub y: T,
}

impl<T> Vector2<T>
where
    T: Number,
{
    pub fn new(x: T, y: T) -> Self {
        Vector2 { x, y }
    }

    pub fn is_negative(&self) -> Vector2<bool> {
        Vector2 {
            x: self.x < T::zero(),
            y: self.y < T::zero(),
        }
    }

    pub fn normalize(v: Self) -> Vector2<Float> {
        v.normalized()
    }

    pub fn normalized(&self) -> Vector2<Float> {
        let mul = to_float(1.0) / self.length();

        Vector2 {
            x: self.x.to_float() * mul,
            y: self.y.to_float() * mul,
        }
    }
}

impl<T> Length<T> for Vector2<T>
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

impl<T: Add<Output = T>> Add for Vector2<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl<T: Add<Output = T> + Copy> Add<T> for Vector2<T> {
    type Output = Self;

    fn add(self, other: T) -> Self {
        Self {
            x: self.x + other,
            y: self.y + other,
        }
    }
}

impl<T: AddAssign + Copy> AddAssign<Self> for Vector2<T> {
    fn add_assign(&mut self, other: Self) {
        self.x += other.x;
        self.y += other.y;
    }
}

impl<T: AddAssign + Copy> AddAssign<T> for Vector2<T> {
    fn add_assign(&mut self, other: T) {
        self.x += other;
        self.y += other;
    }
}

//impl<T: Div<Output=T> + Copy> Div<T> for Vector2<T> {
//    type Output = Self;
//
//    fn div(self, other: T) -> Self {
//        Self {
//            x: self.x / other,
//            y: self.y / other,
//        }
//    }
//}

impl<T> Div<T> for Vector2<T>
where
    T: Num + ToPrimitive + Copy + Div<Output = T> + Mul<Output = T>,
{
    type Output = Self;

    fn div(self, rhs: T) -> Self {
        let inv = (T::one() / rhs) as T;
        Self {
            x: self.x * inv,
            y: self.y * inv,
        }
    }
}

impl<T> Div for Vector2<T>
where
    T: Div<Output = T> + Copy,
{
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        Self {
            x: self.x / rhs.x,
            y: self.y / rhs.y,
        }
    }
}

impl<T: DivAssign + Copy> DivAssign<T> for Vector2<T> {
    fn div_assign(&mut self, other: T) {
        self.x /= other;
        self.y /= other;
    }
}

impl<T> Dot<T> for Vector2<T>
where
    T: Number,
{
    fn dot(&self, other: &Self) -> T {
        self.x * other.x + self.y * other.y
    }

    fn abs_dot(&self, other: &Self) -> T {
        self.dot(other).abs()
    }
}

//impl<I: Integer + Copy, F: Float + Copy> From<Vector2<I>> for Vector2<F> {
//    fn from(v: Vector2<I>) -> Self {
//        Self {
//            x: v.x as F,
//            y: v.y as F,
//        }
//    }
//}

impl<T> Index<i32> for Vector2<T> {
    type Output = T;

    fn index(&self, index: i32) -> &T {
        match index {
            0 => &self.x,
            1 => &self.y,
            _ => panic!("index out of range"),
        }
    }
}

impl<T> IndexMut<i32> for Vector2<T> {
    fn index_mut(&mut self, index: i32) -> &mut T {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            _ => panic!("index out of range"),
        }
    }
}

impl<T: Mul<Output = T> + Copy> Mul for Vector2<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
        }
    }
}

impl<T: Mul<Output = T> + Copy> Mul<T> for Vector2<T> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl<T: Neg<Output = T> + Copy> Neg for Vector2<T> {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: self.x.neg(),
            y: self.y.neg(),
        }
    }
}

impl<T: Num + Copy> Num for Vector2<T> {
    type FromStrRadixErr = ParseVectorError;

    fn from_str_radix(str: &str, radix: u32) -> Result<Self, ParseVectorError> {
        Err(ParseVectorError)
    }
}

impl<T: One + Copy> One for Vector2<T> {
    fn one() -> Self {
        Self {
            x: T::one(),
            y: T::one(),
        }
    }
}

impl<T: MulAssign + Copy> MulAssign<T> for Vector2<T> {
    fn mul_assign(&mut self, rhs: T) {
        self.x *= rhs;
        self.y *= rhs;
    }
}

impl<T: Rem<Output = T> + Copy> Rem for Vector2<T> {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x % rhs.x,
            y: self.y % rhs.y,
        }
    }
}

impl<T: Signed + Copy> Signed for Vector2<T> {
    fn abs(&self) -> Self {
        Self {
            x: self.x.abs(),
            y: self.y.abs(),
        }
    }

    fn abs_sub(&self, other: &Self) -> Self {
        Self {
            x: self.x.abs_sub(&other.x),
            y: self.y.abs_sub(&other.y),
        }
    }

    fn signum(&self) -> Self {
        Self {
            x: self.x.signum(),
            y: self.y.signum(),
        }
    }

    fn is_positive(&self) -> bool {
        self.x.is_positive() && self.y.is_positive()
    }

    fn is_negative(&self) -> bool {
        self.x.is_negative() && self.y.is_negative()
    }
}

impl<T: Sub<Output = T> + Copy> Sub<Self> for Vector2<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl<T: Sub<Output = T> + Copy> Sub<T> for Vector2<T> {
    type Output = Self;

    fn sub(self, rhs: T) -> Self {
        Self {
            x: self.x - rhs,
            y: self.y - rhs,
        }
    }
}

impl<T: SubAssign + Copy> SubAssign<T> for Vector2<T> {
    fn sub_assign(&mut self, rhs: T) {
        self.x -= rhs;
        self.y -= rhs;
    }
}

impl<T: Zero + Copy> Zero for Vector2<T> {
    fn zero() -> Self {
        Self {
            x: T::zero(),
            y: T::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y.is_zero()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn length_squared() {
        let v = Vector2::new(3, 4);
        assert_eq!(v.length_squared(), 25)
    }

    #[test]
    fn length() {
        let v = Vector2::new(3, 4);
        assert_eq!(v.length(), 5.0)
    }

    #[test]
    fn add() {
        assert_eq!(Vector2::new(1, 2) + Vector2::new(1, 2), Vector2::new(2, 4));
        assert_eq!(Vector2::new(1, 2) + 2, Vector2::new(3, 4));
    }

    #[test]
    fn add_assign() {
        let mut v = Vector2::new(1, 2);
        v += 5;
        assert_eq!(v, Vector2::new(6, 7));

        v += Vector2::new(1, 1);
        assert_eq!(v, Vector2::new(7, 8));
    }

    #[test]
    fn div() {
        assert_eq!(Vector2::new(4.0, 8.0) / 2.0, Vector2::new(2.0, 4.0));
    }

    #[test]
    fn div_assign() {
        let mut v = Vector2::new(25, 50);
        v /= 5;
        assert_eq!(v, Vector2::new(5, 10));
    }

    #[test]
    fn add_scalar() {
        assert_eq!(Vector2::new(1, 2) + 5, Vector2::new(6, 7));
    }

    #[test]
    fn index() {
        let v = Vector2::new(1, 2);
        assert_eq!(v[0], 1);
        assert_eq!(v[1], 2);
    }

    #[test]
    fn index_mut() {
        let mut v = Vector2::new(1, 2);
        let mut x = v[0];
        let mut y = v[1];

        x = 5;
        y = 10;
        assert_eq!(v[0], 1);
        assert_eq!(v[1], 2);

        v[0] = 3;
        assert_eq!(v[0], 3);
    }

    #[test]
    fn is_negative() {
        let v = Vector2::new(-1, 1);
        let is_neg = v.is_negative();

        assert_eq!(is_neg.x, true);
        assert_eq!(is_neg.y, false);
    }

    #[test]
    fn mul() {
        let v = Vector2::new(1, 2);
        assert_eq!(v * 5, Vector2::new(5, 10));
    }

    #[test]
    fn mul_assign() {
        let mut v = Vector2::new(1, 2);
        v *= 5;
        assert_eq!(v, Vector2::new(5, 10));
    }

    #[test]
    fn sub() {
        let v = Vector2::new(1, 2);
        assert_eq!(v - 5, Vector2::new(-4, -3));
        assert_eq!(v - Vector2::new(1, 2), Vector2::new(0, 0));
    }

    #[test]
    fn sub_assign() {
        let mut v = Vector2::new(1, 2);
        v -= 5;
        assert_eq!(v, Vector2::new(-4, -3));
    }
}
