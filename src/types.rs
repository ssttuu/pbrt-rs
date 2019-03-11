use crate::math::Sqrt;
use num;
use num::Signed;
use std::cmp;
use std::ops;

pub type Int = i64;
pub type Float = f64;

pub trait ToFloat {
    fn to_float(self) -> Float;
}

macro_rules! impl_to_float {
    ($t:ty) => {
        impl ToFloat for $t {
            fn to_float(self) -> Float {
                self as Float
            }
        }
    };
}

impl_to_float!(i32);
impl_to_float!(i64);
impl_to_float!(f32);
impl_to_float!(f64);

pub fn to_float<T>(v: T) -> Float
where
    T: ToFloat,
{
    v.to_float()
}

pub trait Number:
    num::Num
    + Copy
    + num::ToPrimitive
    + ops::AddAssign
    + ops::DivAssign
    + ops::MulAssign
    + ops::Neg
    + ops::SubAssign
    + PartialOrd
    + Signed
    + Sqrt
    + ToFloat
{
}

impl Number for i32 {}

impl Number for i64 {}

impl Number for f32 {}

impl Number for f64 {}
