use num;
use std::ops;
use crate::math::Sqrt;

pub type Int = i64;
pub type Float = f64;

pub trait ToFloat {
    fn to_float(self) -> Float;
}

macro_rules! impl_to_float {
    ($t:ty) => {
        impl ToFloat for $t {
            fn to_float(self) -> Float { self as Float }
        }
    }
}

impl_to_float!(i32);
impl_to_float!(i64);
impl_to_float!(f32);
impl_to_float!(f64);

pub fn to_float<T>(v: T) -> Float where T: ToFloat {
    v.to_float()
}


pub trait Number:
    num::Num
    + Copy
    + num::ToPrimitive
    + PartialOrd
    + Sqrt
    + ops::MulAssign
    + ops::DivAssign
    + ops::AddAssign
    + ops::SubAssign
    + ToFloat {


}

impl Number for i32 {}

impl Number for i64 {}

impl Number for f32 {}

impl Number for f64 {}
