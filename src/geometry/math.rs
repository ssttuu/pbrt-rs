use num::Num;
use std::f64::NAN;
use core::intrinsics;
use crate::types::{Int, Float};


pub trait Sqrt {
    fn sqrt(self) -> f64;
}

macro_rules! sqrt_impl {
    ( $t:ty, $zero:expr ) => {
        impl Sqrt for $t {
            fn sqrt(self) -> f64 {
                if self < $zero {
                    NAN
                } else {
                    unsafe { intrinsics::sqrtf64(self as f64) }
                }
            }
        }
    }
}


sqrt_impl!(u8, 0);
sqrt_impl!(u16, 0);
sqrt_impl!(u32, 0);
sqrt_impl!(u64, 0);
#[cfg(has_i128)]
sqrt_impl!(u128, 0);
sqrt_impl!(usize, 0);
sqrt_impl!(i8, 0);
sqrt_impl!(i16, 0);
sqrt_impl!(i32, 0);
sqrt_impl!(i64, 0);
#[cfg(has_i128)]
sqrt_impl!(i128, 0);
sqrt_impl!(isize, 0);

sqrt_impl!(f32, 0.0);
sqrt_impl!(f64, 0.0);