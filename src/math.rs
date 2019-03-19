use crate::types::Float;
use crate::types::Int;
use crate::types::Number;
use core::intrinsics;
use num::Num;
use num::ToPrimitive;
use crate::efloat::EFloat;
use std::f64;

pub mod consts {
    use super::next_float_up;
    use super::Float;
    use std::f64;

    pub static INFINITY: Float = f64::INFINITY;
    pub static NAN: Float = f64::NAN;
    pub static PI: Float = f64::consts::PI;
    pub const EPSILON: Float = f64::EPSILON;
}

pub fn acos(v: Float) -> Float {
    v.acos()
}

pub fn atan(v: Float) -> Float {
    v.atan()
}

pub fn atan2(v: Float, other: Float) -> Float {
    v.atan2(other)
}

pub fn ceil(v: Float) -> Float {
    v.ceil()
}

pub fn clamp(v: Float, low: Float, high: Float) -> Float {
    if v < low {
        low
    } else if v > high {
        high
    } else {
        v
    }
}

pub fn cos(v: Float) -> Float {
    v.cos()
}

pub fn degrees(radians: Float) -> Float {
    180.0 / consts::PI * radians
}

pub fn floor(v: Float) -> Float {
    v.floor()
}

pub fn find_interval(size: Int, pred: &Fn(Int) -> bool) -> Int {
    let mut first = 0;
    let mut length = size;
    while length > 0 {
        let half = length >> 1;
        let middle = first + half;
        // bisect range based on value of pred at middle
        if pred(middle) {
            first = middle + 1;
            length -= half + 1;
        } else {
            length = half;
        }
    }

    clamp((first - 1) as Float, 0.0, (size - 2) as Float) as Int
}

pub fn gamma(n: Float) -> Float {
    (n * consts::EPSILON) / (1.0 - n * consts::EPSILON)
}

pub fn is_inf(v: Float) -> bool {
    match v {
        f64::INFINITY => true,
        _ => false,
    }
}

pub fn lerp<T: Number>(t: T, v1: T, v2: T) -> T {
    (T::one() - t) * v1 + t * v2
}

pub fn max<T: Number>(v1: T, v2: T) -> T {
    if v1 >= v2 {
        v1
    } else {
        v2
    }
}

pub fn min<T: Number>(v1: T, v2: T) -> T {
    if v1 <= v2 {
        v1
    } else {
        v2
    }
}

pub fn next_float_up(v: Float) -> Float {
    Float::from_bits(v.to_bits() + 1)
}

pub fn next_float_down(v: Float) -> Float {
    Float::from_bits(v.to_bits() - 1)
}

pub fn radians(degrees: Float) -> Float {
    consts::PI / 180.0 * degrees
}

//pub fn quadratic(a: Float, b: Float, c: Float) -> Option<(EFloat, EFloat)> {
//    let discriminant = b * b - 4. * a * c;
//    if discriminant < 0. {
//        return None;
//    }
//    let root_discriminant = discriminant.sqrt();
//    let float_root_discriminant =
//}

pub fn sin(v: Float) -> Float {
    v.sin()
}

pub fn tan(v: Float) -> Float {
    v.tan()
}

pub trait Sqrt {
    fn sqrt(self) -> f64;
}

macro_rules! sqrt_impl {
    ( $t:ty, $zero:expr ) => {
        impl Sqrt for $t {
            fn sqrt(self) -> f64 {
                if self < $zero {
                    consts::NAN
                } else {
                    unsafe { intrinsics::sqrtf64(self as f64) }
                }
            }
        }
    };
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

pub fn sqrt<T>(v: T) -> f64 where T: Sqrt {
    v.sqrt()
}
