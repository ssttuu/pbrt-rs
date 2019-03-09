use crate::types::Number;
use num::Num;
use num::ToPrimitive;
use std::ops;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Sub;

pub trait Cross {
    fn cross(&self, other: &Self) -> Self;
}

pub trait Dot<T: Num> {
    fn dot(&self, other: &Self) -> T;
    fn abs_dot(&self, other: &Self) -> T;
}

pub trait Length<T>
where
    T: Num + Copy + ToPrimitive + PartialOrd,
{
    fn length_squared(&self) -> T;
    fn length(&self) -> f64;
}

pub fn dot<T, V>(a: &T, b: &T) -> V
where
    T: Dot<V> + Copy,
    V: Num,
{
    a.dot(b)
}

pub fn abs_dot<T, V>(a: &T, b: &T) -> V
where
    T: Dot<V> + Copy,
    V: Num,
{
    a.abs_dot(b)
}

pub fn cross<T>(a: &T, b: &T) -> T
where
    T: Cross,
{
    a.cross(b)
}

pub fn distance_squared<T, V>(a: T, b: T) -> V
where
    T: Sub<Output = T> + Length<V>,
    V: Num + Copy + ToPrimitive + PartialOrd,
{
    (b - a).length_squared()
}

pub fn distance<T, V>(a: T, b: T) -> f64
where
    T: Sub<Output = T> + Length<V>,
    V: Num + Copy + ToPrimitive + PartialOrd,
{
    (b - a).length()
}

pub fn face_forward<T, V>(a: T, b: T) -> T
where
    T: Dot<V> + ops::Neg<Output = T> + Copy,
    V: Number,
{
    if a.dot(&b) < V::zero() {
        -a
    } else {
        a
    }
}
