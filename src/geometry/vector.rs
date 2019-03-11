use crate::types::Float;
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

pub trait Distance {
    fn distance(&self, other: &Self) -> Float;
}

pub trait Dot<T: Number> {
    fn dot(&self, other: &Self) -> T;
    fn abs_dot(&self, other: &Self) -> T;
}

pub trait Length<T>
where
    T: Number,
{
    fn length_squared(&self) -> T;
    fn length(&self) -> Float;
}

pub fn dot<T, V>(a: &T, b: &T) -> V
where
    T: Dot<V> + Copy,
    V: Number,
{
    a.dot(b)
}

pub fn abs_dot<T, V>(a: &T, b: &T) -> V
where
    T: Dot<V> + Copy,
    V: Number,
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
    V: Number,
{
    (b - a).length_squared()
}

pub fn distance<V>(a: V, b: V) -> Float
where
    V: Distance,
{
    a.distance(&b)
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
