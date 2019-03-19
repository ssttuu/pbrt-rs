#![feature(core_intrinsics)]
#![feature(const_fn)]

use crate::geometry::Vector2;
use crate::geometry::Vector3;
use crate::types::Float;
use crate::types::Int;
use std::f32;
use std::f64;

pub mod bounds;
pub mod bssrdf;
pub mod efloat;
pub mod geometry;
pub mod interaction;
pub mod intersect;
pub mod light;
pub mod lights;
pub mod materials;
pub mod math;
pub mod medium;
pub mod primitive;
pub mod quaternion;
pub mod ray;
pub mod reflection;
pub mod sampler;
pub mod sampling;
pub mod scene;
pub mod shape;
pub mod sphere;
pub mod spectrum;
pub mod texture;
pub mod transform;
pub mod types;

pub type Point2<T> = Vector2<T>;
pub type Normal2<T> = Vector2<T>;

pub type Vector2i = Vector2<Int>;
pub type Vector2f = Vector2<Float>;
pub type Point2i = Vector2i;
pub type Point2f = Vector2f;
pub type Normal2i = Vector2i;
pub type Normal2f = Vector2f;

pub type Point3<T> = Vector3<T>;
pub type Normal3<T> = Vector3<T>;

pub type Vector3i = Vector3<Int>;
pub type Vector3f = Vector3<Float>;
pub type Point3i = Vector3i;
pub type Point3f = Vector3f;
pub type Normal3i = Vector3i;
pub type Normal3f = Vector3f;

static SHADOW_EPSILON: Float = 0.0001;
