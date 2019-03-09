use crate::interaction::SurfaceInteraction;

pub trait Texture<T> {
    fn evaluate(&self, si: &SurfaceInteraction) -> T;
}

pub mod constant;

pub use constant::*;
