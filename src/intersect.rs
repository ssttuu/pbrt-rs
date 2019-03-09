use crate::interaction::SurfaceInteraction;
use crate::ray::Ray;
use crate::types::Float;

pub trait Intersect {
    // returns THit
    fn intersect<'ray, 'prim>(
        &'prim self,
        r: &mut Ray,
        si: &'ray mut SurfaceInteraction<'prim>,
    ) -> Option<Float>;
    fn intersect_p(&self, r: &Ray) -> bool;
}
