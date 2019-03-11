use crate::bounds::Bounds3;
use crate::interaction::Interaction;
use crate::interaction::SurfaceInteraction;
use crate::ray::Ray;
use crate::types::Float;
use crate::types::Int;
use crate::Point2f;
use crate::Point3f;
use crate::Vector3f;

pub trait Area {
    fn area(&self) -> Float;
}

pub trait Shape: Area {
    fn intersect(
        &self,
        ray: &mut Ray,
        t_hit: &mut Float,
        si: &SurfaceInteraction,
        test_alpha_texture: bool,
    ) -> bool;
    fn intersect_p(&self, ray: &Ray, test_alpha_texture: bool) -> bool;
    fn object_bound(&self) -> Bounds3;
    fn pdf(&self, int: &Interaction) -> Float;
    fn pdf_wi(&self, int: &Interaction, wi: &mut Vector3f) -> Float;
    fn reverse_orientation(&self) -> bool;
    fn sample(&self, u: &Point2f, pdf: &mut Float) -> Box<Interaction>;
    fn sample_at_interaction(
        &self,
        int: &Interaction,
        u: &Point2f,
        pdf: &mut Float,
    ) -> Box<Interaction>;
    fn solid_angle(&self, p: &Point3f, n_samples: Int) -> Float;
    fn transform_swaps_handedness(&self) -> bool;
    fn world_bound(&self) -> Bounds3;
}
