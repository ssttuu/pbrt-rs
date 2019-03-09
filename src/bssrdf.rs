use crate::interaction::SurfaceInteraction;
use crate::scene::Scene;
use crate::spectrum::Spectrum;
use crate::types::Float;
use crate::Point2f;
use crate::Vector3f;

pub trait BSSRDF {
    fn s(&self, si: &SurfaceInteraction, wi: Vector3f) -> Spectrum;
    fn sample_s(
        &self,
        scene: &Scene,
        u1: Float,
        u2: Point2f,
        si: &SurfaceInteraction,
        pdf: &mut Float,
    ) -> Spectrum;
}

struct Base {
    //    pi: &'a SurfaceInteraction,
}
