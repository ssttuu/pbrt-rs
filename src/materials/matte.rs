use crate::interaction::SurfaceInteraction;
use crate::materials::Material;
use crate::materials::TransportMode;
use crate::spectrum::Spectrum;
use crate::texture::Texture;
use crate::types::Float;
use std::ops::Deref;
use std::rc::Rc;

pub struct MatteMaterial {
    kd: Rc<Texture<Spectrum>>,
    sigma: Rc<Texture<Float>>,
    bump_map: Option<Rc<Texture<Float>>>,
}

impl MatteMaterial {}

impl Material for MatteMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    ) {
        if let Some(bump) = &self.bump_map {
            self.bump(bump.deref(), si)
        }
    }
}
