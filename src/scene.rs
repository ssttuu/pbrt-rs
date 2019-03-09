use crate::bounds::Bounds3f;
use crate::interaction::Interaction;
use crate::interaction::SurfaceInteraction;
use crate::intersect::Intersect;
use crate::light::Light;
use crate::primitive::Primitive;
use crate::ray::Ray;
use crate::sampler::Sampler;
use crate::spectrum::Spectrum;
use crate::types::Float;
use std::rc::Rc;

pub struct Scene {
    pub lights: Vec<Rc<Light>>,
    pub infinite_lights: Vec<Rc<Light>>,

    pub aggregate: Rc<Primitive>,
    pub world_bound: Bounds3f,
}

impl Scene {
    pub fn intersect_tr<'ray, 'prim>(
        &'prim self,
        r: &mut Ray,
        si: &'ray mut SurfaceInteraction<'prim>,
        sampler: &Sampler,
        transmittance: &mut Spectrum,
    ) -> bool {
        transmittance.set_all(1.0);

        loop {
            let t_hit = self.intersect(r, si);

            // TODO
            //            if let Some(medium) = r.medium {
            //                *transmittance *= medium.tr(r, sampler);
            //            }

            if t_hit.is_none() {
                return false;
            }
            if let Some(prim) = si.primitive {
                if prim.get_material().is_none() {
                    return true;
                }
            }
            *r = si.spawn_ray(r.direction);
        }
    }
}

impl Intersect for Scene {
    fn intersect<'ray, 'prim>(
        &'prim self,
        r: &mut Ray,
        si: &'ray mut SurfaceInteraction<'prim>,
    ) -> Option<Float> {
        self.aggregate.intersect(r, si)
    }

    fn intersect_p(&self, r: &Ray) -> bool {
        self.aggregate.intersect_p(r)
    }
}
