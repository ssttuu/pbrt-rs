use crate::interaction::Interaction;
use crate::interaction::SurfaceInteraction;
use crate::intersect::Intersect;
use crate::ray::Ray;
use crate::sampler::Sampler;
use crate::scene::Scene;
use crate::spectrum::Spectrum;
use crate::types::Float;
use crate::Normal3f;
use crate::Point2f;
use crate::Vector3f;

pub enum LightType {
    DeltaPosition = 1 << 1,
    DeltaDirection = 1 << 2,
    Area = 1 << 3,
    Infinite = 1 << 4,
}

pub trait Light {
    fn power(&self) -> Spectrum;
    fn le(&self, ray: &mut Ray) -> Spectrum;
    fn pdf_le(
        &self,
        ray: &mut Ray,
        light_normal: &mut Normal3f,
        pdf_pos: &mut Float,
        pdf_dir: &mut Float,
    );
    fn pdf_li(&self, int: &Interaction, wi: &Vector3f) -> Float;
    fn sample_le(
        &self,
        u1: &Point2f,
        u2: &Point2f,
        time: Float,
        ray: &mut Ray,
        light_normal: &mut Normal3f,
        pdf_pos: &mut Float,
        pdf_dir: &mut Float,
    ) -> Spectrum;
    fn sample_li(
        &self,
        int: &Interaction,
        u: &Point2f,
        wi: &mut Vector3f,
        pdf: &mut Float,
        vis: &VisibilityTester,
    ) -> Spectrum;
}

pub struct VisibilityTester<'ray, 'prim> {
    p0: &'ray Interaction<'prim>,
    p1: &'ray Interaction<'prim>,
}

impl<'ray, 'prim> VisibilityTester<'ray, 'prim> {
    pub fn unoccluded(&self, scene: &'prim Scene) -> bool {
        !scene.intersect_p(&self.p0.spawn_ray_to(self.p1))
    }

    pub fn tr(&self, scene: &'prim Scene, sampler: &'prim Sampler) -> Spectrum {
        let mut ray = self.p0.spawn_ray_to(self.p1);
        let mut tr = Spectrum::new(1.0);

        loop {
            let mut isect = SurfaceInteraction::default();
            let hit_surface = scene.intersect(&mut ray, &mut isect);

            if hit_surface.is_some() {
                if let Some(prim) = isect.primitive {
                    if let Some(material) = prim.get_material() {
                        return Spectrum::new(0.0);
                    }
                }
            }

            //            if let Some(medium) = &ray.medium {
            //                tr *= medium.tr(&mut ray, sampler);
            //            }

            if hit_surface.is_none() {
                break;
            }
            ray = isect.spawn_ray_to(self.p1);
        }

        tr
    }
}
