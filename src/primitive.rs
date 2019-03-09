use crate::bounds::Bounds3;
use crate::bounds::Bounds3f;
use crate::interaction::SurfaceInteraction;
use crate::intersect::Intersect;
use crate::lights::area::AreaLight;
use crate::materials::Material;
use crate::materials::TransportMode;
use crate::medium::MediumInterface;
use crate::ray::Ray;
use crate::shape::Shape;
use crate::types::Float;

pub trait Primitive: Intersect + Material {
    //    fn get_area_light(&self) -> AreaLight;
    fn get_material(&self) -> Option<&Material> {
        None
    }
    fn world_bound(&self) -> Bounds3f;
}

pub struct GeometricPrimitive {
    shape: Box<Shape>,
    material: Box<Material>,
    //    area_light: AreaLight,
    medium_interface: Option<Box<MediumInterface>>,
}

impl GeometricPrimitive {
    pub fn new(shape: Box<Shape>, material: Box<Material>) -> Self {
        Self {
            shape,
            material,
            medium_interface: None,
        }
    }
}

impl Intersect for GeometricPrimitive {
    fn intersect<'ray, 'prim>(
        &'prim self,
        ray: &mut Ray,
        si: &'ray mut SurfaceInteraction<'prim>,
    ) -> Option<Float> {
        let mut t_hit = 0.0;
        if !self.shape.intersect(ray, &mut t_hit, si, true) {
            return None;
        }

        ray.time_max = t_hit;

        si.primitive = Some(self);

        // TODO:
        //        si.medium_interface = match &self.medium_interface {
        //            Some(mi) => {
        //                if mi.is_medium_transition() {
        //                    Some(Box::new(mi))
        //                } else {
        //                    Some(Box::new(&MediumInterface::new()))
        //                }
        //            }
        //            None => None,
        //        };

        Some(t_hit)
    }

    fn intersect_p(&self, r: &Ray) -> bool {
        self.shape.intersect_p(r, true)
    }
}

impl Material for GeometricPrimitive {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    ) {
        self.material
            .compute_scattering_functions(si, mode, allow_multiple_lobes)
    }
}

impl Primitive for GeometricPrimitive {
    fn world_bound(&self) -> Bounds3f {
        self.shape.world_bound()
    }
}
