use crate::bounds::Bounds3;
use crate::interaction::SurfaceInteraction;
use crate::intersect::Intersect;
use crate::lights::area::AreaLight;
use crate::materials::Material;
use crate::materials::TransportMode;
use crate::medium::MediumInterface;
use crate::ray::Ray;
use crate::shape::Shape;
use crate::transform::AnimatedTransform;
use crate::types::Float;

pub trait Primitive: Intersect + Material {
    //    fn get_area_light(&self) -> AreaLight;
    fn get_material(&self) -> Option<&Material> {
        None
    }
    fn world_bound(&self) -> Bounds3;
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
    fn world_bound(&self) -> Bounds3 {
        self.shape.world_bound()
    }
}

pub struct TransformedPrimitive {
    pub primitive: Box<Primitive>,
    pub primitive_to_world: AnimatedTransform,
}

impl TransformedPrimitive {
    pub fn new(primitive: Box<Primitive>, primitive_to_world: AnimatedTransform) -> Self {
        Self {
            primitive,
            primitive_to_world,
        }
    }
}

impl Intersect for TransformedPrimitive {
    fn intersect<'ray, 'prim>(
        &'prim self,
        r: &mut Ray,
        si: &'ray mut SurfaceInteraction<'prim>,
    ) -> Option<Float> {
        let interpolated_prim_to_world = self.primitive_to_world.interpolate(r.time);
        let ray = interpolated_prim_to_world.invert().transform_ray(r);

        let intersects = self.primitive.intersect(r, si);
        if intersects.is_none() {
            return None;
        }
        r.time_max = ray.time_max;

        if !interpolated_prim_to_world.is_identity() {
            *si = interpolated_prim_to_world.transform_surface_interaction(si);
        }

        Some(r.time)
    }

    fn intersect_p(&self, r: &Ray) -> bool {
        let interpolated_prim_to_world = self.primitive_to_world.interpolate(r.time);
        let ray = interpolated_prim_to_world.invert().transform_ray(r);
        self.primitive.intersect_p(&ray)
    }
}

impl Material for TransformedPrimitive {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    ) {
        panic!("should not be called")
    }
}

impl Primitive for TransformedPrimitive {
    fn world_bound(&self) -> Bounds3 {
        self.primitive_to_world
            .motion_bounds(&self.primitive.world_bound())
    }
}
