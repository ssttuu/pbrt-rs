use crate::bounds::Bounds3;
use crate::interaction::Interaction;
use crate::interaction::SurfaceInteraction;
use crate::ray::Ray;
use crate::types::Float;
use crate::types::Int;
use crate::Point2f;
use crate::Point3f;
use crate::Vector3f;
use crate::geometry::vector::distance_squared;
use crate::geometry::vector::abs_dot;
use crate::math::consts;
use crate::math;

pub trait Area {
    fn area(&self) -> Float;
}

pub trait Shape: Area {
    fn intersect<'prim>(
        &'prim self,
        ray: &mut Ray,
        t_hit: &mut Float,
        si: &mut SurfaceInteraction<'prim>,
        test_alpha_texture: bool,
    ) -> bool;
    fn intersect_p(&self, ray: &Ray, test_alpha_texture: bool) -> bool;
    fn object_bound(&self) -> Bounds3;
    fn pdf(&self, int: &Interaction) -> Float;
    fn pdf_wi(&self, int: &Interaction, wi: &mut Vector3f) -> Float;
    fn reverse_orientation(&self) -> bool;
    fn sample<'prim>(&'prim self, u: &Point2f, pdf: &mut Float) -> Box<Interaction<'prim> + 'prim>;
    fn sample_at_interaction<'prim>(
        &'prim self,
        int: &Interaction,
        u: &Point2f,
        pdf: &mut Float,
    ) -> Box<Interaction<'prim> + 'prim>;
    fn solid_angle(&self, p: &Point3f, n_samples: Int) -> Float;
    fn transform_swaps_handedness(&self) -> bool;
    fn world_bound(&self) -> Bounds3;
}

pub fn pdf(shape: &Shape, int: &Interaction) -> Float {
    1. / shape.area()
}

pub fn pdf_wi(shape: &Shape, int: &Interaction, wi: &mut Vector3f) -> Float {
    let mut ray = int.spawn_ray(wi.clone());

    let mut t_hit = Float::default();
    let mut intersect_light = SurfaceInteraction::default();
    let intersects = shape.intersect(&mut ray, &mut t_hit, &mut intersect_light, false);
    if !intersects {
        return 0.;
    }

    let mut pdf = distance_squared(int.get_point(), intersect_light.point) / abs_dot(&intersect_light.normal, &-*wi) * shape.area();
    if math::is_inf(pdf) {
        pdf = 0.;
    }
    pdf
}