use crate::bounds::Bounds3;
use crate::efloat::EFloat;
use crate::geometry::coordinate_system;
use crate::geometry::spherical_direction_xyz;
use crate::geometry::vector::abs_dot;
use crate::geometry::vector::cross;
use crate::geometry::vector::distance;
use crate::geometry::vector::distance_squared;
use crate::geometry::vector::dot;
use crate::geometry::vector::Cross;
use crate::geometry::vector::Length;
use crate::interaction::Interaction;
use crate::interaction::SurfaceInteraction;
use crate::math;
use crate::math::consts;
use crate::math::gamma;
use crate::ray::offset_ray_origin;
use crate::ray::Ray;
use crate::sampling::uniform_cone_pdf;
use crate::sampling::uniform_sample_sphere;
use crate::shape;
use crate::shape::Area;
use crate::shape::Shape;
use crate::transform::Transform;
use crate::types::Float;
use crate::types::Int;
use crate::Point2f;
use crate::Point3f;
use crate::Vector3f;
use num::abs;
use num::Zero;

pub struct Sphere {
    pub object_to_world: Transform,
    pub world_to_object: Transform,
    pub reverse_orientation: bool,
    pub transform_swaps_handedness: bool,

    pub radius: Float,
    pub z_min: Float,
    pub z_max: Float,
    pub theta_min: Float,
    pub theta_max: Float,
    pub phi_max: Float,
}

impl Sphere {
    pub fn new() -> Self {
        Self { ..Self::default() }
    }
}

impl Default for Sphere {
    fn default() -> Self {
        Self {
            object_to_world: Transform::default(),
            world_to_object: Transform::default(),
            reverse_orientation: false,
            transform_swaps_handedness: false,
            radius: 1.,
            z_min: 0.,
            z_max: 0.,
            theta_min: 0.,
            theta_max: 0.,
            phi_max: 360.,
        }
    }
}

impl Area for Sphere {
    fn area(&self) -> Float {
        self.phi_max * self.radius * (self.z_max - self.z_min)
    }
}

impl Shape for Sphere {
    fn intersect<'prim>(
        &'prim self,
        ray: &mut Ray,
        t_hit: &mut Float,
        si: &mut SurfaceInteraction<'prim>,
        test_alpha_texture: bool,
    ) -> bool {
        let ray = self.world_to_object.transform_ray(ray);

        let origin_error = ray.origin_error.unwrap();
        let ox = EFloat::new(ray.origin.x, origin_error.x);
        let oy = EFloat::new(ray.origin.y, origin_error.y);
        let oz = EFloat::new(ray.origin.z, origin_error.z);

        let direction_error = ray.direction_error.unwrap();
        let dx = EFloat::new(ray.direction.x, direction_error.x);
        let dy = EFloat::new(ray.direction.y, direction_error.y);
        let dz = EFloat::new(ray.direction.z, direction_error.z);

        let a = dx * dx + dy * dy + dz * dz;
        let b = dx * ox + dy * oy + dz * oz * 2.0;
        let c = ox * ox + oy * oy + oz * oz
            - EFloat::new(self.radius, 0.) * EFloat::new(self.radius, 0.);

        let intercepts = EFloat::quadratic(a, b, c);
        if intercepts.is_none() {
            return false;
        }

        let (t0, t1) = intercepts.unwrap();
        if t0.high > ray.time_max || t1.low <= 0. {
            return false;
        }

        let mut t_shape_hit = t0;
        if t_shape_hit.low <= 0. {
            t_shape_hit = t1;
            if t_shape_hit.high > ray.time_max {
                return false;
            }
        }

        let mut p_hit = ray.point_at(t_shape_hit.value);
        p_hit *= self.radius / distance(p_hit, Point3f::zero());
        if p_hit.x == 0. && p_hit.y == 0. {
            p_hit.x = 1e-5 * self.radius;
        }

        let mut phi = math::atan2(p_hit.y, p_hit.x);
        if phi < 0. {
            phi *= 2. * consts::PI;
        }

        if (self.z_min > -self.radius && p_hit.z < self.z_min)
            || (self.z_max < self.radius && p_hit.z > self.z_max)
            || phi > self.phi_max
        {
            if t_shape_hit == t1 {
                return false;
            }
            if t1.high > ray.time_max {
                return false;
            }
            t_shape_hit = t1;

            p_hit = ray.point_at(t_shape_hit.value);

            p_hit *= self.radius / distance(p_hit, Point3f::zero());
            if p_hit.x == 0. && p_hit.y == 0. {
                p_hit.x = 1e-5 * self.radius;
            }
            phi = math::atan2(p_hit.y, p_hit.x);
            if phi < 0. {
                phi += 2. * consts::PI;
            }
            if (self.z_min > -self.radius && p_hit.z < self.z_min)
                || (self.z_max < self.radius && p_hit.z > self.z_max)
                || phi > self.phi_max
            {
                return false;
            }
        }

        let u = phi / self.phi_max;
        let theta = math::acos(math::clamp(p_hit.z / self.radius, -1., 1.));
        let v = (theta - self.theta_min) / (self.theta_max - self.theta_min);

        let z_radius = math::sqrt(p_hit.x * p_hit.x + p_hit.y * p_hit.y);
        let inv_z_radius = 1. / z_radius;
        let cos_phi = p_hit.x * inv_z_radius;
        let sin_phi = p_hit.y * inv_z_radius;
        let dpdu = Vector3f::new(-self.phi_max * p_hit.y, self.phi_max * p_hit.x, 0.);
        let dpdv = Vector3f::new(
            p_hit.z * cos_phi,
            p_hit.z * sin_phi,
            -self.radius * math::sin(theta),
        ) * (self.theta_max - self.theta_min);

        let d2pduu = Vector3f::new(p_hit.x, p_hit.y, 0.) * -self.phi_max * self.phi_max;
        let d2pduv = Vector3f::new(-sin_phi, cos_phi, 0.)
            * (self.theta_max - self.theta_min)
            * p_hit.z
            * self.phi_max;
        let d2pdvv = p_hit * -(self.theta_max - self.theta_min) * (self.theta_max - self.theta_min);

        let E = dot(&dpdu, &dpdu);
        let F = dot(&dpdu, &dpdv);
        let G = dot(&dpdv, &dpdv);
        let N = dpdu.cross(&dpdv).normalized();
        let e = dot(&N, &d2pduu);
        let f = dot(&N, &d2pduv);
        let g = dot(&N, &d2pdvv);

        let inv_EGF2 = 1. / (E * G - F * F);
        let dndu = dpdu * (f * F - e * G) * inv_EGF2 + dpdv * (e * F - f * E) * inv_EGF2;
        let dndv = dpdu * (g * F - f * G) * inv_EGF2 + dpdv * (f * F - g * E) * inv_EGF2;

        let p_error = abs(p_hit) * gamma(5.);

        *si = self
            .object_to_world
            .transform_surface_interaction(&SurfaceInteraction::new(
                p_hit,
                p_error,
                Point2f::new(u, v),
                -ray.direction,
                dpdu,
                dpdv,
                dndu,
                dndv,
                ray.time,
                self,
                None,
            ));
        *t_hit = t_shape_hit.value;

        true
    }
    fn intersect_p(&self, ray: &Ray, test_alpha_texture: bool) -> bool {
        let ray = self.world_to_object.transform_ray(ray);

        let origin_error = ray.origin_error.unwrap();
        let ox = EFloat::new(ray.origin.x, origin_error.x);
        let oy = EFloat::new(ray.origin.y, origin_error.y);
        let oz = EFloat::new(ray.origin.z, origin_error.z);

        let direction_error = ray.direction_error.unwrap();
        let dx = EFloat::new(ray.direction.x, direction_error.x);
        let dy = EFloat::new(ray.direction.y, direction_error.y);
        let dz = EFloat::new(ray.direction.z, direction_error.z);

        let a = dx * dx + dy * dy + dz * dz;
        let b = dx * ox + dy * oy + dz * oz * 2.0;
        let c = ox * ox + oy * oy + oz * oz
            - EFloat::new(self.radius, 0.) * EFloat::new(self.radius, 0.);

        let intercepts = EFloat::quadratic(a, b, c);
        if intercepts.is_none() {
            return false;
        }

        let (t0, t1) = intercepts.unwrap();
        if t0.high > ray.time_max || t1.low <= 0. {
            return false;
        }

        let mut t_shape_hit = t0;
        if t_shape_hit.low <= 0. {
            t_shape_hit = t1;
            if t_shape_hit.high > ray.time_max {
                return false;
            }
        }

        let mut p_hit = ray.point_at(t_shape_hit.value);
        p_hit *= self.radius / distance(p_hit, Point3f::zero());
        if p_hit.x == 0. && p_hit.y == 0. {
            p_hit.x = 1e-5 * self.radius;
        }

        let mut phi = math::atan2(p_hit.y, p_hit.x);
        if phi < 0. {
            phi *= 2. * consts::PI;
        }

        if (self.z_min > -self.radius && p_hit.z < self.z_min)
            || (self.z_max < self.radius && p_hit.z > self.z_max)
            || phi > self.phi_max
        {
            if t_shape_hit == t1 {
                return false;
            }
            if t1.high > ray.time_max {
                return false;
            }
            t_shape_hit = t1;

            p_hit = ray.point_at(t_shape_hit.value);

            p_hit *= self.radius / distance(p_hit, Point3f::zero());
            if p_hit.x == 0. && p_hit.y == 0. {
                p_hit.x = 1e-5 * self.radius;
            }
            phi = math::atan2(p_hit.y, p_hit.x);
            if phi < 0. {
                phi += 2. * consts::PI;
            }
            if (self.z_min > -self.radius && p_hit.z < self.z_min)
                || (self.z_max < self.radius && p_hit.z > self.z_max)
                || phi > self.phi_max
            {
                return false;
            }
        }

        true
    }

    fn object_bound(&self) -> Bounds3 {
        Bounds3 {
            min: Point3f::new(-self.radius, -self.radius, -self.z_min),
            max: Point3f::new(self.radius, self.radius, self.z_max),
        }
    }

    fn pdf(&self, int: &Interaction) -> Float {
        Shape::pdf(self, int)
    }

    fn pdf_wi(&self, int: &Interaction, wi: &mut Vector3f) -> Float {
        let (p_center, p_center_error) = self
            .object_to_world
            .transform_point(&Point3f::zero(), &Point3f::zero());
        let p_origin = offset_ray_origin(
            &int.get_point(),
            &int.get_point_error(),
            &int.get_normal(),
            &(p_center - int.get_point()),
        );
        if distance_squared(p_origin, p_center) <= self.radius * self.radius {
            return shape::pdf_wi(self, int, wi);
        }

        let sin_theta_max_2 =
            self.radius * self.radius / distance_squared(int.get_point(), p_center);
        let cos_theta_max = math::sqrt(math::max(0., 1. - sin_theta_max_2));

        uniform_cone_pdf(cos_theta_max)
    }
    fn reverse_orientation(&self) -> bool {
        self.reverse_orientation
    }
    fn sample<'prim>(&'prim self, u: &Point2f, pdf: &mut Float) -> Box<Interaction<'prim> + 'prim> {
        let mut p_obj = uniform_sample_sphere(u) * self.radius;

        let mut it = SurfaceInteraction::default();
        it.normal = self.object_to_world.transform_normal(&p_obj).normalized();
        if self.reverse_orientation {
            it.normal = -it.normal;
        }

        p_obj = p_obj * self.radius / distance(p_obj, Point3f::zero());
        let p_obj_error = abs(p_obj) * gamma(5.);
        let (point, point_error) = self.object_to_world.transform_point(&p_obj, &p_obj_error);

        it.point = point;
        it.point_error = point_error;
        *pdf = 1. / self.area();

        Box::new(it)
    }
    fn sample_at_interaction<'prim>(
        &'prim self,
        int: &Interaction,
        u: &Point2f,
        pdf: &mut Float,
    ) -> Box<Interaction<'prim> + 'prim> {
        let (p_center, p_center_error) = self
            .object_to_world
            .transform_point(&Point3f::zero(), &Point3f::zero());

        let p_origin = offset_ray_origin(
            &int.get_point(),
            &int.get_point_error(),
            &int.get_normal(),
            &(p_center - int.get_point()),
        );
        if distance_squared(p_origin, p_center) <= self.radius * self.radius {
            let mut pdf = Float::default();
            let intr = self.sample(u, &mut pdf);
            let mut wi = intr.get_point() - int.get_point();
            if wi.length_squared() == 0. {
                pdf = 0.;
            } else {
                wi = wi.normalized();
                pdf *= distance_squared(int.get_point(), intr.get_point())
                    / abs_dot(&intr.get_normal(), &-wi);
            }

            if math::is_inf(pdf) {
                pdf = 0.;
            }

            return intr;
        }

        let wc = (p_center - int.get_point()).normalized();
        let (wc_x, wc_y) = coordinate_system(wc);

        let radius_squared = self.radius * self.radius;

        let sin_theta_max_2 = radius_squared / distance_squared(int.get_point(), p_center);
        let cos_theta_max = math::sqrt(math::max(0., 1. - sin_theta_max_2));
        let cos_theta = (1. - u.x) + u.x * cos_theta_max;
        let sin_theta = math::sqrt(math::max(0., 1. - cos_theta * cos_theta));
        let phi = u.y * 2. * consts::PI;

        let dc = distance(int.get_point(), p_center);
        let ds = dc * cos_theta
            - math::sqrt(math::max(
                0.,
                radius_squared - (dc * dc) * (sin_theta * sin_theta),
            ));
        let cos_alpha = (dc * dc + radius_squared - ds * ds) / (2. * dc * self.radius);
        let sin_alpha = math::sqrt(math::max(0., 1. - cos_alpha * cos_alpha));

        let n_world = spherical_direction_xyz(sin_alpha, cos_alpha, phi, &-wc_x, &-wc_y, &-wc);
        let p_world = p_center + n_world * self.radius;

        // TODO: pbrt implementation creates a base class Interaction.  Not sure what to do here.
        // Surface or Medium?
        let mut it = Box::new(SurfaceInteraction::default());
        it.point = p_world;
        it.point_error = abs(p_world) * gamma(5.);
        it.normal = n_world;
        if self.reverse_orientation {
            it.normal *= -1.;
        }

        return it;
    }

    fn solid_angle(&self, p: &Point3f, n_samples: Int) -> Float {
        let (p_center, p_center_error) = self
            .object_to_world
            .transform_point(&Point3f::zero(), &Point3f::zero());
        if distance_squared(p.clone(), p_center.clone()) <= self.radius * self.radius {
            return 4. * consts::PI;
        }

        let sin_theta_2 = self.radius * self.radius / distance_squared(p.clone(), p_center.clone());
        let cos_theta = math::sqrt(math::max(0., 1. - sin_theta_2));

        2. * consts::PI * (1. - cos_theta)
    }

    fn transform_swaps_handedness(&self) -> bool {
        self.transform_swaps_handedness
    }

    fn world_bound(&self) -> Bounds3 {
        self.object_to_world.transform_bounds(&self.object_bound())
    }
}
