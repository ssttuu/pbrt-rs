use crate::bssrdf::BSSRDF;
use crate::geometry::vector::cross;
use crate::geometry::vector::dot;
use crate::geometry::vector::face_forward;
use crate::geometry::vector::Cross;
use crate::medium::Medium;
use crate::medium::MediumInterface;
use crate::primitive::Primitive;
use crate::ray::offset_ray_origin;
use crate::ray::Ray;
use crate::reflection::BSDF;
use crate::shape::Shape;
use crate::types::Float;
use crate::Normal3f;
use crate::Point2f;
use crate::Point3f;
use crate::Vector3f;
use crate::SHADOW_EPSILON;

use crate::math;
use crate::math::consts;
use num::Zero;
use std::rc::Rc;

pub trait Interaction<'prim> {
    fn get_point(&self) -> Point3f;
    fn get_point_error(&self) -> Vector3f;
    fn get_time(&self) -> Float;
    fn get_normal(&self) -> Normal3f;
    fn get_medium(&self, w: Vector3f) -> Option<&'prim Rc<Medium>> {
        match self.get_medium_interface() {
            Some(interface) => {
                let medium = if dot(&w, &self.get_normal()) > 0.0 {
                    &interface.outside
                } else {
                    &interface.inside
                };

                match medium {
                    Some(m) => Some(&m),
                    None => None,
                }
            }
            None => None,
        }
    }
    fn get_medium_interface(&self) -> Option<&'prim MediumInterface>;
    fn set_medium_interface(&mut self, interface: Option<&'prim MediumInterface>);

    fn spawn_ray(&self, direction: Vector3f) -> Ray {
        let origin = offset_ray_origin(
            &self.get_point(),
            &self.get_point_error(),
            &self.get_normal(),
            &direction,
        );
        Ray::new(
            origin,
            direction,
            consts::INFINITY,
            self.get_time(),
            //            self.get_medium(direction),
        )
    }

    fn spawn_ray_to(&self, to: &Interaction<'prim>) -> Ray {
        let origin = offset_ray_origin(
            &self.get_point(),
            &self.get_point_error(),
            &self.get_normal(),
            &(to.get_point() - self.get_point()),
        );
        let target = offset_ray_origin(
            &to.get_point(),
            &to.get_point_error(),
            &to.get_normal(),
            &(origin - self.get_point()),
        );
        let direction = target - origin;
        Ray::new(
            origin,
            direction,
            1.0 - SHADOW_EPSILON,
            self.get_time(),
            //            self.get_medium(direction),
        )
    }
}

#[derive(Clone)]
pub struct Shading {
    pub dpdu: Vector3f,
    pub dpdv: Vector3f,
    pub dndu: Normal3f,
    pub dndv: Normal3f,
    pub normal: Normal3f,
}

#[derive(Clone)]
pub struct SurfaceInteraction<'prim> {
    pub point: Point3f,
    pub point_error: Vector3f,
    pub time: Float,
    pub wo: Vector3f,
    pub normal: Normal3f,
    pub medium_interface: Option<&'prim MediumInterface>,

    pub bsdf: Option<BSDF>,
    pub bssrdf: Option<&'prim BSSRDF>,
    pub dndu: Vector3f,
    pub dndv: Vector3f,
    pub dpdu: Vector3f,
    pub dpdv: Vector3f,
    pub dpdx: Vector3f,
    pub dpdy: Vector3f,
    pub dudx: Float,
    pub dudy: Float,
    pub dvdx: Float,
    pub dvdy: Float,
    pub primitive: Option<&'prim Primitive>,
    pub shading: Shading,
    pub shape: Option<&'prim Shape>,
    pub uv: Point2f,

    pub face_index: Option<usize>,
}

impl<'prim> SurfaceInteraction<'prim> {
    pub fn new(
        point: Point3f,
        point_error: Vector3f,
        uv: Point2f,
        wo: Vector3f,
        dpdu: Vector3f,
        dpdv: Vector3f,
        dndu: Normal3f,
        dndv: Normal3f,
        time: Float,
        shape: &'prim Shape,
        face_index: Option<usize>,
    ) -> Self {
        Self {
            point,
            point_error,
            uv,
            wo,
            dpdu,
            dpdv,
            dndu,
            dndv,
            time,
            shape: Some(shape),
            face_index,
            .. Self::default()
        }
    }

    pub fn set_shading_geometry(
        &mut self,
        dpdu: Vector3f,
        dpdv: Vector3f,
        dndu: Normal3f,
        dndv: Normal3f,
        orientation_is_authoritative: bool,
    ) {
        self.shading.normal = dpdu.cross(&dpdv).normalized();
        if let Some(shape) = self.shape {
            if shape.reverse_orientation() ^ shape.transform_swaps_handedness() {
                self.shading.normal = -self.shading.normal;
            }
        }

        if orientation_is_authoritative {
            self.normal = face_forward(self.normal, self.shading.normal);
        } else {
            self.shading.normal = face_forward(self.shading.normal, self.normal)
        }

        self.shading.dpdu = dpdu;
        self.shading.dpdv = dpdv;
        self.shading.dndu = dndu;
        self.shading.dndv = dndv;
    }
}

impl<'prim> Default for SurfaceInteraction<'prim> {
    fn default() -> Self {
        Self {
            ..Default::default()
        }
        //        Self {
        //            point: Point3f::default(),
        //            point_error: Vector3f::zero(),
        //            time: 0.0,
        //            wo: Vector3f::zero(),
        //            normal: Normal3f::zero(),
        //            medium_interface: None,
        //
        //            bsdf: None,
        //            bssrdf: None,
        //            dndu: Vector3f::zero(),
        //            dndv: Vector3f::zero(),
        //            dpdu: Vector3f::zero(),
        //            dpdv: Vector3f::zero(),
        //            dpdx: Vector3f::zero(),
        //            dpdy: Vector3f::zero(),
        //            dudx: 0.0,
        //            dudy: 0.0,
        //            dvdx: 0.0,
        //            dvdy: 0.0,
        //            primitive: None,
        //            shading: Shading {
        //                dpdu: Vector3f::zero(),
        //                dpdv: Vector3f::zero(),
        //                dndu: Normal3f::zero(),
        //                dndv: Normal3f::zero(),
        //                normal: Normal3f::zero(),
        //            },
        //            shape: None,
        //            uv: Point2f::zero(),
        //            face_index: None,
        //        }
    }
}

impl<'prim> Interaction<'prim> for SurfaceInteraction<'prim> {
    fn get_point(&self) -> Point3f {
        self.point
    }
    fn get_point_error(&self) -> Vector3f {
        self.point_error
    }
    fn get_time(&self) -> Float {
        self.time
    }
    fn get_normal(&self) -> Normal3f {
        self.normal
    }
    fn get_medium_interface(&self) -> Option<&'prim MediumInterface> {
        self.medium_interface
    }
    fn set_medium_interface(&mut self, interface: Option<&'prim MediumInterface>) {
        self.medium_interface = interface;
    }
}

pub struct MediumInteraction<'prim> {
    pub point: Point3f,
    pub point_error: Vector3f,
    pub time: Float,
    pub wo: Vector3f,
    pub normal: Normal3f,
    pub medium_interface: Option<&'prim MediumInterface>,
}

impl<'prim> Interaction<'prim> for MediumInteraction<'prim> {
    fn get_point(&self) -> Point3f {
        self.point
    }
    fn get_point_error(&self) -> Vector3f {
        self.point_error
    }
    fn get_time(&self) -> Float {
        self.time
    }
    fn get_normal(&self) -> Normal3f {
        self.normal
    }
    fn get_medium_interface(&self) -> Option<&'prim MediumInterface> {
        self.medium_interface
    }
    fn set_medium_interface(&mut self, interface: Option<&'prim MediumInterface>) {
        self.medium_interface = interface;
    }
}

//pub enum Interaction<'prim: 'ray> {
//    Surface(SurfaceInteraction<'prim>),
//    Medium(MediumInteraction<'prim>),
//}

fn accept_medium(mi: MediumInteraction) -> bool {
    true
}
