use crate::interaction::SurfaceInteraction;

pub enum TransportMode {
    Radiance,
    Importance,
}

pub trait Material {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    );

    fn bump(&self, d: &Texture<Float>, si: &mut SurfaceInteraction) {
        let mut du = 0.5 * (abs(si.dudx) + abs(si.dudy));
        if du == 0.0 {
            // TODO: I hate this.
            du = 0.00005;
        }
        si.point = si.point + si.shading.dpdu * du;
        si.uv.x += du;
        si.normal = (cross(&si.shading.dpdu, &si.shading.dpdv) + si.dndu * du).normalized();
        let displace_u = d.evaluate(si);

        let mut dv = 0.5 * (abs(si.dvdx) + abs(si.dvdy));
        if dv == 0.0 {
            dv = 0.00005;
        }
        si.point = si.point + si.shading.dpdv * dv;
        si.uv.y += dv;
        si.normal = (cross(&si.shading.dpdu, &si.shading.dpdv) + si.dndv * dv);
        let displace_v = d.evaluate(si);

        let displace = d.evaluate(si);

        let dpdu = si.shading.dpdu
            + Vector3f::fill(displace_u - displace) / si.shading.normal * du
            + Vector3f::fill(displace) * si.shading.dndu;
        let dpdv = si.shading.dpdv
            + Vector3f::fill(displace_v - displace) / si.shading.normal * dv
            + Vector3f::fill(displace) * si.shading.dndv;

        si.set_shading_geometry(dpdu, dpdv, si.shading.dndu, si.shading.dndv, false);
    }
}

pub mod matte;

use crate::geometry::vector::cross;
use crate::texture::Texture;
use crate::types::Float;
use crate::Vector2f;
use crate::Vector3f;
pub use matte::*;
use num::abs;
