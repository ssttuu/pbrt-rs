use crate::types::Float;
use crate::Normal3f;
use crate::Vector3f;
use std::rc::Rc;

#[derive(Clone)]
pub struct BSDF {
    pub eta: Float,
    pub ns: Normal3f,
    pub ng: Normal3f,
    pub ss: Vector3f,
    pub ts: Vector3f,
    pub bxdf_count: i32,
    pub bxdfs: Vec<Rc<BxDF>>,
}

pub enum BxDFType {
    Reflection = 1 << 1,
    Transmission = 1 << 2,
    Diffuse = 1 << 3,
    Glossy = 1 << 4,
    Specular = 1 << 5,
    All = (1 << 6) - 1,
}

pub trait BxDF {
    fn get_type(&self) -> BxDFType;
}
