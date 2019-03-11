use crate::interaction::MediumInteraction;
use crate::ray::Ray;
use crate::sampler::Sampler;
use crate::spectrum::Spectrum;
use crate::types::Float;
use crate::Point2f;
use crate::Vector3f;
use std::rc::Rc;

pub trait PhaseFunction {
    fn p(&self, wo: &mut Vector3f, wi: &mut Vector3f) -> Float;
    fn sample_p(&self, wo: &mut Vector3f, u: &mut Point2f) -> Float;
}

pub trait Medium {
    fn tr(&self, ray: &mut Ray, s: &Sampler) -> Spectrum;
    fn sample(&self, ray: &mut Ray, s: &Sampler, mi: &MediumInteraction) -> Spectrum;
}

pub struct MediumInterface {
    pub inside: Option<Rc<Medium>>,
    pub outside: Option<Rc<Medium>>,
}

impl MediumInterface {
    pub fn new() -> Self {
        Self {
            inside: None,
            outside: None,
        }
    }

    pub fn is_medium_transition(&self) -> bool {
        match (&self.inside, &self.outside) {
            (Some(m1), Some(m2)) => !Rc::ptr_eq(m1, m2),
            _ => false,
        }
    }
}

// TODO:
pub struct MediumImpl {
    pub v: Float,
}

impl PartialEq for MediumImpl {
    fn eq(&self, other: &MediumImpl) -> bool {
        self.v == other.v
    }
}

impl Medium for MediumImpl {
    fn tr(&self, ray: &mut Ray, s: &Sampler) -> Spectrum {
        Spectrum::new(1.0)
    }
    fn sample(&self, ray: &mut Ray, s: &Sampler, mi: &MediumInteraction) -> Spectrum {
        Spectrum::new(1.0)
    }
}
