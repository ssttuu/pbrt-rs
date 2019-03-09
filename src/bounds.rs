use crate::types::Float;
use crate::types::Int;
use crate::types::Number;
use crate::Point2;
use crate::Point3;

pub struct Bounds2<T: Number> {
    pub min: Point2<T>,
    pub max: Point2<T>,
}
pub struct Bounds3<T: Number> {
    pub min: Point3<T>,
    pub max: Point3<T>,
}

pub type Bounds2i = Bounds2<Int>;
pub type Bounds2f = Bounds2<Float>;
pub type Bounds3i = Bounds3<Int>;
pub type Bounds3f = Bounds3<Float>;
