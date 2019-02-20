mod math;
mod vector;
mod vector2;
mod vector3;

pub use self::vector::*;
pub use self::vector2::*;
pub use self::vector3::*;
use crate::types::Int;
use crate::types::Float;

pub type Point2i = Vector2<Int>;
pub type Point2f = Vector2<Float>;

//pub type Normal2<T> = Vector2<T>;
//pub type Point3<T> = Vector3<T>;
//pub type Normal3<T> = Vector3<T>;


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn bob() {
        let p = Point2i::newp();
    }

//    #[test]
//    fn vector2() {
//        let a = Point2::new(1.0, 2.0);
//        let b = Normal2::new(2.0, 4.0);
//        assert_eq!(a + b, Vector2::new(3.0, 6.0));
//        assert_eq!(a * b, Vector2::new(2.0, 8.0));
//        assert_eq!(a - b, Vector2::new(-1.0, -2.0));
//        assert_eq!(a / b, Vector2::new(0.5, 0.5));
//
//        assert_eq!(dot(a, b), 10.0);
//        assert_eq!(abs_dot(a, b), 10.0);
//
//
//    }
}