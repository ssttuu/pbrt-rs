pub struct Int(i64);
pub struct Float(f64);

pub trait ToFloat {
    fn to_float(self) -> Float;
}

impl ToFloat for Int {
    fn to_float(self) -> Float { self as Float }
}

impl ToFloat for Float {
    fn to_float(self) -> Float { self }
}
