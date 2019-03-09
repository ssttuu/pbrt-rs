use crate::interaction::SurfaceInteraction;
use crate::texture::Texture;

pub struct ConstantTexture<T: Copy> {
    value: T,
}

impl<T> ConstantTexture<T>
where
    T: Copy,
{
    pub fn new(v: T) -> Self {
        Self { value: v }
    }
}

impl<T> Texture<T> for ConstantTexture<T>
where
    T: Copy,
{
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        self.value
    }
}
