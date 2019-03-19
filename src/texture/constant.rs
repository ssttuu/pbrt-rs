use crate::interaction::SurfaceInteraction;
use crate::texture::Texture;

pub struct ConstantTexture<T: Clone> {
    value: T,
}

impl<T> ConstantTexture<T>
where
    T: Clone,
{
    pub fn new(v: T) -> Self {
        Self { value: v }
    }
}

impl<T> Texture<T> for ConstantTexture<T>
where
    T: Clone,
{
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        self.value.clone()
    }
}
