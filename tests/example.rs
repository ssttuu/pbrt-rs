extern crate pbrt;

use pbrt::primitive::GeometricPrimitive;
use pbrt::primitive::Primitive;
use pbrt::spectrum::Spectrum;
use pbrt::texture::constant::ConstantTexture;
use pbrt::transform;
use pbrt::types::Float;
use pbrt::Vector3f;

#[test]
fn test_render() {
    let f: Float = 1.0;
    println!("{}", f);

    let mut primitives: Vec<Box<Primitive>> = Vec::new();

    let n = 8;

    for k in 0..n {
        for i in 0..3 {
            let (mut x, mut y, mut z) = (0, 0, 0);
            let color: Spectrum;
            match i {
                0 => {
                    x = k / n * 100;
                    color = Spectrum::new(1.0);
                }
            }

            xform = transform::translate(Vector3f {
                x: x as Float,
                y: y as Float,
                z: z as Float,
            });
            kd = ConstantTexture::new(color);
            sigma = ConstantTexture::new(0.0);
            geoPrim = GeometricPrimitive::new(sphere, MatteMaterial::new(kd, sigma, nil));
            xformPrim =
                TransformedPrimitive::new(geoPrim, AnimatedTransform::new(xform, xform, 0, 1));
            primitives.push(xformPrim);
        }
    }
}
