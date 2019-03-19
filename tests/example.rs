extern crate pbrt;

use pbrt::materials::matte::MatteMaterial;
use pbrt::primitive::GeometricPrimitive;
use pbrt::primitive::Primitive;
use pbrt::spectrum::Spectrum;
use pbrt::texture::constant::ConstantTexture;
use pbrt::transform;
use pbrt::types::Float;
use pbrt::Vector3f;
use pbrt::primitive::TransformedPrimitive;
use pbrt::transform::AnimatedTransform;
use pbrt::sphere::Sphere;

#[test]
fn test_render() {
    let f: Float = 1.0;
    println!("{}", f);

    let mut primitives: Vec<Box<Primitive>> = Vec::new();

    let n = 8;

    for k in 0..n {
        for i in 0..3 {
            let (mut x, mut y, mut z) = (0, 0, 0);
            let mut color = Spectrum::default();
            match i {
                0 => {
                    x = k / n * 100;
                    color = Spectrum { data: vec![1., 0., 0.] };
                },
                1 => {
                    y = k / n * 100;
                    color = Spectrum { data: vec![0., 1., 0.] };
                },
                2 => {
                    z = k / n * 100;
                    color = Spectrum { data: vec![0., 0., 1.] };
                },
                _ => {}
            }

            let sphere = Box::new(Sphere::new());

            let xform = transform::translate(&Vector3f {
                x: x as Float,
                y: y as Float,
                z: z as Float,
            });
            let kd = ConstantTexture::new(color);
            let sigma = ConstantTexture::new(0.0);
            let geo_prim = GeometricPrimitive::new(
                sphere,
                Box::new(MatteMaterial {
                    kd,
                    sigma,
                    bump_map: None,
                }),
            );
            let xform_prim =
                TransformedPrimitive::new(geo_prim, AnimatedTransform::new(xform, xform, 0., 1.));
            primitives.push(xform_prim);
        }
    }
}
