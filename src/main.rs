use std::{f64::consts::PI, fmt::Debug, fs::File, io::Write};

use cgmath::{InnerSpace, Point3, Vector2, Vector3};
use rand::Rng;
use roots::{find_roots_quadratic, Roots};

type DirectionVector = Vector3<f64>;
type WorldCoordinate = Point3<f64>;
type TextureCoordinate = Vector2<f64>;
type Rgb = Vector3<f64>;

trait Mix {
    fn mix(self, other: Self, ratio: f64) -> Self;
}

impl Mix for Rgb {
    fn mix(self, other: Self, ratio: f64) -> Self {
        self * (1. - ratio) + other * ratio
    }
}

trait Clamp {
    fn clamp(self, min: f64, max: f64) -> Self;
}

impl Clamp for Rgb {
    fn clamp(self, min: f64, max: f64) -> Self {
        Rgb::new(
            self[0].min(max).max(min),
            self[1].min(max).max(min),
            self[2].min(max).max(min),
        )
    }
}

#[derive(Debug)]
struct Ray {
    origin: WorldCoordinate,
    direction: DirectionVector,
}

trait Interactable {
    fn intersect(&self, ray: &Ray) -> Option<(f64, usize, TextureCoordinate)>;
    fn surface_data(
        &self,
        point_hit: WorldCoordinate,
        uv: TextureCoordinate,
        index: usize,
    ) -> (DirectionVector, TextureCoordinate);
    fn colour(&self, st: TextureCoordinate) -> Rgb;
    fn albedo(&self) -> f64 {
        0.18
    }
    fn name(&self) -> &str {
        "Interactable"
    }
}

impl Debug for dyn Interactable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Interactable")
            .field("name", &self.name())
            .finish()
    }
}

struct Sphere {
    centre: WorldCoordinate,
    radius: f64,
    base_colour: Rgb,
}

impl Interactable for Sphere {
    fn intersect(&self, ray: &Ray) -> Option<(f64, usize, TextureCoordinate)> {
        // Solve analytically.

        let hypotenuse: DirectionVector = ray.origin - self.centre;
        let a2 = ray.direction.magnitude2();
        let a1 = 2. * ray.direction.dot(hypotenuse);
        let a0 = hypotenuse.magnitude2() - self.radius.powi(2);

        match find_roots_quadratic(a2, a1, a0) {
            Roots::No(_) => None,
            Roots::One([root]) => {
                if root >= 0. {
                    Some((root, 0, TextureCoordinate::new(0., 0.)))
                } else {
                    None
                }
            }
            Roots::Two([t0, t1]) => {
                // t0 is guaranteed to be less than t1.
                // We want to return the closest point of intersection.
                if t0 >= 0. {
                    Some((t0, 0, TextureCoordinate::new(0., 0.)))
                } else if t1 >= 0. {
                    Some((t1, 0, TextureCoordinate::new(0., 0.)))
                } else {
                    None
                }
            }
            _ => unreachable!("Can't have more than two roots to a quadratic equation"),
        }
    }

    fn surface_data(
        &self,
        point_hit: WorldCoordinate,
        _uv: TextureCoordinate,
        _index: usize,
    ) -> (DirectionVector, TextureCoordinate) {
        let normal = (point_hit - self.centre).normalize();
        let tex = TextureCoordinate::new(
            0.5 * (1. + normal[2].atan2(normal[0])) / PI,
            normal[1].acos() / PI,
        );
        (normal, tex)
    }

    fn colour(&self, _st: TextureCoordinate) -> Rgb {
        self.base_colour
    }
}

// struct Plane {
//     normal: DirectionVector,
//     point: WorldCoordinate,
//     base_colour: Rgb,
// }

// impl Interactable for Plane {
//     fn intersect(&self, ray: &Ray) -> Option<f64> {
//         let denominator = self.normal.dot(ray.direction);
//         if denominator.abs() > 1e-6 {
//             let root = (self.point - ray.origin).dot(self.normal) / denominator;
//             if root >= 0. {
//                 return Some(root);
//             }
//         }

//         None
//     }

//     fn surface_data(
//         &self,
//         _point_hit: WorldCoordinate,
//         _uv: TextureCoordinate,
//         _index: usize,
//     ) -> (DirectionVector, TextureCoordinate) {
//         todo!()
//     }

//     fn colour(&self) -> &Rgb {
//         &self.base_colour
//     }
// }

struct TriangularMesh {
    vertices: Vec<WorldCoordinate>,
    triangles: usize,
    vertex_index: Vec<usize>,
    // normals: Vec<DirectionVector>,
    st_coordinates: Vec<TextureCoordinate>,
}

impl Interactable for TriangularMesh {
    fn intersect(&self, ray: &Ray) -> Option<(f64, usize, Vector2<f64>)> {
        let mut nearest_t = f64::MAX;
        let mut pair = None;
        for index in 0..self.triangles {
            let a = self.vertex_index[index * 3];
            let b = self.vertex_index[index * 3 + 1];
            let c = self.vertex_index[index * 3 + 2];

            let e1 = self.vertices[b] - self.vertices[a];
            let e2 = self.vertices[c] - self.vertices[a];
            let p = ray.direction.cross(e2);
            let det = e1.dot(p);

            if det < 1e-6 {
                continue;
            }

            let inv_det = 1. / det;
            let tvec = ray.origin - self.vertices[a];
            let u = tvec.dot(p) * inv_det;
            if u < 0. || u > 1. {
                continue;
            }

            let q = tvec.cross(e1);
            let v = ray.direction.dot(q) * inv_det;
            if v < 0. || u + v > 1. {
                continue;
            }

            let t = e2.dot(q) * inv_det;
            if t < nearest_t {
                nearest_t = t;
                pair = Some((t, index, TextureCoordinate::new(u, v)));
            }
        }
        pair
    }

    fn surface_data(
        &self,
        _point_hit: WorldCoordinate,
        uv: TextureCoordinate,
        index: usize,
    ) -> (DirectionVector, TextureCoordinate) {
        let a = self.vertex_index[index * 3];
        let b = self.vertex_index[index * 3 + 1];
        let c = self.vertex_index[index * 3 + 2];

        let e1 = (self.vertices[b] - self.vertices[a]).normalize();
        let e2 = (self.vertices[c] - self.vertices[b]).normalize();
        let normal = e1.cross(e2).normalize();

        let st0 = self.st_coordinates[self.vertex_index[index * 3]];
        let st1 = self.st_coordinates[self.vertex_index[index * 3 + 1]];
        let st2 = self.st_coordinates[self.vertex_index[index * 3 + 2]];

        (normal, st0 * (1. - uv.x - uv.y) + st1 * uv.x + st2 * uv.y)
    }

    fn colour(&self, st: TextureCoordinate) -> Rgb {
        let scale = 5.;

        let pattern = ((st.x * scale) % 1. > 0.5) ^ ((st.y * scale) % 1. > 0.5);
        // return mix(Vec3f(0.815, 0.235, 0.031), Vec3f(0.937, 0.937, 0.231), pattern);
        Rgb::new(0.815, 0.031, 0.235).mix(Rgb::new(0.937, 0.0, 0.937), pattern as u8 as f64)
    }

    fn name(&self) -> &str {
        "TriangularMesh"
    }
}

trait Light {
    fn intensity(&self, point: WorldCoordinate) -> f64;
    fn colour(&self) -> &Rgb;
    fn direction(&self, point: WorldCoordinate) -> DirectionVector;
}

struct DistantLight {
    intensity: f64,
    colour: Rgb,
    direction: DirectionVector,
}

impl Light for DistantLight {
    fn intensity(&self, _point: WorldCoordinate) -> f64 {
        self.intensity
    }

    fn colour(&self) -> &Rgb {
        &self.colour
    }

    fn direction(&self, _point: WorldCoordinate) -> DirectionVector {
        self.direction
    }
}

// #[derive(Default)]
struct Scene {
    objects: Vec<Box<dyn Interactable>>,
    light: Box<dyn Light>,
}

#[derive(Debug)]
struct RenderOptions {
    width: u32,
    height: u32,
    field_of_view: f64,
    shadow_bias: f64,
}

impl Default for RenderOptions {
    fn default() -> Self {
        RenderOptions {
            width: 800,
            height: 600,
            field_of_view: 90.,
            shadow_bias: 1e-4,
        }
    }
}

struct IntersectData<'a> {
    t: f64,
    object: &'a dyn Interactable,
    index: usize,
    uv: TextureCoordinate,
}

struct RayTraceRenderer {
    options: RenderOptions,
}

impl RayTraceRenderer {
    fn new(options: RenderOptions) -> Self {
        RayTraceRenderer { options }
    }

    fn trace<'a>(&'a self, ray: &Ray, scene: &'a Scene) -> Option<IntersectData> {
        let mut nearest_t = f64::MAX;
        let mut data = None;
        scene.objects.iter().for_each(|object| {
            if let Some((t, index, uv)) = object.intersect(ray) {
                if t < nearest_t {
                    nearest_t = t;
                    data = Some(IntersectData {
                        t,
                        object: &**object,
                        index,
                        uv,
                    });
                };
            };
        });
        data
    }

    fn cast_ray(&self, ray: Ray, scene: &Scene, _recursion_depth: u8) -> Rgb {
        match self.trace(&ray, scene) {
            Some(IntersectData {
                t,
                object,
                index,
                uv,
            }) => {
                let point_hit: WorldCoordinate = ray.origin + ray.direction * t;
                let (normal, st) = object.surface_data(point_hit, uv, index);

                let light_vector = -scene.light.direction(point_hit);
                // compute the color of a diffuse surface illuminated
                // by a single distant light source.

                let shadow_ray = Ray {
                    origin: point_hit + normal * self.options.shadow_bias,
                    direction: light_vector,
                };
                match self.trace(&shadow_ray, scene) {
                    Some(_) => Rgb::new(0., 0., 0.),
                    None => {
                        // object.albedo() / PI
                        //     * scene.light.intensity(point_hit)
                        //     * scene.light.colour()
                        //     * normal.dot(light_vector).max(0.)
                        object.colour(st)
                            * scene.light.intensity(point_hit)
                            * normal.dot(light_vector).max(0.)
                    }
                }

                // let scale = 4.;
                // let pattern = ((tex[0] * scale) % 1. > 0.5) ^ ((tex[1] * scale) % 1. > 0.5);

                // normal.dot(-1. * ray.direction)
                //     * object
                //         .colour()
                //         .mix(&(object.colour() * 0.8), pattern as u8 as f64)
            }
            // None => Rgb::new(0., 0.4, 0.8),
            None => Rgb::new(0.6, 1., 1.), // None => Rgb::new(0., 0., 0.),
        }
    }

    fn render(&self, scene: Scene) {
        println!("Rendering image with options: {:?}", self.options);

        let width_f64 = f64::from(self.options.width);
        let height_f64 = f64::from(self.options.height);

        let scale: f64 = (self.options.field_of_view * PI / 360.).tan();
        let aspect_ratio: f64 = width_f64 / height_f64;

        // Cast rays per pixel.
        let mut framebuffer: Vec<Rgb> = Vec::new();
        for j in 0..self.options.height {
            for i in 0..self.options.width {
                // Translate from pixel coordinates to world coordinate system.

                // Center the pixel
                let mut x = f64::from(i) + 0.5;
                let mut y = f64::from(j) + 0.5;

                // Translate into NDC space (Normalized Device Coordinates)
                // 0 <= x, y < 1
                x /= width_f64;
                y /= height_f64;

                // Centre around the origin. Flip y-axis, since pixels start at the _top_ left.
                x = 2. * x - 1.;
                y = 1. - 2. * y;

                // Scale for field of view. x-axis needs adjusting for aspect ratio too
                x *= aspect_ratio * scale;
                y *= scale;

                let ray_direction: DirectionVector = DirectionVector::new(x, y, -1.).normalize();
                let ray = Ray {
                    origin: WorldCoordinate::new(0., 0., 0.),
                    direction: ray_direction,
                };
                framebuffer.push(self.cast_ray(ray, &scene, 0));
            }
        }
        let mut f = File::create("out.ppm").expect("Unable to create file");
        f.write_all(
            format!("P6\n{} {}\n255\n", self.options.width, self.options.height).as_bytes(),
        )
        .expect("Unable to write data");
        framebuffer.into_iter().for_each(|rgb| {
            let pixelated = (rgb.clamp(0., 1.) * 255.).map(|channel| channel.round() as u8);
            f.write_all(&[pixelated.x, pixelated.y, pixelated.z])
                .expect("Unable to write data");
        });
    }
}

fn main() {
    let mut rng = rand::thread_rng();

    // set up stage
    let mut objects: Vec<Box<dyn Interactable>> = Vec::new();

    // for _ in 0..12 {
    //     // 32 spheres :)
    //     objects.push(Box::new(Sphere {
    //         centre: WorldCoordinate::new(
    //             5. - 10. * rng.gen::<f64>(),
    //             5. - 10. * rng.gen::<f64>(),
    //             -15. - 10. * rng.gen::<f64>(),
    //         ),
    //         radius: 1. + 0.5 * rng.gen::<f64>(),
    //         base_colour: Rgb::new(rng.gen(), rng.gen(), rng.gen()),
    //     }));
    // }
    objects.push(Box::new(Sphere {
        centre: WorldCoordinate::new(-0., -5., -20.),
        radius: 1.5,
        base_colour: Rgb::new(1., 1., 0.),
    }));
    objects.push(Box::new(Sphere {
        centre: WorldCoordinate::new(2.6, 0., -19.2),
        radius: 1.,
        base_colour: Rgb::new(1., 0., 0.),
    }));
    objects.push(Box::new(TriangularMesh {
        vertices: vec![
            WorldCoordinate::new(-10., -8., -10.),
            WorldCoordinate::new(10., -8., -10.),
            WorldCoordinate::new(10., -8., -30.),
            WorldCoordinate::new(-10., -8., -30.),
        ],
        vertex_index: vec![0, 1, 3, 1, 2, 3],
        triangles: 2,
        st_coordinates: vec![
            TextureCoordinate::new(0., 0.),
            TextureCoordinate::new(1., 0.),
            TextureCoordinate::new(1., 1.),
            TextureCoordinate::new(0., 1.),
        ],
    }));
    let stage = Scene {
        objects,
        light: Box::new(DistantLight {
            intensity: 1.,
            colour: Rgb::new(1., 0., 0.),
            direction: DirectionVector::new(-0.3, -1., -0.).normalize(),
        }),
    };
    // Vec3f verts[4] = {{-5,-3,-6}, {5,-3,-6}, {5,-3,-16}, {-5,-3,-16}};
    // uint32_t vertIndex[6] = {0, 1, 3, 1, 2, 3};
    // Vec2f st[4] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    // set up options
    // let options = RenderOptions::default();
    let options = RenderOptions {
        width: 400,
        height: 400,
        field_of_view: 51.52,
        shadow_bias: 1e-4,
    };

    // render
    let raytracer = RayTraceRenderer::new(options);
    raytracer.render(stage);
}
