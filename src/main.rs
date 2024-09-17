use std::{f64::consts::PI, fs::File, io::Write};

use ndarray::{arr1, Array1};
use roots::{find_roots_quadratic, Roots};

type DirectionVector = Array1<f64>;
type WorldCoordinate = Array1<f64>;
type TextureCoordinate = Array1<f64>;
type RGB = Array1<f64>;

trait Normalize {
    fn normalize(self) -> Self;
}

impl Normalize for DirectionVector {
    fn normalize(self) -> Self {
        let magnitude = self.dot(&self).sqrt();
        if magnitude > 0. {
            return self / magnitude;
        }
        self
    }
}

trait Mix {
    fn mix(&self, other: &Self, ratio: f64) -> Self;
}

impl Mix for RGB {
    fn mix(&self, other: &Self, ratio: f64) -> Self {
        return self * (1. - ratio) + other * ratio;
    }
}

struct Ray {
    origin: WorldCoordinate,
    direction: DirectionVector,
}

trait Interactable {
    fn intersect(&self, ray: &Ray) -> Option<f64>;
    fn surface_data(&self, point_hit: WorldCoordinate) -> (DirectionVector, TextureCoordinate);
    fn colour(&self) -> &RGB;
}

struct Sphere {
    centre: WorldCoordinate,
    radius: f64,
    base_colour: RGB,
}

impl Interactable for Sphere {
    fn intersect(&self, ray: &Ray) -> Option<f64> {
        // Solve analytically.

        let hypotenuse: DirectionVector = &ray.origin - &self.centre;
        let a2 = ray.direction.dot(&ray.direction);
        let a1 = 2. * ray.direction.dot(&hypotenuse);
        let a0 = hypotenuse.dot(&hypotenuse) - &self.radius.powi(2);

        return match find_roots_quadratic(a2, a1, a0) {
            Roots::No(_) => None,
            Roots::One([root]) => {
                if root >= 0. {
                    Some(root)
                } else {
                    None
                }
            }
            Roots::Two([t0, t1]) => {
                // t0 is guaranteed to be less than t1.
                // We want to return the closest point of intersection.
                if t0 >= 0. {
                    Some(t0)
                } else if t1 >= 0. {
                    Some(t1)
                } else {
                    None
                }
            }
            _ => unreachable!("Can't have more than two roots to a quadratic equation"),
        };
    }

    fn surface_data(&self, point_hit: WorldCoordinate) -> (DirectionVector, TextureCoordinate) {
        let normal = (point_hit - &self.centre).normalize();
        let tex = arr1(&[
            0.5 * (1. + &normal[2].atan2(normal[0])) / PI,
            &normal[1].acos() / PI,
        ]);
        (normal, tex)
    }

    fn colour(&self) -> &RGB {
        &self.base_colour
    }
}

#[derive(Default)]
struct Scene {
    objects: Vec<Box<dyn Interactable>>,
}

#[derive(Debug)]
struct RenderOptions {
    width: u32,
    height: u32,
    field_of_view: f64,
}

impl Default for RenderOptions {
    fn default() -> Self {
        RenderOptions {
            width: 800,
            height: 600,
            field_of_view: 90.,
        }
    }
}

struct RayTraceRenderer {
    options: RenderOptions,
}

impl RayTraceRenderer {
    fn new(options: RenderOptions) -> Self {
        RayTraceRenderer { options }
    }

    fn trace<'a>(&'a self, ray: &Ray, scene: &'a Scene) -> Option<(f64, &Box<dyn Interactable>)> {
        let mut nearest_t = f64::MAX;
        let mut pair = None;
        scene.objects.iter().for_each(|object| {
            match object.intersect(ray) {
                Some(t) => {
                    if t < nearest_t {
                        nearest_t = t;
                    };
                    pair = Some((t, object));
                }
                None => {}
            };
        });
        pair
    }

    fn cast_ray(&self, ray: Ray, scene: &Scene, recursion_depth: u8) -> RGB {
        match self.trace(&ray, scene) {
            Some((t, object)) => {
                let point_hit: WorldCoordinate = ray.origin + &ray.direction * t;
                let (normal, tex) = object.surface_data(point_hit);
                let scale = 4.;
                let pattern = ((tex[0] * scale) % 1. > 0.5) ^ ((tex[1] * scale) % 1. > 0.5);

                normal.dot(&(-1. * ray.direction))
                    * &object
                        .colour()
                        .mix(&(object.colour() * 0.8), pattern as u8 as f64)
            }
            None => arr1(&[0., 0., 0.]),
        }
    }

    fn render(&self, scene: Scene) {
        println!("Rendering image with options: {:?}", self.options);

        let width_f64 = f64::from(self.options.width);
        let height_f64 = f64::from(self.options.height);

        let scale: f64 = (self.options.field_of_view * PI / 360.).tan();
        let aspect_ratio: f64 = width_f64 / height_f64;

        // let mut rng = rand::thread_rng();

        // Cast rays per pixel.
        let mut framebuffer: Vec<Array1<f64>> = Vec::new();
        for j in 0..self.options.height {
            for i in 0..self.options.width {
                // Translate from pixel coordinates to world coordinate system.

                // Center the pixel
                let mut x = f64::from(i) + 0.5;
                let mut y = f64::from(j) + 0.5;

                // Translate into NDC space (Normalized Device Coordinates)
                // 0 <= x, y < 1
                x = x / width_f64;
                y = y / height_f64;

                // Centre around the origin. Flip y-axis, since pixels start at the _top_ left.
                x = 2. * x - 1.;
                y = 1. - 2. * y;

                // Scale for field of view. x-axis needs adjusting for aspect ratio too
                x *= aspect_ratio * scale;
                y *= scale;

                let ray_direction: DirectionVector = arr1(&[x, y, -1.]).normalize();
                let ray = Ray {
                    origin: arr1(&[0., 0., 0.]),
                    direction: ray_direction,
                };
                framebuffer.push(self.cast_ray(ray, &scene, 0));
                // framebuffer.push([rng.gen(), rng.gen(), rng.gen()]);
            }
        }
        let mut f = File::create("out.ppm").expect("Unable to create file");
        f.write(format!("P6\n{} {}\n255\n", self.options.width, self.options.height).as_bytes())
            .expect("Unable to write data");
        // for (uint32_t i = 0; i < options.height * options.width; ++i) {
        //     char r = (char)(255 * clamp(0, 1, framebuffer[i].x));
        //     char g = (char)(255 * clamp(0, 1, framebuffer[i].y));
        //     char b = (char)(255 * clamp(0, 1, framebuffer[i].z));
        //     ofs << r << g << b;
        // }
        framebuffer.into_iter().for_each(|rgb| {
            f.write_all(
                &(rgb.clamp(0., 1.) * 255.)
                    .iter()
                    .map(|channel| channel.round() as u8)
                    .collect::<Vec<_>>(),
            )
            .expect("Unable to write data");
        });
        // f.write_all(framebuffer.into_iter()..as_bytes())
        // .expect("Unable to write data");
    }
}

fn main() {
    // set up stage
    let mut objects: Vec<Box<dyn Interactable>> = Vec::new();
    objects.push(Box::new(Sphere {
        centre: arr1(&[0., 0., -5.]),
        radius: 1.,
        base_colour: arr1(&[1., 0., 0.]),
    }));
    let stage = Scene { objects };
    // let stage = Scene::default();
    // set up options
    let options = RenderOptions::default();

    // render
    let raytracer = RayTraceRenderer::new(options);
    raytracer.render(stage);
}
