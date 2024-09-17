use std::{f64::consts::PI, fs::File, io::Write};

use ndarray::{arr1, Array1};

struct Scene {}

type FloatVector = Array1<f64>;

trait Normalize {
    fn normalize(self) -> Self;
}

impl Normalize for FloatVector {
    fn normalize(self) -> Self {
        let magnitude = self.dot(&self).sqrt();
        if magnitude > 0. {
            return self / magnitude;
        }
        self
    }
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

    fn cast_ray(&self, direction: FloatVector, recursion_depth: u8) -> FloatVector {
        0.5 * direction + 0.5
    }

    fn render(&self, scene: Scene) {
        println!("Rendering image with options: {:?}", self.options);

        let width_f64 = f64::from(self.options.width);
        let height_f64 = f64::from(self.options.height);

        let scale: f64 = (self.options.field_of_view * PI / 180.).tan();
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

                let ray_direction: FloatVector = arr1(&[x, y, -1.]).normalize();
                framebuffer.push(self.cast_ray(ray_direction, 0));
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
    let stage = Scene {};
    // set up options
    let options = RenderOptions::default();

    // render
    let raytracer = RayTraceRenderer::new(options);
    raytracer.render(stage);
}
