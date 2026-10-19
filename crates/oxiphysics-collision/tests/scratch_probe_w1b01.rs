//! Temporary probe (deleted before finishing): how does `Gjk::query` classify tiny gaps?

use oxiphysics_collision::Gjk;
use oxiphysics_collision::narrowphase::gjk::GjkResult;
use oxiphysics_core::Transform;
use oxiphysics_core::math::Vec3;
use oxiphysics_geometry::{BoxShape, Sphere};

fn at(x: f64, y: f64, z: f64) -> Transform {
    Transform::from_position(Vec3::new(x, y, z))
}

fn describe(result: &GjkResult) -> String {
    match result {
        GjkResult::Intersecting(_) => "Intersecting".to_string(),
        GjkResult::Separated {
            distance,
            point_a,
            point_b,
        } => format!("Separated d={distance:e} a={point_a:?} b={point_b:?}"),
    }
}

#[test]
fn probe_tiny_gaps() {
    let cube = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let sphere = Sphere::new(1.0);
    for gap in [
        1e-1, 1e-2, 1e-3, 1e-4, 3e-5, 1e-5, 1e-6, 1e-7, 1e-8, 1e-10, 0.0, -1e-6,
    ] {
        let r1 = Gjk::query(&cube, &at(0.0, 0.0, 0.0), &cube, &at(2.0 + gap, 0.0, 0.0));
        let r2 = Gjk::query(&cube, &at(0.0, 0.0, 0.0), &cube, &at(2.0 + gap, 0.5, 0.0));
        let r3 = Gjk::query(
            &sphere,
            &at(0.0, 0.0, 0.0),
            &sphere,
            &at(2.0 + gap, 0.0, 0.0),
        );
        println!("gap {gap:e}");
        println!("  box   {}", describe(&r1));
        println!("  boxY  {}", describe(&r2));
        println!("  sph   {}", describe(&r3));
    }
}
