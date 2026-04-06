//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{FractalDimension, Point2};

/// Utility: compute the total length of a set of line segments.
pub fn total_segment_length(segments: &[(Point2, Point2)]) -> f64 {
    segments.iter().map(|(a, b)| a.distance(b)).sum()
}
/// Utility: collect all endpoints of segments into a point set.
pub fn segments_to_points(segments: &[(Point2, Point2)]) -> Vec<Point2> {
    let mut points = Vec::with_capacity(segments.len() * 2);
    for (a, b) in segments {
        points.push(*a);
        points.push(*b);
    }
    points
}
/// Utility: compute the bounding box of a set of segments.
pub fn segments_bounding_box(segments: &[(Point2, Point2)]) -> (f64, f64, f64, f64) {
    let points = segments_to_points(segments);
    FractalDimension::bounding_box(&points)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::fractal_geometry::AffineTransform2D;
    use crate::fractal_geometry::BarnsleyFern;
    use crate::fractal_geometry::DragonCurve;
    use crate::fractal_geometry::HilbertCurve;
    use crate::fractal_geometry::JuliaSet;
    use crate::fractal_geometry::KochCurve;
    use crate::fractal_geometry::LSystem;
    use crate::fractal_geometry::MandelbrotSet;
    use crate::fractal_geometry::OrbitType;
    use crate::fractal_geometry::SierpinskiArrowhead;
    use crate::fractal_geometry::SierpinskiTriangle;
    pub(super) const EPS: f64 = 1e-10;
    #[test]
    fn test_affine_identity() {
        let id = AffineTransform2D::identity();
        let p = Point2::new(3.0, 7.0);
        let q = id.apply(&p);
        assert!((q.x - p.x).abs() < EPS);
        assert!((q.y - p.y).abs() < EPS);
    }
    #[test]
    fn test_affine_translation() {
        let t = AffineTransform2D::translation(1.0, 2.0);
        let p = Point2::new(3.0, 4.0);
        let q = t.apply(&p);
        assert!((q.x - 4.0).abs() < EPS);
        assert!((q.y - 6.0).abs() < EPS);
    }
    #[test]
    fn test_affine_scaling() {
        let s = AffineTransform2D::scaling(2.0, 3.0);
        let p = Point2::new(1.0, 1.0);
        let q = s.apply(&p);
        assert!((q.x - 2.0).abs() < EPS);
        assert!((q.y - 3.0).abs() < EPS);
    }
    #[test]
    fn test_affine_rotation_90() {
        let r = AffineTransform2D::rotation(std::f64::consts::FRAC_PI_2);
        let p = Point2::new(1.0, 0.0);
        let q = r.apply(&p);
        assert!(q.x.abs() < EPS);
        assert!((q.y - 1.0).abs() < EPS);
    }
    #[test]
    fn test_affine_compose() {
        let t = AffineTransform2D::translation(1.0, 0.0);
        let s = AffineTransform2D::scaling(2.0, 2.0);
        let composed = t.compose(&s);
        let p = Point2::new(1.0, 0.0);
        let q = composed.apply(&p);
        assert!((q.x - 4.0).abs() < EPS);
        assert!((q.y - 0.0).abs() < EPS);
    }
    #[test]
    fn test_affine_determinant() {
        let t = AffineTransform2D::new(1.0, 2.0, 3.0, 4.0, 0.0, 0.0);
        assert!((t.determinant() - (-2.0)).abs() < EPS);
    }
    #[test]
    fn test_affine_inverse() {
        let t = AffineTransform2D::new(2.0, 1.0, 1.0, 3.0, 5.0, 7.0);
        let inv = t.inverse().unwrap();
        let p = Point2::new(3.0, 4.0);
        let q = t.apply(&p);
        let r = inv.apply(&q);
        assert!((r.x - p.x).abs() < 1e-9);
        assert!((r.y - p.y).abs() < 1e-9);
    }
    #[test]
    fn test_ifs_num_transforms() {
        let ifs = SierpinskiTriangle::ifs();
        assert_eq!(ifs.num_transforms(), 3);
    }
    #[test]
    fn test_sierpinski_chaos_game_produces_points() {
        let points = SierpinskiTriangle::generate(1000);
        assert_eq!(points.len(), 1000);
    }
    #[test]
    fn test_sierpinski_points_bounded() {
        let points = SierpinskiTriangle::generate(5000);
        for p in &points {
            assert!(p.x >= -0.1 && p.x <= 1.1, "x out of bounds: {}", p.x);
            assert!(p.y >= -0.1 && p.y <= 0.9, "y out of bounds: {}", p.y);
        }
    }
    #[test]
    fn test_sierpinski_theoretical_dimension() {
        let dim = SierpinskiTriangle::theoretical_dimension();
        assert!((dim - 1.585).abs() < 0.001);
    }
    #[test]
    fn test_barnsley_fern_produces_points() {
        let points = BarnsleyFern::generate(500);
        assert_eq!(points.len(), 500);
    }
    #[test]
    fn test_barnsley_fern_ifs_probabilities() {
        let ifs = BarnsleyFern::ifs();
        let sum: f64 = ifs.probabilities.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_deterministic_iteration_growth() {
        let ifs = SierpinskiTriangle::ifs();
        let initial = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.5, 0.866),
        ];
        let gen1 = ifs.deterministic_iterate(&initial, 1);
        assert_eq!(gen1.len(), 9);
        let gen2 = ifs.deterministic_iterate(&initial, 2);
        assert_eq!(gen2.len(), 27);
    }
    #[test]
    fn test_lsystem_koch_generation_0() {
        let ls = KochCurve::lsystem();
        let s = ls.iterate(0);
        assert_eq!(s, "F");
    }
    #[test]
    fn test_lsystem_koch_generation_1() {
        let ls = KochCurve::lsystem();
        let s = ls.iterate(1);
        assert_eq!(s, "F+F--F+F");
    }
    #[test]
    fn test_lsystem_koch_length_grows() {
        let ls = KochCurve::lsystem();
        let l0 = ls.iterate(0).len();
        let l1 = ls.iterate(1).len();
        let l2 = ls.iterate(2).len();
        assert!(l1 > l0);
        assert!(l2 > l1);
    }
    #[test]
    fn test_lsystem_iterate_length() {
        let ls = KochCurve::lsystem();
        for n in 0..5 {
            let actual = ls.iterate(n).len();
            let predicted = ls.iterate_length(n);
            assert_eq!(actual, predicted, "mismatch at generation {n}");
        }
    }
    #[test]
    fn test_koch_segment_count() {
        assert_eq!(KochCurve::segment_count(0), 1);
        assert_eq!(KochCurve::segment_count(1), 4);
        assert_eq!(KochCurve::segment_count(2), 16);
        assert_eq!(KochCurve::segment_count(3), 64);
    }
    #[test]
    fn test_koch_theoretical_dimension() {
        let dim = KochCurve::theoretical_dimension();
        assert!((dim - 1.2619).abs() < 0.001);
    }
    #[test]
    fn test_turtle_graphics_produces_segments() {
        let segments = KochCurve::generate(2, 1.0);
        assert_eq!(segments.len(), 16);
    }
    #[test]
    fn test_mandelbrot_origin_in_set() {
        let m = MandelbrotSet::default_set();
        assert!(m.is_in_set(0.0, 0.0));
    }
    #[test]
    fn test_mandelbrot_large_value_escapes() {
        let m = MandelbrotSet::default_set();
        assert!(!m.is_in_set(10.0, 10.0));
        let esc = m.escape_time(10.0, 10.0);
        assert!(esc.is_some());
        assert!(esc.unwrap() <= 1, "should escape within 1 iteration");
    }
    #[test]
    fn test_mandelbrot_minus_one_in_set() {
        let m = MandelbrotSet::default_set();
        assert!(m.is_in_set(-1.0, 0.0));
    }
    #[test]
    fn test_mandelbrot_boundary() {
        let m = MandelbrotSet::default_set();
        assert!(m.is_boundary(0.26, 0.0, 0.05));
    }
    #[test]
    fn test_mandelbrot_smooth_escape() {
        let m = MandelbrotSet::default_set();
        let s = m.smooth_escape_time(0.5, 0.5);
        assert!(s > 0.0);
        assert!(s < m.max_iter as f64);
    }
    #[test]
    fn test_mandelbrot_render_grid() {
        let m = MandelbrotSet::new(64, 2.0);
        let grid = m.render_grid(-2.0, 1.0, -1.5, 1.5, 10, 10);
        assert_eq!(grid.len(), 10);
        assert_eq!(grid[0].len(), 10);
    }
    #[test]
    fn test_julia_classic_center() {
        let j = JuliaSet::classic();
        let esc = j.escape_time(0.1, 0.1);
        assert!(esc.is_none(), "z0=(0.1,0.1) should be in the Julia set");
    }
    #[test]
    fn test_julia_escape_far_point() {
        let j = JuliaSet::classic();
        assert!(!j.is_in_set(10.0, 10.0));
    }
    #[test]
    fn test_julia_orbit_classification_escaping() {
        let j = JuliaSet::classic();
        let orbit = j.classify_orbit(10.0, 10.0);
        assert_eq!(orbit, OrbitType::Escaping);
    }
    #[test]
    fn test_julia_render_grid() {
        let j = JuliaSet::new(-0.4, 0.6, 64, 2.0);
        let grid = j.render_grid(-1.5, 1.5, -1.5, 1.5, 8, 8);
        assert_eq!(grid.len(), 8);
        assert_eq!(grid[0].len(), 8);
    }
    #[test]
    fn test_box_counting_dimension_sierpinski() {
        let points = SierpinskiTriangle::generate(50000);
        let (dim, _r_sq) = FractalDimension::box_counting(&points, 10);
        assert!(dim > 1.3, "Sierpinski dim too low: {dim}");
        assert!(dim < 1.9, "Sierpinski dim too high: {dim}");
    }
    #[test]
    fn test_box_counting_empty() {
        let (dim, r_sq) = FractalDimension::box_counting(&[], 10);
        assert_eq!(dim, 0.0);
        assert_eq!(r_sq, 0.0);
    }
    #[test]
    fn test_hausdorff_estimate() {
        let points = SierpinskiTriangle::generate(20000);
        let dim = FractalDimension::hausdorff_estimate(&points, 12);
        assert!(dim > 1.0, "Hausdorff estimate too low: {dim}");
        assert!(dim < 2.0, "Hausdorff estimate too high: {dim}");
    }
    #[test]
    fn test_dragon_curve_generates() {
        let segments = DragonCurve::generate(4, 1.0);
        assert!(!segments.is_empty());
        assert_eq!(segments.len(), 16);
    }
    #[test]
    fn test_hilbert_curve_generates() {
        let segments = HilbertCurve::generate(2, 1.0);
        assert!(!segments.is_empty());
    }
    #[test]
    fn test_sierpinski_arrowhead_dimension() {
        let dim = SierpinskiArrowhead::theoretical_dimension();
        assert!((dim - 1.585).abs() < 0.001);
    }
    #[test]
    fn test_total_segment_length() {
        let segments = vec![
            (Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)),
            (Point2::new(1.0, 0.0), Point2::new(1.0, 1.0)),
        ];
        let len = total_segment_length(&segments);
        assert!((len - 2.0).abs() < EPS);
    }
    #[test]
    fn test_point2_distance() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(3.0, 4.0);
        assert!((a.distance(&b) - 5.0).abs() < EPS);
    }
    #[test]
    fn test_lsystem_add_rule() {
        let mut ls = LSystem::new(vec!['A'], vec![], "A");
        ls.add_rule('A', "AB");
        let s = ls.iterate(1);
        assert_eq!(s, "AB");
    }
    #[test]
    fn test_affine_contractivity() {
        let s = AffineTransform2D::scaling(0.5, 0.5);
        let c = s.contractivity();
        assert!((c - (0.5_f64).hypot(0.5)).abs() < EPS);
    }
}
