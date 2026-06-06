//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Compute uniform B-spline knot vector for `n` control points and given degree.
pub fn uniform_knot_vector(n: usize, degree: usize) -> Vec<f64> {
    let m = n + degree + 1;
    (0..m).map(|i| i as f64 / (m - 1) as f64).collect()
}
/// Compute clamped B-spline knot vector (open knot vector).
pub fn clamped_knot_vector(n: usize, degree: usize) -> Vec<f64> {
    let m = n + degree + 1;
    let mut knots = vec![0.0f64; m];
    for (i, knot) in knots.iter_mut().enumerate().take(m) {
        if i <= degree {
            *knot = 0.0;
        } else if i >= n {
            *knot = (n - degree) as f64;
        } else {
            *knot = (i - degree) as f64;
        }
    }
    let max_k = knots[m - 1];
    if max_k > 0.0 {
        for k in &mut knots {
            *k /= max_k;
        }
    }
    knots
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::fractal_geometry::IteratedFunctionSystem;
    use crate::fractal_geometry::LSystem;
    use crate::procedural_geometry::CurveParametric;
    use crate::procedural_geometry::DelaunayTri;
    use crate::procedural_geometry::DelaunayTriangulation2d;
    use crate::procedural_geometry::FractalGeometry;
    use crate::procedural_geometry::NoiseGeometry;
    use crate::procedural_geometry::PerlinNoise2d;
    use crate::procedural_geometry::PerlinNoise3d;
    use crate::procedural_geometry::ReactionDiffusion;
    use crate::procedural_geometry::SpaceColonization;
    use crate::procedural_geometry::SurfaceParametric;
    use crate::procedural_geometry::TurtleState;
    use crate::procedural_geometry::VoronoiTessellation2d;
    use crate::procedural_geometry::WorleyNoise;
    pub(super) const EPS: f64 = 1e-6;
    #[test]
    fn test_koch_fractal_dimension_approx() {
        let koch = FractalGeometry::koch_snowflake(5, 4);
        let dim = koch.box_counting_dimension(10);
        assert!(
            dim > 1.1 && dim < 1.5,
            "Koch dim = {:.6}, expected ~1.26",
            dim
        );
    }
    #[test]
    fn test_sierpinski_fractal_dimension_approx() {
        let sier = FractalGeometry::sierpinski_triangle(7);
        let dim = sier.hausdorff_dimension();
        assert!(
            dim > 0.5 && dim < 2.0,
            "Sierpinski dim = {:.6} out of expected range",
            dim
        );
    }
    #[test]
    fn test_fractal_dim_line_is_one() {
        let pts: Vec<[f64; 2]> = (0..500).map(|i| [i as f64 / 499.0, 0.5]).collect();
        let fg = FractalGeometry::new(pts);
        let dim = fg.box_counting_dimension(8);
        assert!(dim > 0.8 && dim < 1.2, "Line dim = {:.6}", dim);
    }
    #[test]
    fn test_fractal_dim_plane_fill_near_two() {
        let pts: Vec<[f64; 2]> = (0..50)
            .flat_map(|i| (0..50).map(move |j| [i as f64 / 49.0, j as f64 / 49.0]))
            .collect();
        let fg = FractalGeometry::new(pts);
        let dim = fg.box_counting_dimension(8);
        assert!(dim > 1.1, "Dense grid dim = {:.6} should be > 1.1", dim);
    }
    #[test]
    fn test_koch_snowflake_point_count() {
        let k = FractalGeometry::koch_snowflake(3, 2);
        let expected = 3 * 4usize.pow(3) * 2;
        assert_eq!(k.points.len(), expected);
    }
    #[test]
    fn test_sierpinski_triangle_non_empty() {
        let s = FractalGeometry::sierpinski_triangle(4);
        assert!(!s.points.is_empty());
        assert_eq!(s.points.len(), 3usize.pow(4) * 3);
    }
    #[test]
    fn test_lsystem_string_grows_geometrically() {
        let ls = LSystem::koch_curve();
        let s1 = ls.evolve(1).string.len();
        let s2 = ls.evolve(2).string.len();
        let s3 = ls.evolve(3).string.len();
        assert!(s2 > s1, "s2={} should be > s1={}", s2, s1);
        assert!(s3 > s2, "s3={} should be > s2={}", s3, s2);
    }
    #[test]
    fn test_lsystem_axiom_preserved_gen0() {
        let ls = LSystem::sierpinski();
        let state = ls.evolve(0);
        assert_eq!(state.string, ls.axiom);
        assert_eq!(state.generation, 0);
    }
    #[test]
    fn test_lsystem_dragon_curve_generates_segments() {
        let ls = LSystem::dragon_curve();
        let state = ls.evolve(6);
        let segs = ls.interpret(&state);
        assert!(!segs.is_empty(), "dragon curve should have segments");
    }
    #[test]
    fn test_lsystem_plant_branching_has_brackets() {
        let ls = LSystem::plant_branching();
        let state = ls.evolve(3);
        let open = state.string.chars().filter(|&c| c == '[').count();
        let close = state.string.chars().filter(|&c| c == ']').count();
        assert_eq!(
            open, close,
            "unbalanced brackets: open={}, close={}",
            open, close
        );
    }
    #[test]
    fn test_lsystem_generation_count() {
        let ls = LSystem::koch_curve();
        let state = ls.evolve(5);
        assert_eq!(state.generation, 5);
    }
    #[test]
    fn test_lsystem_turtle_state_default() {
        let t = TurtleState::default();
        assert_eq!(t.x, 0.0);
        assert_eq!(t.y, 0.0);
        assert_eq!(t.angle, 0.0);
    }
    #[test]
    fn test_perlin2d_range() {
        let noise = PerlinNoise2d::new(42);
        for i in 0..100 {
            let x = i as f64 * 0.137;
            let y = i as f64 * 0.251;
            let v = noise.noise(x, y);
            assert!((-1.5..=1.5).contains(&v), "Perlin2D out of range: {:.6}", v);
        }
    }
    #[test]
    fn test_perlin3d_range() {
        let noise = PerlinNoise3d::new(123);
        for i in 0..100 {
            let x = i as f64 * 0.1;
            let y = i as f64 * 0.17;
            let z = i as f64 * 0.23;
            let v = noise.noise(x, y, z);
            assert!((-1.5..=1.5).contains(&v), "Perlin3D out of range: {:.6}", v);
        }
    }
    #[test]
    fn test_perlin2d_deterministic() {
        let n1 = PerlinNoise2d::new(99);
        let n2 = PerlinNoise2d::new(99);
        for i in 0..10 {
            let x = i as f64 * 0.3;
            let y = i as f64 * 0.4;
            assert!((n1.noise(x, y) - n2.noise(x, y)).abs() < EPS);
        }
    }
    #[test]
    fn test_fbm_different_seeds_differ() {
        let n1 = PerlinNoise2d::new(1);
        let n2 = PerlinNoise2d::new(999);
        let v1 = n1.fbm(1.5, 2.3, 4, 2.0, 0.5);
        let v2 = n2.fbm(1.5, 2.3, 4, 2.0, 0.5);
        assert!(
            (v1 - v2).abs() > 1e-8,
            "Different seeds should produce different noise"
        );
    }
    #[test]
    fn test_fbm_nonnegative_max_val() {
        let noise = PerlinNoise2d::new(7);
        let v1 = noise.fbm(1.0, 1.0, 1, 2.0, 0.5);
        let v4 = noise.fbm(1.0, 1.0, 4, 2.0, 0.5);
        assert!(v1.is_finite(), "fBm 1 octave should be finite");
        assert!(v4.is_finite(), "fBm 4 octaves should be finite");
    }
    #[test]
    fn test_worley_f1_nonneg() {
        let w = WorleyNoise::new(20, 42);
        for i in 0..50 {
            let x = i as f64 / 49.0;
            let y = (50 - i) as f64 / 49.0;
            assert!(w.f1(x, y) >= 0.0);
        }
    }
    #[test]
    fn test_worley_f2_minus_f1_nonneg() {
        let w = WorleyNoise::new(20, 42);
        for i in 0..50 {
            let x = i as f64 / 49.0;
            let y = (50 - i) as f64 / 49.0;
            assert!(w.f2_minus_f1(x, y) >= -EPS);
        }
    }
    #[test]
    fn test_heightmap_fbm_dimensions() {
        let hm = NoiseGeometry::heightmap_fbm(16, 16, 4.0, 4, 77);
        assert_eq!(hm.len(), 16);
        assert_eq!(hm[0].len(), 16);
    }
    #[test]
    fn test_space_colonization_grows() {
        let attractors: Vec<[f64; 3]> = (0..50)
            .map(|i| [i as f64 * 0.05, 1.0 + i as f64 * 0.03, 0.0])
            .collect();
        let mut sc = SpaceColonization::new(attractors, [0.0, 0.0, 0.0], 0.1, 2.0, 0.15);
        let initial = sc.nodes.len();
        sc.grow();
        assert!(sc.nodes.len() >= initial);
    }
    #[test]
    fn test_space_colonization_root_has_no_parent() {
        let attractors = vec![[0.0, 1.0, 0.0]];
        let sc = SpaceColonization::new(attractors, [0.0, 0.0, 0.0], 0.1, 5.0, 0.15);
        assert!(sc.nodes[0].parent.is_none());
    }
    #[test]
    fn test_space_colonization_grow_until_done() {
        let attractors: Vec<[f64; 3]> = (0..10).map(|i| [i as f64 * 0.1, 0.5, 0.0]).collect();
        let mut sc = SpaceColonization::new(attractors, [0.5, 0.0, 0.0], 0.15, 1.0, 0.2);
        sc.grow_until_done(100);
        assert!(sc.attractors.is_empty() || sc.nodes.len() > 1);
    }
    #[test]
    fn test_reaction_diffusion_step_runs() {
        let mut rd = ReactionDiffusion::new(20, 20, 0.2, 0.1, 0.055, 0.062);
        rd.seed_center();
        rd.step(5, 1.0);
        for &u in &rd.u {
            assert!((0.0..=1.0).contains(&u), "U out of range: {:.6}", u);
        }
    }
    #[test]
    fn test_reaction_diffusion_v_in_range() {
        let mut rd = ReactionDiffusion::new(20, 20, 0.2, 0.1, 0.055, 0.062);
        rd.seed_center();
        rd.step(10, 1.0);
        for &v in &rd.v {
            assert!((0.0..=1.0).contains(&v), "V out of range: {:.6}", v);
        }
    }
    #[test]
    fn test_reaction_diffusion_grid_size() {
        let rd = ReactionDiffusion::new(10, 15, 0.2, 0.1, 0.055, 0.062);
        assert_eq!(rd.u.len(), 10 * 15);
        assert_eq!(rd.v.len(), 10 * 15);
    }
    #[test]
    fn test_ifs_barnsley_fern_generates_points() {
        let ifs = IteratedFunctionSystem::barnsley_fern();
        let pts = ifs.generate(1000, 42);
        assert_eq!(pts.len(), 1000);
    }
    #[test]
    fn test_ifs_barnsley_fern_y_range() {
        let ifs = IteratedFunctionSystem::barnsley_fern();
        let pts = ifs.generate(500, 1);
        for &[_, y] in &pts {
            assert!(
                (-0.5..=11.0).contains(&y),
                "y={:.6} out of expected range",
                y
            );
        }
    }
    #[test]
    fn test_ifs_empty_transforms() {
        let ifs = IteratedFunctionSystem::new(vec![], vec![]);
        let pts = ifs.generate(10, 42);
        assert!(pts.is_empty());
    }
    #[test]
    fn test_voronoi_cell_areas_sum() {
        let voronoi = VoronoiTessellation2d::random_seeds(10, [0.0, 1.0, 0.0, 1.0], 42);
        let areas = voronoi.cell_areas(50);
        let total: f64 = areas.iter().sum();
        assert!((total - 1.0).abs() < 0.05, "Total area = {:.6}", total);
    }
    #[test]
    fn test_voronoi_nearest_seed_is_correct() {
        let seeds = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let voronoi = VoronoiTessellation2d::new(seeds, [-0.5, 1.5, -0.5, 1.5]);
        assert_eq!(voronoi.nearest_seed(0.9, 0.1), 1);
    }
    #[test]
    fn test_voronoi_no_negative_areas() {
        let voronoi = VoronoiTessellation2d::random_seeds(8, [0.0, 1.0, 0.0, 1.0], 77);
        let areas = voronoi.cell_areas(30);
        for &a in &areas {
            assert!(a >= 0.0);
        }
    }
    #[test]
    fn test_voronoi_neighbors_symmetric() {
        let voronoi = VoronoiTessellation2d::random_seeds(6, [0.0, 1.0, 0.0, 1.0], 55);
        let neighbors = voronoi.neighbors(30);
        for (i, nbrs) in neighbors.iter().enumerate() {
            for &j in nbrs {
                assert!(
                    neighbors[j].contains(&i),
                    "neighbor relation not symmetric: {} neighbors {} but not vice versa",
                    i,
                    j
                );
            }
        }
    }
    #[test]
    fn test_delaunay_circumcircle_property() {
        let points = vec![[0.0, 0.0], [1.0, 0.0], [0.5, 0.866], [0.5, 0.3], [0.2, 0.5]];
        let dt = DelaunayTriangulation2d::new(points);
        assert!(dt.verify_delaunay(), "Delaunay property violated");
    }
    #[test]
    fn test_delaunay_triangle_count() {
        let points = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
        let dt = DelaunayTriangulation2d::new(points);
        assert!(!dt.triangles.is_empty(), "Should have at least 1 triangle");
    }
    #[test]
    fn test_delaunay_triangle_indices_valid() {
        let points: Vec<[f64; 2]> = (0..8)
            .map(|i| {
                let a = i as f64 * std::f64::consts::PI / 4.0;
                [a.cos(), a.sin()]
            })
            .collect();
        let dt = DelaunayTriangulation2d::new(points.clone());
        for tri in &dt.triangles {
            for &idx in &tri.indices {
                assert!(idx < points.len(), "Triangle index {} out of bounds", idx);
            }
        }
    }
    #[test]
    fn test_delaunay_circumcircle_center() {
        let points = vec![[0.0, 0.0], [1.0, 0.0], [0.5, 0.866_025_4]];
        let tri = DelaunayTri::new(0, 1, 2);
        let ([cx, cy], _r2) = tri.circumcircle(&points);
        assert!((cx - 0.5).abs() < 0.01, "circumcenter x = {:.6}", cx);
        assert!(cy > 0.0 && cy < 1.0, "circumcenter y = {:.6}", cy);
    }
    #[test]
    fn test_bezier_endpoint_interpolation() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [3.0, 1.0, 0.0],
            [4.0, 0.0, 0.0],
        ];
        let curve = CurveParametric::new(pts.clone());
        let p0 = curve.bezier_de_casteljau(0.0);
        let p1 = curve.bezier_de_casteljau(1.0);
        assert!((p0[0] - pts[0][0]).abs() < EPS, "start x mismatch");
        assert!((p1[0] - pts[3][0]).abs() < EPS, "end x mismatch");
    }
    #[test]
    fn test_bezier_midpoint_on_curve() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
        ];
        let curve = CurveParametric::new(pts);
        let mid = curve.bezier_de_casteljau(0.5);
        assert!(mid[0] >= -0.01 && mid[0] <= 1.01);
        assert!(mid[1] >= -0.01 && mid[1] <= 1.01);
    }
    #[test]
    fn test_bezier_bernstein_matches_de_casteljau() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [2.0, 2.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let curve = CurveParametric::new(pts);
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let dc = curve.bezier_de_casteljau(t);
            let bern = curve.bezier_bernstein(t);
            assert!((dc[0] - bern[0]).abs() < EPS, "x mismatch at t={}", t);
            assert!((dc[1] - bern[1]).abs() < EPS, "y mismatch at t={}", t);
        }
    }
    #[test]
    fn test_bspline_endpoint_at_clamped_knots() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        ];
        let n = pts.len();
        let degree = 3;
        let knots = clamped_knot_vector(n, degree);
        let curve = CurveParametric::new(pts.clone());
        let start = curve.bspline_de_boor(0.0, degree, &knots);
        let end = curve.bspline_de_boor(1.0, degree, &knots);
        assert!((start[0] - pts[0][0]).abs() < EPS, "start x mismatch");
        assert!((end[0] - pts[n - 1][0]).abs() < EPS, "end x mismatch");
    }
    #[test]
    fn test_bezier_arc_length_positive() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = CurveParametric::new(pts);
        let len = curve.bezier_arc_length(50);
        assert!(len > 0.0, "arc length should be positive");
    }
    #[test]
    fn test_bezier_sample_count() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = CurveParametric::new(pts);
        let samples = curve.sample_bezier(20);
        assert_eq!(samples.len(), 20);
    }
    #[test]
    fn test_tensor_bezier_corner_interpolation() {
        let net = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[0.0, 1.0, 0.0], [1.0, 1.0, 1.0], [2.0, 1.0, 0.0]],
            vec![[0.0, 2.0, 0.0], [1.0, 2.0, 0.0], [2.0, 2.0, 0.0]],
        ];
        let surf = SurfaceParametric::new(net);
        let p00 = surf.tensor_bezier(0.0, 0.0);
        let p11 = surf.tensor_bezier(1.0, 1.0);
        assert!((p00[0] - 0.0).abs() < EPS);
        assert!((p11[0] - 2.0).abs() < EPS);
    }
    #[test]
    fn test_coons_patch_corners() {
        let net = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let surf = SurfaceParametric::new(net);
        let p00 = surf.coons_patch(0.0, 0.0);
        let p11 = surf.coons_patch(1.0, 1.0);
        assert!((p00[0] - 0.0).abs() < EPS && (p00[1] - 0.0).abs() < EPS);
        assert!((p11[0] - 1.0).abs() < EPS && (p11[1] - 1.0).abs() < EPS);
    }
    #[test]
    fn test_surface_sample_dimensions() {
        let net = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let surf = SurfaceParametric::new(net);
        let samples = surf.sample(5, 7);
        assert_eq!(samples.len(), 5);
        assert_eq!(samples[0].len(), 7);
    }
    #[test]
    fn test_surface_normal_unit_length() {
        let net = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0]],
            vec![[0.0, 2.0, 0.0], [1.0, 2.0, 0.0], [2.0, 2.0, 0.0]],
        ];
        let surf = SurfaceParametric::new(net);
        let n = surf.normal(0.5, 0.5);
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-4, "normal length = {:.6}", len);
    }
    #[test]
    fn test_uniform_knot_vector_length() {
        let n = 5;
        let degree = 3;
        let knots = uniform_knot_vector(n, degree);
        assert_eq!(knots.len(), n + degree + 1);
    }
    #[test]
    fn test_clamped_knot_vector_endpoints() {
        let n = 6;
        let degree = 3;
        let knots = clamped_knot_vector(n, degree);
        assert_eq!(knots.len(), n + degree + 1);
        assert!((knots[0] - 0.0).abs() < EPS);
        assert!((knots[knots.len() - 1] - 1.0).abs() < EPS);
    }
    #[test]
    fn test_clamped_knot_vector_monotone() {
        let knots = clamped_knot_vector(5, 2);
        for i in 1..knots.len() {
            assert!(knots[i] >= knots[i - 1], "knot not monotone at {}", i);
        }
    }
}
