//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use crate::Color;
    use crate::streamlines::functions::symmetric_eigenvalues_ascending;
    use crate::streamlines::*;
    use crate::tensor_visualization::eig3;
    use oxiphysics_core::math::Vec3;
    #[test]
    fn test_streamline_arc_length_uniform() {
        let pts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
        ];
        let sl = Streamline::new(pts, Color::white());
        assert!((sl.arc_length() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_streamline_is_empty_false() {
        let sl = Streamline::new(vec![Vec3::zeros()], Color::white());
        assert!(!sl.is_empty());
        assert_eq!(sl.len(), 1);
    }
    #[test]
    fn test_streamline_empty() {
        let sl = Streamline::new(vec![], Color::white());
        assert!(sl.is_empty());
    }
    #[test]
    fn test_euler_uniform_field() {
        let field = |_: Vec3| Vec3::new(1.0, 0.0, 0.0);
        let integ = StreamlineIntegrator::EulerForward;
        let pos = integ.step(Vec3::zeros(), 0.1, &field);
        assert!((pos.x - 0.1).abs() < 1e-12);
        assert!(pos.y.abs() < 1e-12);
    }
    #[test]
    fn test_rk4_uniform_field_exact() {
        let c = 2.0_f64;
        let dt = 0.3_f64;
        let field = move |_: Vec3| Vec3::new(c, 0.0, 0.0);
        let integ = StreamlineIntegrator::RungeKutta4;
        let pos = integ.step(Vec3::zeros(), dt, &field);
        assert!((pos.x - c * dt).abs() < 1e-12);
    }
    #[test]
    fn test_rk4_vs_euler_circular() {
        let field = |p: Vec3| Vec3::new(p.y, -p.x, 0.0);
        let seed = Vec3::new(1.0, 0.0, 0.0);
        let dt = 0.05;
        let steps = 100;
        let mut euler_pos = seed;
        let mut rk4_pos = seed;
        for _ in 0..steps {
            euler_pos = StreamlineIntegrator::EulerForward.step(euler_pos, dt, &field);
            rk4_pos = StreamlineIntegrator::RungeKutta4.step(rk4_pos, dt, &field);
        }
        let rk4_err = (rk4_pos.norm() - 1.0).abs();
        let euler_err = (euler_pos.norm() - 1.0).abs();
        assert!(
            rk4_err < euler_err,
            "RK4 err={rk4_err} should < Euler err={euler_err}"
        );
    }
    #[test]
    fn test_tracer_uniform_field() {
        let seeds = vec![[0.0_f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let tracer = StreamlineTracer {
            seed_points: seeds,
            max_steps: 10,
            step_size: 0.1,
            integrator: StreamlineIntegrator::EulerForward,
            stagnation_threshold: 1e-12,
            color: Color::white(),
        };
        let lines = tracer.trace(&|_pos| [1.0, 0.0, 0.0]);
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0].points.len(), 11);
    }
    #[test]
    fn test_tracer_stagnation_stop() {
        let tracer = StreamlineTracer::new(vec![[0.0, 0.0, 0.0]]);
        let lines = tracer.trace(&|_| [0.0_f64, 0.0, 0.0]);
        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0].points.len(), 1);
    }
    #[test]
    fn test_tracer_velocity_magnitudes_filled() {
        let tracer = StreamlineTracer::new(vec![[0.0, 0.0, 0.0]]);
        let lines = tracer.trace(&|_| [2.0_f64, 0.0, 0.0]);
        let sl = &lines[0];
        assert_eq!(sl.velocity_magnitudes.len(), sl.points.len());
        for &m in &sl.velocity_magnitudes {
            assert!((m - 2.0).abs() < 1e-10);
        }
    }
    #[test]
    fn test_pathline_uniform_time_dep() {
        let field = |_t: f64, _p: Vec3| Vec3::new(1.0, 0.0, 0.0);
        let pl = Pathline::trace(Vec3::zeros(), 0.0, 0.1, 5, &field, Color::blue());
        assert_eq!(pl.len(), 6);
        let last = pl.points.last().unwrap();
        assert!((last.x - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_pathline_time_stamps() {
        let field = |_t: f64, _p: Vec3| Vec3::new(1.0, 0.0, 0.0);
        let dt = 0.2;
        let steps = 4;
        let pl = Pathline::trace(Vec3::zeros(), 0.0, dt, steps, &field, Color::red());
        assert_eq!(pl.times.len(), steps + 1);
        assert!((pl.times[0] - 0.0).abs() < 1e-12);
        assert!((pl.times[steps] - dt * steps as f64).abs() < 1e-10);
    }
    #[test]
    fn test_pathline_stagnation() {
        let field = |_t: f64, _p: Vec3| Vec3::zeros();
        let pl = Pathline::trace(
            Vec3::new(1.0, 2.0, 3.0),
            0.0,
            0.1,
            100,
            &field,
            Color::white(),
        );
        assert_eq!(pl.len(), 1);
    }
    #[test]
    fn test_streakline_particle_count() {
        let field = |_p: Vec3| Vec3::new(1.0, 0.0, 0.0);
        let sl = Streakline::compute_steady(Vec3::zeros(), 2.0, 0.5, 0.1, &field, Color::cyan());
        assert_eq!(sl.len(), 5);
    }
    #[test]
    fn test_streakline_zero_interval() {
        let field = |_p: Vec3| Vec3::new(1.0, 0.0, 0.0);
        let sl = Streakline::compute_steady(Vec3::zeros(), 1.0, 0.0, 0.1, &field, Color::white());
        assert!(sl.is_empty());
    }
    #[test]
    fn test_streakline_particles_spread() {
        let field = |_p: Vec3| Vec3::new(1.0, 0.0, 0.0);
        let sl = Streakline::compute_steady(Vec3::zeros(), 1.0, 0.5, 0.05, &field, Color::white());
        assert!(sl.len() >= 2);
        let oldest_x = sl.particle_positions[0].x;
        let newest_x = sl.particle_positions[sl.len() - 1].x;
        assert!(
            oldest_x > newest_x,
            "oldest particle should be furthest: {oldest_x} > {newest_x}"
        );
    }
    #[test]
    fn test_q_criterion_pure_rotation() {
        let j = [[0.0_f64, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let q = q_criterion(&j);
        assert!(
            (q - 1.0).abs() < 1e-10,
            "Q for solid-body rotation should be 1.0, got {q}"
        );
    }
    #[test]
    fn test_q_criterion_pure_strain_negative() {
        let j = [[1.0_f64, 0.0, 0.0], [0.0, -0.5, 0.0], [0.0, 0.0, -0.5]];
        let q = q_criterion(&j);
        assert!(q < 0.0, "Q for pure strain should be negative, got {q}");
    }
    #[test]
    fn test_detect_vortex_cores_circular_field() {
        let field = |p: Vec3| Vec3::new(p.y, -p.x, 0.0);
        let samples: Vec<Vec3> = vec![
            Vec3::new(0.1, 0.1, 0.0),
            Vec3::new(-0.1, 0.2, 0.0),
            Vec3::new(0.3, -0.1, 0.0),
        ];
        let cores = detect_vortex_cores(&samples, &field, 1e-4, 0.5, 0.0);
        assert!(
            !cores.is_empty(),
            "should detect vortex cores in circular flow"
        );
    }
    #[test]
    fn test_detect_vortex_cores_uniform_field_none() {
        let field = |_: Vec3| Vec3::new(1.0, 0.0, 0.0);
        let samples: Vec<Vec3> = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 0.0)];
        let cores = detect_vortex_cores(&samples, &field, 1e-4, 0.1, -0.1);
        assert!(
            cores.is_empty(),
            "uniform field should yield no vortex cores"
        );
    }
    #[test]
    fn test_lic_output_dimensions() {
        let params = LicParams::unit_square(32, 24);
        let field = |_: [f64; 3]| [1.0_f64, 0.0, 0.0];
        let tex = compute_lic(&params, &field);
        assert_eq!(tex.width, 32);
        assert_eq!(tex.height, 24);
        assert_eq!(tex.pixels.len(), 32 * 24);
    }
    #[test]
    fn test_lic_pixels_in_range() {
        let params = LicParams::unit_square(16, 16);
        let field = |pos: [f64; 3]| [-pos[1], pos[0], 0.0];
        let tex = compute_lic(&params, &field);
        for &p in &tex.pixels {
            assert!((0.0..=1.0).contains(&p), "LIC pixel out of [0,1]: {p}");
        }
    }
    #[test]
    fn test_lic_uniform_field_non_uniform_output() {
        let params = LicParams::unit_square(16, 16);
        let field = |_: [f64; 3]| [1.0_f64, 0.0, 0.0];
        let tex = compute_lic(&params, &field);
        let min = tex.pixels.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = tex.pixels.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(
            max - min > 0.01,
            "LIC output should have spatial variation (max-min={:.3})",
            max - min
        );
    }
    #[test]
    fn test_lic_texture_get() {
        let params = LicParams::unit_square(4, 4);
        let field = |_: [f64; 3]| [1.0_f64, 0.0, 0.0];
        let tex = compute_lic(&params, &field);
        let col = 2;
        let row = 3;
        assert_eq!(tex.get(col, row), tex.pixels[row * tex.width + col]);
    }
    #[test]
    fn test_trace_streamline_uniform_field() {
        let field = |_pos: Vec3| Vec3::new(1.0, 0.0, 0.0);
        let sl = trace_streamline(&field, Vec3::zeros(), 0.1, 10);
        assert_eq!(sl.points.len(), 11);
        let last = sl.points.last().unwrap();
        assert!((last.x - 1.0).abs() < 1e-10);
        assert!(last.y.abs() < 1e-10);
        assert!(last.z.abs() < 1e-10);
    }
    #[test]
    fn test_trace_streamlines_grid_count() {
        let field = |_pos: Vec3| Vec3::new(0.0, 1.0, 0.0);
        let seeds = vec![
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
        ];
        let lines = trace_streamlines_grid(&field, &seeds, 0.1, 5);
        assert_eq!(lines.len(), 3);
    }
    #[test]
    fn test_trace_streamline_circular_field_curves() {
        let field = |pos: Vec3| Vec3::new(pos.y, -pos.x, 0.0);
        let seed = Vec3::new(1.0, 0.0, 0.0);
        let sl = trace_streamline(&field, seed, 0.1, 5);
        assert!(sl.points.len() >= 3);
        let p0 = sl.points[0];
        let p1 = sl.points[1];
        let p2 = sl.points[2];
        let d01 = (p1 - p0).normalize();
        let d12 = (p2 - p1).normalize();
        let cross = d01.cross(&d12);
        assert!(
            cross.norm() > 1e-4,
            "streamline should curve in circular flow"
        );
    }
    #[test]
    fn test_stream_surface_n_streamlines_matches_seed_count() {
        let field = |_: [f64; 3]| [1.0_f64, 0.0, 0.0];
        let seeds: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 1.0, 0.0]];
        let surf = StreamSurface::generate(&seeds, &field, 10, 0.1, Color::white());
        assert_eq!(
            surf.n_streamlines(),
            3,
            "should have 3 streamlines for 3 seeds"
        );
    }
    #[test]
    fn test_stream_surface_total_points_nonzero() {
        let field = |_: [f64; 3]| [1.0_f64, 0.0, 0.0];
        let seeds: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let surf = StreamSurface::generate(&seeds, &field, 20, 0.05, Color::cyan());
        assert!(surf.total_points() > 0, "surface should have points");
    }
    #[test]
    fn test_stream_surface_area_positive() {
        let field = |_: [f64; 3]| [1.0_f64, 0.0, 0.0];
        let seeds: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 1.0, 0.0]];
        let surf = StreamSurface::generate(&seeds, &field, 20, 0.1, Color::white());
        let area = surf.surface_area();
        assert!(area > 0.0, "stream surface area should be positive: {area}");
    }
    #[test]
    fn test_stream_ribbon_len_positive() {
        let field = |_: Vec3| Vec3::new(1.0, 0.0, 0.0);
        let ribbon = StreamRibbon::generate(
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
            0.1,
            &field,
            10,
            0.1,
            Color::white(),
        );
        assert!(!ribbon.is_empty(), "ribbon should have points");
        assert_eq!(ribbon.left_edge().len(), ribbon.len());
        assert_eq!(ribbon.right_edge().len(), ribbon.len());
    }
    #[test]
    fn test_stream_ribbon_edges_offset_from_centre() {
        let field = |_: Vec3| Vec3::new(1.0, 0.0, 0.0);
        let hw = 0.2;
        let ribbon = StreamRibbon::generate(
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
            hw,
            &field,
            5,
            0.1,
            Color::white(),
        );
        let left = ribbon.left_edge();
        let right = ribbon.right_edge();
        assert!(!left.is_empty() && !right.is_empty());
        let spread = (left[0] - right[0]).norm();
        assert!(
            spread > 0.0,
            "left and right edges should be separated: spread={spread}"
        );
    }
    #[test]
    fn test_enhanced_lic_output_dimensions() {
        let params = EnhancedLicParams::unit_square(16, 12, 2);
        let field = |_: [f64; 3]| [1.0_f64, 0.0, 0.0];
        let tex = compute_enhanced_lic(&params, &field);
        assert_eq!(tex.width, 16);
        assert_eq!(tex.height, 12);
        assert_eq!(tex.pixels.len(), 16 * 12);
    }
    #[test]
    fn test_enhanced_lic_pixels_in_range() {
        let params = EnhancedLicParams::unit_square(16, 16, 1);
        let field = |p: [f64; 3]| [-p[1], p[0], 0.0];
        let tex = compute_enhanced_lic(&params, &field);
        for &p in &tex.pixels {
            assert!((0.0..=1.0 + 1e-9).contains(&p), "pixel out of range: {p}");
        }
    }
    #[test]
    fn test_topology_stagnation_detected_near_origin() {
        let field = |p: [f64; 3]| [p[0], -p[1], 0.0];
        let points = detect_topology_points([-0.5, 0.5], [-0.5, 0.5], 0.1, 0.0, &field, 0.1);
        assert!(
            !points.is_empty(),
            "should detect topology point near origin"
        );
    }
    #[test]
    fn test_topology_kind_saddle_for_hyperbolic_flow() {
        let field = |p: [f64; 3]| [p[0], -p[1], 0.0];
        let points = detect_topology_points([-0.05, 0.05], [-0.05, 0.05], 0.05, 0.0, &field, 0.1);
        let has_saddle = points.iter().any(|p| p.kind == TopologyPointKind::Saddle);
        assert!(has_saddle, "hyperbolic flow should yield a saddle point");
    }
    #[test]
    fn test_ftle_output_dimensions() {
        let field = |_: [f64; 3]| [1.0_f64, 0.0, 0.0];
        let ftle = compute_ftle_2d([0.0, 1.0], [0.0, 1.0], 0.0, 8, 8, 0.5, 0.1, &field);
        assert_eq!(ftle.width, 8);
        assert_eq!(ftle.height, 8);
        assert_eq!(ftle.values.len(), 64);
    }
    #[test]
    fn test_ftle_uniform_field_near_zero() {
        let field = |_: [f64; 3]| [1.0_f64, 0.0, 0.0];
        let ftle = compute_ftle_2d([0.0, 1.0], [0.0, 1.0], 0.0, 8, 8, 0.5, 0.1, &field);
        for &v in &ftle.values {
            assert!(v < 1e-6, "uniform field FTLE should be ~0, got {v}");
        }
    }
    #[test]
    fn test_ftle_shear_flow_nonzero() {
        let field = |p: [f64; 3]| [p[1], 0.0_f64, 0.0];
        let ftle = compute_ftle_2d([0.0, 1.0], [0.0, 1.0], 0.0, 8, 8, 1.0, 0.1, &field);
        let max = ftle.max_value();
        assert!(
            max > 0.0,
            "shear flow should produce positive FTLE: max={max}"
        );
    }
    #[test]
    fn test_ftle_threshold_returns_correct_indices() {
        let field = |p: [f64; 3]| [p[1], 0.0_f64, 0.0];
        let ftle = compute_ftle_2d([0.0, 1.0], [0.0, 1.0], 0.0, 8, 8, 1.0, 0.1, &field);
        let thresh = ftle.max_value() * 0.5;
        let above = ftle.threshold(thresh);
        for (col, row) in above {
            assert!(
                ftle.get(col, row) >= thresh,
                "threshold cell has value below threshold"
            );
        }
    }
    #[test]
    fn test_seed_uniform_grid_count() {
        let seeds = generate_seeds(
            [0.0; 3],
            [1.0; 3],
            SeedStrategy::UniformGrid {
                nx: 3,
                ny: 2,
                nz: 4,
            },
        );
        assert_eq!(seeds.len(), 3 * 2 * 4);
    }
    #[test]
    fn test_seed_uniform_grid_single_cell() {
        let seeds = generate_seeds(
            [0.0; 3],
            [1.0; 3],
            SeedStrategy::UniformGrid {
                nx: 1,
                ny: 1,
                nz: 1,
            },
        );
        assert_eq!(seeds.len(), 1);
        assert!((seeds[0].x - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_seed_random_count() {
        let seeds = generate_seeds(
            [-1.0; 3],
            [1.0; 3],
            SeedStrategy::Random {
                count: 20,
                seed: 42,
            },
        );
        assert_eq!(seeds.len(), 20);
    }
    #[test]
    fn test_seed_random_within_bounds() {
        let min = [-2.0, 0.0, 1.0];
        let max = [2.0, 3.0, 5.0];
        let seeds = generate_seeds(min, max, SeedStrategy::Random { count: 50, seed: 7 });
        for s in &seeds {
            assert!(s.x >= min[0] && s.x <= max[0]);
            assert!(s.y >= min[1] && s.y <= max[1]);
            assert!(s.z >= min[2] && s.z <= max[2]);
        }
    }
    #[test]
    fn test_velocity_color_zero_is_blue() {
        let c = velocity_color(0.0);
        assert!(
            c.b > c.r && c.b > c.g,
            "zero velocity should be bluish: {:?}",
            c
        );
    }
    #[test]
    fn test_velocity_color_one_is_red() {
        let c = velocity_color(1.0);
        assert!(
            c.r > c.b && c.r > c.g,
            "max velocity should be reddish: {:?}",
            c
        );
    }
    #[test]
    fn test_velocity_color_clamped() {
        let c_neg = velocity_color(-1.0);
        let c_pos = velocity_color(2.0);
        let c0 = velocity_color(0.0);
        let c1 = velocity_color(1.0);
        assert_eq!(c_neg.r, c0.r);
        assert_eq!(c_neg.g, c0.g);
        assert_eq!(c_neg.b, c0.b);
        assert_eq!(c_pos.r, c1.r);
    }
    #[test]
    fn test_color_by_velocity_length_matches_points() {
        let pts = vec![
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
        ];
        let mags = vec![0.5, 1.5, 3.0];
        let sl = Streamline::with_magnitudes(pts, mags, Color::white());
        let colors = color_by_velocity(&sl, 0.0, 3.0);
        assert_eq!(colors.len(), sl.points.len());
    }
    #[test]
    fn test_color_by_velocity_empty_mags_returns_white() {
        let pts = vec![Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0)];
        let sl = Streamline::new(pts, Color::white());
        let colors = color_by_velocity(&sl, 0.0, 1.0);
        for c in &colors {
            assert!((c.r - 1.0).abs() < 1e-6);
            assert!((c.g - 1.0).abs() < 1e-6);
            assert!((c.b - 1.0).abs() < 1e-6);
        }
    }
    #[test]
    fn test_arrow_glyphs_empty_when_short_line() {
        let sl = Streamline::new(vec![Vec3::zeros()], Color::white());
        let glyphs = arrow_glyphs_along_streamline(&sl, 1.0, 0.1, Color::red());
        assert!(glyphs.is_empty());
    }
    #[test]
    fn test_arrow_glyphs_count_for_long_line() {
        let pts: Vec<Vec3> = (0..=20)
            .map(|i| Vec3::new(i as f64 * 0.1, 0.0, 0.0))
            .collect();
        let sl = Streamline::new(pts, Color::white());
        let glyphs = arrow_glyphs_along_streamline(&sl, 0.5, 0.1, Color::red());
        assert!(!glyphs.is_empty(), "should place at least one arrow");
        assert!(
            glyphs.len() >= 3,
            "expected multiple arrows, got {}",
            glyphs.len()
        );
    }
    #[test]
    fn test_arrow_glyphs_direction_along_x() {
        let pts: Vec<Vec3> = (0..=10)
            .map(|i| Vec3::new(i as f64 * 0.1, 0.0, 0.0))
            .collect();
        let sl = Streamline::new(pts, Color::white());
        let glyphs = arrow_glyphs_along_streamline(&sl, 0.3, 0.05, Color::green());
        for g in &glyphs {
            assert!(
                g.direction.x > 0.9,
                "arrow should point along +x: {:?}",
                g.direction
            );
        }
    }
    #[test]
    fn test_velocity_ribbon_len_matches_streamline() {
        let pts: Vec<Vec3> = (0..=5).map(|i| Vec3::new(i as f64, 0.0, 0.0)).collect();
        let mags = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let sl = Streamline::with_magnitudes(pts, mags, Color::white());
        let ribbon = StreamRibbonVelocity::from_streamline(&sl, 0.1, Color::cyan());
        assert_eq!(ribbon.len(), sl.points.len());
    }
    #[test]
    fn test_velocity_ribbon_wider_at_higher_speed() {
        let pts: Vec<Vec3> = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
        ];
        let mags = vec![1.0, 2.0, 1.5];
        let sl = Streamline::with_magnitudes(pts, mags, Color::white());
        let ribbon = StreamRibbonVelocity::from_streamline(&sl, 0.1, Color::white());
        assert!(
            ribbon.half_widths[1] > ribbon.half_widths[0],
            "higher velocity → wider ribbon"
        );
    }
    #[test]
    fn test_velocity_ribbon_edges_offset() {
        let pts: Vec<Vec3> = (0..=3).map(|i| Vec3::new(i as f64, 0.0, 0.0)).collect();
        let sl = Streamline::new(pts, Color::white());
        let ribbon = StreamRibbonVelocity::from_streamline(&sl, 0.2, Color::white());
        let left = ribbon.left_edge();
        let right = ribbon.right_edge();
        assert_eq!(left.len(), ribbon.len());
        assert_eq!(right.len(), ribbon.len());
        let spread = (left[0] - right[0]).norm();
        assert!(spread > 0.0, "ribbon edges should be separated");
    }
    #[test]
    fn test_jobard_lefer_produces_at_least_one_line() {
        let params = JobardLeferParams::unit_square(0.2);
        let field = |_: [f64; 3]| [1.0_f64, 0.0, 0.0];
        let lines = compute_jobard_lefer(&params, &field);
        assert!(
            !lines.is_empty(),
            "Jobard–Lefer should produce at least one streamline"
        );
    }
    #[test]
    fn test_jobard_lefer_lines_have_points() {
        let params = JobardLeferParams::unit_square(0.25);
        let field = |_: [f64; 3]| [1.0_f64, 0.0, 0.0];
        let lines = compute_jobard_lefer(&params, &field);
        for sl in &lines {
            assert!(
                sl.points.len() >= 2,
                "each streamline should have ≥2 points"
            );
        }
    }
    #[test]
    fn test_jobard_lefer_circular_field_multiple_lines() {
        let params = JobardLeferParams {
            d_sep: 0.15,
            d_test: 0.075,
            step_size: 0.02,
            max_steps: 200,
            initial_seed: [0.8, 0.5],
            z_plane: 0.0,
            domain: [0.0, 1.0, 0.0, 1.0],
        };
        let field = |p: [f64; 3]| [-(p[1] - 0.5), p[0] - 0.5, 0.0];
        let lines = compute_jobard_lefer(&params, &field);
        assert!(
            !lines.is_empty(),
            "circular field should produce at least one streamline"
        );
    }

    // ---- λ₂ criterion: eigenvalues of S² + W² ----

    type Mat3 = [[f64; 3]; 3];

    fn mat_mul(a: &Mat3, b: &Mat3) -> Mat3 {
        std::array::from_fn(|i| {
            std::array::from_fn(|k| (0..3).map(|m| a[i][m] * b[m][k]).sum::<f64>())
        })
    }
    fn transpose(a: &Mat3) -> Mat3 {
        std::array::from_fn(|i| std::array::from_fn(|k| a[k][i]))
    }
    /// `A = S² + W²` built straight from the definition.
    fn s2_plus_w2(j: &Mat3) -> Mat3 {
        let jt = transpose(j);
        let s: Mat3 = std::array::from_fn(|i| std::array::from_fn(|k| 0.5 * (j[i][k] + jt[i][k])));
        let w: Mat3 = std::array::from_fn(|i| std::array::from_fn(|k| 0.5 * (j[i][k] - jt[i][k])));
        let (ss, ww) = (mat_mul(&s, &s), mat_mul(&w, &w));
        std::array::from_fn(|i| std::array::from_fn(|k| ss[i][k] + ww[i][k]))
    }
    /// Reference eigenvalues of a symmetric matrix from Jacobi rotations, ascending.
    fn jacobi_eigenvalues_ascending(a: &Mat3) -> [f64; 3] {
        let e = eig3(&[a[0][0], a[1][1], a[2][2], a[0][1], a[1][2], a[0][2]]);
        [e.values[2], e.values[1], e.values[0]]
    }
    /// Uniform pseudo-random number in `[0, 1)` from a 64-bit LCG.
    fn next_unit(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (*state >> 11) as f64 / (1_u64 << 53) as f64
    }
    /// Rotation matrix about the unit `axis` by `angle` (Rodrigues' formula).
    fn rotation(axis: [f64; 3], angle: f64) -> Mat3 {
        let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        let [x, y, z] = [axis[0] / norm, axis[1] / norm, axis[2] / norm];
        let k = [[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]];
        let k2 = mat_mul(&k, &k);
        let (s, c) = angle.sin_cos();
        std::array::from_fn(|i| {
            std::array::from_fn(|m| {
                let id = if i == m { 1.0 } else { 0.0 };
                id + s * k[i][m] + (1.0 - c) * k2[i][m]
            })
        })
    }
    fn linear_velocity(j: Mat3) -> impl Fn(Vec3) -> Vec3 {
        move |p: Vec3| {
            Vec3::new(
                j[0][0] * p.x + j[0][1] * p.y + j[0][2] * p.z,
                j[1][0] * p.x + j[1][1] * p.y + j[1][2] * p.z,
                j[2][0] * p.x + j[2][1] * p.y + j[2][2] * p.z,
            )
        }
    }
    #[test]
    fn test_lambda2_eigenvalues_general_gradient() {
        // The old Cardano variant returned [-0.706, -0.706, -0.068] here.
        let j = [[0.3, 1.0, 0.2], [-0.8, 0.1, 0.5], [0.4, -0.3, -0.4]];
        let eigs = lambda2_eigenvalues(&j);
        let expected = [-0.975479730706, -0.742856172094, 0.238335902800];
        for k in 0..3 {
            assert!(
                (eigs[k] - expected[k]).abs() < 1e-9,
                "eigenvalue {k}: {} vs {}",
                eigs[k],
                expected[k]
            );
        }
        let reference = jacobi_eigenvalues_ascending(&s2_plus_w2(&j));
        for k in 0..3 {
            assert!(
                (eigs[k] - reference[k]).abs() < 1e-9,
                "{eigs:?} vs {reference:?}"
            );
        }
    }
    #[test]
    fn test_lambda2_eigenvalues_match_jacobi_on_random_gradients() {
        let mut state = 0x2545_F491_4F6C_DD1D_u64;
        for case in 0..300 {
            let j: Mat3 =
                std::array::from_fn(|_| std::array::from_fn(|_| 4.0 * next_unit(&mut state) - 2.0));
            let a = s2_plus_w2(&j);
            let eigs = lambda2_eigenvalues(&j);
            let reference = jacobi_eigenvalues_ascending(&a);
            for k in 0..3 {
                assert!(
                    (eigs[k] - reference[k]).abs() < 1e-9,
                    "case {case}, eigenvalue {k}: {eigs:?} vs {reference:?}"
                );
            }
            assert!(
                eigs[0] <= eigs[1] && eigs[1] <= eigs[2],
                "{eigs:?} not ascending"
            );
            // tr(S² + W²) = |S|² − |W|² = −2Q.
            let sum: f64 = eigs.iter().sum();
            assert!(
                (sum + 2.0 * q_criterion(&j)).abs() < 1e-9,
                "case {case}: trace {sum} vs -2Q {}",
                -2.0 * q_criterion(&j)
            );
        }
    }
    #[test]
    fn test_lambda2_eigenvalues_diagonal_with_tiny_coupling() {
        // J = J^T, so A = J²; the 1e-6 coupling shifts the eigenvalues of
        // diag(1, 2, 3) by far less than 1e-6.
        let j = [
            [1.0, 1e-6, 0.0],
            [1e-6, 2.0_f64.sqrt(), 0.0],
            [0.0, 0.0, 3.0_f64.sqrt()],
        ];
        let eigs = lambda2_eigenvalues(&j);
        for (got, want) in eigs.iter().zip([1.0, 2.0, 3.0]) {
            assert!((got - want).abs() < 1e-6, "{eigs:?}");
        }
    }
    #[test]
    fn test_lambda2_eigenvalues_canonical_flows() {
        // Solid-body rotation: A = W² = diag(-1, -1, 0).
        let rotation_flow = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        assert_eq!(lambda2_eigenvalues(&rotation_flow), [-1.0, -1.0, 0.0]);
        // Simple shear u = (y, 0, 0): S² and W² cancel exactly, A = 0.
        let shear = [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        assert_eq!(lambda2_eigenvalues(&shear), [0.0; 3]);
        // Pure strain diag(1, -0.5, -0.5): A = diag(1, 0.25, 0.25).
        let strain = [[1.0, 0.0, 0.0], [0.0, -0.5, 0.0], [0.0, 0.0, -0.5]];
        assert_eq!(lambda2_eigenvalues(&strain), [0.25, 0.25, 1.0]);
    }
    #[test]
    fn test_lambda2_eigenvalues_rotated_symmetric_gradient() {
        // J = R·diag(d)·Rᵀ is symmetric, so A = J² has eigenvalues d².
        let d = [0.7, -1.3, 2.1];
        let r = rotation([1.0, 2.0, 2.0], 0.9);
        let diag: Mat3 = [[d[0], 0.0, 0.0], [0.0, d[1], 0.0], [0.0, 0.0, d[2]]];
        let j = mat_mul(&mat_mul(&r, &diag), &transpose(&r));
        let eigs = lambda2_eigenvalues(&j);
        let mut expected = [d[0] * d[0], d[1] * d[1], d[2] * d[2]];
        expected.sort_by(f64::total_cmp);
        for k in 0..3 {
            assert!(
                (eigs[k] - expected[k]).abs() < 1e-12,
                "{eigs:?} vs {expected:?}"
            );
        }
    }
    #[test]
    fn test_lambda2_eigenvalues_rotated_stretched_vortex() {
        // u = (-Ω·y - a·x/2, Ω·x - a·y/2, a·z) has A = diag(a²/4 - Ω², a²/4 - Ω², a²);
        // the similarity transform R·J·Rᵀ leaves the eigenvalues unchanged but
        // fills A.  The double eigenvalue limits the closed form to about √ε.
        let r = rotation([0.3, -1.0, 0.5], 1.1);
        for (omega, a) in [(2.0_f64, 1.0_f64), (0.2, 1.0)] {
            let j0 = [
                [-0.5 * a, -omega, 0.0],
                [omega, -0.5 * a, 0.0],
                [0.0, 0.0, a],
            ];
            let j = mat_mul(&mat_mul(&r, &j0), &transpose(&r));
            let eigs = lambda2_eigenvalues(&j);
            let double = 0.25 * a * a - omega * omega;
            let mut expected = [double, double, a * a];
            expected.sort_by(f64::total_cmp);
            for k in 0..3 {
                assert!(
                    (eigs[k] - expected[k]).abs() < 1e-6,
                    "omega {omega}: {eigs:?} vs {expected:?}"
                );
            }
            let field = linear_velocity(j);
            let cores = detect_vortex_cores(&[Vec3::new(0.2, -0.1, 0.4)], &field, 1e-3, 0.0, 0.0);
            assert_eq!(cores.len(), usize::from(omega > 0.5 * a), "omega {omega}");
        }
    }
    #[test]
    fn test_detect_vortex_cores_uses_the_true_lambda2() {
        let probe = [Vec3::new(0.3, -0.2, 0.1)];
        // λ₂ = +0.1555 > 0 although Q = 0.41 > 0.1: not a vortex.  The old
        // formula returned λ₂ = -0.144 and flagged it.
        let strain_dominated = [[-0.2, 0.5, 0.9], [-0.6, 0.0, 0.9], [-1.0, 0.5, 0.8]];
        let field = linear_velocity(strain_dominated);
        assert!(detect_vortex_cores(&probe, &field, 1e-3, 0.1, 0.0).is_empty());
        // λ₂ = -0.04275 < 0 with Q = 0.18: a vortex the old formula missed
        // (it returned λ₂ = +0.050).
        let weak_vortex = [[0.4, 0.6, 0.9], [0.7, 1.0, 0.8], [-0.6, -0.9, -0.4]];
        let field = linear_velocity(weak_vortex);
        let cores = detect_vortex_cores(&probe, &field, 1e-3, 0.1, 0.0);
        assert_eq!(cores.len(), 1);
        assert!(
            (cores[0].lambda2 + 0.042_750_055).abs() < 1e-6,
            "{}",
            cores[0].lambda2
        );
        assert!(
            (cores[0].q_criterion - 0.18).abs() < 1e-6,
            "{}",
            cores[0].q_criterion
        );
    }
    #[test]
    fn test_symmetric_eigenvalues_degenerate_and_special_matrices() {
        // I + 𝟙𝟙ᵀ has eigenvalues 1, 1, 4.
        let degenerate = [[2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]];
        let eigs = symmetric_eigenvalues_ascending(&degenerate);
        for (got, want) in eigs.iter().zip([1.0, 1.0, 4.0]) {
            assert!((got - want).abs() < 1e-6, "{eigs:?}");
        }
        // 2I - 𝟙𝟙ᵀ has eigenvalues -1, 2, 2 (the double one on top).
        let negative = [[1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]];
        let eigs = symmetric_eigenvalues_ascending(&negative);
        for (got, want) in eigs.iter().zip([-1.0, 2.0, 2.0]) {
            assert!((got - want).abs() < 1e-6, "{eigs:?}");
        }
        let scalar = [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]];
        assert_eq!(symmetric_eigenvalues_ascending(&scalar), [5.0; 3]);
        assert_eq!(symmetric_eigenvalues_ascending(&[[0.0; 3]; 3]), [0.0; 3]);
        // Entries far below 1 still resolve: the diagonal test is relative.
        let tiny = [[2e-20, 1e-20, 0.0], [1e-20, 2e-20, 0.0], [0.0, 0.0, 5e-20]];
        let eigs = symmetric_eigenvalues_ascending(&tiny);
        for (got, want) in eigs.iter().zip([1e-20, 3e-20, 5e-20]) {
            assert!((got - want).abs() < 1e-26, "{eigs:?}");
        }
        let poisoned = [[1.0, f64::NAN, 0.0], [f64::NAN, 2.0, 0.0], [0.0, 0.0, 3.0]];
        assert!(
            symmetric_eigenvalues_ascending(&poisoned)
                .iter()
                .all(|v| v.is_nan())
        );
        let nan_gradient = [[f64::NAN, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!(
            lambda2_eigenvalues(&nan_gradient)
                .iter()
                .all(|v| v.is_nan())
        );
        let field = linear_velocity(nan_gradient);
        assert!(detect_vortex_cores(&[Vec3::zeros()], &field, 1e-3, 0.0, 0.0).is_empty());
    }
    // ---- 2-D flow topology ----

    /// Classify the single critical point of `field` at the origin: one cell
    /// of side 1 centred there.
    fn classify_at_origin(field: &dyn Fn([f64; 3]) -> [f64; 3]) -> Vec<FlowTopologyPoint> {
        detect_topology_points([-0.5, 0.5], [-0.5, 0.5], 1.0, 0.0, field, 0.1)
    }
    fn kind_at_origin(field: &dyn Fn([f64; 3]) -> [f64; 3]) -> TopologyPointKind {
        let points = classify_at_origin(field);
        assert_eq!(points.len(), 1, "{points:?}");
        assert_eq!(points[0].position, [0.0, 0.0]);
        points[0].kind
    }
    #[test]
    fn test_topology_classifies_every_linear_critical_point() {
        // (y, x): J = [[0, 1], [1, 0]], det -1.  The old det = ∂u/∂x·∂v/∂y
        // omitted the cross term and reported a stagnation point.
        assert_eq!(
            kind_at_origin(&|p| [p[1], p[0], 0.0]),
            TopologyPointKind::Saddle
        );
        assert_eq!(
            kind_at_origin(&|p| [p[0], -p[1], 0.0]),
            TopologyPointKind::Saddle
        );
        // Rotation: complex eigenvalues ±i.
        assert_eq!(
            kind_at_origin(&|p| [-p[1], p[0], 0.0]),
            TopologyPointKind::VortexCentre
        );
        // Spiral foci (tr = ±2, det = 2): still circulating flow.
        assert_eq!(
            kind_at_origin(&|p| [p[0] - p[1], p[0] + p[1], 0.0]),
            TopologyPointKind::VortexCentre
        );
        assert_eq!(
            kind_at_origin(&|p| [-p[0] - p[1], p[0] - p[1], 0.0]),
            TopologyPointKind::VortexCentre
        );
        // Radial source / sink, including a star node (disc = 0) and a proper node.
        assert_eq!(
            kind_at_origin(&|p| [p[0], p[1], 0.0]),
            TopologyPointKind::Source
        );
        assert_eq!(
            kind_at_origin(&|p| [-p[0], -p[1], 0.0]),
            TopologyPointKind::Sink
        );
        assert_eq!(
            kind_at_origin(&|p| [2.0 * p[0], 3.0 * p[1], 0.0]),
            TopologyPointKind::Source
        );
        assert_eq!(
            kind_at_origin(&|p| [-2.0 * p[0], -3.0 * p[1], 0.0]),
            TopologyPointKind::Sink
        );
    }
    #[test]
    fn test_topology_singular_gradient_is_a_degenerate_stagnation_point() {
        assert_eq!(
            kind_at_origin(&|_| [0.0, 0.0, 0.0]),
            TopologyPointKind::Stagnation
        );
        // Shear (y, 0): J = [[0, 1], [0, 0]] has det 0.
        assert_eq!(
            kind_at_origin(&|p| [p[1], 0.0, 0.0]),
            TopologyPointKind::Stagnation
        );
        // Rank-one gradient (x, 0).
        assert_eq!(
            kind_at_origin(&|p| [p[0], 0.0, 0.0]),
            TopologyPointKind::Stagnation
        );
    }
    #[test]
    fn test_topology_classification_is_scale_free() {
        // A very slow saddle: det = -1e-18 would pass an absolute
        // degeneracy threshold.
        let tiny = 1e-9;
        let saddle = move |p: [f64; 3]| [tiny * p[1], tiny * p[0], 0.0];
        let points = detect_topology_points([-0.5, 0.5], [-0.5, 0.5], 1.0, 0.0, &saddle, 1.0);
        assert_eq!(points.len(), 1);
        assert_eq!(points[0].kind, TopologyPointKind::Saddle);
        // The same flow at huge amplitude.
        let huge = 1e9;
        let vortex = move |p: [f64; 3]| [-huge * p[1], huge * p[0], 0.0];
        let points = detect_topology_points([-0.5, 0.5], [-0.5, 0.5], 1.0, 0.0, &vortex, 1e12);
        assert_eq!(points.len(), 1);
        assert_eq!(points[0].kind, TopologyPointKind::VortexCentre);
    }
    #[test]
    fn test_topology_rejects_ill_defined_scans_instead_of_hanging() {
        let field = |_: [f64; 3]| [0.0_f64, 0.0, 0.0];
        let scan = |x: [f64; 2], y: [f64; 2], cell: f64| {
            detect_topology_points(x, y, cell, 0.0, &field, 1.0)
        };
        // A non-positive cell size never advanced the old `x += cell_size` loop.
        assert!(scan([0.0, 1.0], [0.0, 1.0], 0.0).is_empty());
        assert!(scan([0.0, 1.0], [0.0, 1.0], -0.1).is_empty());
        assert!(scan([0.0, 1.0], [0.0, 1.0], f64::NAN).is_empty());
        assert!(scan([0.0, 1.0], [0.0, 1.0], f64::INFINITY).is_empty());
        // An infinite range, or one needing more cells than `usize` holds.
        assert!(scan([0.0, f64::INFINITY], [0.0, 1.0], 0.5).is_empty());
        assert!(scan([0.0, 1.0], [f64::NEG_INFINITY, 1.0], 0.5).is_empty());
        assert!(scan([0.0, f64::NAN], [0.0, 1.0], 0.5).is_empty());
        assert!(scan([0.0, 1e300], [0.0, 1.0], 1e-10).is_empty());
        // Empty and reversed ranges contain no cells.
        assert!(scan([1.0, 1.0], [0.0, 1.0], 0.5).is_empty());
        assert!(scan([0.0, 1.0], [1.0, 0.0], 0.5).is_empty());
    }
    #[test]
    fn test_topology_scans_every_cell_centre_once() {
        let field = |_: [f64; 3]| [0.0_f64, 0.0, 0.0];
        let points = detect_topology_points([0.0, 1.0], [0.0, 0.5], 0.25, 0.0, &field, 1.0);
        assert_eq!(points.len(), 4 * 2);
        for (k, p) in points.iter().enumerate() {
            let (ix, iy) = (k / 2, k % 2);
            assert!(
                (p.position[0] - (0.125 + 0.25 * ix as f64)).abs() < 1e-12,
                "{p:?}"
            );
            assert!(
                (p.position[1] - (0.125 + 0.25 * iy as f64)).abs() < 1e-12,
                "{p:?}"
            );
            assert_eq!(p.kind, TopologyPointKind::Stagnation);
            assert_eq!(p.velocity_magnitude, 0.0);
        }
        // The last cell may overhang a range that is not a multiple of the cell size.
        let points = detect_topology_points([0.0, 1.0], [0.0, 0.3], 0.3, 0.0, &field, 1.0);
        assert_eq!(points.len(), 4);
        assert!((points[3].position[0] - 1.05).abs() < 1e-12);
    }
    #[test]
    fn test_topology_skips_fast_and_non_finite_cells() {
        // Only the cell centred on the saddle of (x - 0.5, 0.5 - y) is slow.
        let field = |p: [f64; 3]| [p[0] - 0.5, 0.5 - p[1], 0.0];
        let points = detect_topology_points([0.0, 1.0], [0.0, 1.0], 0.5, 0.0, &field, 0.1);
        assert!(
            points.is_empty(),
            "no cell centre coincides with the saddle: {points:?}"
        );
        let points = detect_topology_points([0.0, 1.0], [0.0, 1.0], 1.0, 0.0, &field, 0.1);
        assert_eq!(points.len(), 1);
        assert_eq!(points[0].kind, TopologyPointKind::Saddle);
        assert_eq!(points[0].position, [0.5, 0.5]);
        // NaN velocities never form a candidate, and a non-finite stencil is skipped.
        let nan_field = |_: [f64; 3]| [f64::NAN, 0.0, 0.0];
        assert!(classify_at_origin(&nan_field).is_empty());
        let half_nan = |p: [f64; 3]| {
            if p[0] > 0.0 {
                [f64::NAN, 0.0, 0.0]
            } else {
                [0.0, 0.0, 0.0]
            }
        };
        assert!(classify_at_origin(&half_nan).is_empty());
    }
}
