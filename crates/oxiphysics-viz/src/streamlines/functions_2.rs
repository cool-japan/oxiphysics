//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use crate::Color;
    use crate::streamlines::*;
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
}
