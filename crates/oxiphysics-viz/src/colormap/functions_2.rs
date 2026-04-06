//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::types::*;
#[cfg(test)]
mod tests_new_colormap {
    use super::*;
    use crate::Color;
    #[test]
    fn test_builder_empty_returns_black() {
        let builder = CustomColormapBuilder::new();
        let c = builder.sample(0.5);
        assert!(c.r < 0.01 && c.g < 0.01 && c.b < 0.01);
    }
    #[test]
    fn test_builder_single_stop() {
        let color = Color {
            r: 0.5,
            g: 0.0,
            b: 0.8,
            a: 1.0,
        };
        let builder = CustomColormapBuilder::new().add_stop(0.5, color);
        let c = builder.sample(0.0);
        assert!((c.r - 0.5).abs() < 1e-4);
    }
    #[test]
    fn test_builder_two_stops_midpoint() {
        let black = Color {
            r: 0.0,
            g: 0.0,
            b: 0.0,
            a: 1.0,
        };
        let white = Color {
            r: 1.0,
            g: 1.0,
            b: 1.0,
            a: 1.0,
        };
        let builder = CustomColormapBuilder::new()
            .add_stop(0.0, black)
            .add_stop(1.0, white);
        let mid = builder.sample(0.5);
        assert!(
            (mid.r - 0.5).abs() < 1e-4,
            "r midpoint should be 0.5, got {}",
            mid.r
        );
    }
    #[test]
    fn test_builder_clamped_to_endpoints() {
        let red = Color {
            r: 1.0,
            g: 0.0,
            b: 0.0,
            a: 1.0,
        };
        let blue = Color {
            r: 0.0,
            g: 0.0,
            b: 1.0,
            a: 1.0,
        };
        let builder = CustomColormapBuilder::new()
            .add_stop(0.0, red)
            .add_stop(1.0, blue);
        let c0 = builder.sample(-0.5);
        assert!((c0.r - 1.0).abs() < 1e-4, "should clamp to first stop");
        let c1 = builder.sample(1.5);
        assert!((c1.b - 1.0).abs() < 1e-4, "should clamp to last stop");
    }
    #[test]
    fn test_builder_smooth_step_different_from_linear() {
        let black = Color {
            r: 0.0,
            g: 0.0,
            b: 0.0,
            a: 1.0,
        };
        let white = Color {
            r: 1.0,
            g: 1.0,
            b: 1.0,
            a: 1.0,
        };
        let linear = CustomColormapBuilder::new()
            .add_stop(0.0, black)
            .add_stop(1.0, white)
            .with_interpolation(LinearInterp::Linear);
        let smooth = CustomColormapBuilder::new()
            .add_stop(0.0, black)
            .add_stop(1.0, white)
            .with_interpolation(LinearInterp::SmoothStep);
        let t = 0.25_f64;
        let lv = linear.sample(t).r;
        let sv = smooth.sample(t).r;
        assert!(
            (lv - sv).abs() > 0.05,
            "smooth vs linear should differ at t=0.25"
        );
    }
    #[test]
    fn test_builder_build_length() {
        let builder = CustomColormapBuilder::new()
            .add_stop(
                0.0,
                Color {
                    r: 0.0,
                    g: 0.0,
                    b: 0.0,
                    a: 1.0,
                },
            )
            .add_stop(
                1.0,
                Color {
                    r: 1.0,
                    g: 1.0,
                    b: 1.0,
                    a: 1.0,
                },
            );
        let palette = builder.build(16);
        assert_eq!(palette.len(), 16);
    }
    #[test]
    fn test_builder_build_endpoints() {
        let red = Color {
            r: 1.0,
            g: 0.0,
            b: 0.0,
            a: 1.0,
        };
        let blue = Color {
            r: 0.0,
            g: 0.0,
            b: 1.0,
            a: 1.0,
        };
        let builder = CustomColormapBuilder::new()
            .add_stop(0.0, red)
            .add_stop(1.0, blue);
        let pal = builder.build(5);
        assert!((pal[0].r - 1.0).abs() < 1e-4, "first color should be red");
        assert!((pal[4].b - 1.0).abs() < 1e-4, "last color should be blue");
    }
    #[test]
    fn test_spline_empty_returns_black() {
        let sp = SplineColormap::new(std::iter::empty());
        let c = sp.sample(0.5);
        assert!(c.r < 0.01 && c.g < 0.01 && c.b < 0.01);
    }
    #[test]
    fn test_spline_single_stop() {
        let c = Color {
            r: 0.3,
            g: 0.7,
            b: 0.1,
            a: 1.0,
        };
        let sp = SplineColormap::new([(0.5, c)]);
        let sampled = sp.sample(0.0);
        assert!((sampled.r - 0.3).abs() < 1e-4);
    }
    #[test]
    fn test_spline_two_stops_midpoint() {
        let black = Color {
            r: 0.0,
            g: 0.0,
            b: 0.0,
            a: 1.0,
        };
        let white = Color {
            r: 1.0,
            g: 1.0,
            b: 1.0,
            a: 1.0,
        };
        let sp = SplineColormap::new([(0.0, black), (1.0, white)]);
        let mid = sp.sample(0.5);
        assert!((mid.r - 0.5).abs() < 0.05, "spline 2-stop midpoint ≈ 0.5");
    }
    #[test]
    fn test_spline_four_stops_sorted() {
        let stops = vec![
            (
                0.0,
                Color {
                    r: 0.0,
                    g: 0.0,
                    b: 1.0,
                    a: 1.0,
                },
            ),
            (
                0.33,
                Color {
                    r: 0.0,
                    g: 1.0,
                    b: 1.0,
                    a: 1.0,
                },
            ),
            (
                0.66,
                Color {
                    r: 1.0,
                    g: 1.0,
                    b: 0.0,
                    a: 1.0,
                },
            ),
            (
                1.0,
                Color {
                    r: 1.0,
                    g: 0.0,
                    b: 0.0,
                    a: 1.0,
                },
            ),
        ];
        let sp = SplineColormap::new(stops);
        let c = sp.sample(0.0);
        assert!(c.b > 0.5, "t=0 should be bluish");
    }
    #[test]
    fn test_spline_build_length() {
        let sp = SplineColormap::new([
            (
                0.0,
                Color {
                    r: 0.0,
                    g: 0.0,
                    b: 0.0,
                    a: 1.0,
                },
            ),
            (
                1.0,
                Color {
                    r: 1.0,
                    g: 1.0,
                    b: 1.0,
                    a: 1.0,
                },
            ),
        ]);
        assert_eq!(sp.build(20).len(), 20);
    }
    #[test]
    fn test_spline_in_range() {
        let sp = SplineColormap::new([
            (
                0.0,
                Color {
                    r: 1.0,
                    g: 0.0,
                    b: 0.0,
                    a: 1.0,
                },
            ),
            (
                0.5,
                Color {
                    r: 0.0,
                    g: 1.0,
                    b: 0.0,
                    a: 1.0,
                },
            ),
            (
                1.0,
                Color {
                    r: 0.0,
                    g: 0.0,
                    b: 1.0,
                    a: 1.0,
                },
            ),
        ]);
        for i in 0..=20 {
            let c = sp.sample(i as f64 / 20.0);
            assert!(c.r >= 0.0 && c.r <= 1.0, "r out of range: {}", c.r);
            assert!(c.g >= 0.0 && c.g <= 1.0, "g out of range: {}", c.g);
            assert!(c.b >= 0.0 && c.b <= 1.0, "b out of range: {}", c.b);
        }
    }
    #[test]
    fn test_cyclic_hsv_wraps() {
        let cm = CyclicColormap::new(CyclicStyle::HsvWheel);
        let c0 = cm.sample(0.0);
        let c1 = cm.sample(1.0);
        let diff = (c0.r - c1.r).abs() + (c0.g - c1.g).abs() + (c0.b - c1.b).abs();
        assert!(diff < 0.1, "HSV wheel should wrap: diff={diff}");
    }
    #[test]
    fn test_cyclic_twilight_in_range() {
        let cm = CyclicColormap::new(CyclicStyle::Twilight);
        for i in 0..=20 {
            let c = cm.sample(i as f64 / 20.0);
            assert!(c.r >= 0.0 && c.r <= 1.0);
            assert!(c.g >= 0.0 && c.g <= 1.0);
            assert!(c.b >= 0.0 && c.b <= 1.0);
        }
    }
    #[test]
    fn test_cyclic_phase_in_range() {
        let cm = CyclicColormap::new(CyclicStyle::Phase);
        for i in 0..=20 {
            let c = cm.sample(i as f64 / 20.0);
            assert!(c.r >= 0.0 && c.r <= 1.0);
            assert!(c.g >= 0.0 && c.g <= 1.0);
            assert!(c.b >= 0.0 && c.b <= 1.0);
        }
    }
    #[test]
    fn test_cyclic_multi_cycle_wraps() {
        let cm = CyclicColormap::new(CyclicStyle::HsvWheel).with_cycles(3.0);
        let c0 = cm.sample(0.0);
        let c3 = cm.sample(1.0 / 3.0);
        let diff = (c0.r - c3.r).abs() + (c0.g - c3.g).abs() + (c0.b - c3.b).abs();
        assert!(diff < 0.1, "Multi-cycle wrap failed: diff={diff}");
    }
    #[test]
    fn test_cyclic_build_length() {
        let cm = CyclicColormap::new(CyclicStyle::Twilight);
        assert_eq!(cm.build(12).len(), 12);
    }
    #[test]
    fn test_hsv_rainbow_range() {
        let cm = HsvColormap::rainbow();
        for i in 0..=20 {
            let c = cm.sample(i as f64 / 20.0);
            assert!(c.r >= 0.0 && c.r <= 1.0);
            assert!(c.g >= 0.0 && c.g <= 1.0);
            assert!(c.b >= 0.0 && c.b <= 1.0);
        }
    }
    #[test]
    fn test_hsv_rainbow_red_at_t0() {
        let cm = HsvColormap::rainbow();
        let c = cm.sample(0.0);
        assert!(
            c.r > 0.9 && c.g < 0.1,
            "t=0 should be red: r={}, g={}",
            c.r,
            c.g
        );
    }
    #[test]
    fn test_hsv_custom_build_length() {
        let cm = HsvColormap::new(0.0, 1.0, 0.8, 0.9);
        assert_eq!(cm.build(8).len(), 8);
    }
    #[test]
    fn test_hsv_saturation_zero_is_grey() {
        let cm = HsvColormap::new(0.3, 0.7, 0.0, 0.7);
        let c = cm.sample(0.5);
        assert!(
            (c.r - c.g).abs() < 1e-4 && (c.g - c.b).abs() < 1e-4,
            "S=0 should produce grey: r={} g={} b={}",
            c.r,
            c.g,
            c.b
        );
    }
    #[test]
    fn test_turbo_in_range() {
        for i in 0..=20 {
            let c = sample_turbo(i as f32 / 20.0);
            assert!(c.r >= 0.0 && c.r <= 1.0, "r out of range at t={}", i);
            assert!(c.g >= 0.0 && c.g <= 1.0, "g out of range at t={}", i);
            assert!(c.b >= 0.0 && c.b <= 1.0, "b out of range at t={}", i);
        }
    }
    #[test]
    fn test_turbo_endpoints_differ() {
        let c0 = sample_turbo(0.0);
        let c1 = sample_turbo(1.0);
        let diff = (c0.r - c1.r).abs() + (c0.g - c1.g).abs() + (c0.b - c1.b).abs();
        assert!(
            diff > 0.3,
            "turbo endpoints should differ significantly: diff={diff}"
        );
    }
    #[test]
    fn test_jet_in_range() {
        for i in 0..=20 {
            let c = sample_jet(i as f32 / 20.0);
            assert!(c.r >= 0.0 && c.r <= 1.0);
            assert!(c.g >= 0.0 && c.g <= 1.0);
            assert!(c.b >= 0.0 && c.b <= 1.0);
        }
    }
    #[test]
    fn test_jet_blue_at_zero() {
        let c = sample_jet(0.0);
        assert!(c.b > 0.4, "jet at t=0 should be blue: b={}", c.b);
        assert!(c.r < 0.1, "jet at t=0 should have little red: r={}", c.r);
    }
    #[test]
    fn test_jet_red_at_one() {
        let c = sample_jet(1.0);
        assert!(c.r > 0.4, "jet at t=1 should be reddish: r={}", c.r);
        assert!(c.b < 0.1, "jet at t=1 should have little blue: b={}", c.b);
    }
    #[test]
    fn test_custom_uniformity_score_constant_is_zero() {
        let score = custom_colormap_uniformity_score(
            |_| Color {
                r: 0.5,
                g: 0.5,
                b: 0.5,
                a: 1.0,
            },
            20,
        );
        assert!(
            score < 1e-6,
            "constant colormap should have zero uniformity score"
        );
    }
    #[test]
    fn test_lightness_white() {
        let white = Color {
            r: 1.0,
            g: 1.0,
            b: 1.0,
            a: 1.0,
        };
        let l = lightness(white);
        assert!(l > 95.0, "white should have L* near 100, got {l}");
    }
    #[test]
    fn test_lightness_black() {
        let black = Color {
            r: 0.0,
            g: 0.0,
            b: 0.0,
            a: 1.0,
        };
        let l = lightness(black);
        assert!(l < 5.0, "black should have L* near 0, got {l}");
    }
    #[test]
    fn test_is_monotone_lightness_linear_grey() {
        let result = is_monotone_lightness(
            |t| Color {
                r: t as f32,
                g: t as f32,
                b: t as f32,
                a: 1.0,
            },
            20,
        );
        assert!(result, "linear grey ramp should be monotone in lightness");
    }
    #[test]
    fn test_adjust_brightness_doubles() {
        let c = Color {
            r: 0.3,
            g: 0.4,
            b: 0.2,
            a: 1.0,
        };
        let brighter = adjust_brightness(c, 2.0);
        assert!((brighter.r - 0.6).abs() < 1e-4, "r should double");
        assert!((brighter.g - 0.8).abs() < 1e-4, "g should double");
    }
    #[test]
    fn test_adjust_brightness_clamps() {
        let c = Color {
            r: 0.8,
            g: 0.9,
            b: 1.0,
            a: 1.0,
        };
        let bright = adjust_brightness(c, 10.0);
        assert!(
            bright.r <= 1.0 && bright.g <= 1.0 && bright.b <= 1.0,
            "should clamp to 1.0"
        );
    }
    #[test]
    fn test_adjust_saturation_zero_is_grey() {
        let c = Color {
            r: 1.0,
            g: 0.0,
            b: 0.0,
            a: 1.0,
        };
        let grey = adjust_saturation(c, 0.0);
        assert!(
            (grey.r - grey.g).abs() < 0.05 && (grey.g - grey.b).abs() < 0.05,
            "zero saturation should be grey"
        );
    }
    #[test]
    fn test_map_scalar_custom_midpoint() {
        let c = map_scalar_custom(5.0, 0.0, 10.0, |t| Color {
            r: t as f32,
            g: 0.0,
            b: 0.0,
            a: 1.0,
        });
        assert!((c.r - 0.5).abs() < 1e-4, "midpoint should map to r=0.5");
    }
    #[test]
    fn test_apply_custom_colormap_length() {
        let values = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let result = apply_custom_colormap(&values, 0.0, 5.0, |t| Color {
            r: t as f32,
            g: t as f32,
            b: t as f32,
            a: 1.0,
        });
        assert_eq!(result.len(), 6);
    }
    #[test]
    fn test_apply_custom_colormap_alpha_255() {
        let values = vec![0.5];
        let result = apply_custom_colormap(&values, 0.0, 1.0, |_| Color {
            r: 0.0,
            g: 0.0,
            b: 0.0,
            a: 1.0,
        });
        assert_eq!(result[0][3], 255);
    }
    #[test]
    fn test_twilight_in_range() {
        for i in 0..=20 {
            let c = sample_twilight(i as f32 / 20.0);
            assert!(c.r >= 0.0 && c.r <= 1.0);
            assert!(c.g >= 0.0 && c.g <= 1.0);
            assert!(c.b >= 0.0 && c.b <= 1.0);
        }
    }
    #[test]
    fn test_phase_in_range() {
        for i in 0..=20 {
            let c = sample_phase(i as f32 / 20.0);
            assert!(c.r >= 0.0 && c.r <= 1.0);
            assert!(c.g >= 0.0 && c.g <= 1.0);
            assert!(c.b >= 0.0 && c.b <= 1.0);
        }
    }
}
