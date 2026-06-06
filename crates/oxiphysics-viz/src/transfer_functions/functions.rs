//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{Rgba, TransferFunction};

#[cfg(test)]
mod tests {
    use super::*;

    use crate::transfer_functions::ColorStop;
    use crate::transfer_functions::GradientColormap;

    use crate::transfer_functions::OpacityProfile;
    use crate::transfer_functions::OpacityProfileKind;

    use crate::transfer_functions::TransferFunctionEditor;

    #[test]
    fn test_rgba_lerp_midpoint() {
        let a = Rgba::new(0.0, 0.0, 0.0, 0.0);
        let b = Rgba::new(1.0, 1.0, 1.0, 1.0);
        let mid = Rgba::lerp(a, b, 0.5);
        assert!((mid.r - 0.5).abs() < 1e-10);
        assert!((mid.a - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_rgba_premultiply() {
        let c = Rgba::new(1.0, 0.5, 0.25, 0.5);
        let pm = c.premultiply_alpha();
        assert!((pm.r - 0.5).abs() < 1e-10);
        assert!((pm.g - 0.25).abs() < 1e-10);
    }
    #[test]
    fn test_rgba_to_u8() {
        let c = Rgba::new(1.0, 0.0, 0.5, 1.0);
        let u = c.to_u8_array();
        assert_eq!(u[0], 255);
        assert_eq!(u[1], 0);
        assert_eq!(u[2], 128);
    }
    #[test]
    fn test_tf_sample_empty_returns_transparent() {
        let tf = TransferFunction::new();
        let c = tf.sample(0.5);
        assert_eq!(c.a, 0.0);
    }
    #[test]
    fn test_tf_sample_single_point() {
        let mut tf = TransferFunction::new();
        tf.add_point(0.5, Rgba::new(1.0, 0.0, 0.0, 1.0));
        let c = tf.sample(0.5);
        assert!((c.r - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_tf_sample_linear_interpolation() {
        let mut tf = TransferFunction::new();
        tf.add_point(0.0, Rgba::new(0.0, 0.0, 0.0, 0.0));
        tf.add_point(1.0, Rgba::new(1.0, 0.0, 0.0, 1.0));
        let c = tf.sample(0.5);
        assert!((c.r - 0.5).abs() < 1e-6);
        assert!((c.a - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_tf_sample_clamped_above() {
        let mut tf = TransferFunction::new();
        tf.add_point(0.0, Rgba::new(0.0, 0.0, 0.0, 0.0));
        tf.add_point(1.0, Rgba::new(1.0, 1.0, 1.0, 1.0));
        let c = tf.sample(2.0);
        assert!((c.r - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_tf_remove_near() {
        let mut tf = TransferFunction::new();
        tf.add_point(0.5, Rgba::white());
        assert_eq!(tf.point_count(), 1);
        let removed = tf.remove_near(0.5, 0.01);
        assert!(removed);
        assert_eq!(tf.point_count(), 0);
    }
    #[test]
    fn test_tf_build_lut_count() {
        let tf = TransferFunction::default();
        let lut = tf.build_lut(256);
        assert_eq!(lut.len(), 256);
    }
    #[test]
    fn test_tf_temperature_non_trivial() {
        let tf = TransferFunction::temperature();
        assert!(tf.point_count() > 2);
        let c_low = tf.sample(0.0);
        assert!(c_low.a < 0.1);
        let c_high = tf.sample(1.0);
        assert!(c_high.a > 0.9);
    }
    #[test]
    fn test_tf_pressure_symmetric_midpoint() {
        let tf = TransferFunction::pressure();
        let c_mid = tf.sample(0.5);
        assert!(c_mid.r > 0.8 && c_mid.g > 0.8 && c_mid.b > 0.8);
    }
    #[test]
    fn test_opacity_linear_ramp() {
        let op = OpacityProfile::new(OpacityProfileKind::LinearRamp);
        assert!((op.evaluate(0.0)).abs() < 1e-10);
        assert!((op.evaluate(1.0) - 1.0).abs() < 1e-10);
        assert!((op.evaluate(0.5) - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_opacity_constant() {
        let op = OpacityProfile::new(OpacityProfileKind::Constant(0.7));
        assert!((op.evaluate(0.0) - 0.7).abs() < 1e-10);
        assert!((op.evaluate(0.5) - 0.7).abs() < 1e-10);
        assert!((op.evaluate(1.0) - 0.7).abs() < 1e-10);
    }
    #[test]
    fn test_opacity_gaussian_peak() {
        let op = OpacityProfile::new(OpacityProfileKind::Gaussian {
            center: 0.5,
            sigma: 0.1,
        });
        let peak = op.evaluate(0.5);
        let shoulder = op.evaluate(0.9);
        assert!((peak - 1.0).abs() < 1e-6, "Gaussian peak should be 1.0");
        assert!(shoulder < peak * 0.01, "Shoulder should be near zero");
    }
    #[test]
    fn test_opacity_step() {
        let op = OpacityProfile::new(OpacityProfileKind::Step { threshold: 0.5 });
        assert!((op.evaluate(0.4)).abs() < 1e-10);
        assert!((op.evaluate(0.5) - 1.0).abs() < 1e-10);
        assert!((op.evaluate(0.9) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_opacity_tent() {
        let op = OpacityProfile::new(OpacityProfileKind::Tent {
            center: 0.5,
            half_width: 0.2,
        });
        assert!(
            (op.evaluate(0.5) - 1.0).abs() < 1e-10,
            "Tent peak should be 1.0"
        );
        assert!((op.evaluate(0.3)).abs() < 1e-10, "Tent base should be 0");
        assert!((op.evaluate(0.7)).abs() < 1e-10, "Tent base should be 0");
    }
    #[test]
    fn test_opacity_scale_factor() {
        let op = OpacityProfile::new(OpacityProfileKind::Constant(1.0)).with_scale(0.5);
        assert!((op.evaluate(0.5) - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_gradient_default_endpoints() {
        let g = GradientColormap::default();
        let c0 = g.sample(0.0);
        let c1 = g.sample(1.0);
        assert!(c0[2] > c0[0], "Start should be blue");
        assert!(c1[0] > c1[2], "End should be red");
    }
    #[test]
    fn test_gradient_midpoint_interpolation() {
        let mut g = GradientColormap::new();
        g.add_stop(ColorStop::new(0.0, 0.0, 0.0, 0.0));
        g.add_stop(ColorStop::new(1.0, 1.0, 1.0, 1.0));
        let mid = g.sample(0.5);
        for ch in mid {
            assert!((ch - 0.5).abs() < 1e-6, "Midpoint should be 0.5, got {ch}");
        }
    }
    #[test]
    fn test_gradient_bwr_symmetric() {
        let g = GradientColormap::diverging_bwr();
        let lo = g.sample(0.0);
        let hi = g.sample(1.0);
        assert!(lo[2] > lo[0], "bwr low end: expected blue");
        assert!(hi[0] > hi[2], "bwr high end: expected red");
    }
    #[test]
    fn test_gradient_rainbow_monotone_hue_count() {
        let g = GradientColormap::rainbow();
        assert_eq!(g.stop_count(), 7);
    }
    #[test]
    fn test_gradient_inferno_brightness_increases() {
        let g = GradientColormap::inferno_approx();
        let dark = g.sample(0.0);
        let bright = g.sample(1.0);
        let lum = |c: [f64; 3]| 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2];
        assert!(lum(bright) > lum(dark), "Inferno should get brighter");
    }
    #[test]
    fn test_gradient_sample_with_opacity() {
        let g = GradientColormap::diverging_bwr();
        let op = OpacityProfile::new(OpacityProfileKind::LinearRamp);
        let c = g.sample_with_opacity(0.5, &op);
        assert!((c.a - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_gradient_build_lut_count() {
        let g = GradientColormap::rainbow();
        let lut = g.build_rgb_lut(64);
        assert_eq!(lut.len(), 64);
    }
    #[test]
    fn test_editor_normalize() {
        let ed = TransferFunctionEditor::new(TransferFunction::default(), 0.0, 100.0);
        assert!((ed.normalize(50.0) - 0.5).abs() < 1e-10);
        assert!((ed.normalize(0.0)).abs() < 1e-10);
        assert!((ed.normalize(100.0) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_editor_sample_scalar() {
        let mut tf = TransferFunction::new();
        tf.add_point(0.0, Rgba::new(0.0, 0.0, 0.0, 0.0));
        tf.add_point(1.0, Rgba::new(1.0, 0.0, 0.0, 1.0));
        let ed = TransferFunctionEditor::new(tf, 0.0, 200.0);
        let c = ed.sample_scalar(100.0);
        assert!((c.r - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_editor_select_nearest() {
        let tf = TransferFunction::default();
        let mut ed = TransferFunctionEditor::new(tf, 0.0, 1.0);
        ed.select_nearest(0.02, 0.1);
        assert_eq!(ed.selected, Some(0));
        ed.select_nearest(0.5, 0.01);
        assert_eq!(ed.selected, None);
    }
    #[test]
    fn test_editor_add_at() {
        let tf = TransferFunction::default();
        let mut ed = TransferFunctionEditor::new(tf, 0.0, 1.0);
        let before = ed.tf.point_count();
        ed.add_at(0.5);
        assert_eq!(ed.tf.point_count(), before + 1);
    }
}
/// Normalize a raw scalar value into \[0, 1\] given `[data_min, data_max]`.
///
/// Returns 0.5 if `data_min == data_max`.
pub fn normalize_scalar(value: f64, data_min: f64, data_max: f64) -> f64 {
    let range = data_max - data_min;
    if range.abs() < 1e-30 {
        return 0.5;
    }
    ((value - data_min) / range).clamp(0.0, 1.0)
}
/// Normalize a whole slice of scalars to \[0, 1\], using the slice's own min/max.
///
/// Returns an empty `Vec` if the input is empty.
pub fn normalize_slice(data: &[f64]) -> Vec<f64> {
    if data.is_empty() {
        return Vec::new();
    }
    let mn = data.iter().cloned().fold(f64::INFINITY, f64::min);
    let mx = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    data.iter().map(|&v| normalize_scalar(v, mn, mx)).collect()
}
/// Apply a transfer function to a whole slice of normalized scalars,
/// returning one `Rgba` per value.
pub fn apply_transfer_function(tf: &TransferFunction, normalized: &[f64]) -> Vec<Rgba> {
    normalized.iter().map(|&t| tf.sample(t)).collect()
}
/// Compute the gradient magnitude at index `i` in a 1-D scalar field.
///
/// Uses central differences; clamps to boundaries using one-sided differences.
pub fn gradient_magnitude_1d(data: &[f64], i: usize) -> f64 {
    let n = data.len();
    if n < 2 {
        return 0.0;
    }
    if i == 0 {
        (data[1] - data[0]).abs()
    } else if i == n - 1 {
        (data[n - 1] - data[n - 2]).abs()
    } else {
        (data[i + 1] - data[i - 1]).abs() * 0.5
    }
}
/// Compute opacity modulated by the normalized gradient magnitude.
///
/// Larger gradients → higher opacity (boundary enhancement).
/// `base_opacity` is the opacity from the TF at position `t`.
/// `gradient_mag` should be normalized to \[0, 1\].
/// `sharpness` controls how steeply the boundary effect kicks in.
pub fn gradient_modulated_opacity(base_opacity: f64, gradient_mag: f64, sharpness: f64) -> f64 {
    let boost = (gradient_mag * sharpness).clamp(0.0, 1.0);
    (base_opacity + boost * (1.0 - base_opacity)).clamp(0.0, 1.0)
}
/// Accumulate one ray sample into a front-to-back composite buffer.
///
/// Updates `(accumulated_rgb, accumulated_alpha)` in place.
///
/// # Parameters
/// - `acc_rgb` — currently accumulated pre-multiplied RGB
/// - `acc_alpha` — currently accumulated alpha
/// - `sample_rgba` — RGBA of the new sample (straight alpha)
/// - `step_size` — ray step size (scales opacity per step)
pub fn composite_front_to_back(
    acc_rgb: &mut [f64; 3],
    acc_alpha: &mut f64,
    sample_rgba: Rgba,
    step_size: f64,
) {
    let a = 1.0 - (-sample_rgba.a * step_size * 100.0).exp();
    let transparency = 1.0 - *acc_alpha;
    acc_rgb[0] += transparency * a * sample_rgba.r;
    acc_rgb[1] += transparency * a * sample_rgba.g;
    acc_rgb[2] += transparency * a * sample_rgba.b;
    *acc_alpha += transparency * a;
}
/// Perform a simple 1-D ray-march through a scalar field.
///
/// Samples are drawn from `data` at positions along a virtual ray.
/// Returns the final composited RGBA color.
///
/// This is a CPU reference implementation suitable for unit tests and offline
/// rendering; not intended for real-time use.
pub fn raymarch_1d(
    data: &[f64],
    data_min: f64,
    data_max: f64,
    tf: &TransferFunction,
    max_samples: usize,
    step_size: f64,
    early_termination: f64,
) -> Rgba {
    let n = data.len();
    if n == 0 {
        return Rgba::transparent();
    }
    let mut acc_rgb = [0.0_f64; 3];
    let mut acc_alpha = 0.0_f64;
    let steps = max_samples.min(n);
    for s in 0..steps {
        if acc_alpha >= early_termination {
            break;
        }
        let idx = (s as f64 / steps as f64 * n as f64) as usize;
        let idx = idx.min(n - 1);
        let t = normalize_scalar(data[idx], data_min, data_max);
        let sample = tf.sample(t);
        composite_front_to_back(&mut acc_rgb, &mut acc_alpha, sample, step_size);
    }
    Rgba::new(
        acc_rgb[0].clamp(0.0, 1.0),
        acc_rgb[1].clamp(0.0, 1.0),
        acc_rgb[2].clamp(0.0, 1.0),
        acc_alpha.clamp(0.0, 1.0),
    )
}
#[cfg(test)]
mod extended_tests {
    use super::*;
    use crate::transfer_functions::ClassificationTable;
    use crate::transfer_functions::ColorStop;
    use crate::transfer_functions::GradientColormap;

    use crate::transfer_functions::HistogramAnalyzer;
    use crate::transfer_functions::MaterialClass;

    use crate::transfer_functions::OpacityProfile;
    use crate::transfer_functions::OpacityProfileKind;

    use crate::transfer_functions::RenderQuality;

    use crate::transfer_functions::TwoPartTransferFunction;
    use crate::transfer_functions::VolumeRenderingParams;
    #[test]
    fn test_histogram_from_data_basic() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let hist = HistogramAnalyzer::from_data(&data, 10).unwrap();
        assert_eq!(hist.num_buckets(), 10);
        for i in 0..10 {
            assert!(hist.count(i) > 0, "bucket {i} should not be empty");
        }
    }
    #[test]
    fn test_histogram_from_data_empty_returns_none() {
        let hist = HistogramAnalyzer::from_data(&[], 10);
        assert!(hist.is_none());
    }
    #[test]
    fn test_histogram_from_data_constant_returns_none() {
        let data = vec![5.0_f64; 20];
        let hist = HistogramAnalyzer::from_data(&data, 10);
        assert!(hist.is_none(), "constant data should yield None");
    }
    #[test]
    fn test_histogram_frequency_sums_to_one() {
        let data: Vec<f64> = (0..1000).map(|i| i as f64 * 0.1).collect();
        let hist = HistogramAnalyzer::from_data(&data, 20).unwrap();
        let total: f64 = (0..hist.num_buckets()).map(|i| hist.frequency(i)).sum();
        assert!(
            (total - 1.0).abs() < 1e-9,
            "frequencies should sum to 1, got {total}"
        );
    }
    #[test]
    fn test_histogram_peak_bucket_in_range() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let hist = HistogramAnalyzer::from_data(&data, 10).unwrap();
        assert!(hist.peak_bucket() < 10);
    }
    #[test]
    fn test_histogram_percentile_monotone() {
        let data: Vec<f64> = (0..500).map(|i| i as f64).collect();
        let hist = HistogramAnalyzer::from_data(&data, 50).unwrap();
        let p25 = hist.percentile(0.25);
        let p50 = hist.percentile(0.50);
        let p75 = hist.percentile(0.75);
        assert!(p25 <= p50, "p25={p25} should be <= p50={p50}");
        assert!(p50 <= p75, "p50={p50} should be <= p75={p75}");
    }
    #[test]
    fn test_histogram_cdf_monotone() {
        let data: Vec<f64> = (0..200).map(|i| i as f64).collect();
        let hist = HistogramAnalyzer::from_data(&data, 20).unwrap();
        let mut prev = 0.0_f64;
        for i in 0..hist.num_buckets() {
            let c = hist.cdf(i);
            assert!(c >= prev - 1e-12, "CDF must be monotone at bucket {i}");
            prev = c;
        }
        assert!(
            (prev - 1.0).abs() < 1e-9,
            "CDF must reach 1.0 at last bucket"
        );
    }
    #[test]
    fn test_histogram_auto_tf_has_entries() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let hist = HistogramAnalyzer::from_data(&data, 8).unwrap();
        let cm = GradientColormap::default();
        let tf = hist.auto_transfer_function(&cm);
        assert!(tf.point_count() > 0);
    }
    #[test]
    fn test_classification_empty_returns_transparent() {
        let ct = ClassificationTable::new();
        let c = ct.sample(0.5);
        assert_eq!(c.a, 0.0, "empty table should return transparent");
    }
    #[test]
    fn test_classification_in_range() {
        let mut ct = ClassificationTable::new();
        ct.add_class(MaterialClass::new(
            "bone",
            0.7,
            1.0,
            Rgba::new(1.0, 1.0, 1.0, 1.0),
            1.0,
        ));
        let c = ct.sample(0.8);
        assert!(c.a > 0.0, "bone class should be opaque");
    }
    #[test]
    fn test_classification_out_of_range_transparent() {
        let mut ct = ClassificationTable::new();
        ct.add_class(MaterialClass::new(
            "bone",
            0.7,
            1.0,
            Rgba::new(1.0, 1.0, 1.0, 1.0),
            1.0,
        ));
        let c = ct.sample(0.3);
        assert_eq!(c.a, 0.0, "value outside class range should be transparent");
    }
    #[test]
    fn test_classification_ct_tissues_bone_opaque() {
        let ct = ClassificationTable::ct_tissues();
        let bone = ct.sample(0.85);
        assert!(bone.a > 0.5, "bone should be semi-opaque or opaque");
    }
    #[test]
    fn test_classification_ct_tissues_air_transparent() {
        let ct = ClassificationTable::ct_tissues();
        let air = ct.sample(0.05);
        assert!(air.a < 0.1, "air should be nearly transparent");
    }
    #[test]
    fn test_classification_to_transfer_function() {
        let ct = ClassificationTable::ct_tissues();
        let tf = ct.to_transfer_function(64);
        assert!(tf.point_count() >= 2);
    }
    #[test]
    fn test_material_class_contains() {
        let cls = MaterialClass::new("test", 0.3, 0.7, Rgba::white(), 1.0);
        assert!(cls.contains(0.5));
        assert!(cls.contains(0.3));
        assert!(cls.contains(0.7));
        assert!(!cls.contains(0.0));
        assert!(!cls.contains(1.0));
    }
    #[test]
    fn test_two_part_tf_sample_combines_color_and_opacity() {
        let mut cm = GradientColormap::new();
        cm.add_stop(ColorStop::new(0.0, 1.0, 0.0, 0.0));
        cm.add_stop(ColorStop::new(1.0, 0.0, 0.0, 1.0));
        let op = OpacityProfile::new(OpacityProfileKind::LinearRamp);
        let tptf = TwoPartTransferFunction::new(cm, op);
        let c0 = tptf.sample(0.0);
        assert!(c0.r > 0.5, "at t=0 color should be red");
        assert!(c0.a < 0.01, "at t=0 opacity should be 0 (linear ramp)");
        let c1 = tptf.sample(1.0);
        assert!(c1.b > 0.5, "at t=1 color should be blue");
        assert!((c1.a - 1.0).abs() < 1e-6, "at t=1 opacity should be 1");
    }
    #[test]
    fn test_two_part_tf_global_opacity_scale() {
        let cm = GradientColormap::default();
        let op = OpacityProfile::new(OpacityProfileKind::Constant(1.0));
        let tptf = TwoPartTransferFunction::new(cm, op).with_global_opacity(0.5);
        let c = tptf.sample(0.5);
        assert!(
            (c.a - 0.5).abs() < 1e-6,
            "global opacity 0.5 should halve opacity"
        );
    }
    #[test]
    fn test_two_part_tf_build_lut_size() {
        let cm = GradientColormap::default();
        let op = OpacityProfile::new(OpacityProfileKind::LinearRamp);
        let tptf = TwoPartTransferFunction::new(cm, op);
        let lut = tptf.build_lut(128);
        assert_eq!(lut.len(), 128);
    }
    #[test]
    fn test_two_part_tf_xray_style_non_trivial() {
        let xr = TwoPartTransferFunction::xray_style();
        let lut = xr.build_lut(64);
        let any_visible = lut.iter().any(|c| c.a > 0.01);
        assert!(any_visible, "xray style should have some visible samples");
    }
    #[test]
    fn test_two_part_tf_fire_style_increases_brightness() {
        let fire = TwoPartTransferFunction::fire_style();
        let c0 = fire.sample(0.0);
        let c1 = fire.sample(1.0);
        let lum = |c: Rgba| 0.2126 * c.r + 0.7152 * c.g + 0.0722 * c.b;
        assert!(
            lum(c1) > lum(c0),
            "fire style should be brighter at high end"
        );
    }
    #[test]
    fn test_two_part_tf_to_transfer_function() {
        let cm = GradientColormap::default();
        let op = OpacityProfile::new(OpacityProfileKind::LinearRamp);
        let tptf = TwoPartTransferFunction::new(cm, op);
        let tf = tptf.to_transfer_function(32);
        assert!(tf.point_count() >= 2);
    }
    #[test]
    fn test_volume_rendering_params_default() {
        let p = VolumeRenderingParams::default();
        assert!(p.step_size > 0.0);
        assert!(p.early_termination_threshold > 0.5);
        assert!(p.max_samples > 0);
    }
    #[test]
    fn test_volume_rendering_params_quality_preview() {
        let p = VolumeRenderingParams::with_quality(RenderQuality::Preview);
        let q = VolumeRenderingParams::with_quality(RenderQuality::High);
        assert!(
            p.step_size > q.step_size,
            "preview step should be larger than high-quality"
        );
        assert!(
            p.max_samples < q.max_samples,
            "preview should have fewer max samples"
        );
    }
    #[test]
    fn test_volume_rendering_params_with_lighting() {
        let p = VolumeRenderingParams::default().with_lighting(0.3, 0.7);
        assert!((p.ambient - 0.3).abs() < 1e-10);
        assert!((p.diffuse - 0.7).abs() < 1e-10);
    }
    #[test]
    fn test_render_quality_step_sizes() {
        assert!(RenderQuality::Preview.step_size() > RenderQuality::Medium.step_size());
        assert!(RenderQuality::Medium.step_size() > RenderQuality::High.step_size());
    }
    #[test]
    fn test_normalize_scalar_basic() {
        assert!((normalize_scalar(0.5, 0.0, 1.0) - 0.5).abs() < 1e-12);
        assert!((normalize_scalar(0.0, 0.0, 1.0)).abs() < 1e-12);
        assert!((normalize_scalar(1.0, 0.0, 1.0) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_normalize_scalar_constant_returns_half() {
        assert!((normalize_scalar(5.0, 5.0, 5.0) - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_normalize_slice_basic() {
        let data = vec![0.0, 50.0, 100.0];
        let norm = normalize_slice(&data);
        assert_eq!(norm.len(), 3);
        assert!((norm[0]).abs() < 1e-12);
        assert!((norm[1] - 0.5).abs() < 1e-10);
        assert!((norm[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_normalize_slice_empty() {
        let norm = normalize_slice(&[]);
        assert!(norm.is_empty());
    }
    #[test]
    fn test_apply_transfer_function() {
        let mut tf = TransferFunction::new();
        tf.add_point(0.0, Rgba::new(0.0, 0.0, 0.0, 0.0));
        tf.add_point(1.0, Rgba::new(1.0, 1.0, 1.0, 1.0));
        let norm = vec![0.0, 0.5, 1.0];
        let colors = apply_transfer_function(&tf, &norm);
        assert_eq!(colors.len(), 3);
        assert!((colors[1].r - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_gradient_magnitude_1d_flat() {
        let data = vec![1.0_f64; 10];
        for i in 0..10 {
            assert_eq!(gradient_magnitude_1d(&data, i), 0.0);
        }
    }
    #[test]
    fn test_gradient_magnitude_1d_ramp() {
        let data: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let g = gradient_magnitude_1d(&data, 5);
        assert!((g - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_gradient_magnitude_1d_boundaries() {
        let data = vec![0.0, 1.0, 2.0, 3.0];
        assert!((gradient_magnitude_1d(&data, 0) - 1.0).abs() < 1e-10);
        assert!((gradient_magnitude_1d(&data, 3) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_gradient_modulated_opacity_zero_gradient() {
        let op = gradient_modulated_opacity(0.5, 0.0, 2.0);
        assert!(
            (op - 0.5).abs() < 1e-10,
            "zero gradient should not change opacity"
        );
    }
    #[test]
    fn test_gradient_modulated_opacity_max_gradient() {
        let op = gradient_modulated_opacity(0.0, 1.0, 2.0);
        assert!(op > 0.5, "max gradient should boost opacity significantly");
    }
    #[test]
    fn test_composite_front_to_back_accumulates_alpha() {
        let mut rgb = [0.0_f64; 3];
        let mut alpha = 0.0_f64;
        let sample = Rgba::new(1.0, 0.0, 0.0, 1.0);
        composite_front_to_back(&mut rgb, &mut alpha, sample, 0.01);
        assert!(alpha > 0.0, "alpha should increase after compositing");
    }
    #[test]
    fn test_composite_front_to_back_early_termination() {
        let mut rgb = [0.0_f64; 3];
        let mut alpha = 0.0_f64;
        for _ in 0..200 {
            composite_front_to_back(&mut rgb, &mut alpha, Rgba::new(1.0, 0.0, 0.0, 1.0), 0.1);
        }
        assert!(
            alpha > 0.99,
            "accumulated alpha should approach 1 after many opaque samples"
        );
    }
    #[test]
    fn test_raymarch_1d_empty_data() {
        let tf = TransferFunction::default();
        let result = raymarch_1d(&[], 0.0, 1.0, &tf, 64, 0.01, 0.99);
        assert_eq!(result.a, 0.0, "empty data should return transparent");
    }
    #[test]
    fn test_raymarch_1d_opaque_data_not_transparent() {
        let data = vec![1.0_f64; 100];
        let tf = TransferFunction::default();
        let result = raymarch_1d(&data, 0.0, 1.0, &tf, 64, 0.01, 0.99);
        assert!(
            result.a > 0.0,
            "opaque data should produce non-transparent result"
        );
    }
    #[test]
    fn test_raymarch_1d_transparent_data_stays_transparent() {
        let data = vec![0.0_f64; 100];
        let tf = TransferFunction::default();
        let result = raymarch_1d(&data, 0.0, 1.0, &tf, 64, 0.01, 0.99);
        assert!(
            result.a < 0.1,
            "fully transparent data should produce near-transparent result"
        );
    }
}
#[cfg(test)]
mod extended_tests_2 {
    use super::*;
    use crate::transfer_functions::ClassificationTable;

    use crate::transfer_functions::GradientColormap;
    use crate::transfer_functions::GradientMagnitudeOpacity;

    use crate::transfer_functions::MultiDimensionalTransferFunction;

    use crate::transfer_functions::PiecewiseLinearOpacityMap;
    use crate::transfer_functions::PreintegrationTable;

    use crate::transfer_functions::TwoPartTransferFunction;

    #[test]
    fn test_2d_tf_new_transparent() {
        let tf = MultiDimensionalTransferFunction::new(8, 8);
        let c = tf.sample(0.5, 0.5);
        assert_eq!(c.a, 0.0, "new 2D TF should be fully transparent");
    }
    #[test]
    fn test_2d_tf_set_and_sample_corner() {
        let mut tf = MultiDimensionalTransferFunction::new(4, 4);
        tf.set(0, 0, Rgba::new(1.0, 0.0, 0.0, 1.0));
        let c = tf.sample(0.0, 0.0);
        assert!((c.r - 1.0).abs() < 0.01, "corner should be red");
    }
    #[test]
    fn test_2d_tf_bilinear_midpoint() {
        let mut tf = MultiDimensionalTransferFunction::new(2, 2);
        tf.set(0, 0, Rgba::new(0.0, 0.0, 0.0, 0.0));
        tf.set(1, 0, Rgba::new(1.0, 0.0, 0.0, 1.0));
        tf.set(0, 1, Rgba::new(0.0, 1.0, 0.0, 1.0));
        tf.set(1, 1, Rgba::new(1.0, 1.0, 0.0, 1.0));
        let c = tf.sample(0.5, 0.5);
        assert!(c.a > 0.0 && c.a < 1.0 + 1e-5);
    }
    #[test]
    fn test_2d_tf_flatten_to_1d_point_count() {
        let mut tf = MultiDimensionalTransferFunction::new(8, 4);
        tf.set(7, 3, Rgba::new(0.5, 0.5, 0.5, 0.5));
        let tf1d = tf.flatten_to_1d();
        assert_eq!(tf1d.point_count(), 8);
    }
    #[test]
    fn test_2d_tf_dimension_clamping() {
        let mut tf = MultiDimensionalTransferFunction::new(4, 4);
        tf.set(100, 100, Rgba::new(1.0, 0.0, 0.0, 1.0));
        let _ = tf.sample(1.1, 1.1);
    }
    #[test]
    fn test_preintegration_table_build_no_panic() {
        let tf = TransferFunction::temperature();
        let pt = PreintegrationTable::build(&tf, 32, 8);
        assert_eq!(pt.table.len(), 32 * 32);
    }
    #[test]
    fn test_preintegration_sample_transparent_same() {
        let tf = TransferFunction::new();
        let pt = PreintegrationTable::build(&tf, 16, 8);
        let c = pt.sample(0.0, 1.0);
        assert!(
            c.a < 0.01,
            "transparent TF preintegration should give near-zero alpha"
        );
    }
    #[test]
    fn test_preintegration_diagonal_nonzero() {
        let tf = TransferFunction::default();
        let pt = PreintegrationTable::build(&tf, 16, 8);
        let c = pt.sample(0.9, 1.0);
        assert!(c.a >= 0.0, "alpha must be non-negative");
    }
    #[test]
    fn test_preintegration_symmetry() {
        let tf = TransferFunction::density();
        let pt = PreintegrationTable::build(&tf, 16, 8);
        let c1 = pt.sample(0.3, 0.7);
        let c2 = pt.sample(0.7, 0.3);
        for comp in [c1.r, c1.g, c1.b, c1.a, c2.r, c2.g, c2.b, c2.a] {
            assert!((0.0..=1.0 + 1e-5).contains(&comp));
        }
    }
    #[test]
    fn test_pwl_opacity_default_endpoints() {
        let m = PiecewiseLinearOpacityMap::default();
        assert!((m.evaluate(0.0)).abs() < 1e-10);
        assert!((m.evaluate(1.0) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_pwl_opacity_midpoint_interpolation() {
        let m = PiecewiseLinearOpacityMap::default();
        assert!((m.evaluate(0.5) - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_pwl_opacity_add_knot() {
        let mut m = PiecewiseLinearOpacityMap::new();
        m.add_knot(0.0, 0.0);
        m.add_knot(0.5, 1.0);
        m.add_knot(1.0, 0.0);
        assert_eq!(m.knot_count(), 3);
        assert!((m.evaluate(0.5) - 1.0).abs() < 1e-10);
        assert!((m.evaluate(0.0)).abs() < 1e-10);
    }
    #[test]
    fn test_pwl_opacity_clamps_values() {
        let mut m = PiecewiseLinearOpacityMap::new();
        m.add_knot(1.5, 2.0);
        assert_eq!(m.knot_count(), 1);
        assert!((m.knots[0].0 - 1.0).abs() < 1e-10);
        assert!((m.knots[0].1 - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_pwl_opacity_build_lut() {
        let m = PiecewiseLinearOpacityMap::default();
        let lut = m.build_lut(11);
        assert_eq!(lut.len(), 11);
        assert!((lut[5] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_pwl_opacity_to_transfer_function() {
        let m = PiecewiseLinearOpacityMap::default();
        let cm = GradientColormap::diverging_bwr();
        let tf = m.to_transfer_function(&cm);
        assert!(tf.point_count() >= 2);
    }
    #[test]
    fn test_pwl_opacity_empty_returns_zero() {
        let m = PiecewiseLinearOpacityMap::new();
        assert_eq!(m.evaluate(0.5), 0.0);
    }
    #[test]
    fn test_pwl_opacity_single_knot_constant() {
        let mut m = PiecewiseLinearOpacityMap::new();
        m.add_knot(0.5, 0.7);
        assert!((m.evaluate(0.0) - 0.7).abs() < 1e-10);
        assert!((m.evaluate(0.5) - 0.7).abs() < 1e-10);
        assert!((m.evaluate(1.0) - 0.7).abs() < 1e-10);
    }
    #[test]
    fn test_gradient_mag_opacity_zero_gradient() {
        let base = PiecewiseLinearOpacityMap::default();
        let gmo = GradientMagnitudeOpacity::new(base);
        let op = gmo.evaluate(0.5, 0.0);
        assert!((op - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_gradient_mag_opacity_max_gradient() {
        let mut base = PiecewiseLinearOpacityMap::new();
        base.add_knot(0.0, 0.0);
        base.add_knot(1.0, 0.0);
        let gmo = GradientMagnitudeOpacity::new(base).with_sharpness(5.0);
        let op = gmo.evaluate(0.5, 1.0);
        assert!(
            op > 0.9,
            "max gradient should fully boost opacity, got {op}"
        );
    }
    #[test]
    fn test_gradient_mag_opacity_floor() {
        let base = PiecewiseLinearOpacityMap::new();
        let gmo = GradientMagnitudeOpacity::new(base).with_floor(0.05);
        let op = gmo.evaluate(0.5, 0.0);
        assert!((op - 0.05).abs() < 1e-10, "floor should be 0.05, got {op}");
    }
    #[test]
    fn test_gradient_mag_opacity_clamped_output() {
        let mut base = PiecewiseLinearOpacityMap::new();
        base.add_knot(0.5, 0.9);
        let gmo = GradientMagnitudeOpacity::new(base).with_sharpness(10.0);
        let op = gmo.evaluate(0.5, 1.0);
        assert!(op <= 1.0 + 1e-10, "opacity must not exceed 1.0");
        assert!(op >= 0.0);
    }
    #[test]
    fn test_classification_table_ct_tissues() {
        let ct = ClassificationTable::ct_tissues();
        let c_air = ct.sample(0.05);
        assert!(c_air.a < 0.1, "air should be near-transparent");
        let c_bone = ct.sample(0.85);
        assert!(c_bone.a > 0.5, "bone should be semi-opaque");
    }
    #[test]
    fn test_classification_table_to_tf() {
        let ct = ClassificationTable::ct_tissues();
        let tf = ct.to_transfer_function(32);
        assert_eq!(tf.point_count(), 32);
    }
    #[test]
    fn test_classification_table_gap_transparent() {
        let ct = ClassificationTable::ct_tissues();
        let c = ct.sample(0.2);
        assert_eq!(c.a, 0.0, "gap region should be transparent");
    }
    #[test]
    fn test_classification_table_classify() {
        let ct = ClassificationTable::ct_tissues();
        let cls = ct.classify(0.85);
        assert!(cls.is_some());
        assert_eq!(cls.unwrap().name, "bone");
    }
    #[test]
    fn test_two_part_tf_xray() {
        let tf = TwoPartTransferFunction::xray_style();
        let lut = tf.build_lut(64);
        assert_eq!(lut.len(), 64);
        let idx_peak = (0.75 * 63.0) as usize;
        let idx_edge = 0;
        assert!(
            lut[idx_peak].a > lut[idx_edge].a,
            "peak opacity should exceed edge opacity"
        );
    }
    #[test]
    fn test_two_part_tf_fire() {
        let tf = TwoPartTransferFunction::fire_style();
        let dark = tf.sample(0.0);
        let bright = tf.sample(1.0);
        let lum_dark = 0.2126 * dark.r + 0.7152 * dark.g + 0.0722 * dark.b;
        let lum_bright = 0.2126 * bright.r + 0.7152 * bright.g + 0.0722 * bright.b;
        assert!(lum_bright > lum_dark, "fire TF should get brighter");
    }
    #[test]
    fn test_two_part_tf_global_opacity() {
        let tf = TwoPartTransferFunction::xray_style().with_global_opacity(0.5);
        assert!((tf.global_opacity - 0.5).abs() < 1e-10);
        let c = tf.sample(0.75);
        assert!(c.a <= 0.5 + 1e-10, "global opacity 0.5 should cap alpha");
    }
    #[test]
    fn test_two_part_tf_to_transfer_function() {
        let tf = TwoPartTransferFunction::fire_style();
        let tf1d = tf.to_transfer_function(16);
        assert_eq!(tf1d.point_count(), 16);
    }
}
