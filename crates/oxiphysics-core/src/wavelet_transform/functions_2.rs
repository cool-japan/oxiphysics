//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[cfg(test)]
mod tests {
    use super::*;
    use crate::wavelet_transform::DaubechiesWavelet;
    use crate::wavelet_transform::HaarTransform;
    use crate::wavelet_transform::LiftingHaar;
    use crate::wavelet_transform::MotherWavelet;
    use crate::wavelet_transform::ThresholdMode;
    use crate::wavelet_transform::WaveletFamily;
    use std::f64::consts::PI;
    pub(super) const TOL: f64 = 1e-6;
    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }
    #[test]
    fn test_haar_forward_inverse() {
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let (a, d) = HaarTransform::forward(&signal);
        let reconstructed = HaarTransform::inverse(&a, &d);
        for (i, (&orig, &rec)) in signal.iter().zip(reconstructed.iter()).enumerate() {
            assert!(
                approx_eq(orig, rec, TOL),
                "Haar mismatch at {}: {} vs {}",
                i,
                orig,
                rec
            );
        }
    }
    #[test]
    fn test_haar_constant_signal() {
        let signal = vec![5.0; 8];
        let (a, d) = HaarTransform::forward(&signal);
        for &di in &d {
            assert!(di.abs() < TOL, "Haar detail should be zero for constant");
        }
        assert_eq!(a.len(), 4);
    }
    #[test]
    fn test_haar_multilevel() {
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let result = HaarTransform::forward_multilevel(&signal, 3);
        assert_eq!(result.len(), 4);
        assert_eq!(result[0].len(), 1);
    }
    #[test]
    fn test_dwt_haar_single_level() {
        let signal = vec![1.0, 3.0, 5.0, 7.0];
        let level = dwt_single(&signal, WaveletFamily::Haar);
        assert_eq!(level.approx.len(), 2);
        assert_eq!(level.detail.len(), 2);
    }
    #[test]
    fn test_dwt_multilevel_haar() {
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let decomp = dwt(&signal, WaveletFamily::Haar, 3);
        assert_eq!(decomp.details.len(), 3);
    }
    #[test]
    fn test_dwt_idwt_haar_roundtrip() {
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let decomp = dwt(&signal, WaveletFamily::Haar, 2);
        let rec = idwt(&decomp);
        assert_eq!(rec.len(), signal.len());
        for (i, (&o, &r)) in signal.iter().zip(rec.iter()).enumerate() {
            assert!(
                approx_eq(o, r, 1e-4),
                "Haar roundtrip mismatch at {}: {} vs {}",
                i,
                o,
                r
            );
        }
    }
    #[test]
    fn test_daubechies_db2_forward() {
        let db = DaubechiesWavelet::new(2);
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let level = db.forward(&signal);
        assert_eq!(level.approx.len(), 5);
        assert_eq!(level.detail.len(), 5);
    }
    #[test]
    fn test_daubechies_db3_forward() {
        let db = DaubechiesWavelet::new(3);
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let level = db.forward(&signal);
        assert!(!level.approx.is_empty());
        assert!(!level.detail.is_empty());
    }
    #[test]
    fn test_daubechies_db4_forward() {
        let db = DaubechiesWavelet::new(4);
        let signal: Vec<f64> = (0..16).map(|i| (i as f64).sin()).collect();
        let level = db.forward(&signal);
        assert!(!level.approx.is_empty());
    }
    #[test]
    fn test_daubechies_db5_forward() {
        let db = DaubechiesWavelet::new(5);
        let signal: Vec<f64> = (0..16).map(|i| (i as f64 * 0.5).cos()).collect();
        let level = db.forward(&signal);
        assert!(!level.approx.is_empty());
    }
    #[test]
    fn test_daubechies_db6_forward() {
        let db = DaubechiesWavelet::new(6);
        let signal: Vec<f64> = (0..32).map(|i| (i as f64 * 0.3).sin()).collect();
        let level = db.forward(&signal);
        assert!(!level.approx.is_empty());
    }
    #[test]
    fn test_multiresolution_analysis() {
        let signal: Vec<f64> = (0..32).map(|i| (i as f64 * PI / 8.0).sin()).collect();
        let mra = multiresolution_analysis(&signal, WaveletFamily::Haar, 3);
        assert_eq!(mra.approximations.len(), 4);
        assert_eq!(mra.detail_contributions.len(), 3);
    }
    #[test]
    fn test_wavelet_packet_decompose() {
        let signal: Vec<f64> = (0..16).map(|i| (i as f64).sin()).collect();
        let tree = wavelet_packet_decompose(&signal, WaveletFamily::Haar, 2);
        assert!(tree.nodes.len() >= 2);
        assert_eq!(tree.nodes[0].len(), 1);
    }
    #[test]
    fn test_best_basis_selection() {
        let signal: Vec<f64> = (0..16).map(|i| (i as f64 * 0.5).sin()).collect();
        let tree = wavelet_packet_decompose(&signal, WaveletFamily::Haar, 2);
        let basis = best_basis_selection(&tree);
        assert!(!basis.is_empty());
    }
    #[test]
    fn test_cwt_morlet() {
        let signal: Vec<f64> = (0..64)
            .map(|i| (2.0 * PI * i as f64 / 16.0).sin())
            .collect();
        let scales = log_scales(1.0, 8, 0.5);
        let result = cwt(&signal, &scales, MotherWavelet::morlet(), 1.0);
        assert_eq!(result.coefficients.len(), scales.len());
        assert_eq!(result.coefficients[0].len(), signal.len());
    }
    #[test]
    fn test_cwt_mexican_hat() {
        let signal: Vec<f64> = (0..32).map(|i| (2.0 * PI * i as f64 / 8.0).sin()).collect();
        let scales = linear_scales(1.0, 5.0, 5);
        let result = cwt(&signal, &scales, MotherWavelet::mexican_hat(), 1.0);
        assert_eq!(result.coefficients.len(), 5);
    }
    #[test]
    fn test_morlet_evaluation() {
        let morlet = MotherWavelet::morlet();
        let v = morlet.evaluate(0.0);
        assert!(v > 0.9, "Morlet at t=0 should be close to 1");
    }
    #[test]
    fn test_mexican_hat_evaluation() {
        let mh = MotherWavelet::mexican_hat();
        let v0 = mh.evaluate(0.0);
        assert!(v0 > 0.0, "Mexican hat at t=0 should be positive");
        let v5 = mh.evaluate(5.0);
        assert!(v5.abs() < 0.01, "Mexican hat should decay at large t");
    }
    #[test]
    fn test_scalogram() {
        let signal: Vec<f64> = (0..32).map(|i| (i as f64 * 0.5).sin()).collect();
        let scales = log_scales(1.0, 4, 0.5);
        let cwt_result = cwt(&signal, &scales, MotherWavelet::morlet(), 1.0);
        let scalo = scalogram(&cwt_result);
        assert_eq!(scalo.energy.len(), 4);
        assert!(scalo.scale_energy.iter().all(|&e| e >= 0.0));
    }
    #[test]
    fn test_global_wavelet_spectrum() {
        let signal: Vec<f64> = (0..32).map(|i| (i as f64 * 0.3).sin()).collect();
        let scales = log_scales(1.0, 6, 0.5);
        let cwt_result = cwt(&signal, &scales, MotherWavelet::morlet(), 1.0);
        let scalo = scalogram(&cwt_result);
        let gws = global_wavelet_spectrum(&scalo);
        assert_eq!(gws.len(), 6);
        assert!(gws.iter().all(|&e| e >= 0.0));
    }
    #[test]
    fn test_hard_thresholding() {
        let coeffs = vec![0.1, 0.5, -0.3, 1.0, -0.05, 2.0];
        let result = apply_threshold(&coeffs, 0.4, ThresholdMode::Hard);
        assert_eq!(result[0], 0.0);
        assert_eq!(result[1], 0.5);
        assert_eq!(result[4], 0.0);
    }
    #[test]
    fn test_soft_thresholding() {
        let coeffs = vec![1.0, -2.0, 0.3];
        let result = apply_threshold(&coeffs, 0.5, ThresholdMode::Soft);
        assert!(approx_eq(result[0], 0.5, TOL));
        assert!(approx_eq(result[1], -1.5, TOL));
        assert_eq!(result[2], 0.0);
    }
    #[test]
    fn test_universal_threshold() {
        let thresh = universal_threshold(1.0, 100);
        assert!(thresh > 0.0);
        assert!(thresh < 5.0);
    }
    #[test]
    fn test_estimate_noise_sigma() {
        let noise: Vec<f64> = (0..100).map(|i| (i as f64 * 0.7).sin() * 0.1).collect();
        let sigma = estimate_noise_sigma(&noise);
        assert!(sigma > 0.0);
    }
    #[test]
    fn test_wavelet_denoise_soft() {
        let clean: Vec<f64> = (0..64)
            .map(|i| (2.0 * PI * i as f64 / 16.0).sin())
            .collect();
        let noisy: Vec<f64> = clean.iter().map(|&x| x + 0.1).collect();
        let denoised = wavelet_denoise(&noisy, WaveletFamily::Haar, 3, ThresholdMode::Soft, None);
        assert_eq!(denoised.len(), noisy.len());
    }
    #[test]
    fn test_wavelet_denoise_hard() {
        let signal: Vec<f64> = (0..32).map(|i| (i as f64 * 0.4).cos()).collect();
        let denoised = wavelet_denoise(
            &signal,
            WaveletFamily::Db2,
            2,
            ThresholdMode::Hard,
            Some(0.1),
        );
        assert_eq!(denoised.len(), signal.len());
    }
    #[test]
    fn test_level_energy() {
        let coeffs = vec![1.0, 2.0, 3.0];
        assert!(approx_eq(level_energy(&coeffs), 14.0, TOL));
    }
    #[test]
    fn test_energy_distribution() {
        let signal: Vec<f64> = (0..32).map(|i| (i as f64 * 0.2).sin()).collect();
        let decomp = dwt(&signal, WaveletFamily::Haar, 3);
        let (detail_e, approx_e) = energy_distribution(&decomp);
        assert_eq!(detail_e.len(), 3);
        assert!(approx_e >= 0.0);
    }
    #[test]
    fn test_relative_energy() {
        let signal: Vec<f64> = (0..16).map(|i| (i as f64 * 0.5).sin()).collect();
        let decomp = dwt(&signal, WaveletFamily::Haar, 2);
        let rel = relative_energy(&decomp);
        let total: f64 = rel.iter().sum();
        assert!(
            approx_eq(total, 1.0, 0.01),
            "Relative energy should sum to ~1.0"
        );
    }
    #[test]
    fn test_wavelet_entropy() {
        let signal: Vec<f64> = (0..32).map(|i| (i as f64).sin()).collect();
        let decomp = dwt(&signal, WaveletFamily::Haar, 3);
        let ent = wavelet_entropy(&decomp);
        assert!(ent >= 0.0);
    }
    #[test]
    fn test_swt() {
        let signal: Vec<f64> = (0..16).map(|i| (i as f64 * 0.3).cos()).collect();
        let result = swt(&signal, WaveletFamily::Haar, 2);
        assert_eq!(result.details.len(), 2);
        assert_eq!(result.details[0].len(), 16);
        assert_eq!(result.approx.len(), 16);
    }
    #[test]
    fn test_wavelet_cross_spectrum() {
        let x: Vec<f64> = (0..32).map(|i| (i as f64 * 0.3).sin()).collect();
        let y: Vec<f64> = (0..32).map(|i| (i as f64 * 0.3).cos()).collect();
        let scales = log_scales(1.0, 3, 1.0);
        let cwt_x = cwt(&x, &scales, MotherWavelet::morlet(), 1.0);
        let cwt_y = cwt(&y, &scales, MotherWavelet::morlet(), 1.0);
        let cross = wavelet_cross_spectrum(&cwt_x, &cwt_y);
        assert_eq!(cross.len(), 3);
    }
    #[test]
    fn test_wavelet_coherence() {
        let x: Vec<f64> = (0..16).map(|i| (i as f64 * 0.5).sin()).collect();
        let y = x.clone();
        let scales = log_scales(1.0, 3, 1.0);
        let cwt_x = cwt(&x, &scales, MotherWavelet::morlet(), 1.0);
        let cwt_y = cwt(&y, &scales, MotherWavelet::morlet(), 1.0);
        let coh = wavelet_coherence(&cwt_x, &cwt_y, 2);
        assert_eq!(coh.len(), 3);
        for row in &coh {
            for &c in row {
                assert!((0.0..=1.0 + TOL).contains(&c));
            }
        }
    }
    #[test]
    fn test_log_scales() {
        let scales = log_scales(1.0, 5, 0.5);
        assert_eq!(scales.len(), 5);
        assert!(approx_eq(scales[0], 1.0, TOL));
        for i in 1..scales.len() {
            assert!(scales[i] > scales[i - 1]);
        }
    }
    #[test]
    fn test_linear_scales() {
        let scales = linear_scales(1.0, 10.0, 10);
        assert_eq!(scales.len(), 10);
        assert!(approx_eq(scales[0], 1.0, TOL));
        assert!(approx_eq(scales[9], 10.0, TOL));
    }
    #[test]
    fn test_detect_ridges() {
        let signal: Vec<f64> = (0..64)
            .map(|i| (2.0 * PI * i as f64 / 16.0).sin())
            .collect();
        let scales = log_scales(1.0, 10, 0.25);
        let cwt_result = cwt(&signal, &scales, MotherWavelet::morlet(), 1.0);
        let scalo = scalogram(&cwt_result);
        let ridges = detect_ridges(&scalo);
        assert!(
            ridges
                .iter()
                .all(|&(s, t)| s < scales.len() && t < signal.len())
        );
    }
    #[test]
    fn test_wavelet_compress() {
        let signal: Vec<f64> = (0..32).map(|i| (i as f64 * 0.3).sin()).collect();
        let compressed = wavelet_compress(&signal, WaveletFamily::Haar, 3, 0.5);
        assert_eq!(compressed.len(), signal.len());
    }
    #[test]
    fn test_compression_ratio() {
        let signal: Vec<f64> = (0..32).map(|i| (i as f64 * 0.2).sin()).collect();
        let decomp = dwt(&signal, WaveletFamily::Haar, 3);
        let ratio = compression_ratio(&decomp);
        assert!((0.0..=1.0).contains(&ratio));
    }
    #[test]
    fn test_modwt() {
        let signal: Vec<f64> = (0..32).map(|i| (i as f64 * 0.4).sin()).collect();
        let result = modwt(&signal, WaveletFamily::Haar, 3);
        assert_eq!(result.details.len(), 3);
        assert_eq!(result.details[0].len(), 32);
        assert_eq!(result.approx.len(), 32);
    }
    #[test]
    fn test_lifting_haar_roundtrip() {
        let original = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let mut data = original.clone();
        LiftingHaar::forward(&mut data);
        LiftingHaar::inverse(&mut data);
        for (i, (&o, &r)) in original.iter().zip(data.iter()).enumerate() {
            assert!(
                approx_eq(o, r, TOL),
                "Lifting roundtrip mismatch at {}: {} vs {}",
                i,
                o,
                r
            );
        }
    }
    #[test]
    fn test_scale_frequency_roundtrip() {
        let wavelet = MotherWavelet::morlet();
        let dt = 0.01;
        let freq = 10.0;
        let scale = frequency_to_scale(freq, dt, wavelet);
        let freq_back = scale_to_frequency(scale, dt, wavelet);
        assert!(
            approx_eq(freq, freq_back, 0.01),
            "Scale-freq roundtrip: {} vs {}",
            freq,
            freq_back
        );
    }
    #[test]
    fn test_reconstruction_error() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![1.0, 2.0, 3.0];
        assert!(approx_eq(reconstruction_error(&a, &b), 0.0, TOL));
    }
    #[test]
    fn test_reconstruction_snr_perfect() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![1.0, 2.0, 3.0];
        assert_eq!(reconstruction_snr(&a, &b), f64::INFINITY);
    }
    #[test]
    fn test_bayes_shrink_denoise() {
        let signal: Vec<f64> = (0..32).map(|i| (i as f64 * 0.4).sin()).collect();
        let denoised = bayes_shrink_denoise(&signal, WaveletFamily::Haar, 3, ThresholdMode::Soft);
        assert_eq!(denoised.len(), signal.len());
    }
    #[test]
    fn test_wavelet_features() {
        let signal: Vec<f64> = (0..32).map(|i| (i as f64 * 0.2).cos()).collect();
        let decomp = dwt(&signal, WaveletFamily::Haar, 3);
        let features = wavelet_features(&decomp);
        assert_eq!(features.len(), 4);
    }
    #[test]
    fn test_cone_of_influence() {
        let scales = vec![1.0, 2.0, 4.0];
        let coi = cone_of_influence(&scales, MotherWavelet::morlet());
        assert_eq!(coi.len(), 3);
        assert!(coi[0] < coi[1]);
        assert!(coi[1] < coi[2]);
    }
    #[test]
    fn test_wavelet_variance() {
        let signal: Vec<f64> = (0..32).map(|i| (i as f64 * 0.5).sin()).collect();
        let scales = log_scales(1.0, 4, 0.5);
        let cwt_result = cwt(&signal, &scales, MotherWavelet::morlet(), 1.0);
        let var = wavelet_variance(&cwt_result);
        assert_eq!(var.len(), 4);
        assert!(var.iter().all(|&v| v >= 0.0));
    }
    #[test]
    fn test_wavelet_power() {
        let signal: Vec<f64> = (0..16).map(|i| (i as f64 * 0.3).sin()).collect();
        let scales = log_scales(1.0, 3, 1.0);
        let cwt_result = cwt(&signal, &scales, MotherWavelet::morlet(), 1.0);
        let power = wavelet_power(&cwt_result);
        assert_eq!(power.len(), 3);
        for row in &power {
            assert!(row.iter().all(|&p| p >= 0.0));
        }
    }
}
