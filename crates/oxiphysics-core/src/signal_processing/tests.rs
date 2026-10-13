//! Unit tests for the signal-processing module (FFT, windows, filters, PSD, wavelets).

use super::*;

const EPS: f64 = 1e-9;

fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

// --- Complex64 ---

#[test]
fn test_complex_basic_ops() {
    let a = Complex64::new(3.0, 4.0);
    let b = Complex64::new(1.0, -2.0);
    let sum = a + b;
    assert_eq!(sum.re, 4.0);
    assert_eq!(sum.im, 2.0);
    let prod = a * b;
    // (3+4i)(1-2i) = 3 - 6i + 4i + 8 = 11 - 2i
    assert!(approx_eq(prod.re, 11.0, EPS));
    assert!(approx_eq(prod.im, -2.0, EPS));
}

#[test]
fn test_complex_abs_and_arg() {
    let c = Complex64::new(3.0, 4.0);
    assert!(approx_eq(c.abs(), 5.0, EPS));
}

#[test]
fn test_complex_conj() {
    let c = Complex64::new(1.0, -2.0);
    let cc = c.conj();
    assert_eq!(cc.re, 1.0);
    assert_eq!(cc.im, 2.0);
}

#[test]
fn test_complex_exp_i() {
    let c = Complex64::exp_i(0.0);
    assert!(approx_eq(c.re, 1.0, EPS));
    assert!(approx_eq(c.im, 0.0, EPS));
    let c90 = Complex64::exp_i(PI / 2.0);
    assert!(approx_eq(c90.re, 0.0, 1e-14));
    assert!(approx_eq(c90.im, 1.0, 1e-14));
}

// --- FFT ---

#[test]
fn test_fft_dc() {
    // DC signal: all ones -> FFT has only X[0] = N, rest 0
    let signal = vec![1.0f64; 8];
    let spec = fft(&signal);
    assert!(approx_eq(spec[0].re, 8.0, 1e-10));
    assert!(approx_eq(spec[0].im, 0.0, 1e-10));
    for s in spec.iter().skip(1) {
        assert!(approx_eq(s.re, 0.0, 1e-10));
        assert!(approx_eq(s.im, 0.0, 1e-10));
    }
}

#[test]
fn test_fft_single_tone() {
    // cos(2π*k*n/N) -> peaks at bin k and N-k
    let n = 16usize;
    let k = 2usize;
    let signal: Vec<f64> = (0..n)
        .map(|i| (2.0 * PI * k as f64 * i as f64 / n as f64).cos())
        .collect();
    let spec = fft(&signal);
    // Bin k should have magnitude ≈ N/2
    assert!(approx_eq(spec[k].abs(), n as f64 / 2.0, 1e-8));
}

#[test]
fn test_fft_ifft_roundtrip() {
    let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let spec = fft(&signal);
    let recovered = ifft_real(&spec);
    for (a, b) in signal.iter().zip(recovered.iter()) {
        assert!(approx_eq(*a, *b, 1e-10));
    }
}

#[test]
fn test_fft_parseval() {
    // Parseval's theorem: sum |x[n]|^2 = (1/N) * sum |X[k]|^2
    let signal = vec![1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 0.5, -0.5];
    let n = signal.len() as f64;
    let time_energy: f64 = signal.iter().map(|x| x * x).sum();
    let spec = fft(&signal);
    let freq_energy: f64 = spec.iter().map(|c| c.norm_sq()).sum::<f64>() / n;
    assert!(approx_eq(time_energy, freq_energy, 1e-8));
}

// --- Windows ---

#[test]
fn test_hann_window_endpoints() {
    let w = hann_window(8);
    assert!(approx_eq(w[0], 0.0, 1e-10));
}

#[test]
fn test_hamming_window_range() {
    let w = hamming_window(64);
    for v in &w {
        assert!(*v >= 0.0 && *v <= 1.0);
    }
}

#[test]
fn test_blackman_window_endpoints() {
    let w = blackman_window(16);
    assert!(w[0].abs() < 1e-10);
}

#[test]
fn test_kaiser_window_unity_center() {
    let w = kaiser_window(9, 4.0);
    // Center should be 1.0
    assert!(approx_eq(w[4], 1.0, 1e-10));
}

#[test]
fn test_bessel_i0() {
    // I0(0) = 1
    assert!(approx_eq(bessel_i0(0.0), 1.0, EPS));
    // I0(1) ≈ 1.26606
    assert!(approx_eq(bessel_i0(1.0), 1.2660658, 1e-6));
}

// --- Convolution ---

#[test]
fn test_convolve_identity() {
    let x = vec![1.0, 2.0, 3.0];
    let h = vec![1.0, 0.0, 0.0];
    let y = convolve(&x, &h);
    assert!(approx_eq(y[0], 1.0, 1e-10));
    assert!(approx_eq(y[1], 2.0, 1e-10));
    assert!(approx_eq(y[2], 3.0, 1e-10));
}

#[test]
fn test_convolve_box() {
    // Convolve [1,1,1] with [1,1] -> [1,2,2,1]
    let x = vec![1.0, 1.0, 1.0];
    let h = vec![1.0, 1.0];
    let y = convolve(&x, &h);
    assert!(approx_eq(y[0], 1.0, 1e-10));
    assert!(approx_eq(y[1], 2.0, 1e-10));
    assert!(approx_eq(y[2], 2.0, 1e-10));
    assert!(approx_eq(y[3], 1.0, 1e-10));
}

#[test]
fn test_auto_correlate_peak() {
    let signal = vec![1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0];
    let ac = auto_correlate(&signal);
    // Peak should be at lag = signal.len() - 1
    let mid = signal.len() - 1;
    let peak_idx = ac
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap();
    assert_eq!(peak_idx, mid);
}

// --- FIR filter ---

#[test]
fn test_fir_lowpass_dc_gain() {
    let fir = FirFilter::lowpass(31, 0.5);
    // DC gain should be ≈ 1.0
    let dc_gain: f64 = fir.coeffs.iter().sum();
    assert!(approx_eq(dc_gain, 1.0, 0.05));
}

#[test]
fn test_fir_highpass_dc_reject() {
    let fir = FirFilter::highpass(31, 0.5);
    let dc_gain: f64 = fir.coeffs.iter().sum();
    assert!(dc_gain.abs() < 0.1);
}

#[test]
fn test_fir_apply_length() {
    let fir = FirFilter::lowpass(15, 0.3);
    let x = vec![1.0; 64];
    let y = fir.apply(&x);
    assert_eq!(y.len(), 64 + 15 - 1);
}

// --- Biquad IIR ---

#[test]
fn test_biquad_passthrough() {
    // Unity gain: b=[1,0,0], a=[0,0]
    let mut bq = Biquad::new([1.0, 0.0, 0.0], [0.0, 0.0]);
    let x = vec![1.0, 2.0, 3.0, 4.0];
    let y = bq.process(&x);
    assert!(approx_eq(y[0], 1.0, EPS));
    assert!(approx_eq(y[1], 2.0, EPS));
}

#[test]
fn test_biquad_reset() {
    let mut bq = Biquad::new([0.5, 0.5, 0.0], [0.0, 0.0]);
    let _ = bq.process(&[1.0, 1.0, 1.0]);
    bq.reset();
    let y = bq.process(&[1.0]);
    assert!(approx_eq(y[0], 0.5, EPS));
}

// --- Hilbert transform ---

#[test]
fn test_hilbert_of_cos_is_sin() {
    let n = 256;
    let fs = 256.0f64;
    let f = 10.0f64;
    let cos_sig: Vec<f64> = (0..n)
        .map(|i| (2.0 * PI * f * i as f64 / fs).cos())
        .collect();
    let sin_sig: Vec<f64> = (0..n)
        .map(|i| (2.0 * PI * f * i as f64 / fs).sin())
        .collect();
    let ht = hilbert_transform(&cos_sig);
    // Compare middle section to avoid edge effects
    for i in 16..n - 16 {
        assert!(approx_eq(ht[i], sin_sig[i], 0.02));
    }
}

#[test]
fn test_instantaneous_amplitude() {
    let n = 256;
    let fs = 256.0f64;
    let f = 10.0f64;
    let amp = 3.0f64;
    let sig: Vec<f64> = (0..n)
        .map(|i| amp * (2.0 * PI * f * i as f64 / fs).cos())
        .collect();
    let env = instantaneous_amplitude(&sig);
    // Middle section should be ≈ amp
    for (_, &e) in env.iter().enumerate().take(n - 16).skip(16) {
        assert!(approx_eq(e, amp, 0.1));
    }
}

// --- Haar DWT ---

#[test]
fn test_haar_decompose_reconstruct() {
    let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let (approx, detail) = haar_decompose(&signal);
    let recon = haar_reconstruct(&approx, &detail);
    for (a, b) in signal.iter().zip(recon.iter()) {
        assert!(approx_eq(*a, *b, 1e-10));
    }
}

#[test]
fn test_haar_dwt_energy_preservation() {
    let signal: Vec<f64> = (0..8).map(|i| i as f64 + 1.0).collect();
    let coeffs = haar_dwt(&signal, 3);
    let orig_energy: f64 = signal.iter().map(|x| x * x).sum();
    let coeff_energy: f64 = coeffs.iter().map(|x| x * x).sum();
    assert!(approx_eq(orig_energy, coeff_energy, 1e-8));
}

// --- Spectrogram ---

#[test]
fn test_spectrogram_shape() {
    let n = 256;
    let fs = 1000.0f64;
    let signal: Vec<f64> = (0..n)
        .map(|i| (2.0 * PI * 100.0 * i as f64 / fs).sin())
        .collect();
    let sg = spectrogram(&signal, fs, 32, 16, WindowType::Hann);
    assert!(!sg.times.is_empty());
    assert!(!sg.freqs.is_empty());
    for frame in &sg.magnitude {
        assert_eq!(frame.len(), sg.freqs.len());
    }
}

// --- PSD ---

#[test]
fn test_welch_psd_bins_count() {
    let signal: Vec<f64> = (0..512).map(|i| (i as f64 * 0.1).sin()).collect();
    let result = welch_psd(&signal, 1000.0, 128, 64, WindowType::Hann).unwrap();
    assert!(!result.psd.is_empty());
    assert_eq!(result.freqs.len(), result.psd.len());
    assert_eq!(result.psd.len(), 128 / 2 + 1);
}

#[test]
fn test_periodogram_length() {
    let signal = vec![0.0f64; 64];
    let psd = periodogram_psd(&signal, 100.0, WindowType::Hamming).unwrap();
    assert!(!psd.psd.is_empty());
    assert_eq!(psd.psd.len(), 33);
}

fn assert_invalid_input<T: std::fmt::Debug>(result: Result<T, PhysicsError>, field: &str) {
    match result {
        Err(PhysicsError::InvalidInput { field: f, .. }) => assert_eq!(f, field),
        other => panic!("expected InvalidInput for `{field}`, got {other:?}"),
    }
}

#[test]
fn test_periodogram_of_empty_signal_is_an_error_not_a_hang() {
    // Regression: segment_len = 0 and overlap = 0 made the segment loop spin forever.
    assert_invalid_input(periodogram_psd(&[], 1.0, WindowType::Hann), "segment_len");
}

#[test]
fn test_welch_rejects_overlap_not_smaller_than_segment() {
    let signal: Vec<f64> = (0..256).map(|i| (i as f64 * 0.3).sin()).collect();
    // overlap == segment_len used to give a zero hop (infinite loop).
    assert_invalid_input(welch_psd(&signal, 1.0, 64, 64, WindowType::Hann), "overlap");
    // overlap > segment_len used to underflow `segment_len - overlap`.
    assert_invalid_input(welch_psd(&signal, 1.0, 64, 80, WindowType::Hann), "overlap");
}

#[test]
fn test_welch_rejects_degenerate_arguments() {
    let signal: Vec<f64> = (0..100).map(|i| i as f64).collect();
    assert_invalid_input(
        welch_psd(&signal, 1.0, 0, 0, WindowType::Hann),
        "segment_len",
    );
    assert_invalid_input(welch_psd(&signal, 1.0, 128, 0, WindowType::Hann), "signal");
    assert_invalid_input(welch_psd(&[], 1.0, 8, 0, WindowType::Hann), "signal");
    for fs in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_invalid_input(welch_psd(&signal, fs, 16, 8, WindowType::Hann), "fs");
    }
    // A length-2 Hann or Blackman window is identically zero.
    assert_invalid_input(
        welch_psd(&signal, 1.0, 2, 0, WindowType::Hann),
        "window_type",
    );
    assert_invalid_input(
        welch_psd(&signal, 1.0, 2, 0, WindowType::Blackman),
        "window_type",
    );
}

#[test]
fn test_bartlett_psd_validates_its_arguments() {
    let signal = vec![1.0f64; 64];
    assert_invalid_input(bartlett_psd(&[], 1.0, 16), "signal");
    assert_invalid_input(bartlett_psd(&signal, 1.0, 0), "segment_len");
    assert_invalid_input(bartlett_psd(&signal, 1.0, 65), "signal");
    assert!(bartlett_psd(&signal, 1.0, 16).is_ok());
}

#[test]
fn test_welch_with_single_sample_segments_is_finite() {
    // hann_window(1) used to be NaN (0/0); the 1-sample window is 1.0.
    let signal = vec![1.0, -2.0, 3.0, 0.5];
    let result = welch_psd(&signal, 2.0, 1, 0, WindowType::Hann).unwrap();
    assert_eq!(result.psd.len(), 1);
    assert!(result.psd[0].is_finite());
    // PSD = |x|^2 / (fs * sum w^2) averaged over the samples.
    let expected = signal.iter().map(|v| v * v).sum::<f64>() / signal.len() as f64 / 2.0;
    assert!(approx_eq(result.psd[0], expected, 1e-12));
}

#[test]
fn test_bartlett_psd_uses_a_rectangular_window() {
    // For a constant signal the DC bin of a rectangular-window periodogram is
    // L^2 / (fs * L) = L / fs; a Hann window would give only (2/3) L / fs.
    let signal = vec![1.0f64; 256];
    let bartlett = bartlett_psd(&signal, 100.0, 64).unwrap();
    assert!(approx_eq(bartlett.psd[0], 64.0 / 100.0, 1e-9));
    let rectangular = welch_psd(&signal, 100.0, 64, 0, WindowType::Rectangular).unwrap();
    assert_eq!(bartlett.psd, rectangular.psd);
    let hann = welch_psd(&signal, 100.0, 64, 0, WindowType::Hann).unwrap();
    assert!(hann.psd[0] < 0.7 * bartlett.psd[0]);
}

#[test]
fn test_welch_psd_integrates_to_the_signal_power_for_every_window() {
    // A bin-centred sine of amplitude 2 has power A^2 / 2 = 2.
    let fs = 1000.0;
    let n = 2048;
    let signal: Vec<f64> = (0..n)
        .map(|i| 2.0 * (2.0 * PI * 125.0 * i as f64 / fs).sin())
        .collect();
    let windows = [
        WindowType::Rectangular,
        WindowType::Hann,
        WindowType::Hamming,
        WindowType::Blackman,
        WindowType::BlackmanHarris,
        WindowType::FlatTop,
        WindowType::Bartlett,
        WindowType::Kaiser,
        WindowType::Nuttall,
    ];
    for window in windows {
        let psd = welch_psd(&signal, fs, 256, 128, window).unwrap();
        let df = psd.freqs[1] - psd.freqs[0];
        let power: f64 = psd.psd.iter().sum::<f64>() * df;
        assert!(
            (power - 2.0).abs() < 0.06,
            "{window:?}: integrated power {power}"
        );
        // The spectral peak sits at the sine frequency.
        let (peak_bin, _) = psd
            .psd
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .unwrap();
        assert!(
            approx_eq(psd.freqs[peak_bin], 125.0, 1e-9),
            "{window:?}: peak at {}",
            psd.freqs[peak_bin]
        );
    }
}

#[test]
fn test_window_functions_handle_lengths_zero_and_one() {
    type WindowFn = fn(usize) -> Vec<f64>;
    let windows: [(&str, WindowFn); 7] = [
        ("hann", hann_window),
        ("hamming", hamming_window),
        ("blackman", blackman_window),
        ("blackman_harris", blackman_harris_window),
        ("nuttall", nuttall_window),
        ("bartlett", bartlett_window),
        ("flat_top", flat_top_window),
    ];
    for (name, window) in windows {
        assert!(window(0).is_empty(), "{name}(0)");
        assert_eq!(window(1), vec![1.0], "{name}(1)");
        let two = window(2);
        assert_eq!(two.len(), 2, "{name}(2)");
        assert!(two.iter().all(|v| v.is_finite()), "{name}(2)");
    }
    assert!(kaiser_window(0, 5.0).is_empty());
    assert_eq!(kaiser_window(1, 5.0), vec![1.0]);
}

#[test]
fn test_window_type_generate_matches_the_window_functions() {
    let n = 33;
    assert_eq!(WindowType::Rectangular.generate(n), vec![1.0; n]);
    assert_eq!(WindowType::Hann.generate(n), hann_window(n));
    assert_eq!(WindowType::Hamming.generate(n), hamming_window(n));
    assert_eq!(WindowType::Blackman.generate(n), blackman_window(n));
    assert_eq!(
        WindowType::BlackmanHarris.generate(n),
        blackman_harris_window(n)
    );
    assert_eq!(WindowType::FlatTop.generate(n), flat_top_window(n));
    assert_eq!(WindowType::Bartlett.generate(n), bartlett_window(n));
    assert_eq!(WindowType::Nuttall.generate(n), nuttall_window(n));
    assert_eq!(
        WindowType::Kaiser.generate(n),
        kaiser_window(n, DEFAULT_KAISER_BETA)
    );
}

#[test]
fn test_spectrogram_and_stft_with_zero_hop_or_window_return_no_frames() {
    let signal = vec![1.0f64; 64];
    let sg = spectrogram(&signal, 100.0, 16, 0, WindowType::Hann);
    assert!(sg.times.is_empty() && sg.magnitude.is_empty());
    let sg = spectrogram(&signal, 100.0, 0, 8, WindowType::Hann);
    assert!(sg.times.is_empty() && sg.magnitude.is_empty());
    assert!(stft(&signal, 16, 0, WindowType::Hann).is_empty());
    assert!(stft(&signal, 0, 8, WindowType::Hann).is_empty());
}

#[test]
fn test_spectrogram_honours_the_rectangular_window() {
    // DC bin of a constant frame = sum of the window: 8 for a rectangular
    // window of length 8, only 3.5 for a Hann window.
    let signal = vec![1.0f64; 8];
    let rect = spectrogram(&signal, 8.0, 8, 4, WindowType::Rectangular);
    assert!(approx_eq(rect.magnitude[0][0], 8.0, 1e-9));
    let hann = spectrogram(&signal, 8.0, 8, 4, WindowType::Hann);
    assert!(approx_eq(hann.magnitude[0][0], 3.5, 1e-9));
}

// --- Resampling ---

#[test]
fn test_resample_linear_length() {
    let signal = vec![0.0, 1.0, 2.0, 3.0];
    let out = resample_linear(&signal, 8);
    assert_eq!(out.len(), 8);
}

#[test]
fn test_upsample_downsample() {
    let signal = vec![1.0, 2.0, 3.0, 4.0];
    let up = upsample(&signal, 2);
    assert_eq!(up.len(), 8);
    assert_eq!(up[0], 1.0);
    assert_eq!(up[2], 2.0);
    let dn = downsample(&up, 2);
    assert_eq!(dn, signal);
}

// --- Utilities ---

#[test]
fn test_rms() {
    let signal = vec![1.0, -1.0, 1.0, -1.0];
    assert!(approx_eq(rms(&signal), 1.0, EPS));
}

#[test]
fn test_signal_energy() {
    let signal = vec![3.0, 4.0];
    assert!(approx_eq(signal_energy(&signal), 25.0, EPS));
}

#[test]
fn test_normalize_energy() {
    let signal = vec![3.0, 4.0];
    let normed = normalize_energy(&signal);
    assert!(approx_eq(signal_energy(&normed), 1.0, 1e-10));
}

#[test]
fn test_find_peaks() {
    let signal = vec![0.0, 1.0, 0.5, 2.0, 0.0, 1.5, 0.0];
    let peaks = find_peaks(&signal);
    assert!(peaks.contains(&1));
    assert!(peaks.contains(&3));
    assert!(peaks.contains(&5));
}

#[test]
fn test_moving_average() {
    let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let ma = moving_average(&signal, 3);
    assert_eq!(ma.len(), signal.len());
    assert!(approx_eq(ma[4], 4.0, EPS));
}

#[test]
fn test_snr_db() {
    let snr = snr_db(100.0, 1.0);
    assert!(approx_eq(snr, 20.0, 1e-6));
}

#[test]
fn test_generate_sine() {
    let sine = generate_sine(100.0, 1000.0, 100, 1.0, 0.0);
    assert_eq!(sine.len(), 100);
    // First sample should be ≈ 0
    assert!(sine[0].abs() < 1e-10);
}

#[test]
fn test_zero_pad() {
    let sig = vec![1.0, 2.0, 3.0];
    let padded = zero_pad(&sig, 8);
    assert_eq!(padded.len(), 8);
    assert_eq!(padded[3], 0.0);
}

#[test]
fn test_fft_frequencies() {
    let freqs = fft_frequencies(8, 8.0);
    assert_eq!(freqs.len(), 8);
    assert!(approx_eq(freqs[1], 1.0, EPS));
    assert!(approx_eq(freqs[4], 4.0, EPS));
}

#[test]
fn test_rfft_frequencies() {
    let freqs = rfft_frequencies(8, 8.0);
    assert_eq!(freqs.len(), 5); // 0..=4
    assert!(approx_eq(freqs[4], 4.0, EPS));
}

#[test]
fn test_real_cepstrum_length() {
    let signal: Vec<f64> = (0..16).map(|i| (i as f64).sin()).collect();
    let cep = real_cepstrum(&signal);
    assert_eq!(cep.len(), 16);
}

#[test]
fn test_bartlett_window_peak() {
    let w = bartlett_window(9);
    // Peak at center index 4
    assert!(approx_eq(w[4], 1.0, 1e-10));
}

#[test]
fn test_blackman_harris_range() {
    let w = blackman_harris_window(64);
    for &v in &w {
        assert!((-0.01..=1.01).contains(&v));
    }
}

#[test]
fn test_db2_decompose_reconstruct() {
    let signal = vec![1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0];
    let (approx, detail) = db2_decompose(&signal);
    assert!(!approx.is_empty());
    assert!(!detail.is_empty());
    let recon = db2_reconstruct(&approx, &detail);
    assert!(!recon.is_empty());
}

#[test]
fn test_stft_istft_roundtrip_shape() {
    let n = 128;
    let signal: Vec<f64> = (0..n).map(|i| (i as f64 * 0.1).sin()).collect();
    let frames = stft(&signal, 32, 16, WindowType::Hann);
    assert!(!frames.is_empty());
    for frame in &frames {
        assert_eq!(frame.len(), 17); // 32/2 + 1
    }
}

#[test]
fn test_butterworth_lowpass_creates_filter() {
    let filt = butterworth_lowpass(4, 0.5);
    assert_eq!(filt.sections.len(), 2);
}

#[test]
fn test_chebyshev1_lowpass_creates_filter() {
    let filt = chebyshev1_lowpass(4, 0.4, 1.0);
    assert!(!filt.sections.is_empty());
}
