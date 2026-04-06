//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
use std::f64::consts::PI;

/// Mother wavelet type for CWT.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MotherWavelet {
    /// Morlet wavelet: exp(-t^2/2) * cos(omega_0 * t).
    Morlet {
        /// Central frequency parameter (default: 6.0).
        omega0: f64,
    },
    /// Mexican hat (Ricker) wavelet: (1 - t^2) * exp(-t^2/2).
    MexicanHat,
}
impl MotherWavelet {
    /// Evaluate the mother wavelet at point t.
    pub fn evaluate(&self, t: f64) -> f64 {
        match self {
            MotherWavelet::Morlet { omega0 } => (-t * t / 2.0).exp() * (omega0 * t).cos(),
            MotherWavelet::MexicanHat => {
                let t2 = t * t;
                (2.0 / (3.0_f64.sqrt() * PI.powf(0.25))) * (1.0 - t2) * (-t2 / 2.0).exp()
            }
        }
    }
    /// Create a default Morlet wavelet (omega0 = 6.0).
    pub fn morlet() -> Self {
        MotherWavelet::Morlet { omega0: 6.0 }
    }
    /// Create a Mexican hat wavelet.
    pub fn mexican_hat() -> Self {
        MotherWavelet::MexicanHat
    }
}
/// Haar wavelet transform (standalone implementation).
///
/// Provides a simple O(n) forward and inverse Haar wavelet transform.
pub struct HaarTransform;
impl HaarTransform {
    /// Forward Haar transform (single level).
    ///
    /// Returns `(approx, detail)` coefficients.
    pub fn forward(signal: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = signal.len() / 2;
        let inv_sqrt2 = std::f64::consts::FRAC_1_SQRT_2;
        let mut approx = Vec::with_capacity(n);
        let mut detail = Vec::with_capacity(n);
        for i in 0..n {
            let a = signal[2 * i];
            let b = signal[2 * i + 1];
            approx.push((a + b) * inv_sqrt2);
            detail.push((a - b) * inv_sqrt2);
        }
        (approx, detail)
    }
    /// Inverse Haar transform (single level).
    pub fn inverse(approx: &[f64], detail: &[f64]) -> Vec<f64> {
        let n = approx.len().min(detail.len());
        let inv_sqrt2 = std::f64::consts::FRAC_1_SQRT_2;
        let mut out = vec![0.0; 2 * n];
        for i in 0..n {
            out[2 * i] = (approx[i] + detail[i]) * inv_sqrt2;
            out[2 * i + 1] = (approx[i] - detail[i]) * inv_sqrt2;
        }
        out
    }
    /// Multi-level forward Haar transform (in-place style, returns all coefficients).
    pub fn forward_multilevel(signal: &[f64], levels: usize) -> Vec<Vec<f64>> {
        let mut current = signal.to_vec();
        let mut all_details = Vec::with_capacity(levels);
        for _ in 0..levels {
            if current.len() < 2 {
                break;
            }
            let (a, d) = Self::forward(&current);
            all_details.push(d);
            current = a;
        }
        let mut result = vec![current];
        result.extend(all_details.into_iter().rev());
        result
    }
}
/// Node in a wavelet packet tree.
#[derive(Debug, Clone)]
pub struct WaveletPacketNode {
    /// Coefficients at this node.
    pub coefficients: Vec<f64>,
    /// Level in the tree (0 = root).
    pub level: usize,
    /// Position in the level (0-indexed).
    pub position: usize,
    /// Shannon entropy of coefficients (for best-basis selection).
    pub entropy: f64,
}
/// Lifting scheme implementation for the Haar wavelet.
///
/// The lifting scheme provides a computationally efficient, in-place
/// implementation of the wavelet transform.
pub struct LiftingHaar;
impl LiftingHaar {
    /// Forward lifting transform (in-place).
    ///
    /// After transform, even indices contain approximation and
    /// odd indices contain detail coefficients.
    pub fn forward(data: &mut [f64]) {
        let n = data.len();
        if n < 2 {
            return;
        }
        let half = n / 2;
        let mut even: Vec<f64> = (0..half).map(|i| data[2 * i]).collect();
        let mut odd: Vec<f64> = (0..half).map(|i| data[2 * i + 1]).collect();
        for i in 0..half {
            odd[i] -= even[i];
        }
        for i in 0..half {
            even[i] += odd[i] / 2.0;
        }
        data[..half].copy_from_slice(&even);
        data[half..2 * half].copy_from_slice(&odd);
    }
    /// Inverse lifting transform (in-place).
    pub fn inverse(data: &mut [f64]) {
        let n = data.len();
        if n < 2 {
            return;
        }
        let half = n / 2;
        let mut even: Vec<f64> = data[..half].to_vec();
        let mut odd: Vec<f64> = data[half..half * 2].to_vec();
        for i in 0..half {
            even[i] -= odd[i] / 2.0;
        }
        for i in 0..half {
            odd[i] += even[i];
        }
        for i in 0..half {
            data[2 * i] = even[i];
            data[2 * i + 1] = odd[i];
        }
    }
}
/// Thresholding mode for wavelet denoising.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThresholdMode {
    /// Hard thresholding: zero out coefficients below threshold.
    Hard,
    /// Soft thresholding: shrink coefficients toward zero.
    Soft,
}
/// Multi-level DWT decomposition result.
#[derive(Debug, Clone)]
pub struct DwtDecomposition {
    /// Detail coefficients at each level (from finest to coarsest).
    pub details: Vec<Vec<f64>>,
    /// Final approximation coefficients (coarsest level).
    pub approx: Vec<f64>,
    /// Wavelet family used for decomposition.
    pub wavelet: WaveletFamily,
    /// Original signal lengths at each level (for reconstruction).
    pub lengths: Vec<usize>,
}
/// Daubechies wavelet family (db2 through db6).
///
/// Provides methods for forward/inverse transforms at various vanishing moments.
pub struct DaubechiesWavelet {
    /// The wavelet family variant.
    pub family: WaveletFamily,
}
impl DaubechiesWavelet {
    /// Create a new Daubechies wavelet with the given order (2..=6).
    ///
    /// # Panics
    /// Panics if order is not in 2..=6.
    pub fn new(order: usize) -> Self {
        let family = match order {
            2 => WaveletFamily::Db2,
            3 => WaveletFamily::Db3,
            4 => WaveletFamily::Db4,
            5 => WaveletFamily::Db5,
            6 => WaveletFamily::Db6,
            _ => panic!("Daubechies order must be 2..=6, got {}", order),
        };
        Self { family }
    }
    /// Single-level forward DWT using this wavelet.
    pub fn forward(&self, signal: &[f64]) -> DwtLevel {
        dwt_single(signal, self.family)
    }
    /// Single-level inverse DWT using this wavelet.
    pub fn inverse(&self, level: &DwtLevel, target_len: usize) -> Vec<f64> {
        idwt_single(level, self.family, target_len)
    }
    /// Multi-level forward DWT.
    pub fn decompose(&self, signal: &[f64], levels: usize) -> DwtDecomposition {
        dwt(signal, self.family, levels)
    }
    /// Multi-level inverse DWT.
    pub fn reconstruct(&self, decomp: &DwtDecomposition) -> Vec<f64> {
        idwt(decomp)
    }
}
/// Multiresolution analysis: reconstructs the signal at each resolution level.
///
/// Returns a vector of approximations from coarsest to finest.
pub struct MultiresolutionAnalysis {
    /// Approximation at each level (index 0 = coarsest).
    pub approximations: Vec<Vec<f64>>,
    /// Detail contributions at each level.
    pub detail_contributions: Vec<Vec<f64>>,
}
/// Scalogram: squared magnitude of CWT coefficients.
///
/// Represents the wavelet energy spectrum in the time-scale plane.
#[derive(Debug, Clone)]
pub struct Scalogram {
    /// Energy values: `energy[scale_idx][position]`.
    pub energy: Vec<Vec<f64>>,
    /// Scales.
    pub scales: Vec<f64>,
    /// Total energy at each scale.
    pub scale_energy: Vec<f64>,
}
/// Wavelet packet decomposition tree.
pub struct WaveletPacketTree {
    /// All nodes organized by level.
    pub nodes: Vec<Vec<WaveletPacketNode>>,
    /// Wavelet family used.
    pub wavelet: WaveletFamily,
    /// Maximum decomposition depth.
    pub max_level: usize,
}
/// Result of a continuous wavelet transform.
#[derive(Debug, Clone)]
pub struct CwtResult {
    /// CWT coefficient matrix: `coefficients[scale_idx][position]`.
    pub coefficients: Vec<Vec<f64>>,
    /// Scales used for the transform.
    pub scales: Vec<f64>,
    /// Mother wavelet used.
    pub wavelet: MotherWavelet,
}
/// Stationary (undecimated) wavelet transform: no downsampling, shift-invariant.
///
/// Returns detail and approximation coefficients at each level,
/// all with the same length as the input signal.
pub struct SwtDecomposition {
    /// Detail coefficients at each level.
    pub details: Vec<Vec<f64>>,
    /// Approximation at each level (last is the coarsest).
    pub approx: Vec<f64>,
    /// Wavelet used.
    pub wavelet: WaveletFamily,
}
/// Maximal Overlap DWT (MODWT) - shift-invariant variant of the DWT.
///
/// All coefficient vectors have the same length as the input signal.
pub struct ModwtDecomposition {
    /// Detail coefficients at each level.
    pub details: Vec<Vec<f64>>,
    /// Final approximation coefficients.
    pub approx: Vec<f64>,
    /// Wavelet used.
    pub wavelet: WaveletFamily,
}
/// Result of a single-level DWT decomposition.
#[derive(Debug, Clone)]
pub struct DwtLevel {
    /// Approximation coefficients (low-pass).
    pub approx: Vec<f64>,
    /// Detail coefficients (high-pass).
    pub detail: Vec<f64>,
}
/// Supported wavelet families.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WaveletFamily {
    /// Haar (Daubechies-1) wavelet.
    Haar,
    /// Daubechies-2 wavelet.
    Db2,
    /// Daubechies-3 wavelet.
    Db3,
    /// Daubechies-4 wavelet.
    Db4,
    /// Daubechies-5 wavelet.
    Db5,
    /// Daubechies-6 wavelet.
    Db6,
}
impl WaveletFamily {
    /// Get the low-pass decomposition filter for the wavelet.
    pub fn lo_dec(&self) -> &'static [f64] {
        match self {
            WaveletFamily::Haar => &HAAR_LO,
            WaveletFamily::Db2 => &DB2_LO,
            WaveletFamily::Db3 => &DB3_LO,
            WaveletFamily::Db4 => &DB4_LO,
            WaveletFamily::Db5 => &DB5_LO,
            WaveletFamily::Db6 => &DB6_LO,
        }
    }
    /// Get the high-pass decomposition filter for the wavelet.
    pub fn hi_dec(&self) -> &'static [f64] {
        match self {
            WaveletFamily::Haar => &HAAR_HI,
            WaveletFamily::Db2 => &DB2_HI,
            WaveletFamily::Db3 => &DB3_HI,
            WaveletFamily::Db4 => &DB4_HI,
            WaveletFamily::Db5 => &DB5_HI,
            WaveletFamily::Db6 => &DB6_HI,
        }
    }
    /// Get the low-pass reconstruction filter (reverse of decomposition).
    pub fn lo_rec(&self) -> Vec<f64> {
        self.lo_dec().iter().rev().copied().collect()
    }
    /// Get the high-pass reconstruction filter (reverse of decomposition).
    pub fn hi_rec(&self) -> Vec<f64> {
        self.hi_dec().iter().rev().copied().collect()
    }
    /// Filter length for this wavelet family.
    pub fn filter_length(&self) -> usize {
        self.lo_dec().len()
    }
}
