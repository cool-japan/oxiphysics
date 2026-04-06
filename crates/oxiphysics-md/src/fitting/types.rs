//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
/// Generic parametric potential: V = Σ_k c_k * f_k(x)
/// where f_k are basis functions and c_k are the coefficients.
///
/// This struct holds the basis functions and provides linear least-squares
/// fitting to energies.
#[allow(dead_code)]
pub struct ParametricPotential {
    /// Basis functions f_k(x). Must all be the same length when evaluated.
    pub basis_fns: Vec<Box<dyn Fn(f64) -> f64 + Send + Sync>>,
}
impl ParametricPotential {
    /// Create a new parametric potential from basis functions.
    #[allow(dead_code)]
    pub fn new(basis_fns: Vec<Box<dyn Fn(f64) -> f64 + Send + Sync>>) -> Self {
        Self { basis_fns }
    }
    /// Evaluate all basis functions at `x`.
    #[allow(dead_code)]
    pub fn features(&self, x: f64) -> Vec<f64> {
        self.basis_fns.iter().map(|f| f(x)).collect()
    }
    /// Fit coefficients to (x, energy) data using ordinary least squares.
    ///
    /// Solves the normal equations: (Φᵀ Φ) c = Φᵀ y.
    #[allow(dead_code)]
    pub fn fit(&self, xs: &[f64], energies: &[f64]) -> Vec<f64> {
        let n = xs.len();
        let k = self.basis_fns.len();
        assert_eq!(n, energies.len());
        assert!(k > 0);
        let phi: Vec<Vec<f64>> = xs.iter().map(|&x| self.features(x)).collect();
        let mut ata = vec![vec![0.0f64; k]; k];
        let mut aty = vec![0.0f64; k];
        for i in 0..n {
            for a in 0..k {
                aty[a] += phi[i][a] * energies[i];
                for b in 0..k {
                    ata[a][b] += phi[i][a] * phi[i][b];
                }
            }
        }
        for a in 0..k {
            ata[a][a] += 1e-12;
        }
        solve_linear_system(&ata, &aty)
    }
    /// Predict energy at `x` using fitted coefficients.
    #[allow(dead_code)]
    pub fn predict(&self, x: f64, coefficients: &[f64]) -> f64 {
        self.features(x)
            .iter()
            .zip(coefficients.iter())
            .map(|(f, c)| f * c)
            .sum()
    }
}
/// Result of a RESP (Restrained ElectroStatic Potential) charge fit.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RespFitResult {
    /// Fitted partial charges (one per atom).
    pub charges: Vec<f64>,
    /// RMSE of ESP fit.
    pub rmse: f64,
    /// Net charge constraint satisfaction error.
    pub net_charge_error: f64,
}
/// Result of a Lennard-Jones parameter fit.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct LjFitResult {
    /// Potential-well depth ε.
    pub epsilon: f64,
    /// Zero-crossing distance σ.
    pub sigma: f64,
    /// Root-mean-square error of the fit.
    pub rmse: f64,
}
/// Result of a harmonic bond fit.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct HarmonicBondFit {
    /// Equilibrium bond length.
    pub r0: f64,
    /// Spring constant k.
    pub k: f64,
}
/// Bootstrap resampling for fitting uncertainty.
///
/// Draws `n_bootstrap` samples with replacement, fits each, and returns
/// the standard deviation of the k and r0 estimates.
#[allow(dead_code)]
pub struct BootstrapResult {
    /// Standard deviation of k estimates.
    pub std_k: f64,
    /// Standard deviation of r0 estimates.
    pub std_r0: f64,
    /// Number of bootstrap samples.
    pub n_bootstrap: usize,
}
/// Result of a harmonic angle fit.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct HarmonicAngleFit {
    /// Equilibrium angle θ₀ (radians).
    pub theta0: f64,
    /// Spring constant k.
    pub k: f64,
}
/// Result of a single-term dihedral (OPLS-style) fit.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct DihedralFit {
    /// Torsion barrier height k.
    pub k: f64,
    /// Periodicity n.
    pub n: u32,
    /// Phase offset δ (radians).
    pub delta: f64,
}
