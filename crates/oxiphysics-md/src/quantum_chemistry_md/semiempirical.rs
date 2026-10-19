// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Semiempirical-flavoured model pair potentials for Born-Oppenheimer MD.
//!
//! [`SemiempiricalMd`] evaluates a *model* core-core energy for each
//! [`SemiempiricalMethod`]: a screened Coulomb repulsion between the cores plus
//! the Gaussian modifiers that characterise the AM1/PM3 core-core functions,
//!
//! ```text
//! E = Σ_{i<j} f(r_ij),   f(r) = Z_i Z_j e^{−α r} / r + Σ_k a_k e^{−w_k (r − c_k)²}.
//! ```
//!
//! The parameters (α, a_k, w_k, c_k) are fixed per method and are **not** the
//! element-specific NDDO parameter sets of AM1, PM3, GFN2-xTB or DFTB: no
//! electronic Hamiltonian is solved.  The energy and its analytic gradient are,
//! however, exactly consistent with each other, so the driver produces
//! energy-conserving dynamics on the model surface.

use super::{SemiempiricalMethod, dist};

/// Distance below which a pair is treated as coincident (no defined direction).
const MIN_PAIR_DISTANCE: f64 = 1e-6;

/// Parameters of one model pair potential.
#[derive(Debug, Clone, Copy)]
struct PairModel {
    /// Screening exponent α of the core-core Coulomb term \[1/Å\].
    alpha: f64,
    /// Gaussian modifiers `(a, w, c)`: a·exp(−w (r − c)²).
    gaussians: &'static [(f64, f64, f64)],
}

impl PairModel {
    fn for_method(method: &SemiempiricalMethod) -> Self {
        match method {
            SemiempiricalMethod::Am1 => Self {
                alpha: 1.2,
                gaussians: &[(0.05, 0.5, 1.8)],
            },
            SemiempiricalMethod::Pm3 => Self {
                alpha: 1.4,
                gaussians: &[(0.03, 0.6, 1.5), (0.02, 0.6, 2.1)],
            },
            SemiempiricalMethod::Xtb => Self {
                alpha: 1.0,
                gaussians: &[],
            },
            SemiempiricalMethod::Dftb => Self {
                alpha: 1.1,
                gaussians: &[],
            },
        }
    }

    /// Pair energy f(r) in eV for the core charge product `zz`.
    fn energy(&self, zz: f64, r: f64) -> f64 {
        let screened = zz * (-self.alpha * r).exp() / r;
        let modifiers: f64 = self
            .gaussians
            .iter()
            .map(|&(a, w, c)| a * (-w * (r - c).powi(2)).exp())
            .sum();
        screened + modifiers
    }

    /// Radial derivative f'(r) in eV/Å.
    fn derivative(&self, zz: f64, r: f64) -> f64 {
        let screened = -zz * (-self.alpha * r).exp() * (self.alpha * r + 1.0) / (r * r);
        let modifiers: f64 = self
            .gaussians
            .iter()
            .map(|&(a, w, c)| -2.0 * w * a * (r - c) * (-w * (r - c).powi(2)).exp())
            .sum();
        screened + modifiers
    }
}

/// Semiempirical-flavoured model MD driver.
///
/// Holds the nuclear configuration and evaluates the model pair potential of
/// the selected [`SemiempiricalMethod`] (see the module documentation for the
/// functional form and its limits).  Energies are in eV, lengths in Å.
#[derive(Debug, Clone)]
pub struct SemiempiricalMd {
    /// Model pair potential to use.
    pub method: SemiempiricalMethod,
    /// Atomic numbers (used as core charges of the screened Coulomb term).
    pub atomic_numbers: Vec<u32>,
    /// Cartesian positions \[Angstrom\].
    pub positions: Vec<[f64; 3]>,
    /// Velocities \[Angstrom/fs\].
    pub velocities: Vec<[f64; 3]>,
    /// Masses \[amu\].
    pub masses: Vec<f64>,
}

impl SemiempiricalMd {
    /// Create a new semiempirical MD driver.
    ///
    /// # Arguments
    /// * `method` – Model pair potential type.
    /// * `atomic_numbers` – Atomic numbers Z for each atom.
    /// * `positions` – Initial Cartesian positions in Angstrom.
    /// * `masses` – Atomic masses in amu.
    pub fn new(
        method: SemiempiricalMethod,
        atomic_numbers: Vec<u32>,
        positions: Vec<[f64; 3]>,
        masses: Vec<f64>,
    ) -> Self {
        let n = positions.len();
        Self {
            method,
            atomic_numbers,
            positions,
            velocities: vec![[0.0; 3]; n],
            masses,
        }
    }

    /// Total model energy in eV for the configured method.
    pub fn energy(&self) -> f64 {
        self.pair_energy(&PairModel::for_method(&self.method))
    }

    /// Analytic energy gradient ∂E/∂R_i in eV/Å (the negative of the forces).
    ///
    /// For every pair, g_i += f'(r) (R_i − R_j)/r and g_j −= the same, where
    /// f is the pair potential of the configured method, so the gradient is the
    /// exact derivative of [`energy`](Self::energy).  Coincident atoms
    /// (r < 1e-6 Å) contribute nothing because the direction is undefined.
    pub fn gradient(&self) -> Vec<[f64; 3]> {
        let model = PairModel::for_method(&self.method);
        let n = self.positions.len();
        let mut g = vec![[0.0f64; 3]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist(&self.positions[i], &self.positions[j]);
                if r < MIN_PAIR_DISTANCE {
                    continue;
                }
                let zz = self.core_charge(i) * self.core_charge(j);
                let scale = model.derivative(zz, r) / r;
                for k in 0..3 {
                    let component = scale * (self.positions[i][k] - self.positions[j][k]);
                    g[i][k] += component;
                    g[j][k] -= component;
                }
            }
        }
        g
    }

    /// Forces −∂E/∂R_i in eV/Å.
    pub fn forces(&self) -> Vec<[f64; 3]> {
        self.gradient()
            .into_iter()
            .map(|g| [-g[0], -g[1], -g[2]])
            .collect()
    }

    /// AM1-flavoured model energy in eV.
    ///
    /// f(r) = Z_i Z_j e^{−1.2 r}/r + 0.05 e^{−0.5 (r − 1.8)²}.
    pub fn am1_energy(&self) -> f64 {
        self.pair_energy(&PairModel::for_method(&SemiempiricalMethod::Am1))
    }

    /// PM3-flavoured model energy in eV.
    ///
    /// f(r) = Z_i Z_j e^{−1.4 r}/r + 0.03 e^{−0.6 (r − 1.5)²} + 0.02 e^{−0.6 (r − 2.1)²}.
    pub fn pm3_energy(&self) -> f64 {
        self.pair_energy(&PairModel::for_method(&SemiempiricalMethod::Pm3))
    }

    /// xTB-flavoured model energy in eV: f(r) = Z_i Z_j e^{−r}/r.
    pub fn xtb_energy(&self) -> f64 {
        self.pair_energy(&PairModel::for_method(&SemiempiricalMethod::Xtb))
    }

    /// DFTB-flavoured model energy in eV: f(r) = Z_i Z_j e^{−1.1 r}/r.
    pub fn dftb_energy(&self) -> f64 {
        self.pair_energy(&PairModel::for_method(&SemiempiricalMethod::Dftb))
    }

    fn core_charge(&self, i: usize) -> f64 {
        self.atomic_numbers.get(i).copied().unwrap_or(0) as f64
    }

    fn pair_energy(&self, model: &PairModel) -> f64 {
        let n = self.positions.len();
        let mut e = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist(&self.positions[i], &self.positions[j]).max(MIN_PAIR_DISTANCE);
                let zz = self.core_charge(i) * self.core_charge(j);
                e += model.energy(zz, r);
            }
        }
        e
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ALL_METHODS: [SemiempiricalMethod; 4] = [
        SemiempiricalMethod::Am1,
        SemiempiricalMethod::Pm3,
        SemiempiricalMethod::Xtb,
        SemiempiricalMethod::Dftb,
    ];

    fn triatomic(method: SemiempiricalMethod) -> SemiempiricalMd {
        // Non-collinear geometry so that every Cartesian component is exercised.
        SemiempiricalMd::new(
            method,
            vec![8, 1, 6],
            vec![[0.0, 0.0, 0.0], [0.96, 0.31, -0.12], [-0.4, 1.35, 0.8]],
            vec![16.0, 1.0, 12.0],
        )
    }

    #[test]
    fn test_semi_am1_energy_positive() {
        let md = SemiempiricalMd::new(
            SemiempiricalMethod::Am1,
            vec![1, 6],
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![1.0, 12.0],
        );
        assert!(md.am1_energy() > 0.0);
    }

    #[test]
    fn test_semi_pm3_energy_positive() {
        let md = SemiempiricalMd::new(
            SemiempiricalMethod::Pm3,
            vec![1, 6],
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![1.0, 12.0],
        );
        assert!(md.pm3_energy() > 0.0);
    }

    #[test]
    fn test_semi_energy_dispatch_matches_method_routines() {
        let base = triatomic(SemiempiricalMethod::Am1);
        let expected = [
            base.am1_energy(),
            base.pm3_energy(),
            base.xtb_energy(),
            base.dftb_energy(),
        ];
        for (method, e_ref) in ALL_METHODS.into_iter().zip(expected) {
            let md = triatomic(method);
            assert!((md.energy() - e_ref).abs() < 1e-12);
            assert!(md.energy() > 0.0);
        }
    }

    #[test]
    fn test_semi_am1_pair_energy_closed_form() {
        let r: f64 = 1.3;
        let md = SemiempiricalMd::new(
            SemiempiricalMethod::Am1,
            vec![6, 8],
            vec![[0.0, 0.0, 0.0], [r, 0.0, 0.0]],
            vec![12.0, 16.0],
        );
        let expected = 48.0 * (-1.2 * r).exp() / r + 0.05 * (-0.5 * (r - 1.8).powi(2)).exp();
        assert!((md.energy() - expected).abs() < 1e-12);
    }

    #[test]
    fn test_semi_gradient_matches_finite_difference_all_methods() {
        let h = 1e-5;
        for method in ALL_METHODS {
            let md = triatomic(method.clone());
            let g = md.gradient();
            for atom in 0..3 {
                for k in 0..3 {
                    let mut plus = md.clone();
                    plus.positions[atom][k] += h;
                    let mut minus = md.clone();
                    minus.positions[atom][k] -= h;
                    let fd = (plus.energy() - minus.energy()) / (2.0 * h);
                    assert!(
                        (fd - g[atom][k]).abs() < 1e-6,
                        "{method:?} atom {atom} k {k}: fd {fd} vs analytic {}",
                        g[atom][k]
                    );
                }
            }
        }
    }

    #[test]
    fn test_semi_gradient_is_repulsive_derivative_at_short_range() {
        // At short range the screened Coulomb term dominates: E decreases as the
        // atoms separate, so the gradient on atom 1 (at +x) points along −x and
        // the force pushes the atoms apart.
        let md = SemiempiricalMd::new(
            SemiempiricalMethod::Xtb,
            vec![6, 6],
            vec![[0.0, 0.0, 0.0], [0.8, 0.0, 0.0]],
            vec![12.0, 12.0],
        );
        let g = md.gradient();
        let f = md.forces();
        assert!(g[1][0] < 0.0);
        assert!(f[1][0] > 0.0 && f[0][0] < 0.0);
        assert!((g[0][0] + g[1][0]).abs() < 1e-12);
    }

    #[test]
    fn test_semi_gradient_shape() {
        let md = SemiempiricalMd::new(
            SemiempiricalMethod::Pm3,
            vec![1, 1],
            vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
            vec![1.0, 1.0],
        );
        assert_eq!(md.gradient().len(), 2);
    }

    #[test]
    fn test_semi_gradient_sums_to_zero() {
        for method in ALL_METHODS {
            let g = triatomic(method).gradient();
            for k in 0..3 {
                let total: f64 = g.iter().map(|gi| gi[k]).sum();
                assert!(total.abs() < 1e-12);
            }
        }
    }
}
