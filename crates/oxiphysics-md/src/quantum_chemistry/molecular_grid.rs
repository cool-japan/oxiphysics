// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Becke-partitioned molecular integration grid for Kohn–Sham DFT.
//!
//! Every atom carries a spherical product grid:
//!
//! * **radial** — Gauss–Chebyshev quadrature of the second kind with Becke's
//!   mapping `r = r_m (1 + x)/(1 − x)`, where `r_m` is half the Bragg–Slater
//!   radius of the element (the full radius for hydrogen), A. D. Becke,
//!   J. Chem. Phys. 88, 2547 (1988);
//! * **angular** — Gauss–Legendre quadrature in `cos θ` times the uniform
//!   (trapezoidal) rule in `φ`.  With `n_θ` polar and `2 n_θ` azimuthal points
//!   the rule integrates every spherical harmonic of degree `≤ 2 n_θ − 1`
//!   exactly; it plays the role of a Lebedev grid of comparable degree but
//!   needs no tabulated coefficients.
//!
//! The atomic grids are fused with Becke's fuzzy-cell weights (three iterations
//! of the cutoff polynomial) including the heteronuclear atomic-size
//! adjustment based on Bragg–Slater radii.

use std::f64::consts::PI;

use super::error::QuantumChemistryError;

/// Ångström → Bohr.
const ANGSTROM_TO_BOHR: f64 = 1.889_726_124_565_062;

/// Resolution of the molecular grid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MolecularGridSettings {
    /// Radial points per atom.
    pub radial_points: usize,
    /// Polar (Gauss–Legendre in cos θ) points per radial shell; the azimuthal
    /// direction uses twice as many points.
    pub polar_points: usize,
}

impl Default for MolecularGridSettings {
    /// 75 radial × (14 × 28) angular points per atom: total energies of
    /// first-row molecules converged to about 10⁻⁶ Hartree.
    fn default() -> Self {
        Self {
            radial_points: 75,
            polar_points: 14,
        }
    }
}

impl MolecularGridSettings {
    /// A dense grid (150 radial × 24 × 48 angular points per atom) for
    /// reference-quality energies.
    pub fn fine() -> Self {
        Self {
            radial_points: 150,
            polar_points: 24,
        }
    }
}

/// Integration points and weights covering all of space.
#[derive(Debug, Clone, PartialEq)]
pub struct MolecularGrid {
    /// Grid points \[Bohr\].
    pub points: Vec<[f64; 3]>,
    /// Quadrature weights \[Bohr³\] (Becke partition weight included).
    pub weights: Vec<f64>,
}

/// Bragg–Slater radius \[Å\] (Slater 1964; Becke's 0.35 Å for hydrogen).
fn bragg_slater_radius(z: u32) -> f64 {
    match z {
        1 | 2 => 0.35,
        3 => 1.45,
        4 => 1.05,
        5 => 0.85,
        6 => 0.70,
        7 => 0.65,
        8 => 0.60,
        9 => 0.50,
        10 => 0.45,
        11 => 1.80,
        12 => 1.50,
        13 => 1.25,
        14 => 1.10,
        15..=17 => 1.00,
        18 => 0.95,
        _ => 1.30,
    }
}

/// Gauss–Legendre nodes and weights on \[−1, 1\] (Newton iteration on P_n).
pub(crate) fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = vec![0.0_f64; n];
    let mut weights = vec![0.0_f64; n];
    let nf = n as f64;
    for i in 0..n.div_ceil(2) {
        // Tricomi's initial guess for the i-th root.
        let mut x = (PI * (i as f64 + 0.75) / (nf + 0.5)).cos();
        let mut derivative = 1.0;
        for _ in 0..100 {
            let mut p0 = 1.0;
            let mut p1 = x;
            for k in 2..=n {
                let kf = k as f64;
                let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
                p0 = p1;
                p1 = p2;
            }
            let pn = if n == 0 { 1.0 } else { p1 };
            let pn_minus_1 = if n <= 1 { 1.0 } else { p0 };
            derivative = nf * (x * pn - pn_minus_1) / (x * x - 1.0);
            let dx = pn / derivative;
            x -= dx;
            if dx.abs() < 1e-15 {
                break;
            }
        }
        let w = 2.0 / ((1.0 - x * x) * derivative * derivative);
        nodes[i] = -x;
        weights[i] = w;
        nodes[n - 1 - i] = x;
        weights[n - 1 - i] = w;
    }
    (nodes, weights)
}

/// Unit vectors and weights (summing to 4π) of the angular product rule.
pub(crate) fn angular_rule(polar_points: usize) -> Vec<([f64; 3], f64)> {
    let (cos_theta, w_theta) = gauss_legendre(polar_points);
    let n_phi = 2 * polar_points;
    let w_phi = 2.0 * PI / n_phi as f64;
    let mut out = Vec::with_capacity(polar_points * n_phi);
    for (&ct, &wt) in cos_theta.iter().zip(w_theta.iter()) {
        let st = (1.0 - ct * ct).max(0.0).sqrt();
        for k in 0..n_phi {
            let phi = (k as f64 + 0.5) * w_phi;
            out.push(([st * phi.cos(), st * phi.sin(), ct], wt * w_phi));
        }
    }
    out
}

/// Radial points and weights (including r² and the Jacobian) for `∫₀^∞ f(r) r² dr`.
pub(crate) fn radial_rule(n: usize, r_m: f64) -> Vec<(f64, f64)> {
    let h = PI / (n as f64 + 1.0);
    (1..=n)
        .map(|i| {
            let theta = i as f64 * h;
            let x = theta.cos();
            let r = r_m * (1.0 + x) / (1.0 - x);
            // ∫₋₁¹ g(x) dx with g = f r² dr/dx, dr/dx = 2 r_m/(1 − x)²;
            // Chebyshev-II weight h sin²θ divided by √(1 − x²) = sin θ.
            let w = h * theta.sin() * 2.0 * r_m / ((1.0 - x) * (1.0 - x)) * r * r;
            (r, w)
        })
        .collect()
}

/// Becke's cutoff profile `s(ν) = ½(1 − f(f(f(ν))))`, `f(x) = 1.5x − 0.5x³`.
fn becke_cutoff(nu: f64) -> f64 {
    let mut f = nu;
    for _ in 0..3 {
        f = 1.5 * f - 0.5 * f * f * f;
    }
    0.5 * (1.0 - f)
}

impl MolecularGrid {
    /// Build the grid for nuclei with charges `nuclear_charges` at
    /// `positions` \[Bohr\].
    ///
    /// # Errors
    /// [`QuantumChemistryError::DimensionMismatch`] if the arrays differ in
    /// length, [`QuantumChemistryError::MissingData`] for no atoms or a zero
    /// grid resolution, [`QuantumChemistryError::CoincidentNuclei`] if two
    /// nuclei share a position.
    pub fn new(
        nuclear_charges: &[f64],
        positions: &[[f64; 3]],
        settings: MolecularGridSettings,
    ) -> Result<Self, QuantumChemistryError> {
        if nuclear_charges.len() != positions.len() {
            return Err(QuantumChemistryError::DimensionMismatch {
                what: "nuclear positions",
                expected: nuclear_charges.len(),
                found: positions.len(),
            });
        }
        if positions.is_empty() {
            return Err(QuantumChemistryError::MissingData("atoms"));
        }
        if settings.radial_points == 0 || settings.polar_points == 0 {
            return Err(QuantumChemistryError::MissingData("grid resolution"));
        }
        let n_atoms = positions.len();
        let radii: Vec<f64> = nuclear_charges
            .iter()
            .map(|&q| {
                let z = q.round().clamp(0.0, 200.0) as u32;
                bragg_slater_radius(z) * ANGSTROM_TO_BOHR
            })
            .collect();
        // Inter-nuclear distances and size adjustments a_ij.
        let mut inv_dist = vec![vec![0.0_f64; n_atoms]; n_atoms];
        let mut adjust = vec![vec![0.0_f64; n_atoms]; n_atoms];
        for i in 0..n_atoms {
            for j in 0..n_atoms {
                if i == j {
                    continue;
                }
                let d = distance(positions[i], positions[j]);
                if d < 1e-8 {
                    return Err(QuantumChemistryError::CoincidentNuclei {
                        first: i.min(j),
                        second: i.max(j),
                    });
                }
                inv_dist[i][j] = 1.0 / d;
                let chi = radii[i] / radii[j];
                let u = (chi - 1.0) / (chi + 1.0);
                let a = if u.abs() < 1e-14 {
                    0.0
                } else {
                    u / (u * u - 1.0)
                };
                adjust[i][j] = a.clamp(-0.5, 0.5);
            }
        }
        let angular = angular_rule(settings.polar_points);
        let mut points = Vec::new();
        let mut weights = Vec::new();
        let mut cell = vec![0.0_f64; n_atoms];
        let mut dist = vec![0.0_f64; n_atoms];
        for (a, &centre) in positions.iter().enumerate() {
            let r_m = if nuclear_charges[a].round() as i64 == 1 {
                radii[a]
            } else {
                0.5 * radii[a]
            };
            for (r, w_r) in radial_rule(settings.radial_points, r_m) {
                for &(dir, w_a) in &angular {
                    let p = [
                        centre[0] + r * dir[0],
                        centre[1] + r * dir[1],
                        centre[2] + r * dir[2],
                    ];
                    let partition = if n_atoms == 1 {
                        1.0
                    } else {
                        for (dk, pos) in dist.iter_mut().zip(positions.iter()) {
                            *dk = distance(p, *pos);
                        }
                        for i in 0..n_atoms {
                            let mut prod = 1.0;
                            for j in 0..n_atoms {
                                if i == j {
                                    continue;
                                }
                                let mu = (dist[i] - dist[j]) * inv_dist[i][j];
                                let nu = mu + adjust[i][j] * (1.0 - mu * mu);
                                prod *= becke_cutoff(nu);
                                if prod == 0.0 {
                                    break;
                                }
                            }
                            cell[i] = prod;
                        }
                        let total: f64 = cell.iter().sum();
                        if total > 0.0 { cell[a] / total } else { 0.0 }
                    };
                    let w = w_r * w_a * partition;
                    if w > 1e-15 {
                        points.push(p);
                        weights.push(w);
                    }
                }
            }
        }
        Ok(Self { points, weights })
    }

    /// Number of grid points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Whether the grid has no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Integrate a scalar function over the grid.
    pub fn integrate(&self, f: impl Fn([f64; 3]) -> f64) -> f64 {
        self.points
            .iter()
            .zip(self.weights.iter())
            .map(|(&p, &w)| w * f(p))
            .sum()
    }
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gauss_legendre_integrates_polynomials_exactly() {
        let (x, w) = gauss_legendre(7);
        assert!((w.iter().sum::<f64>() - 2.0).abs() < 1e-14);
        // ∫ x^12 = 2/13 (degree ≤ 13 is exact for 7 points)
        let v: f64 = x.iter().zip(w.iter()).map(|(x, w)| w * x.powi(12)).sum();
        assert!((v - 2.0 / 13.0).abs() < 1e-14);
    }

    #[test]
    fn angular_rule_integrates_spherical_harmonics() {
        let rule = angular_rule(10);
        let total: f64 = rule.iter().map(|(_, w)| w).sum();
        assert!((total - 4.0 * PI).abs() < 1e-12);
        // ∫ x² dΩ = 4π/3, ∫ x²y²z² dΩ = 4π/105, ∫ x y dΩ = 0
        let x2: f64 = rule.iter().map(|(d, w)| w * d[0] * d[0]).sum();
        assert!((x2 - 4.0 * PI / 3.0).abs() < 1e-12);
        let xyz2: f64 = rule
            .iter()
            .map(|(d, w)| w * (d[0] * d[1] * d[2]).powi(2))
            .sum();
        assert!((xyz2 - 4.0 * PI / 105.0).abs() < 1e-12);
        let xy: f64 = rule.iter().map(|(d, w)| w * d[0] * d[1]).sum();
        assert!(xy.abs() < 1e-13);
    }

    #[test]
    fn radial_rule_integrates_slater_density() {
        // ∫₀^∞ e^{−2r} r² dr = 1/4
        let v: f64 = radial_rule(80, 0.35 * ANGSTROM_TO_BOHR)
            .iter()
            .map(|(r, w)| w * (-2.0 * r).exp())
            .sum();
        assert!((v - 0.25).abs() < 1e-10, "{v}");
    }

    #[test]
    fn molecular_grid_integrates_atom_centred_gaussians() {
        // Normalised Gaussians on three centres: each integrates to 1.
        let pos = [
            [0.0, -0.143, 0.0],
            [1.638, 1.137, 0.0],
            [-1.638, 1.137, 0.0],
        ];
        for (settings, tol) in [
            (MolecularGridSettings::default(), 5e-5),
            (
                MolecularGridSettings {
                    radial_points: 100,
                    polar_points: 18,
                },
                2e-6,
            ),
        ] {
            let grid = MolecularGrid::new(&[8.0, 1.0, 1.0], &pos, settings).expect("grid");
            for (alpha, c) in [(130.0, pos[0]), (0.17, pos[1]), (1.2, pos[2])] {
                let norm = (alpha / PI).powf(1.5);
                let v = grid.integrate(|p| {
                    let r2 = (p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2) + (p[2] - c[2]).powi(2);
                    norm * (-alpha * r2).exp()
                });
                assert!((v - 1.0).abs() < tol, "{settings:?} alpha {alpha}: {v}");
            }
        }
    }

    #[test]
    fn coincident_nuclei_are_rejected() {
        let err = MolecularGrid::new(
            &[1.0, 1.0],
            &[[0.0; 3], [0.0; 3]],
            MolecularGridSettings::default(),
        );
        assert_eq!(
            err,
            Err(QuantumChemistryError::CoincidentNuclei {
                first: 0,
                second: 1
            })
        );
    }
}
