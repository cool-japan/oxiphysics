// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Uniform real-space grid, finite-difference operators and a free-space
//! Poisson solver for the grid Kohn-Sham DFT driver.
//!
//! Grid functions are stored row-major with the z index fastest:
//! `index = (i · n_y + j) · n_z + k`, matching the layout of
//! [`oxifft::rfft3d`].  All finite-difference operators use fourth-order
//! central stencils with homogeneous Dirichlet boundaries (the function is zero
//! outside the box), so the gradient is exactly antisymmetric and the
//! Laplacian exactly symmetric as matrices.

use oxifft::{Complex, irfft3d, rfft3d};

use super::error::QuantumChemistryMdError;

/// ∫ 1/|r| d³r over the unit cube centred on the origin:
/// 3 ln((√3 + 1)/(√3 − 1)) − π/2.
const UNIT_CUBE_COULOMB_SELF_INTEGRAL: f64 = 2.380_077_363_979_554;

/// Fourth-order first-derivative stencil weights for offsets 1 and 2 (× 1/h).
const D1_W1: f64 = 8.0 / 12.0;
const D1_W2: f64 = -1.0 / 12.0;
/// Fourth-order second-derivative stencil weights for offsets 0, 1, 2 (× 1/h²).
const D2_W0: f64 = -30.0 / 12.0;
const D2_W1: f64 = 16.0 / 12.0;
const D2_W2: f64 = -1.0 / 12.0;

/// A uniform orthogonal grid of cubic voxels.
#[derive(Debug, Clone, PartialEq)]
pub struct UniformGrid {
    /// Coordinates of grid point (0, 0, 0) \[Bohr\].
    pub origin: [f64; 3],
    /// Grid spacing h along every axis \[Bohr\].
    pub spacing: f64,
    /// Number of grid points along x, y and z.
    pub shape: [usize; 3],
}

impl UniformGrid {
    /// Create a grid from its origin, spacing and shape.
    ///
    /// # Errors
    /// [`QuantumChemistryMdError::InvalidParameter`] for a non-positive or
    /// non-finite spacing, a non-finite origin, or an axis with fewer than
    /// 5 points (the fourth-order stencils need at least 5).
    pub fn new(
        origin: [f64; 3],
        spacing: f64,
        shape: [usize; 3],
    ) -> Result<Self, QuantumChemistryMdError> {
        if !(spacing.is_finite() && spacing > 0.0) {
            return Err(QuantumChemistryMdError::InvalidParameter {
                name: "grid spacing",
                value: spacing,
            });
        }
        if let Some(&x) = origin.iter().find(|x| !x.is_finite()) {
            return Err(QuantumChemistryMdError::InvalidParameter {
                name: "grid origin",
                value: x,
            });
        }
        if let Some(&n) = shape.iter().find(|&&n| n < 5) {
            return Err(QuantumChemistryMdError::InvalidParameter {
                name: "grid points per axis",
                value: n as f64,
            });
        }
        Ok(Self {
            origin,
            spacing,
            shape,
        })
    }

    /// Create a grid whose centre is `center`.
    ///
    /// # Errors
    /// As for [`new`](Self::new).
    pub fn centered(
        center: [f64; 3],
        spacing: f64,
        shape: [usize; 3],
    ) -> Result<Self, QuantumChemistryMdError> {
        let origin = [
            center[0] - 0.5 * spacing * (shape[0].saturating_sub(1)) as f64,
            center[1] - 0.5 * spacing * (shape[1].saturating_sub(1)) as f64,
            center[2] - 0.5 * spacing * (shape[2].saturating_sub(1)) as f64,
        ];
        Self::new(origin, spacing, shape)
    }

    /// Total number of grid points.
    pub fn n_points(&self) -> usize {
        self.shape[0] * self.shape[1] * self.shape[2]
    }

    /// Volume element dV = h³ \[Bohr³\].
    pub fn dv(&self) -> f64 {
        self.spacing * self.spacing * self.spacing
    }

    /// Flat index of grid point (i, j, k).
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        (i * self.shape[1] + j) * self.shape[2] + k
    }

    /// Cartesian coordinates of the grid point with flat index `idx`.
    pub fn point(&self, idx: usize) -> [f64; 3] {
        let k = idx % self.shape[2];
        let j = (idx / self.shape[2]) % self.shape[1];
        let i = idx / (self.shape[1] * self.shape[2]);
        [
            self.origin[0] + self.spacing * i as f64,
            self.origin[1] + self.spacing * j as f64,
            self.origin[2] + self.spacing * k as f64,
        ]
    }

    /// Coordinates of all grid points in storage order.
    pub fn points(&self) -> Vec<[f64; 3]> {
        (0..self.n_points()).map(|idx| self.point(idx)).collect()
    }

    /// Stride of the flat index along `axis`.
    fn stride(&self, axis: usize) -> usize {
        match axis {
            0 => self.shape[1] * self.shape[2],
            1 => self.shape[2],
            _ => 1,
        }
    }

    /// Integer coordinate of flat index `idx` along `axis`.
    fn coordinate(&self, idx: usize, axis: usize) -> usize {
        (idx / self.stride(axis)) % self.shape[axis]
    }

    /// Fourth-order Laplacian ∇²f with zero Dirichlet boundaries.
    pub(super) fn laplacian(&self, f: &[f64]) -> Vec<f64> {
        let inv_h2 = 1.0 / (self.spacing * self.spacing);
        let mut out: Vec<f64> = f.iter().map(|&v| 3.0 * D2_W0 * v).collect();
        for axis in 0..3 {
            let stride = self.stride(axis);
            let n_axis = self.shape[axis];
            for (idx, o) in out.iter_mut().enumerate() {
                let c = self.coordinate(idx, axis);
                if c >= 1 {
                    *o += D2_W1 * f[idx - stride];
                }
                if c >= 2 {
                    *o += D2_W2 * f[idx - 2 * stride];
                }
                if c + 1 < n_axis {
                    *o += D2_W1 * f[idx + stride];
                }
                if c + 2 < n_axis {
                    *o += D2_W2 * f[idx + 2 * stride];
                }
            }
        }
        for o in &mut out {
            *o *= inv_h2;
        }
        out
    }

    /// Fourth-order central derivative ∂f/∂x_axis with zero Dirichlet boundaries.
    pub(super) fn derivative(&self, f: &[f64], axis: usize) -> Vec<f64> {
        let inv_h = 1.0 / self.spacing;
        let stride = self.stride(axis);
        let n_axis = self.shape[axis];
        (0..f.len())
            .map(|idx| {
                let c = self.coordinate(idx, axis);
                let mut d = 0.0;
                if c + 1 < n_axis {
                    d += D1_W1 * f[idx + stride];
                }
                if c + 2 < n_axis {
                    d += D1_W2 * f[idx + 2 * stride];
                }
                if c >= 1 {
                    d -= D1_W1 * f[idx - stride];
                }
                if c >= 2 {
                    d -= D1_W2 * f[idx - 2 * stride];
                }
                d * inv_h
            })
            .collect()
    }

    /// Symbol of −∇² (fourth order) for the periodic Fourier mode `m` of an
    /// axis with `n` points.
    fn minus_laplacian_symbol(&self, m: usize, n: usize) -> f64 {
        let theta = 2.0 * std::f64::consts::PI * m as f64 / n as f64;
        (30.0 - 32.0 * theta.cos() + 2.0 * (2.0 * theta).cos()) / (12.0 * self.spacing.powi(2))
    }
}

/// Error function with ~1e-15 relative accuracy.
///
/// Maclaurin series for |x| < 2, and the Lentz continued fraction of erfc for
/// |x| ≥ 2 (erf(x) = 1 to double precision beyond |x| ≈ 6).
pub(super) fn erf(x: f64) -> f64 {
    let ax = x.abs();
    let value = if ax < 2.0 {
        // erf x = 2/√π Σ (−1)^n x^{2n+1} / (n! (2n+1)).
        let x2 = ax * ax;
        let mut term = ax;
        let mut sum = ax;
        for n in 1..80 {
            let n = n as f64;
            term *= -x2 / n;
            let contribution = term / (2.0 * n + 1.0);
            sum += contribution;
            if contribution.abs() <= 1e-17 * sum.abs() {
                break;
            }
        }
        2.0 / std::f64::consts::PI.sqrt() * sum
    } else if ax < 6.5 {
        1.0 - erfc_continued_fraction(ax)
    } else {
        1.0
    };
    if x < 0.0 { -value } else { value }
}

/// erfc(x) for x ≥ 2 from the continued fraction
/// erfc x = e^{−x²}/√π · 1/(x + (1/2)/(x + 1/(x + (3/2)/(x + …)))).
fn erfc_continued_fraction(x: f64) -> f64 {
    let tiny = 1e-300;
    let mut f = x;
    let mut c = x;
    let mut d = 0.0;
    for n in 1..300 {
        let a = 0.5 * n as f64;
        d = x + a * d;
        if d.abs() < tiny {
            d = tiny;
        }
        c = x + a / c;
        if c.abs() < tiny {
            c = tiny;
        }
        d = 1.0 / d;
        let delta = c * d;
        f *= delta;
        if (delta - 1.0).abs() < 1e-16 {
            break;
        }
    }
    (-x * x).exp() / (std::f64::consts::PI.sqrt() * f)
}

/// Free-space Poisson solver on a [`UniformGrid`] (Hockney zero-padding).
///
/// Computes V(r_g) = Σ_{g'} n(r_{g'}) dV / |r_g − r_{g'}| — the potential of an
/// isolated charge distribution with no periodic images — by a linear
/// convolution evaluated with real FFTs on a grid doubled along every axis.
/// The singular self term uses the exact Coulomb integral over the voxel,
/// 2.380077·h², so a smooth density reproduces ∇²V = −4πn.
#[derive(Debug, Clone)]
pub(super) struct FreeSpacePoisson {
    grid: UniformGrid,
    padded: [usize; 3],
    kernel_hat: Vec<Complex<f64>>,
}

impl FreeSpacePoisson {
    /// Precompute the Fourier transform of the Coulomb kernel for `grid`.
    pub(super) fn new(grid: &UniformGrid) -> Self {
        let padded = [2 * grid.shape[0], 2 * grid.shape[1], 2 * grid.shape[2]];
        let h = grid.spacing;
        let dv = grid.dv();
        let mut kernel = vec![0.0f64; padded[0] * padded[1] * padded[2]];
        let wrap = |i: usize, n: usize| -> f64 {
            if i <= n / 2 {
                i as f64
            } else {
                i as f64 - n as f64
            }
        };
        for i in 0..padded[0] {
            let x = wrap(i, padded[0]);
            for j in 0..padded[1] {
                let y = wrap(j, padded[1]);
                for k in 0..padded[2] {
                    let z = wrap(k, padded[2]);
                    let r = h * (x * x + y * y + z * z).sqrt();
                    let idx = (i * padded[1] + j) * padded[2] + k;
                    kernel[idx] = if r > 0.0 {
                        dv / r
                    } else {
                        UNIT_CUBE_COULOMB_SELF_INTEGRAL * h * h
                    };
                }
            }
        }
        let kernel_hat = rfft3d(&kernel, padded[0], padded[1], padded[2]);
        Self {
            grid: grid.clone(),
            padded,
            kernel_hat,
        }
    }

    /// The grid this solver was built for.
    pub(super) fn grid(&self) -> &UniformGrid {
        &self.grid
    }

    /// Electrostatic potential of the charge density `n` (same layout as the grid).
    pub(super) fn potential(&self, n: &[f64]) -> Vec<f64> {
        let [nx, ny, nz] = self.grid.shape;
        let [px, py, pz] = self.padded;
        let mut padded = vec![0.0f64; px * py * pz];
        for i in 0..nx {
            for j in 0..ny {
                let src = (i * ny + j) * nz;
                let dst = (i * py + j) * pz;
                padded[dst..dst + nz].copy_from_slice(&n[src..src + nz]);
            }
        }
        let mut n_hat = rfft3d(&padded, px, py, pz);
        for (a, b) in n_hat.iter_mut().zip(self.kernel_hat.iter()) {
            *a = Complex::new(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
        }
        let conv = irfft3d(&n_hat, px, py, pz);
        let mut v = vec![0.0f64; nx * ny * nz];
        for i in 0..nx {
            for j in 0..ny {
                let src = (i * py + j) * pz;
                let dst = (i * ny + j) * nz;
                v[dst..dst + nz].copy_from_slice(&conv[src..src + nz]);
            }
        }
        v
    }
}

/// Kinetic-energy preconditioner (−½∇²_periodic + s)⁻¹ applied with FFTs.
///
/// Uses the Fourier symbol of the same fourth-order Laplacian as the
/// Hamiltonian, so high-frequency residual components — which dominate the
/// condition number of the grid Hamiltonian — are damped by 1/(k²/2).
#[derive(Debug, Clone)]
pub(super) struct KineticPreconditioner {
    grid: UniformGrid,
    kinetic_symbol: Vec<f64>,
}

impl KineticPreconditioner {
    pub(super) fn new(grid: &UniformGrid) -> Self {
        let [nx, ny, nz] = grid.shape;
        let nz_half = nz / 2 + 1;
        let sx: Vec<f64> = (0..nx)
            .map(|m| grid.minus_laplacian_symbol(m, nx))
            .collect();
        let sy: Vec<f64> = (0..ny)
            .map(|m| grid.minus_laplacian_symbol(m, ny))
            .collect();
        let sz: Vec<f64> = (0..nz_half)
            .map(|m| grid.minus_laplacian_symbol(m, nz))
            .collect();
        let mut kinetic_symbol = Vec::with_capacity(nx * ny * nz_half);
        for &a in &sx {
            for &b in &sy {
                for &c in &sz {
                    kinetic_symbol.push(0.5 * (a + b + c));
                }
            }
        }
        Self {
            grid: grid.clone(),
            kinetic_symbol,
        }
    }

    pub(super) fn grid(&self) -> &UniformGrid {
        &self.grid
    }

    /// Apply (T + shift)⁻¹ to `r`; `shift` must be positive.
    pub(super) fn apply(&self, r: &[f64], shift: f64) -> Vec<f64> {
        let [nx, ny, nz] = self.grid.shape;
        let mut r_hat = rfft3d(r, nx, ny, nz);
        for (c, &t) in r_hat.iter_mut().zip(self.kinetic_symbol.iter()) {
            let scale = 1.0 / (t + shift);
            *c = Complex::new(c.re * scale, c.im * scale);
        }
        irfft3d(&r_hat, nx, ny, nz)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_erf_reference_values() {
        // Abramowitz & Stegun Table 7.1 / high-precision references.
        let cases = [
            (0.0, 0.0),
            (0.1, 0.112_462_916_018_284_9),
            (0.5, 0.520_499_877_813_046_5),
            (1.0, 0.842_700_792_949_714_9),
            (1.9, 0.992_790_429_235_257_5),
            (2.0, 0.995_322_265_018_952_7),
            (3.0, 0.999_977_909_503_001_4),
            (-1.5, -0.966_105_146_475_310_7),
        ];
        for (x, expected) in cases {
            assert!(
                (erf(x) - expected).abs() < 2e-15,
                "erf({x}) = {} vs {expected}",
                erf(x)
            );
        }
        assert_eq!(erf(10.0), 1.0);
    }

    #[test]
    fn test_grid_point_indexing_roundtrip() {
        let grid = UniformGrid::new([1.0, -2.0, 0.5], 0.25, [5, 6, 7]).expect("grid");
        let idx = grid.index(3, 4, 5);
        let p = grid.point(idx);
        assert!((p[0] - 1.75).abs() < 1e-15);
        assert!((p[1] - (-1.0)).abs() < 1e-15);
        assert!((p[2] - 1.75).abs() < 1e-15);
        assert_eq!(grid.points().len(), 210);
        assert!(UniformGrid::new([0.0; 3], 0.0, [8, 8, 8]).is_err());
        assert!(UniformGrid::new([0.0; 3], 0.2, [8, 4, 8]).is_err());
    }

    #[test]
    fn test_laplacian_of_gaussian_and_symmetry() {
        let grid = UniformGrid::centered([0.0; 3], 0.15, [24, 24, 24]).expect("grid");
        let alpha = 1.5;
        let f: Vec<f64> = grid
            .points()
            .iter()
            .map(|p| (-alpha * (p[0] * p[0] + p[1] * p[1] + p[2] * p[2])).exp())
            .collect();
        let lap = grid.laplacian(&f);
        // Analytic: (4α²r² − 6α) e^{−αr²}; check near the centre.
        for (idx, p) in grid.points().iter().enumerate() {
            let r2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
            if r2 < 1.0 {
                let exact = (4.0 * alpha * alpha * r2 - 6.0 * alpha) * (-alpha * r2).exp();
                assert!(
                    (lap[idx] - exact).abs() < 1e-3 * 6.0 * alpha,
                    "{} vs {exact}",
                    lap[idx]
                );
            }
        }
        // ⟨g|∇²f⟩ = ⟨f|∇²g⟩ and ⟨g|∂f⟩ = −⟨f|∂g⟩ for arbitrary grid functions.
        let g: Vec<f64> = (0..grid.n_points())
            .map(|i| ((i * 7919) % 101) as f64 / 101.0 - 0.5)
            .collect();
        let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
        let lg = grid.laplacian(&g);
        assert!((dot(&g, &lap) - dot(&f, &lg)).abs() < 1e-9 * dot(&g, &lap).abs().max(1.0));
        for axis in 0..3 {
            let df = grid.derivative(&f, axis);
            let dg = grid.derivative(&g, axis);
            assert!((dot(&g, &df) + dot(&f, &dg)).abs() < 1e-9);
        }
    }

    #[test]
    fn test_poisson_gaussian_charge_matches_erf_potential() {
        // Normalised Gaussian charge ρ = (α/π)^{3/2} e^{−αr²} has the potential
        // erf(√α r)/r.
        let grid = UniformGrid::centered([0.0; 3], 0.2, [28, 28, 28]).expect("grid");
        let alpha: f64 = 1.2;
        let norm = (alpha / std::f64::consts::PI).powf(1.5);
        let points = grid.points();
        let rho: Vec<f64> = points
            .iter()
            .map(|p| norm * (-alpha * (p[0] * p[0] + p[1] * p[1] + p[2] * p[2])).exp())
            .collect();
        let poisson = FreeSpacePoisson::new(&grid);
        let v = poisson.potential(&rho);
        let mut max_err: f64 = 0.0;
        for (idx, p) in points.iter().enumerate() {
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            let exact = if r > 1e-12 {
                erf(alpha.sqrt() * r) / r
            } else {
                2.0 * (alpha / std::f64::consts::PI).sqrt()
            };
            max_err = max_err.max((v[idx] - exact).abs());
        }
        // Midpoint-rule discretisation of the singular kernel: < 0.5 % of the
        // peak potential 2√(α/π) ≈ 1.24.
        assert!(max_err < 6e-3, "max error {max_err}");
    }

    #[test]
    fn test_preconditioner_inverts_kinetic_plus_shift_on_plane_wave() {
        let grid = UniformGrid::new([0.0; 3], 0.3, [8, 10, 12]).expect("grid");
        let pre = KineticPreconditioner::new(&grid);
        let shift = 0.7;
        // A periodic plane wave along z is an eigenfunction of the periodic
        // operator with eigenvalue ½ symbol(m) + shift.
        let m = 2;
        let r: Vec<f64> = (0..grid.n_points())
            .map(|idx| {
                let k = idx % 12;
                (2.0 * std::f64::consts::PI * (m * k) as f64 / 12.0).cos()
            })
            .collect();
        let out = pre.apply(&r, shift);
        let lambda = 0.5 * grid.minus_laplacian_symbol(m, 12) + shift;
        for (o, x) in out.iter().zip(r.iter()) {
            assert!((o * lambda - x).abs() < 1e-10);
        }
    }
}
