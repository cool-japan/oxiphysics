// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Periodic voxel-FE homogenization of an [`RveGeometry`].
//!
//! The cubic unit cell is meshed with `n³` trilinear hex8 voxels (2×2×2 Gauss
//! quadrature).  Each voxel takes the phase (inclusion or matrix) found at its
//! centre.  For each of the six unit macro strains `E_j` the periodic
//! fluctuation `ũ` solves
//!
//! ```text
//!     K ũ = −∫ Bᵀ C E_j dV,      ũ periodic, one node pinned,
//! ```
//!
//! and the effective stiffness column is the volume average
//! `C_eff[:, j] = ⟨C (E_j + B ũ)⟩`.
//!
//! Voigt order (crate convention): `[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]` with
//! engineering shear strains.

use super::{RveGeometry, initial_isotropic_stiffness};
use crate::error::{Error, Result};
use crate::piezo::{hex8_b_matrix, hex8_shape_derivatives};
use crate::solvers::pcg_solve;
use crate::sparse::CsrMatrix;

/// Natural coordinates of the hex8 nodes (same ordering as the piezo hex8 helpers).
const HEX8_NODES: [[f64; 3]; 8] = [
    [-1.0, -1.0, -1.0],
    [1.0, -1.0, -1.0],
    [1.0, 1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [1.0, -1.0, 1.0],
    [1.0, 1.0, 1.0],
    [-1.0, 1.0, 1.0],
];
/// 2-point Gauss abscissa (weights are 1).
const GAUSS_2: [f64; 2] = [-0.577_350_269_189_625_8, 0.577_350_269_189_625_8];
/// Row permutation from the piezo B ordering `[xx, yy, zz, yz, xz, xy]` to the
/// crate Voigt ordering `[xx, yy, zz, xy, yz, xz]`.
const VOIGT_FROM_PIEZO: [usize; 6] = [0, 1, 2, 5, 3, 4];
/// Relative residual tolerance of the fluctuation solves.
const PCG_RTOL: f64 = 1e-12;

/// Result of a periodic RVE homogenization.
#[derive(Debug, Clone)]
pub struct RveHomogenization {
    /// Effective 6×6 stiffness (crate Voigt order, engineering shear).
    pub effective_stiffness: Vec<Vec<f64>>,
    /// Periodic fluctuation fields for the six unit macro strains; entry `j`
    /// holds `3·n³` nodal values (node index `ix + n·(iy + n·iz)`).
    pub fluctuations: Vec<Vec<f64>>,
    /// Inclusion volume fraction of the voxelised microstructure.
    pub inclusion_voxel_fraction: f64,
    /// Voxels per side used for the solve.
    pub n_per_side: usize,
    /// PCG iterations needed for each of the six load cases.
    pub solver_iterations: Vec<usize>,
}

impl RveHomogenization {
    /// Homogenized macro stress Σ = C_eff E.
    pub fn average_stress(&self, strain: &[f64; 6]) -> [f64; 6] {
        let mut s = [0.0; 6];
        for (i, si) in s.iter_mut().enumerate() {
            *si = (0..6)
                .map(|j| self.effective_stiffness[i][j] * strain[j])
                .sum();
        }
        s
    }

    /// Periodic fluctuation field for an arbitrary macro strain (superposition
    /// of the six unit-strain solutions).
    pub fn fluctuation(&self, strain: &[f64; 6]) -> Vec<f64> {
        let n = self.fluctuations.first().map_or(0, Vec::len);
        let mut u = vec![0.0; n];
        for (field, &e) in self.fluctuations.iter().zip(strain.iter()) {
            for (dst, &v) in u.iter_mut().zip(field.iter()) {
                *dst += e * v;
            }
        }
        u
    }
}

/// Micro-scale fields of a solved RVE at all Gauss points.
#[derive(Debug, Clone)]
pub struct RveMicroFields {
    /// Micro stress at every Gauss point (Voigt).
    pub stresses: Vec<[f64; 6]>,
    /// Micro strain `E + B ũ` at every Gauss point (Voigt).
    pub strains: Vec<[f64; 6]>,
    /// Integration weights (`det J`, Gauss weights are 1).
    pub weights: Vec<f64>,
}

/// Precomputed voxel kinematics (identical for every voxel of the grid).
struct VoxelKinematics {
    /// B matrices at the 8 Gauss points (crate Voigt order, 6 × 24).
    b: Vec<[[f64; 24]; 6]>,
    /// Gauss weight × det J.
    weight: f64,
}

impl VoxelKinematics {
    fn new(h: f64) -> Self {
        let mut b = Vec::with_capacity(8);
        for &zeta in &GAUSS_2 {
            for &eta in &GAUSS_2 {
                for &xi in &GAUSS_2 {
                    let dn_nat = hex8_shape_derivatives(xi, eta, zeta);
                    // Cube of side h: J = (h/2)·I.
                    let dndx = dn_nat.map(|d| d.map(|v| v * 2.0 / h));
                    let b_piezo = hex8_b_matrix(&dndx);
                    let mut b_voigt = [[0.0; 24]; 6];
                    for (row, &src) in b_voigt.iter_mut().zip(VOIGT_FROM_PIEZO.iter()) {
                        *row = b_piezo[src];
                    }
                    b.push(b_voigt);
                }
            }
        }
        Self {
            b,
            weight: h * h * h / 8.0,
        }
    }

    /// Element stiffness Σ Bᵀ C B w for a constant C.
    fn stiffness(&self, c: &[Vec<f64>]) -> Vec<[f64; 24]> {
        let mut ke = vec![[0.0; 24]; 24];
        for b in &self.b {
            let mut cb = [[0.0; 24]; 6];
            for r in 0..6 {
                for col in 0..24 {
                    cb[r][col] = (0..6).map(|k| c[r][k] * b[k][col]).sum();
                }
            }
            for (p, row) in ke.iter_mut().enumerate() {
                for (q, v) in row.iter_mut().enumerate() {
                    *v += (0..6).map(|k| b[k][p] * cb[k][q]).sum::<f64>() * self.weight;
                }
            }
        }
        ke
    }

    /// Σ_g Bᵀ C w (24 × 6) for a constant C.
    fn coupling(&self, c: &[Vec<f64>]) -> Vec<[f64; 6]> {
        let mut g = vec![[0.0; 6]; 24];
        for b in &self.b {
            for (p, row) in g.iter_mut().enumerate() {
                for (j, v) in row.iter_mut().enumerate() {
                    *v += (0..6).map(|k| b[k][p] * c[k][j]).sum::<f64>() * self.weight;
                }
            }
        }
        g
    }

    /// Strain B ũ_e at Gauss point `g`.
    fn strain(&self, g: usize, u_e: &[f64; 24]) -> [f64; 6] {
        let mut eps = [0.0; 6];
        for (r, e) in eps.iter_mut().enumerate() {
            *e = (0..24).map(|c| self.b[g][r][c] * u_e[c]).sum();
        }
        eps
    }
}

impl RveGeometry {
    /// Whether the point lies in the spherical inclusion of radius
    /// r = L·(3φ/4π)^{1/3} centred in the cell.
    fn in_inclusion(&self, x: f64, y: f64, z: f64) -> bool {
        let cx = self.side_length / 2.0;
        let r_incl = self.side_length
            * (3.0 * self.inclusion_fraction / (4.0 * std::f64::consts::PI)).powf(1.0 / 3.0);
        let dist = ((x - cx).powi(2) + (y - cx).powi(2) + (z - cx).powi(2)).sqrt();
        self.inclusion_fraction > 0.0 && dist <= r_incl
    }

    /// Local elastic phase `(E, ν)` at a point: inclusion or matrix.
    pub fn local_phase(&self, x: f64, y: f64, z: f64) -> (f64, f64) {
        if self.in_inclusion(x, y, z) {
            (self.e_inclusion, self.nu_inclusion)
        } else {
            (self.e_matrix, self.nu_matrix)
        }
    }

    /// Whether voxel (ix, iy, iz) belongs to the inclusion (phase at voxel centre).
    fn voxel_is_inclusion(&self, ix: usize, iy: usize, iz: usize) -> bool {
        let h = self.element_size();
        self.in_inclusion(
            (ix as f64 + 0.5) * h,
            (iy as f64 + 0.5) * h,
            (iz as f64 + 0.5) * h,
        )
    }

    /// Inclusion volume fraction of the voxelised microstructure.
    pub fn inclusion_voxel_fraction(&self) -> f64 {
        let n = self.n_per_side;
        if n == 0 {
            return 0.0;
        }
        let mut count = 0usize;
        for iz in 0..n {
            for iy in 0..n {
                for ix in 0..n {
                    if self.voxel_is_inclusion(ix, iy, iz) {
                        count += 1;
                    }
                }
            }
        }
        count as f64 / (n * n * n) as f64
    }

    /// Reuss (uniform-stress) lower bound of the effective stiffness for the
    /// nominal inclusion fraction: `C_R = (φ S_i + (1−φ) S_m)⁻¹`.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if the phase constants are not physical.
    pub fn reuss_effective_stiffness(&self) -> Result<Vec<Vec<f64>>> {
        self.validate()?;
        let phi = self.inclusion_fraction;
        let s_i = isotropic_compliance(self.e_inclusion, self.nu_inclusion);
        let s_m = isotropic_compliance(self.e_matrix, self.nu_matrix);
        // Normal block (3×3) and diagonal shear block.
        let mut s_nn = [[0.0; 3]; 3];
        for (i, row) in s_nn.iter_mut().enumerate() {
            for (j, v) in row.iter_mut().enumerate() {
                *v = phi * s_i.0[i][j] + (1.0 - phi) * s_m.0[i][j];
            }
        }
        let s_shear = phi * s_i.1 + (1.0 - phi) * s_m.1;
        let c_nn = invert_3x3(&s_nn).ok_or(Error::SingularSystem)?;
        let mut c = vec![vec![0.0; 6]; 6];
        for i in 0..3 {
            for j in 0..3 {
                c[i][j] = c_nn[i][j];
            }
            c[3 + i][3 + i] = 1.0 / s_shear;
        }
        Ok(c)
    }

    fn validate(&self) -> Result<()> {
        let phase_ok = |e: f64, nu: f64| e.is_finite() && e > 0.0 && nu > -1.0 && nu < 0.5;
        if !(self.side_length.is_finite() && self.side_length > 0.0) {
            return Err(Error::InvalidInput(format!(
                "RVE side length must be positive, got {}",
                self.side_length
            )));
        }
        if !(0.0..=1.0).contains(&self.inclusion_fraction) {
            return Err(Error::InvalidInput(format!(
                "inclusion fraction must lie in [0, 1], got {}",
                self.inclusion_fraction
            )));
        }
        if !phase_ok(self.e_inclusion, self.nu_inclusion)
            || !phase_ok(self.e_matrix, self.nu_matrix)
        {
            return Err(Error::InvalidInput(
                "RVE phase moduli must satisfy E > 0 and −1 < ν < 0.5".into(),
            ));
        }
        Ok(())
    }

    /// Global periodic node indices of voxel (ix, iy, iz) in hex8 order.
    fn voxel_nodes(&self, ix: usize, iy: usize, iz: usize) -> [usize; 8] {
        let n = self.n_per_side;
        let mut nodes = [0usize; 8];
        for (a, c) in HEX8_NODES.iter().enumerate() {
            let gx = (ix + usize::from(c[0] > 0.0)) % n;
            let gy = (iy + usize::from(c[1] > 0.0)) % n;
            let gz = (iz + usize::from(c[2] > 0.0)) % n;
            nodes[a] = gx + n * (gy + n * gz);
        }
        nodes
    }

    /// Homogenize the RVE with periodic voxel finite elements.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] for `n_per_side == 0` or non-physical constants;
    /// [`Error::General`] if a fluctuation solve does not converge.
    pub fn homogenize_periodic(&self) -> Result<RveHomogenization> {
        self.validate()?;
        let n = self.n_per_side;
        if n == 0 {
            return Err(Error::InvalidInput(
                "RVE needs at least one voxel per side".into(),
            ));
        }
        let n_nodes = n * n * n;
        let n_dof = 3 * n_nodes;
        let kin = VoxelKinematics::new(self.element_size());
        let c_inc = initial_isotropic_stiffness(self.e_inclusion, self.nu_inclusion);
        let c_mat = initial_isotropic_stiffness(self.e_matrix, self.nu_matrix);
        let ke_inc = kin.stiffness(&c_inc);
        let ke_mat = kin.stiffness(&c_mat);
        let g_inc = kin.coupling(&c_inc);
        let g_mat = kin.coupling(&c_mat);

        // Assemble the periodic stiffness and the six load vectors.  DOFs 0..3
        // (node 0) are pinned to remove the rigid translations.
        let n_free = n_dof - 3;
        let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(n_nodes * 576);
        let mut loads = vec![vec![0.0; n_free]; 6];
        let mut load_abs = vec![vec![0.0; n_free]; 6];
        let mut n_inclusion = 0usize;
        for iz in 0..n {
            for iy in 0..n {
                for ix in 0..n {
                    let inclusion = self.voxel_is_inclusion(ix, iy, iz);
                    if inclusion {
                        n_inclusion += 1;
                    }
                    let (ke, g) = if inclusion {
                        (&ke_inc, &g_inc)
                    } else {
                        (&ke_mat, &g_mat)
                    };
                    let dofs = element_dofs(&self.voxel_nodes(ix, iy, iz));
                    for p in 0..24 {
                        let Some(fp) = dofs[p].checked_sub(3) else {
                            continue;
                        };
                        for q in 0..24 {
                            if let Some(fq) = dofs[q].checked_sub(3) {
                                triplets.push((fp, fq, ke[p][q]));
                            }
                        }
                        for j in 0..6 {
                            loads[j][fp] -= g[p][j];
                            load_abs[j][fp] += g[p][j].abs();
                        }
                    }
                }
            }
        }

        let mut fluctuations = vec![vec![0.0; n_dof]; 6];
        let mut solver_iterations = vec![0usize; 6];
        if n_free > 0 {
            let k = CsrMatrix::from_triplets(n_free, n_free, &triplets);
            let mut precond = vec![0.0; n_free];
            for (i, m) in precond.iter_mut().enumerate() {
                let d = k.diagonal(i);
                if !(d.is_finite() && d > 0.0) {
                    return Err(Error::SingularSystem);
                }
                *m = 1.0 / d;
            }
            for j in 0..6 {
                let b_norm = norm(&loads[j]);
                if b_norm <= 1e-13 * norm(&load_abs[j]) {
                    // Homogeneous cell: the load cancels exactly, ũ = 0.
                    continue;
                }
                let (x, stats) = pcg_solve(&k, &loads[j], &precond, 20 * n_free + 100, PCG_RTOL);
                if !stats.converged || x.iter().any(|v| !v.is_finite()) {
                    return Err(Error::General(format!(
                        "RVE fluctuation solve {j} did not converge (residual {:e})",
                        stats.residual
                    )));
                }
                solver_iterations[j] = stats.iterations;
                fluctuations[j][3..].copy_from_slice(&x);
            }
        }

        // Effective stiffness: volume average of C (E_j + B ũ_j).
        let volume = self.volume();
        let mut c_eff = vec![vec![0.0; 6]; 6];
        for iz in 0..n {
            for iy in 0..n {
                for ix in 0..n {
                    let c = if self.voxel_is_inclusion(ix, iy, iz) {
                        &c_inc
                    } else {
                        &c_mat
                    };
                    let dofs = element_dofs(&self.voxel_nodes(ix, iy, iz));
                    for (j, field) in fluctuations.iter().enumerate() {
                        let u_e: [f64; 24] = dofs.map(|d| field[d]);
                        for g in 0..kin.b.len() {
                            let mut eps = kin.strain(g, &u_e);
                            eps[j] += 1.0;
                            for i in 0..6 {
                                c_eff[i][j] +=
                                    (0..6).map(|k| c[i][k] * eps[k]).sum::<f64>() * kin.weight;
                            }
                        }
                    }
                }
            }
        }
        for row in &mut c_eff {
            for v in row.iter_mut() {
                *v /= volume;
            }
        }

        Ok(RveHomogenization {
            effective_stiffness: c_eff,
            fluctuations,
            inclusion_voxel_fraction: n_inclusion as f64 / n_nodes as f64,
            n_per_side: n,
            solver_iterations,
        })
    }

    /// Micro stress/strain fields at every Gauss point for macro strain `strain`
    /// using a solved homogenization of this RVE.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `hom` was computed with a different resolution.
    pub fn micro_fields(
        &self,
        hom: &RveHomogenization,
        strain: &[f64; 6],
    ) -> Result<RveMicroFields> {
        let n = self.n_per_side;
        if hom.n_per_side != n || hom.fluctuations.len() != 6 {
            return Err(Error::InvalidInput(
                "homogenization result does not match this RVE".into(),
            ));
        }
        let kin = VoxelKinematics::new(self.element_size());
        let u = hom.fluctuation(strain);
        let c_inc = initial_isotropic_stiffness(self.e_inclusion, self.nu_inclusion);
        let c_mat = initial_isotropic_stiffness(self.e_matrix, self.nu_matrix);
        let n_gp = n * n * n * kin.b.len();
        let mut fields = RveMicroFields {
            stresses: Vec::with_capacity(n_gp),
            strains: Vec::with_capacity(n_gp),
            weights: Vec::with_capacity(n_gp),
        };
        for iz in 0..n {
            for iy in 0..n {
                for ix in 0..n {
                    let c = if self.voxel_is_inclusion(ix, iy, iz) {
                        &c_inc
                    } else {
                        &c_mat
                    };
                    let dofs = element_dofs(&self.voxel_nodes(ix, iy, iz));
                    let u_e: [f64; 24] = dofs.map(|d| u[d]);
                    for g in 0..kin.b.len() {
                        let mut eps = kin.strain(g, &u_e);
                        for (e, s) in eps.iter_mut().zip(strain.iter()) {
                            *e += s;
                        }
                        let mut sigma = [0.0; 6];
                        for (i, s) in sigma.iter_mut().enumerate() {
                            *s = (0..6).map(|k| c[i][k] * eps[k]).sum();
                        }
                        fields.stresses.push(sigma);
                        fields.strains.push(eps);
                        fields.weights.push(kin.weight);
                    }
                }
            }
        }
        Ok(fields)
    }
}

/// Isotropic compliance: (normal 3×3 block, shear compliance 1/G).
fn isotropic_compliance(e: f64, nu: f64) -> ([[f64; 3]; 3], f64) {
    let mut s = [[-nu / e; 3]; 3];
    for (i, row) in s.iter_mut().enumerate() {
        row[i] = 1.0 / e;
    }
    (s, 2.0 * (1.0 + nu) / e)
}

fn invert_3x3(m: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    let scale = m.iter().flatten().fold(0.0_f64, |a, v| a.max(v.abs()));
    if !det.is_finite() || det.abs() <= 1e-14 * scale.powi(3) {
        return None;
    }
    let inv_det = 1.0 / det;
    Some([
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ])
}

fn element_dofs(nodes: &[usize; 8]) -> [usize; 24] {
    let mut d = [0usize; 24];
    for (a, &node) in nodes.iter().enumerate() {
        for c in 0..3 {
            d[3 * a + c] = 3 * node + c;
        }
    }
    d
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

#[cfg(test)]
mod tests {
    use super::super::{hill_mandel_residual, volume_average_stress};
    use super::*;

    fn assert_matrix_close(a: &[Vec<f64>], b: &[Vec<f64>], tol: f64) {
        assert_eq!(a.len(), b.len());
        for (i, (ra, rb)) in a.iter().zip(b.iter()).enumerate() {
            for (j, (x, y)) in ra.iter().zip(rb.iter()).enumerate() {
                assert!((x - y).abs() < tol, "entry ({i},{j}): {x} vs {y}");
            }
        }
    }

    fn max_abs(c: &[Vec<f64>]) -> f64 {
        c.iter().flatten().fold(0.0_f64, |a, v| a.max(v.abs()))
    }

    /// Smallest eigenvalue of a symmetric 6×6 matrix (via the crate Jacobi solver).
    fn min_eigenvalue(a: &[Vec<f64>]) -> f64 {
        let (vals, _) = crate::eigen::jacobi_eigen_dense(a, 6);
        vals.first().copied().unwrap_or(f64::NAN)
    }

    fn difference(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
        a.iter()
            .zip(b.iter())
            .map(|(ra, rb)| ra.iter().zip(rb.iter()).map(|(x, y)| x - y).collect())
            .collect()
    }

    #[test]
    fn zero_inclusion_reproduces_matrix_stiffness() {
        let rve = RveGeometry::new(1.0, 4, 0.0, 200.0, 70.0, 0.3, 0.33);
        let hom = rve.homogenize_periodic().expect("homogenization");
        let c_m = initial_isotropic_stiffness(70.0, 0.33);
        let scale = max_abs(&c_m);
        assert_matrix_close(&hom.effective_stiffness, &c_m, 1e-8 * scale);
        assert_eq!(hom.inclusion_voxel_fraction, 0.0);
    }

    #[test]
    fn identical_phases_give_homogeneous_stiffness() {
        let rve = RveGeometry::new(2.0, 4, 0.3, 70.0, 70.0, 0.33, 0.33);
        let hom = rve.homogenize_periodic().expect("homogenization");
        let c_m = initial_isotropic_stiffness(70.0, 0.33);
        assert_matrix_close(&hom.effective_stiffness, &c_m, 1e-8 * 140.0);
    }

    #[test]
    fn effective_stiffness_is_symmetric_and_bounded() {
        let rve = RveGeometry::new(1.0, 6, 0.3, 400.0, 70.0, 0.25, 0.33);
        let hom = rve.homogenize_periodic().expect("homogenization");
        let c = &hom.effective_stiffness;
        let scale = max_abs(c);
        let c_t: Vec<Vec<f64>> = (0..6)
            .map(|j| c.iter().map(|row| row[j]).collect())
            .collect();
        assert_matrix_close(c, &c_t, 1e-8 * scale);
        // Bounds for the *voxelised* volume fraction.
        let mut discrete = rve.clone();
        discrete.inclusion_fraction = hom.inclusion_voxel_fraction;
        assert!(discrete.inclusion_fraction > 0.1 && discrete.inclusion_fraction < 0.5);
        let c_voigt = discrete.voigt_effective_stiffness();
        let c_reuss = discrete.reuss_effective_stiffness().expect("reuss");
        let sym = |m: &[Vec<f64>]| -> Vec<Vec<f64>> {
            (0..6)
                .map(|i| (0..6).map(|j| 0.5 * (m[i][j] + m[j][i])).collect())
                .collect()
        };
        let c_sym = sym(c);
        assert!(min_eigenvalue(&difference(&c_voigt, &c_sym)) > -1e-8 * scale);
        assert!(min_eigenvalue(&difference(&c_sym, &c_reuss)) > -1e-8 * scale);
        // Strictly inside the bounds for a genuine two-phase cell.
        assert!(c[0][0] < c_voigt[0][0] - 1e-3 * scale);
        assert!(c[0][0] > c_reuss[0][0] + 1e-3 * scale);
        assert!(c[3][3] < c_voigt[3][3] && c[3][3] > c_reuss[3][3]);
        // Cubic symmetry of a centred sphere: equal axial and shear moduli.
        assert!((c[0][0] - c[1][1]).abs() < 1e-8 * scale);
        assert!((c[1][1] - c[2][2]).abs() < 1e-8 * scale);
        assert!((c[3][3] - c[5][5]).abs() < 1e-8 * scale);
        assert!(hom.solver_iterations.iter().all(|&it| it > 0));
    }

    #[test]
    fn micro_fields_satisfy_averaging_and_hill_mandel() {
        let rve = RveGeometry::new(1.0, 4, 0.2, 300.0, 70.0, 0.3, 0.33);
        let hom = rve.homogenize_periodic().expect("homogenization");
        let strain = [1e-3, -2e-4, 5e-4, 3e-4, -1e-4, 2e-4];
        let fields = rve.micro_fields(&hom, &strain).expect("fields");
        let sigma_avg = volume_average_stress(&fields.stresses, &fields.weights);
        let sigma_macro = hom.average_stress(&strain);
        for i in 0..6 {
            assert!((sigma_avg[i] - sigma_macro[i]).abs() < 1e-9 * 300.0 * 1e-3);
        }
        // Periodic fluctuations average to zero strain.
        let mut eps_avg = [0.0; 6];
        let w_tot: f64 = fields.weights.iter().sum();
        for (e, &w) in fields.strains.iter().zip(fields.weights.iter()) {
            for i in 0..6 {
                eps_avg[i] += e[i] * w / w_tot;
            }
        }
        for i in 0..6 {
            assert!((eps_avg[i] - strain[i]).abs() < 1e-12);
        }
        let hm = hill_mandel_residual(
            &sigma_macro,
            &strain,
            &fields.stresses,
            &fields.strains,
            &fields.weights,
        );
        let work: f64 = (0..6).map(|i| sigma_macro[i] * strain[i]).sum();
        assert!(hm < 1e-8 * work.abs());
    }

    #[test]
    fn single_voxel_cell_and_invalid_input() {
        let one = RveGeometry::new(1.0, 1, 0.5, 200.0, 70.0, 0.3, 0.33);
        let hom = one.homogenize_periodic().expect("single voxel");
        // The single voxel centre lies in the inclusion.
        let c_i = initial_isotropic_stiffness(200.0, 0.3);
        assert!((hom.effective_stiffness[0][0] - c_i[0][0]).abs() < 1e-9 * c_i[0][0]);
        let zero = RveGeometry::new(1.0, 0, 0.3, 200.0, 70.0, 0.3, 0.33);
        assert!(zero.homogenize_periodic().is_err());
        let bad = RveGeometry::new(1.0, 2, 1.5, 200.0, 70.0, 0.3, 0.33);
        assert!(bad.homogenize_periodic().is_err());
    }

    #[test]
    fn local_phase_distinguishes_inclusion_and_matrix() {
        let rve = RveGeometry::new(1.0, 4, 0.2, 200.0, 70.0, 0.25, 0.33);
        assert_eq!(rve.local_phase(0.5, 0.5, 0.5), (200.0, 0.25));
        assert_eq!(rve.local_phase(0.0, 0.0, 0.0), (70.0, 0.33));
        let empty = RveGeometry::new(1.0, 4, 0.0, 200.0, 70.0, 0.25, 0.33);
        assert_eq!(empty.local_phase(0.5, 0.5, 0.5), (70.0, 0.33));
    }
}
