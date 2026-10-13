//! Tetrahedral element kinematics, nodal-force assembly and explicit-stability
//! estimates shared by [`super::TissueModel`] and [`super::SoftTissueModel`].
//!
//! All matrices are 3×3 and stored row-major in `[f64; 9]`.  A tetrahedron
//! with nodes `x0..x3` has the edge matrix `Ds = [x1 - x0, x2 - x0, x3 - x0]`
//! (edges as columns) and the deformation gradient `F = Ds * Dm⁻¹`, where `Dm`
//! is the edge matrix of the rest configuration.  The gradient of the linear
//! shape function of node `k` (`k = 1..3`) with respect to the rest position
//! is the `k`-th row of `Dm⁻¹`; node 0's gradient is minus their sum.
//!
//! For a strain-energy density `W(F)` with first Piola-Kirchhoff stress
//! `P = ∂W/∂F` the nodal forces of an element of rest volume `V0` are
//! `f_k = -V0 * P * ∇N_k`, i.e. the columns of `H = -V0 * P * Dm⁻ᵀ`.

use super::functions::{mat3_det, mat3_mul, norm3, scale3, sub3};

/// A rest tetrahedron whose `|det Dm|` is below this fraction of
/// `|e0| |e1| |e2|` (the volume of the box spanned by its edges) is treated as
/// flat: its edge matrix is not invertible in any useful sense.
const DEGENERATE_TOLERANCE: f64 = 1.0e-12;

/// Fraction of the critical explicit time step actually used when sub-stepping.
///
/// The critical step is derived from a conservative frequency bound of the
/// *rest* configuration; the margin covers the stiffening of nonlinear
/// materials and the discretisation of the damping terms.
pub(super) const STABILITY_SAFETY: f64 = 0.5;

/// Upper bound on the number of sub-steps one `step` call may take to stay
/// below the critical time step.
///
/// Tissue meshes with very stiff, tiny elements would otherwise need an
/// unbounded amount of work per frame.  When the bound is hit the step is no
/// longer guaranteed stable and the caller must lower `dt` or the stiffness
/// (see the models' `critical_time_step`).
pub(super) const MAX_STABILITY_SUBSTEPS: usize = 4096;

/// Rest-configuration data of one tetrahedral element.
#[derive(Debug, Clone, Copy)]
pub(super) struct TetRest {
    /// Rest volume `V0 = |det Dm| / 6` \[m³\].
    pub(super) volume: f64,
    /// Inverse of the rest edge matrix `Dm⁻¹` (row-major).
    pub(super) dm_inv: [f64; 9],
}

impl TetRest {
    /// Build the rest data from the four rest positions.
    ///
    /// Returns `None` for flat, degenerate or non-finite rest tetrahedra.
    pub(super) fn from_positions(x: [[f64; 3]; 4]) -> Option<Self> {
        let dm = edge_matrix(&x);
        let det = mat3_det(dm);
        let e0 = sub3(x[1], x[0]);
        let e1 = sub3(x[2], x[0]);
        let e2 = sub3(x[3], x[0]);
        let box_volume = norm3(e0) * norm3(e1) * norm3(e2);
        if !det.is_finite()
            || !box_volume.is_finite()
            || det.abs() <= DEGENERATE_TOLERANCE * box_volume
        {
            return None;
        }
        let inv_det = 1.0 / det;
        let adjugate = mat3_adjugate(dm);
        let mut dm_inv = [0.0; 9];
        for (inv, adj) in dm_inv.iter_mut().zip(adjugate.iter()) {
            *inv = adj * inv_det;
        }
        Some(Self {
            volume: det.abs() / 6.0,
            dm_inv,
        })
    }

    /// Deformation gradient `Ds * Dm⁻¹` of the vector field given at the four
    /// nodes.
    ///
    /// For the nodal positions this is the deformation gradient `F`; for nodal
    /// displacements it is the displacement gradient `∇u = F - I`; for nodal
    /// velocities it is the rate `Ḟ`.
    pub(super) fn gradient(&self, field: &[[f64; 3]; 4]) -> [f64; 9] {
        mat3_mul(edge_matrix(field), self.dm_inv)
    }

    /// Nodal forces `-V0 * P * Dm⁻ᵀ` for the stress `p` (first Piola-Kirchhoff,
    /// or the Cauchy stress in the small-strain limit).
    ///
    /// The four forces sum to zero (no net force on a stressed element).
    pub(super) fn nodal_forces(&self, p: &[f64; 9]) -> [[f64; 3]; 4] {
        let mut forces = [[0.0; 3]; 4];
        for k in 0..3 {
            let g = [
                self.dm_inv[3 * k],
                self.dm_inv[3 * k + 1],
                self.dm_inv[3 * k + 2],
            ];
            let pg = [
                p[0] * g[0] + p[1] * g[1] + p[2] * g[2],
                p[3] * g[0] + p[4] * g[1] + p[5] * g[2],
                p[6] * g[0] + p[7] * g[1] + p[8] * g[2],
            ];
            let f = scale3(pg, -self.volume);
            forces[k + 1] = f;
            forces[0] = sub3(forces[0], f);
        }
        forces
    }

    /// Gradients `∇N_k` of the four linear shape functions with respect to the
    /// rest position: row `k - 1` of `Dm⁻¹` for the nodes `k = 1..3`, and
    /// minus their sum for node 0.
    fn shape_gradients(&self) -> [[f64; 3]; 4] {
        let mut gradients = [[0.0; 3]; 4];
        for k in 0..3 {
            let row = [
                self.dm_inv[3 * k],
                self.dm_inv[3 * k + 1],
                self.dm_inv[3 * k + 2],
            ];
            gradients[k + 1] = row;
            gradients[0] = sub3(gradients[0], row);
        }
        gradients
    }

    /// Block row sums `Σ_b ‖K_ab‖_F` of the stiffness matrix of this element for
    /// an isotropic linear-elastic material with Lamé parameters `λ` and `μ`,
    /// one entry per node `a`.
    ///
    /// The 3×3 block coupling the displacements of nodes `a` and `b` is
    /// `K_ab = V0 [λ ∇N_a ∇N_bᵀ + μ ∇N_b ∇N_aᵀ + μ (∇N_a·∇N_b) I]`, and the
    /// Frobenius norm bounds the spectral norm.
    ///
    /// These sums bound the natural frequencies of a mesh. For nodal masses
    /// `m_i` the matrix `A = M⁻¹K` has blocks `A_ij = K_ij / m_i`, and
    /// `max_i Σ_j ‖A_ij‖₂` is the operator norm induced by the vector norm
    /// `max_i ‖x_i‖₂`, so it bounds the spectral radius of `A`, i.e. the
    /// largest squared natural frequency. The triangle inequality over the
    /// elements that share a node pair lets the element sums be added up per
    /// node without assembling `K`:
    /// `ω² <= max_i (Σ_{e ∋ i} Σ_{b ∈ e} ‖K^e_ib‖_F) / m_i`.
    /// Constrained nodes only remove rows and columns, which cannot raise the
    /// bound.
    ///
    /// The viscous damping matrix of a Newtonian viscosity `η` has the same
    /// structure with `λ = 0` and `μ = η / 2`.
    pub(super) fn stiffness_row_sums(&self, lame_lambda: f64, shear_modulus: f64) -> [f64; 4] {
        let gradients = self.shape_gradients();
        let mut sums = [0.0; 4];
        for (sum, ga) in sums.iter_mut().zip(gradients.iter()) {
            for gb in &gradients {
                let alignment = ga[0] * gb[0] + ga[1] * gb[1] + ga[2] * gb[2];
                let block: [f64; 9] = std::array::from_fn(|n| {
                    let (r, c) = (n / 3, n % 3);
                    let mut entry = lame_lambda * ga[r] * gb[c] + shear_modulus * gb[r] * ga[c];
                    if r == c {
                        entry += shear_modulus * alignment;
                    }
                    entry
                });
                let frobenius_sq: f64 = block.iter().map(|entry| entry * entry).sum();
                *sum += self.volume * frobenius_sq.sqrt();
            }
        }
        sums
    }
}

/// Edge matrix `[x1 - x0, x2 - x0, x3 - x0]` (edges as columns, row-major).
pub(super) fn edge_matrix(x: &[[f64; 3]; 4]) -> [f64; 9] {
    let e0 = sub3(x[1], x[0]);
    let e1 = sub3(x[2], x[0]);
    let e2 = sub3(x[3], x[0]);
    [
        e0[0], e1[0], e2[0], e0[1], e1[1], e2[1], e0[2], e1[2], e2[2],
    ]
}

/// Adjugate (transposed cofactor matrix) of a 3×3 row-major matrix, so that
/// `m * adj(m) = det(m) * I`.
pub(super) fn mat3_adjugate(m: [f64; 9]) -> [f64; 9] {
    [
        m[4] * m[8] - m[5] * m[7],
        m[2] * m[7] - m[1] * m[8],
        m[1] * m[5] - m[2] * m[4],
        m[5] * m[6] - m[3] * m[8],
        m[0] * m[8] - m[2] * m[6],
        m[2] * m[3] - m[0] * m[5],
        m[3] * m[7] - m[4] * m[6],
        m[1] * m[6] - m[0] * m[7],
        m[0] * m[4] - m[1] * m[3],
    ]
}

/// Transpose of a 3×3 row-major matrix.
pub(super) fn mat3_transpose(m: [f64; 9]) -> [f64; 9] {
    [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]]
}

/// Symmetric part `(m + mᵀ) / 2` of a 3×3 row-major matrix.
pub(super) fn mat3_symmetric_part(m: [f64; 9]) -> [f64; 9] {
    [
        m[0],
        0.5 * (m[1] + m[3]),
        0.5 * (m[2] + m[6]),
        0.5 * (m[3] + m[1]),
        m[4],
        0.5 * (m[5] + m[7]),
        0.5 * (m[6] + m[2]),
        0.5 * (m[7] + m[5]),
        m[8],
    ]
}

/// Linear-elastic stress `σ = λ tr(ε) I + 2 μ ε + σ₀ I` with `ε = sym(h)` for
/// the displacement gradient `h`, Lamé parameters `λ`, `μ` and the isotropic
/// pre-stress `σ₀`.
pub(super) fn linear_elastic_stress(
    h: &[f64; 9],
    lame_lambda: f64,
    shear_modulus: f64,
    pre_stress: f64,
) -> [f64; 9] {
    let eps = mat3_symmetric_part(*h);
    let trace = eps[0] + eps[4] + eps[8];
    let isotropic = lame_lambda * trace + pre_stress;
    let mut sigma = [0.0; 9];
    for (s, e) in sigma.iter_mut().zip(eps.iter()) {
        *s = 2.0 * shear_modulus * e;
    }
    sigma[0] += isotropic;
    sigma[4] += isotropic;
    sigma[8] += isotropic;
    sigma
}

/// Linear-elastic strain energy density `λ/2 tr(ε)² + μ tr(ε²)` of the small
/// strain `ε` (symmetric, row-major).
pub(super) fn linear_elastic_energy_density(
    strain: &[f64; 9],
    lame_lambda: f64,
    shear_modulus: f64,
) -> f64 {
    let trace = strain[0] + strain[4] + strain[8];
    let frobenius_sq: f64 = strain.iter().map(|e| e * e).sum();
    0.5 * lame_lambda * trace * trace + shear_modulus * frobenius_sq
}

/// Largest time step for which the symplectic (semi-implicit) Euler scheme
/// stays stable for a system whose highest angular frequency is `omega_max`
/// \[rad/s\] and whose damping rate at that frequency is `damping_rate` \[1/s\].
///
/// For a mode with frequency `ω` and damping rate `γ` the scheme
/// `v += h(-ω² x - γ v); x += h v` is stable iff
/// `h < 2 / (sqrt(ω² + γ²/4) + γ/2)` (the classical damped central-difference
/// limit `2 (sqrt(1 + ξ²) - ξ) / ω` with `ξ = γ / (2ω)`).  Without stiffness
/// the limit reduces to `2 / γ`; without any damping or stiffness it is
/// infinite.
pub(super) fn explicit_stability_limit(omega_max: f64, damping_rate: f64) -> f64 {
    let omega = omega_max.max(0.0);
    let gamma = damping_rate.max(0.0);
    let denominator = (omega * omega + 0.25 * gamma * gamma).sqrt() + 0.5 * gamma;
    if denominator > 0.0 {
        2.0 / denominator
    } else {
        f64::INFINITY
    }
}

/// Number of equal sub-steps needed to cover `dt` with steps no longer than
/// [`STABILITY_SAFETY`] times the critical step `critical_dt`.
///
/// Always at least 1 and at most [`MAX_STABILITY_SUBSTEPS`]; non-positive or
/// non-finite `dt` yields a single step.
pub(super) fn stability_substeps(dt: f64, critical_dt: f64) -> usize {
    if !dt.is_finite() || dt <= 0.0 {
        return 1;
    }
    let allowed = STABILITY_SAFETY * critical_dt;
    if allowed.is_nan() || allowed <= 0.0 {
        return MAX_STABILITY_SUBSTEPS;
    }
    let needed = (dt / allowed).ceil();
    if needed >= 1.0 {
        (needed as usize).min(MAX_STABILITY_SUBSTEPS)
    } else {
        1
    }
}
