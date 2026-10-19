// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Sparse-grid MLS-MPM kernel for dry, cohesionless sand.
//!
//! The scheme is the moving-least-squares MPM of Hu et al. (2018) with
//! quadratic B-spline weights and an APIC affine velocity per particle, driven
//! by the Hencky elasticity and the Drucker-Prager return mapping of
//! Klár et al. (2016, "Drucker-Prager elastoplasticity for sand animation").
//!
//! One sub-step (`h` = cell size, `dt` = time step):
//!
//! 1. **Particle to grid.** Every particle with elastic deformation gradient
//!    `F = U Σ Vᵀ` has the Hencky strain `ε = ln Σ` and the Kirchhoff stress
//!    `τ = U (2 μ ε + λ tr ε) Uᵀ`. With `w_ip` the quadratic B-spline weight and
//!    `d_ip = x_i − x_p`, node `i` receives
//!    `m_i += w_ip m_p` and
//!    `(m v)_i += w_ip (m_p v_p + (m_p C_p − 4 dt V⁰_p τ_p / h²) d_ip)`.
//!    The factor `4 / h²` is `D⁻¹`, the inverse of the second moment
//!    `Σ_i w_ip d_ip d_ipᵀ = h²/4 · I` of the quadratic B-spline.
//! 2. **Grid update.** `v_i = (m v)_i / m_i + dt g`, then the contact with the
//!    floor plane `y = 0` (nodes with `y_i ≤ 0`): the normal velocity is
//!    clamped to be non-negative and Coulomb friction acts on the tangential
//!    velocity with the impulse that removed the normal motion.
//! 3. **Grid to particle.** `v_p = Σ_i w_ip v_i`,
//!    `C_p = 4/h² Σ_i w_ip v_i ⊗ d_ip`, `F ← (I + dt C_p) F`, followed by the
//!    return mapping below, and finally `x_p ← x_p + dt v_p`.
//!
//! **Return mapping** (`ε = ln Σ`, `tr ε = Σ ε_k`, `ε̂ = ε − tr ε / 3`):
//!
//! * `tr ε ≥ 0` (tension): `F` is set to its rotation `U Vᵀ`; the material is
//!   stress free. The dilation `tr ε` is accounted as plastic volume change.
//! * otherwise `δγ = ‖ε̂‖ + (3λ + 2μ)/(2μ) · α · tr ε` with
//!   `α = √(2/3) · 2 sin φ / (3 − sin φ)`; `δγ ≤ 0` means elastic, otherwise
//!   `H = ε − δγ ε̂ / ‖ε̂‖` and `F = U exp(H) Vᵀ`.
//!
//! The singular value decomposition is never formed: only `F Fᵀ = U Σ² Uᵀ` is
//! eigen-decomposed (cyclic Jacobi), and since
//! `U exp(H) Vᵀ = [U diag(e^H / σ) Uᵀ] F` the projected `F` follows without `V`.
//!
//! The grid is **sparse**: only the nodes touched by particles are stored (in a
//! hash map with a deterministic hasher; all loops run over the insertion-ordered
//! node vector, so results do not depend on hash iteration order), so the domain
//! is unbounded and memory is proportional to the number of particles.

use super::SandParticle;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

/// 3-vector.
pub(super) type Vec3 = [f64; 3];
/// Row-major 3×3 matrix, `m[row][column]`.
pub(super) type Mat3 = [[f64; 3]; 3];

/// 3×3 identity.
pub(super) const IDENTITY: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
/// 3×3 zero matrix.
pub(super) const ZERO: Mat3 = [[0.0; 3]; 3];

/// Offset that maps cell indices to the non-negative 21-bit key components.
const KEY_BIAS: i64 = 1 << 20;
/// Largest `|cell index|` of a particle; the 3-cell stencil stays inside the key range.
const MAX_CELL_INDEX: f64 = (KEY_BIAS - 4) as f64;
/// Smallest principal stretch used by the constitutive law (guards `ln 0`).
const SIGMA_MIN: f64 = 1.0e-3;

// ---------------------------------------------------------------------------
// Small matrix helpers
// ---------------------------------------------------------------------------

#[inline]
pub(super) fn mat_mul(a: &Mat3, b: &Mat3) -> Mat3 {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j])
    })
}

#[inline]
pub(super) fn mat_vec(a: &Mat3, v: Vec3) -> Vec3 {
    std::array::from_fn(|i| a[i][0] * v[0] + a[i][1] * v[1] + a[i][2] * v[2])
}

#[inline]
pub(super) fn mat_transpose(a: &Mat3) -> Mat3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[j][i]))
}

#[inline]
fn mat_add(a: &Mat3, b: &Mat3) -> Mat3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][j] + b[i][j]))
}

#[inline]
fn mat_scale(a: &Mat3, s: f64) -> Mat3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][j] * s))
}

#[cfg(test)]
#[inline]
pub(super) fn mat_det(a: &Mat3) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

/// `true` if every entry of the matrix is finite.
#[inline]
fn mat_is_finite(a: &Mat3) -> bool {
    a.iter().all(|row| row.iter().all(|x| x.is_finite()))
}

/// Matrix from a column-major flat array (the layout of [`SandParticle`]).
#[inline]
pub(super) fn mat_from_col_major(a: &[f64; 9]) -> Mat3 {
    [[a[0], a[3], a[6]], [a[1], a[4], a[7]], [a[2], a[5], a[8]]]
}

/// Column-major flat array of a matrix.
#[inline]
pub(super) fn mat_to_col_major(m: &Mat3) -> [f64; 9] {
    [
        m[0][0], m[1][0], m[2][0], m[0][1], m[1][1], m[2][1], m[0][2], m[1][2], m[2][2],
    ]
}

// ---------------------------------------------------------------------------
// Symmetric eigen-decomposition
// ---------------------------------------------------------------------------

/// Eigen-decomposition `A = V diag(λ) Vᵀ` of a symmetric 3×3 matrix.
///
/// Cyclic Jacobi rotations; the eigenvectors are the **columns** of the returned
/// matrix and are orthonormal to machine precision. The eigenvalues come in no
/// particular order. Only the upper triangle of `a` is read.
pub(super) fn sym_eigen3(a: &Mat3) -> ([f64; 3], Mat3) {
    let mut m = [
        [a[0][0], a[0][1], a[0][2]],
        [a[0][1], a[1][1], a[1][2]],
        [a[0][2], a[1][2], a[2][2]],
    ];
    let mut v = IDENTITY;
    for _ in 0..32 {
        let off = m[0][1] * m[0][1] + m[0][2] * m[0][2] + m[1][2] * m[1][2];
        let diag = m[0][0] * m[0][0] + m[1][1] * m[1][1] + m[2][2] * m[2][2];
        if off <= 1.0e-30 * (diag + off) {
            break;
        }
        for (p, q) in [(0usize, 1usize), (0, 2), (1, 2)] {
            let apq = m[p][q];
            if apq == 0.0 {
                continue;
            }
            // tan of the rotation angle: the smaller root of t² + 2 t θ − 1 = 0.
            let theta = (m[q][q] - m[p][p]) / (2.0 * apq);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            let r = 3 - p - q;
            let (app, aqq) = (m[p][p], m[q][q]);
            m[p][p] = app - t * apq;
            m[q][q] = aqq + t * apq;
            m[p][q] = 0.0;
            m[q][p] = 0.0;
            let (arp, arq) = (m[r][p], m[r][q]);
            m[r][p] = c * arp - s * arq;
            m[r][q] = s * arp + c * arq;
            m[p][r] = m[r][p];
            m[q][r] = m[r][q];
            for row in v.iter_mut() {
                let (vp, vq) = (row[p], row[q]);
                row[p] = c * vp - s * vq;
                row[q] = s * vp + c * vq;
            }
        }
    }
    ([m[0][0], m[1][1], m[2][2]], v)
}

// ---------------------------------------------------------------------------
// Material
// ---------------------------------------------------------------------------

/// Elastic and plastic constants of the sand, derived once per step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Material {
    /// Lamé shear modulus μ \[Pa\].
    pub mu: f64,
    /// Lamé constant λ \[Pa\].
    pub lambda: f64,
    /// Drucker-Prager cone constant `α = √(2/3) · 2 sin φ / (3 − sin φ)`.
    pub alpha: f64,
}

impl Material {
    /// Builds the constants from Young's modulus \[Pa\], Poisson's ratio and the
    /// friction angle \[degrees\]; `None` if any value is outside its range
    /// (`E > 0`, `0 ≤ ν < 0.5`, `0 ≤ φ < 90°`) or not finite.
    pub(super) fn new(
        youngs_modulus_pa: f64,
        poisson_ratio: f64,
        friction_angle_deg: f64,
    ) -> Option<Self> {
        let valid = youngs_modulus_pa.is_finite()
            && youngs_modulus_pa > 0.0
            && (0.0..0.5).contains(&poisson_ratio)
            && friction_angle_deg.is_finite()
            && (0.0..90.0).contains(&friction_angle_deg);
        if !valid {
            return None;
        }
        let mu = youngs_modulus_pa / (2.0 * (1.0 + poisson_ratio));
        let lambda = youngs_modulus_pa * poisson_ratio
            / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
        let sin_phi = friction_angle_deg.to_radians().sin();
        let alpha = (2.0_f64 / 3.0).sqrt() * 2.0 * sin_phi / (3.0 - sin_phi);
        Some(Self { mu, lambda, alpha })
    }

    /// Constrained (P-wave) modulus `λ + 2μ` \[Pa\].
    pub(super) fn p_wave_modulus(&self) -> f64 {
        self.lambda + 2.0 * self.mu
    }
}

/// Principal decomposition of a deformation gradient: `F Fᵀ = U diag(σ²) Uᵀ`.
#[derive(Debug, Clone, Copy)]
pub(super) struct Decomposition {
    /// Left singular vectors (columns).
    pub u: Mat3,
    /// Principal stretches σ (clamped below by `SIGMA_MIN`).
    pub sigma: [f64; 3],
    /// Hencky strain `ε = ln σ`.
    pub eps: [f64; 3],
}

/// Principal decomposition of `F` via the eigenvalues of `F Fᵀ`.
pub(super) fn decompose(f: &Mat3) -> Decomposition {
    let b = mat_mul(f, &mat_transpose(f));
    let (lam, u) = sym_eigen3(&b);
    let floor = SIGMA_MIN * SIGMA_MIN;
    let sigma = lam.map(|l| l.max(floor).sqrt());
    let eps = sigma.map(f64::ln);
    Decomposition { u, sigma, eps }
}

/// `U diag(d) Uᵀ`.
fn rotate_diagonal(u: &Mat3, d: [f64; 3]) -> Mat3 {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            u[i][0] * d[0] * u[j][0] + u[i][1] * d[1] * u[j][1] + u[i][2] * d[2] * u[j][2]
        })
    })
}

/// Kirchhoff stress `τ = U (2 μ ε + λ tr ε) Uᵀ` of the Hencky model.
pub(super) fn kirchhoff_stress(f: &Mat3, mat: &Material) -> Mat3 {
    let dec = decompose(f);
    let tr = dec.eps[0] + dec.eps[1] + dec.eps[2];
    let d = dec.eps.map(|e| 2.0 * mat.mu * e + mat.lambda * tr);
    rotate_diagonal(&dec.u, d)
}

/// Result of the Drucker-Prager return mapping in principal log-strain space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct StrainProjection {
    /// Projected Hencky strain `H`.
    pub strain: [f64; 3],
    /// Volumetric plastic strain increment `tr ε − tr H` (≥ 0).
    pub plastic_volume: f64,
    /// `false` if the trial state was inside the yield surface (`H = ε`).
    pub yielded: bool,
}

/// Klár et al. (2016) return mapping of the trial Hencky strain `ε` onto the
/// Drucker-Prager cone (cohesionless).
pub(super) fn project_log_strain(eps: [f64; 3], mat: &Material) -> StrainProjection {
    let tr = eps[0] + eps[1] + eps[2];
    if tr >= 0.0 {
        return StrainProjection {
            strain: [0.0; 3],
            plastic_volume: tr,
            yielded: true,
        };
    }
    let mean = tr / 3.0;
    let dev = eps.map(|e| e - mean);
    let dev_norm = (dev[0] * dev[0] + dev[1] * dev[1] + dev[2] * dev[2]).sqrt();
    let elastic = StrainProjection {
        strain: eps,
        plastic_volume: 0.0,
        yielded: false,
    };
    if dev_norm <= 1.0e-14 {
        return elastic;
    }
    let delta_gamma =
        dev_norm + (3.0 * mat.lambda + 2.0 * mat.mu) / (2.0 * mat.mu) * mat.alpha * tr;
    if delta_gamma <= 0.0 {
        return elastic;
    }
    StrainProjection {
        strain: std::array::from_fn(|k| eps[k] - delta_gamma * dev[k] / dev_norm),
        plastic_volume: 0.0,
        yielded: true,
    }
}

/// Applies the return mapping to a trial deformation gradient.
///
/// Returns the corrected gradient and the volumetric plastic strain
/// increment. The gradient is returned bit-for-bit unchanged while the state is
/// elastic; a non-finite result falls back to the stress-free identity.
pub(super) fn project_deformation(f: &Mat3, mat: &Material) -> (Mat3, f64) {
    let dec = decompose(f);
    let proj = project_log_strain(dec.eps, mat);
    if !proj.yielded {
        return (*f, 0.0);
    }
    // F' = U exp(H) Vᵀ = [U diag(exp(H)/σ) Uᵀ] F.
    let scale: [f64; 3] = std::array::from_fn(|k| proj.strain[k].exp() / dec.sigma[k]);
    let corrected = mat_mul(&rotate_diagonal(&dec.u, scale), f);
    if mat_is_finite(&corrected) {
        (corrected, proj.plastic_volume)
    } else {
        (IDENTITY, proj.plastic_volume)
    }
}

// ---------------------------------------------------------------------------
// Quadratic B-spline stencil
// ---------------------------------------------------------------------------

/// Quadratic B-spline weights of the three nodes `base, base+1, base+2` for the
/// fractional position `fx ∈ [0.5, 1.5)` (in cells, measured from `base`).
///
/// They are non-negative, sum to one, and `Σ_k w_k (k − fx) = 0`.
#[inline]
pub(super) fn quadratic_weights(fx: f64) -> [f64; 3] {
    [
        0.5 * (1.5 - fx) * (1.5 - fx),
        0.75 - (fx - 1.0) * (fx - 1.0),
        0.5 * (fx - 0.5) * (fx - 0.5),
    ]
}

/// The 3×3×3 node block a particle interacts with.
#[derive(Debug, Clone, Copy)]
pub(super) struct Stencil {
    /// Cell index of the lower corner node; the block is `base + (a, b, c)`.
    pub base: [i64; 3],
    /// Fractional position in cells relative to `base`, in `[0.5, 1.5)`.
    pub fx: [f64; 3],
    /// Per-axis B-spline weights.
    pub w: [[f64; 3]; 3],
}

impl Stencil {
    /// Stencil of a particle at `pos` on a grid with inverse spacing `inv_h`;
    /// `None` for a non-finite position or one outside the supported index range.
    pub(super) fn new(pos: Vec3, inv_h: f64) -> Option<Self> {
        let x = pos.map(|p| p * inv_h);
        if !x.iter().all(|v| v.is_finite() && v.abs() < MAX_CELL_INDEX) {
            return None;
        }
        let lower = x.map(|v| (v - 0.5).floor());
        let fx: [f64; 3] = std::array::from_fn(|d| x[d] - lower[d]);
        Some(Self {
            base: lower.map(|v| v as i64),
            fx,
            w: fx.map(quadratic_weights),
        })
    }

    /// Weight and offset `x_i − x_p` (in units of `h`) of block node `(a, b, c)`.
    #[inline]
    pub(super) fn node(&self, a: usize, b: usize, c: usize) -> (f64, Vec3) {
        let w = self.w[0][a] * self.w[1][b] * self.w[2][c];
        (
            w,
            [
                a as f64 - self.fx[0],
                b as f64 - self.fx[1],
                c as f64 - self.fx[2],
            ],
        )
    }
}

// ---------------------------------------------------------------------------
// Sparse grid
// ---------------------------------------------------------------------------

/// Deterministic mixing hasher for packed node keys (murmur3 finalizer).
#[derive(Default)]
struct KeyHasher(u64);

impl Hasher for KeyHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0 ^ u64::from(b)).wrapping_mul(0x0000_0100_0000_01b3);
        }
    }

    fn write_u64(&mut self, key: u64) {
        let mut x = key;
        x ^= x >> 33;
        x = x.wrapping_mul(0xff51_afd7_ed55_8ccd);
        x ^= x >> 33;
        x = x.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
        x ^= x >> 33;
        self.0 = x;
    }
}

type KeyMap = HashMap<u64, usize, BuildHasherDefault<KeyHasher>>;

/// Packs a cell index triple into a 63-bit key.
#[inline]
fn node_key(cell: [i64; 3]) -> u64 {
    // Each component lies in [0, 2^21): |index| < 2^20 is guaranteed by `MAX_CELL_INDEX`.
    (((cell[0] + KEY_BIAS) as u64) << 42)
        | (((cell[1] + KEY_BIAS) as u64) << 21)
        | ((cell[2] + KEY_BIAS) as u64)
}

/// Sparse background grid: mass and velocity (momentum before the grid update)
/// of the nodes touched by particles, stored in insertion order.
pub(super) struct SparseGrid {
    index: KeyMap,
    /// Cell index of every node.
    pub cells: Vec<[i64; 3]>,
    /// Node mass.
    pub mass: Vec<f64>,
    /// Node momentum during P2G, node velocity after [`SparseGrid::update`].
    pub vel: Vec<Vec3>,
    /// Base cell and node indices of the most recently registered block;
    /// neighbouring particles often share it.
    last_block: Option<([i64; 3], [usize; 27])>,
}

impl SparseGrid {
    /// Empty grid with room for roughly `capacity` nodes.
    pub(super) fn with_capacity(capacity: usize) -> Self {
        Self {
            index: KeyMap::with_capacity_and_hasher(capacity, BuildHasherDefault::default()),
            cells: Vec::with_capacity(capacity),
            mass: Vec::with_capacity(capacity),
            vel: Vec::with_capacity(capacity),
            last_block: None,
        }
    }

    /// Index of the node at `cell`, created on first use.
    fn node(&mut self, cell: [i64; 3]) -> usize {
        let next = self.cells.len();
        let id = *self.index.entry(node_key(cell)).or_insert(next);
        if id == next {
            self.cells.push(cell);
            self.mass.push(0.0);
            self.vel.push([0.0; 3]);
        }
        id
    }

    /// Registers the 27 nodes of `stencil` and returns their indices in the
    /// order `(a, b, c)` = `(0,0,0), (0,0,1), …, (2,2,2)`.
    pub(super) fn register(&mut self, stencil: &Stencil) -> [usize; 27] {
        if let Some((base, ids)) = &self.last_block
            && *base == stencil.base
        {
            return *ids;
        }
        let mut ids = [0usize; 27];
        let mut slot = 0;
        for a in 0..3i64 {
            for b in 0..3i64 {
                for c in 0..3i64 {
                    ids[slot] = self.node([
                        stencil.base[0] + a,
                        stencil.base[1] + b,
                        stencil.base[2] + c,
                    ]);
                    slot += 1;
                }
            }
        }
        self.last_block = Some((stencil.base, ids));
        ids
    }

    /// Number of nodes.
    #[cfg(test)]
    pub(super) fn len(&self) -> usize {
        self.cells.len()
    }

    /// Converts momentum to velocity, adds gravity and applies the boundary
    /// conditions of `c`: the floor contact (`y_i ≤ 0`, Coulomb friction) and,
    /// if set, the symmetry planes of the plane-strain slab (zero `z` velocity
    /// on and beyond them).
    pub(super) fn update(&mut self, c: &StepConstants) {
        let tolerance = 1.0e-9 * c.h;
        for ((mass, vel), cell) in self.mass.iter().zip(self.vel.iter_mut()).zip(&self.cells) {
            if *mass > 0.0 {
                let inv_m = 1.0 / *mass;
                *vel = [
                    vel[0] * inv_m + c.dt * c.gravity[0],
                    vel[1] * inv_m + c.dt * c.gravity[1],
                    vel[2] * inv_m + c.dt * c.gravity[2],
                ];
                if cell[1] <= 0 {
                    floor_contact(vel, c.floor_friction);
                }
                if let Some([z_min, z_max]) = c.plane_strain_z {
                    let z = cell[2] as f64 * c.h;
                    if z <= z_min + tolerance || z >= z_max - tolerance {
                        vel[2] = 0.0;
                    }
                }
            } else {
                *vel = [0.0; 3];
            }
        }
    }
}

/// Floor contact of a node velocity: the normal (y) velocity may not point into
/// the floor and Coulomb friction `mu` (the impulse that removed the normal
/// motion times `mu`) slows the tangential velocity, down to rest. `mu = ∞` is
/// a sticky floor.
pub(super) fn floor_contact(v: &mut Vec3, mu: f64) {
    if v[1] >= 0.0 {
        return;
    }
    let impulse = -v[1];
    v[1] = 0.0;
    let tangential = v[0].hypot(v[2]);
    if tangential > 0.0 {
        let removable = mu * impulse;
        if removable >= tangential {
            v[0] = 0.0;
            v[2] = 0.0;
        } else {
            let keep = 1.0 - removable / tangential;
            v[0] *= keep;
            v[2] *= keep;
        }
    }
}

// ---------------------------------------------------------------------------
// One sub-step
// ---------------------------------------------------------------------------

/// Per-sub-step constants.
#[derive(Debug, Clone, Copy)]
pub(super) struct StepConstants {
    /// Grid cell size \[m\].
    pub h: f64,
    /// Time step \[s\].
    pub dt: f64,
    /// Gravity \[m/s²\].
    pub gravity: Vec3,
    /// Constitutive constants.
    pub material: Material,
    /// Coulomb friction coefficient of the floor plane `y = 0`.
    pub floor_friction: f64,
    /// Optional plane-strain slab `[z_min, z_max]` (symmetry planes).
    pub plane_strain_z: Option<[f64; 2]>,
}

/// `true` if the particle state is finite with a positive mass and volume.
fn state_is_valid(p: &SandParticle) -> bool {
    p.mass.is_finite()
        && p.mass > 0.0
        && p.volume.is_finite()
        && p.volume > 0.0
        && p.position.iter().all(|v| v.is_finite())
        && p.velocity.iter().all(|v| v.is_finite())
        && p.deformation_gradient.iter().all(|v| v.is_finite())
        && p.affine_velocity.iter().all(|v| v.is_finite())
        && p.plastic_strain_vol.is_finite()
}

/// `true` if the particle can take part in the simulation: finite state,
/// positive mass and volume, and a position inside the supported index range.
pub(super) fn is_active(p: &SandParticle, inv_h: f64) -> bool {
    state_is_valid(p) && Stencil::new(p.position, inv_h).is_some()
}

/// Stencil and node indices of an active particle.
struct Prepared {
    stencil: Stencil,
    nodes: [usize; 27],
}

/// Advances all active particles by one explicit MLS-MPM sub-step.
///
/// Inactive particles (see [`is_active`]) are left untouched and do not
/// influence the others.
pub(super) fn substep(particles: &mut [SandParticle], c: &StepConstants) {
    let inv_h = 1.0 / c.h;
    let inv_h2 = inv_h * inv_h;
    let mut grid = SparseGrid::with_capacity(particles.len().saturating_mul(4));

    // Stencils and node indices of the active particles.
    let prepared: Vec<Option<Prepared>> = particles
        .iter()
        .map(|p| {
            if !state_is_valid(p) {
                return None;
            }
            let stencil = Stencil::new(p.position, inv_h)?;
            let nodes = grid.register(&stencil);
            Some(Prepared { stencil, nodes })
        })
        .collect();

    // Particle to grid.
    for (p, prep) in particles.iter().zip(&prepared) {
        let Some(Prepared { stencil, nodes }) = prep else {
            continue;
        };
        let tau = kirchhoff_stress(&mat_from_col_major(&p.deformation_gradient), &c.material);
        let affine = mat_add(
            &mat_scale(&mat_from_col_major(&p.affine_velocity), p.mass),
            &mat_scale(&tau, -4.0 * c.dt * p.volume * inv_h2),
        );
        let momentum = [
            p.mass * p.velocity[0],
            p.mass * p.velocity[1],
            p.mass * p.velocity[2],
        ];
        let mut slot = 0;
        for a in 0..3 {
            for b in 0..3 {
                for k in 0..3 {
                    let (w, offset) = stencil.node(a, b, k);
                    let node = nodes[slot];
                    slot += 1;
                    let a_d = mat_vec(&affine, offset.map(|o| o * c.h));
                    grid.mass[node] += w * p.mass;
                    for ((acc, m), ad) in grid.vel[node].iter_mut().zip(momentum).zip(a_d) {
                        *acc += w * (m + ad);
                    }
                }
            }
        }
    }

    grid.update(c);

    // Grid to particle, deformation update, plasticity, advection.
    for (p, prep) in particles.iter_mut().zip(&prepared) {
        let Some(Prepared { stencil, nodes }) = prep else {
            continue;
        };
        let mut velocity = [0.0; 3];
        let mut c_matrix = ZERO;
        let mut slot = 0;
        for a in 0..3 {
            for b in 0..3 {
                for k in 0..3 {
                    let (w, offset) = stencil.node(a, b, k);
                    let vi = grid.vel[nodes[slot]];
                    slot += 1;
                    let wv = vi.map(|x| w * x);
                    for (acc, x) in velocity.iter_mut().zip(wv) {
                        *acc += x;
                    }
                    for (row, wr) in c_matrix.iter_mut().zip(wv) {
                        for (entry, o) in row.iter_mut().zip(offset) {
                            *entry += wr * o * c.h;
                        }
                    }
                }
            }
        }
        let c_matrix = mat_scale(&c_matrix, 4.0 * inv_h2);

        let f_old = mat_from_col_major(&p.deformation_gradient);
        let f_trial = mat_mul(&mat_add(&IDENTITY, &mat_scale(&c_matrix, c.dt)), &f_old);
        let (f_new, plastic_volume) = if mat_is_finite(&f_trial) {
            project_deformation(&f_trial, &c.material)
        } else {
            (IDENTITY, 0.0)
        };

        p.velocity = velocity;
        p.affine_velocity = mat_to_col_major(&c_matrix);
        p.deformation_gradient = mat_to_col_major(&f_new);
        p.plastic_strain_vol += plastic_volume;
        for (x, v) in p.position.iter_mut().zip(velocity) {
            *x += c.dt * v;
        }
        confine_to_domain(p, c.plane_strain_z);
    }
}

/// Keeps a particle on or above the floor plane and inside the optional
/// plane-strain slab, removing the velocity component that points out of it.
pub(super) fn confine_to_domain(p: &mut SandParticle, plane_strain_z: Option<[f64; 2]>) {
    if p.position[1] < 0.0 {
        p.position[1] = 0.0;
        p.velocity[1] = p.velocity[1].max(0.0);
    }
    if let Some([z_min, z_max]) = plane_strain_z {
        if p.position[2] < z_min {
            p.position[2] = z_min;
            p.velocity[2] = p.velocity[2].max(0.0);
        } else if p.position[2] > z_max {
            p.position[2] = z_max;
            p.velocity[2] = p.velocity[2].min(0.0);
        }
    }
}

/// Largest stable sub-step `cfl · h / (c + v_max)` for the active particles;
/// `f64::INFINITY` if there is none.
///
/// `c = √((λ + 2μ)/ρ_min)` is the P-wave speed of the least dense particle and
/// `v_max` the largest particle speed.
pub(super) fn stable_time_step(
    particles: &[SandParticle],
    h: f64,
    cfl: f64,
    material: &Material,
) -> f64 {
    let inv_h = 1.0 / h;
    let mut v_max: f64 = 0.0;
    let mut rho_min = f64::INFINITY;
    for p in particles.iter().filter(|p| is_active(p, inv_h)) {
        let speed = (p.velocity[0] * p.velocity[0]
            + p.velocity[1] * p.velocity[1]
            + p.velocity[2] * p.velocity[2])
            .sqrt();
        v_max = v_max.max(speed);
        rho_min = rho_min.min(p.mass / p.volume);
    }
    if !rho_min.is_finite() {
        return f64::INFINITY;
    }
    let wave_speed = (material.p_wave_modulus() / rho_min).sqrt();
    cfl * h / (wave_speed + v_max)
}
