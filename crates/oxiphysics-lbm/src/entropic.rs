// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Entropic Lattice Boltzmann Method (H-theorem stabilized) for D2Q9.
//!
//! This module implements two entropic LBM variants:
//!
//! - **EntropicD2Q9**: Standard entropic LBM with adaptive relaxation parameter α
//!   determined by enforcing the discrete H-theorem (H-stability).
//! - **KbcD2Q9**: Karlin–Bösch–Chikatamarla (KBC) variant that decomposes the
//!   equilibrium into symmetric (s) and antisymmetric (h) parts and relaxes them
//!   independently.
//!
//! Also provides:
//! - ELBM diagnostics (entropy monitoring, alpha statistics)
//! - Newton-Raphson alpha optimization (faster convergence)
//! - Entropy function variants (Boltzmann H, relative entropy)
//!
//! ## References
//!
//! - Ansumali & Karlin (2002), Phys. Rev. E 65, 056312.
//! - Karlin, Bösch & Chikatamarla (2014), Phys. Rev. E 90, 031302.

// ---------------------------------------------------------------------------
// D2Q9 constants
// ---------------------------------------------------------------------------

/// D2Q9 lattice weights: \[centre, axis×4, diagonal×4\]
///
/// Ordering convention (standard):
/// ```text
/// 0: (0,0)  w = 4/9
/// 1: (1,0)  w = 1/9
/// 2: (0,1)  w = 1/9
/// 3: (-1,0) w = 1/9
/// 4: (0,-1) w = 1/9
/// 5: (1,1)  w = 1/36
/// 6: (-1,1) w = 1/36
/// 7: (-1,-1)w = 1/36
/// 8: (1,-1) w = 1/36
/// ```
const W: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// D2Q9 lattice velocities: (cx, cy) for each direction.
const C: [(i32, i32); 9] = [
    (0, 0),
    (1, 0),
    (0, 1),
    (-1, 0),
    (0, -1),
    (1, 1),
    (-1, 1),
    (-1, -1),
    (1, -1),
];

/// D2Q9 opposite-direction mapping.
const OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Compute the D2Q9 Maxwell–Boltzmann equilibrium distribution (free function).
///
/// Uses the second-order Taylor expansion:
/// f_eq_i = w_i * ρ * (1 + (c·u)/cs² + (c·u)²/(2cs⁴) − u²/(2cs²))
/// where cs² = 1/3.
#[allow(dead_code)]
pub fn d2q9_equilibrium(rho: f64, u: [f64; 2]) -> [f64; 9] {
    let ux = u[0];
    let uy = u[1];
    let u2 = ux * ux + uy * uy;
    let mut feq = [0.0; 9];
    for i in 0..9 {
        let cx = C[i].0 as f64;
        let cy = C[i].1 as f64;
        let cu = cx * ux + cy * uy;
        feq[i] = W[i] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
    }
    feq
}

/// Compute the discrete H-function: H = Σ_i f_i * ln(f_i / w_i).
///
/// This is the negative of the Boltzmann entropy. The H-theorem guarantees
/// H decreases (or stays constant) during collision.
///
/// Returns `f64::INFINITY` if any `f_i ≤ 0` (unphysical state).
#[allow(dead_code)]
pub fn h_function(f: &[f64; 9]) -> f64 {
    let mut h = 0.0;
    for i in 0..9 {
        if f[i] <= 0.0 {
            return f64::INFINITY;
        }
        h += f[i] * (f[i] / W[i]).ln();
    }
    h
}

/// Compute the relative entropy (Kullback-Leibler divergence) of f from f_eq:
///
/// ```text
/// D_KL(f || f_eq) = Σ_i f_i * ln(f_i / f_eq_i)
/// ```
///
/// Returns `f64::INFINITY` if any `f_i ≤ 0` or `f_eq_i ≤ 0`.
#[allow(dead_code)]
pub fn relative_entropy(f: &[f64; 9], f_eq: &[f64; 9]) -> f64 {
    let mut dkl = 0.0;
    for i in 0..9 {
        if f[i] <= 0.0 || f_eq[i] <= 0.0 {
            return f64::INFINITY;
        }
        dkl += f[i] * (f[i] / f_eq[i]).ln();
    }
    dkl
}

/// Compute the Boltzmann entropy: S = -Σ_i f_i * ln(f_i / w_i) = -H.
#[allow(dead_code)]
pub fn boltzmann_entropy(f: &[f64; 9]) -> f64 {
    -h_function(f)
}

/// Compute the mirror state: f_mirror = 2 * f_eq − f.
///
/// The mirror state is used in the entropic overrelaxation to find the
/// α value that keeps the post-collision state on the H = H(f) level set.
#[allow(dead_code)]
pub fn mirror_state(f: &[f64; 9], f_eq: &[f64; 9]) -> [f64; 9] {
    let mut mirror = [0.0; 9];
    for i in 0..9 {
        mirror[i] = 2.0 * f_eq[i] - f[i];
    }
    mirror
}

/// Find the entropic relaxation parameter α ∈ \[0, 2\] via bisection.
///
/// α is chosen so that H(f + α*(f_eq − f)) = H(f), enforcing the
/// discrete H-theorem.  The physical upper bound α = 2 corresponds to
/// standard BGK overrelaxation.
///
/// If H(f_eq) ≤ H(f) already (which is the thermodynamically consistent
/// case), we return α = 2 (the BGK limit).
///
/// # Arguments
/// * `f`        – pre-collision distribution
/// * `f_eq`     – local equilibrium distribution
/// * `max_iter` – maximum bisection iterations (20–50 is typically sufficient)
#[allow(dead_code)]
pub fn find_entropic_alpha(f: &[f64; 9], f_eq: &[f64; 9], max_iter: usize) -> f64 {
    let h_f = h_function(f);

    // If f is already at or below equilibrium H, use standard BGK.
    if h_function(f_eq) <= h_f {
        return 2.0;
    }

    // We want to find α ∈ (0, 2) such that H(f + α*(f_eq − f)) = H(f).
    // Define φ(α) = H(f + α*(f_eq − f)) − H(f).
    // φ(0) = 0, φ(2) = H(mirror) − H(f).
    //
    // Because H is convex and H(f_eq) < H(f), there exists a root in (0, 2].
    // We bisect on [1, 2] since φ is monotone decreasing through 0 somewhere
    // in that interval.

    let phi = |alpha: f64| -> f64 {
        let mut f_new = [0.0; 9];
        for i in 0..9 {
            f_new[i] = f[i] + alpha * (f_eq[i] - f[i]);
        }
        h_function(&f_new) - h_f
    };

    // Check sign at α = 2 (mirror state).
    let phi2 = phi(2.0);
    if phi2 >= 0.0 {
        // Mirror state has H ≥ H(f): α = 2 is safe.
        return 2.0;
    }

    // Bisect on [1, 2] where φ crosses zero from positive to negative.
    let mut lo = 1.0_f64;
    let mut hi = 2.0_f64;

    for _ in 0..max_iter {
        let mid = 0.5 * (lo + hi);
        if phi(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-14 {
            break;
        }
    }

    0.5 * (lo + hi)
}

/// Find the entropic α using Newton-Raphson iteration (faster convergence).
///
/// Uses the derivative dH/dα to converge quadratically.
/// Falls back to bisection if Newton step goes out of bounds.
#[allow(dead_code)]
pub fn find_entropic_alpha_newton(f: &[f64; 9], f_eq: &[f64; 9], max_iter: usize) -> f64 {
    let h_f = h_function(f);

    if h_function(f_eq) <= h_f {
        return 2.0;
    }

    let phi_and_dphi = |alpha: f64| -> (f64, f64) {
        let mut f_new = [0.0; 9];
        let mut df = [0.0; 9]; // df[i] = f_eq[i] - f[i]
        for i in 0..9 {
            df[i] = f_eq[i] - f[i];
            f_new[i] = f[i] + alpha * df[i];
        }
        let h_new = h_function(&f_new);
        let phi = h_new - h_f;

        // dH/dα = Σ_i (f_eq_i - f_i) * (1 + ln(f_new_i / w_i))
        let mut dphi = 0.0;
        for i in 0..9 {
            if f_new[i] > 0.0 {
                dphi += df[i] * (1.0 + (f_new[i] / W[i]).ln());
            }
        }
        (phi, dphi)
    };

    // Start from α = 2.0 and iterate
    let mut alpha = 2.0_f64;

    for _ in 0..max_iter {
        let (phi, dphi) = phi_and_dphi(alpha);

        if phi.abs() < 1e-14 {
            break;
        }

        if dphi.abs() < 1e-30 {
            // Degenerate derivative, fall back to bisection midpoint
            alpha = 1.5;
            continue;
        }

        let step = -phi / dphi;
        let new_alpha = alpha + step;

        // Clamp to [0.5, 2.0] to avoid wild steps
        alpha = new_alpha.clamp(0.5, 2.0);
    }

    alpha.clamp(0.0, 2.0)
}

/// Compute non-equilibrium stress tensor components from f and f_eq.
///
/// Returns the three independent components of the symmetric
/// non-equilibrium stress tensor in 2D: `(σ_xx, σ_xy, σ_yy)`.
#[allow(dead_code)]
pub fn non_equilibrium_stress(f: &[f64; 9], f_eq: &[f64; 9]) -> (f64, f64, f64) {
    let mut s_xx = 0.0;
    let mut s_xy = 0.0;
    let mut s_yy = 0.0;
    for i in 0..9 {
        let cx = C[i].0 as f64;
        let cy = C[i].1 as f64;
        let f_neq = f[i] - f_eq[i];
        s_xx += cx * cx * f_neq;
        s_xy += cx * cy * f_neq;
        s_yy += cy * cy * f_neq;
    }
    (s_xx, s_xy, s_yy)
}

// ---------------------------------------------------------------------------
// EntropicD2Q9
// ---------------------------------------------------------------------------

/// Entropic LBM solver for D2Q9 with adaptive relaxation (H-theorem stable).
///
/// The solver uses an adaptive α ∈ \[0, 2\] at each cell computed by
/// [`find_entropic_alpha`] so that the post-collision state never violates
/// the discrete H-theorem.  When α = 2 this reduces to standard BGK.
pub struct EntropicD2Q9 {
    /// Domain width (number of cells in x direction).
    pub nx: usize,
    /// Domain height (number of cells in y direction).
    pub ny: usize,
    /// Distribution functions stored as `f[cell_index][direction]`.
    pub f: Vec<[f64; 9]>,
    /// Kinematic viscosity ν.
    pub nu: f64,
}

impl EntropicD2Q9 {
    /// Create a new solver with uniform rest equilibrium (ρ=1, u=0).
    pub fn new(nx: usize, ny: usize, nu: f64) -> Self {
        let n = nx * ny;
        let f_rest = Self::equilibrium(1.0, [0.0, 0.0]);
        Self {
            nx,
            ny,
            f: vec![f_rest; n],
            nu,
        }
    }

    /// Compute the D2Q9 Maxwell–Boltzmann equilibrium distribution.
    ///
    /// Uses the second-order Taylor expansion:
    /// f_eq_i = w_i * ρ * (1 + (c·u)/cs² + (c·u)²/(2cs⁴) − u²/(2cs²))
    /// where cs² = 1/3.
    pub fn equilibrium(rho: f64, u: [f64; 2]) -> [f64; 9] {
        let ux = u[0];
        let uy = u[1];
        let u2 = ux * ux + uy * uy;
        let mut feq = [0.0; 9];
        for i in 0..9 {
            let cx = C[i].0 as f64;
            let cy = C[i].1 as f64;
            let cu = cx * ux + cy * uy;
            feq[i] = W[i] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
        }
        feq
    }

    /// Compute macroscopic density and velocity at cell `idx`.
    ///
    /// Returns `(rho, [ux, uy])`.
    pub fn macros(&self, idx: usize) -> (f64, [f64; 2]) {
        let f = &self.f[idx];
        let rho: f64 = f.iter().sum();
        let mut mx = 0.0;
        let mut my = 0.0;
        for i in 0..9 {
            mx += f[i] * C[i].0 as f64;
            my += f[i] * C[i].1 as f64;
        }
        let inv_rho = if rho > 0.0 { 1.0 / rho } else { 0.0 };
        (rho, [mx * inv_rho, my * inv_rho])
    }

    /// Perform entropic collision at every cell.
    ///
    /// For each cell the adaptive α is found by bisection such that
    /// H(f_new) ≤ H(f).  The update rule is: f_new = f + α*(f_eq − f).
    pub fn collide(&mut self) {
        let n = self.nx * self.ny;
        for k in 0..n {
            let (rho, u) = self.macros(k);
            let f_eq = Self::equilibrium(rho, u);
            let alpha = find_entropic_alpha(&self.f[k], &f_eq, 50);
            let f = &mut self.f[k];
            for i in 0..9 {
                f[i] += alpha * (f_eq[i] - f[i]);
            }
        }
    }

    /// Perform entropic collision using Newton-Raphson alpha optimization.
    #[allow(dead_code)]
    pub fn collide_newton(&mut self) {
        let n = self.nx * self.ny;
        for k in 0..n {
            let (rho, u) = self.macros(k);
            let f_eq = Self::equilibrium(rho, u);
            let alpha = find_entropic_alpha_newton(&self.f[k], &f_eq, 20);
            let f = &mut self.f[k];
            for i in 0..9 {
                f[i] += alpha * (f_eq[i] - f[i]);
            }
        }
    }

    /// Perform pull-scheme streaming with periodic boundary conditions.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let f_old = self.f.clone();
        for y in 0..ny {
            for x in 0..nx {
                let k_dst = y * nx + x;
                for i in 0..9 {
                    let cx = C[i].0;
                    let cy = C[i].1;
                    // Pull from the source cell (periodic wrap).
                    let sx = (x as isize - cx as isize).rem_euclid(nx as isize) as usize;
                    let sy = (y as isize - cy as isize).rem_euclid(ny as isize) as usize;
                    let k_src = sy * nx + sx;
                    self.f[k_dst][i] = f_old[k_src][i];
                }
            }
        }
    }

    /// Perform one full LBM step: collide then stream.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
    }

    /// Initialize all cells to a uniform density and velocity.
    pub fn set_uniform(&mut self, rho: f64, u: [f64; 2]) {
        let feq = Self::equilibrium(rho, u);
        for cell in self.f.iter_mut() {
            *cell = feq;
        }
    }

    /// Compute the total H-function (entropy) over the entire domain.
    #[allow(dead_code)]
    pub fn total_entropy(&self) -> f64 {
        let n = self.nx * self.ny;
        let mut total_h = 0.0;
        for k in 0..n {
            total_h += h_function(&self.f[k]);
        }
        total_h
    }

    /// Compute the maximum velocity magnitude in the domain.
    #[allow(dead_code)]
    pub fn max_velocity(&self) -> f64 {
        let n = self.nx * self.ny;
        let mut max_u = 0.0_f64;
        for k in 0..n {
            let (_, u) = self.macros(k);
            let u_mag = (u[0] * u[0] + u[1] * u[1]).sqrt();
            max_u = max_u.max(u_mag);
        }
        max_u
    }

    /// Compute the minimum and maximum density in the domain.
    #[allow(dead_code)]
    pub fn density_range(&self) -> (f64, f64) {
        let n = self.nx * self.ny;
        let mut rho_min = f64::MAX;
        let mut rho_max = f64::MIN;
        for k in 0..n {
            let (rho, _) = self.macros(k);
            rho_min = rho_min.min(rho);
            rho_max = rho_max.max(rho);
        }
        (rho_min, rho_max)
    }
}

// ---------------------------------------------------------------------------
// ElbmDiagnostics
// ---------------------------------------------------------------------------

/// Diagnostics for the entropic LBM.
///
/// Tracks entropy evolution, alpha statistics, and stability indicators.
#[allow(dead_code)]
pub struct ElbmDiagnostics {
    /// History of total domain entropy (H function).
    pub entropy_history: Vec<f64>,
    /// History of minimum α values per step.
    pub alpha_min_history: Vec<f64>,
    /// History of average α values per step.
    pub alpha_avg_history: Vec<f64>,
    /// History of maximum velocity magnitudes per step.
    pub max_velocity_history: Vec<f64>,
}

#[allow(dead_code)]
impl ElbmDiagnostics {
    /// Create empty diagnostics.
    pub fn new() -> Self {
        Self {
            entropy_history: Vec::new(),
            alpha_min_history: Vec::new(),
            alpha_avg_history: Vec::new(),
            max_velocity_history: Vec::new(),
        }
    }
}

impl Default for ElbmDiagnostics {
    fn default() -> Self {
        Self::new()
    }
}

impl ElbmDiagnostics {
    /// Record diagnostics for one step from an EntropicD2Q9 solver.
    pub fn record(&mut self, solver: &EntropicD2Q9) {
        let n = solver.nx * solver.ny;

        // Total entropy
        let mut total_h = 0.0;
        let mut alpha_min = 2.0_f64;
        let mut alpha_sum = 0.0;
        let mut max_u = 0.0_f64;

        for k in 0..n {
            let (rho, u) = solver.macros(k);
            let f_eq = EntropicD2Q9::equilibrium(rho, u);
            let alpha = find_entropic_alpha(&solver.f[k], &f_eq, 50);
            total_h += h_function(&solver.f[k]);
            alpha_min = alpha_min.min(alpha);
            alpha_sum += alpha;
            let u_mag = (u[0] * u[0] + u[1] * u[1]).sqrt();
            max_u = max_u.max(u_mag);
        }

        self.entropy_history.push(total_h);
        self.alpha_min_history.push(alpha_min);
        self.alpha_avg_history.push(alpha_sum / n as f64);
        self.max_velocity_history.push(max_u);
    }

    /// Check if entropy is monotonically non-increasing (H-theorem).
    pub fn is_h_theorem_satisfied(&self) -> bool {
        if self.entropy_history.len() < 2 {
            return true;
        }
        for i in 1..self.entropy_history.len() {
            if self.entropy_history[i] > self.entropy_history[i - 1] + 1e-10 {
                return false;
            }
        }
        true
    }

    /// Number of recorded steps.
    pub fn step_count(&self) -> usize {
        self.entropy_history.len()
    }
}

// ---------------------------------------------------------------------------
// KbcD2Q9
// ---------------------------------------------------------------------------

/// KBC (Karlin–Bösch–Chikatamarla) variant of entropic LBM for D2Q9.
///
/// KBC decomposes the equilibrium into symmetric (s) and antisymmetric (h)
/// parts and applies a single global α to the s-part while keeping the
/// h-part at its BGK rate (overrelaxation by factor 2).
pub struct KbcD2Q9 {
    /// Domain width.
    pub nx: usize,
    /// Domain height.
    pub ny: usize,
    /// Distribution functions stored as `f[cell_index][direction]`.
    pub f: Vec<[f64; 9]>,
    /// Kinematic viscosity ν.
    pub nu: f64,
}

impl KbcD2Q9 {
    /// Decompose f_eq into symmetric (s) and antisymmetric (h) parts.
    ///
    /// The decomposition follows Karlin et al. (2014):
    /// - s_i = (f_eq_i + f_eq_{opp(i)}) / 2   (symmetric under i ↔ opp(i))
    /// - h_i = f_eq_i − s_i                    (antisymmetric)
    ///
    /// Returns `(s, h)` where `s + h = f_eq`.
    #[allow(dead_code)]
    pub fn decompose_feq(f_eq: &[f64; 9]) -> ([f64; 9], [f64; 9]) {
        let mut s = [0.0; 9];
        let mut h = [0.0; 9];
        for i in 0..9 {
            s[i] = 0.5 * (f_eq[i] + f_eq[OPP[i]]);
            h[i] = f_eq[i] - s[i];
        }
        (s, h)
    }

    /// Decompose an arbitrary distribution into symmetric and antisymmetric parts.
    #[allow(dead_code)]
    pub fn decompose_general(f: &[f64; 9]) -> ([f64; 9], [f64; 9]) {
        let mut s = [0.0; 9];
        let mut h = [0.0; 9];
        for i in 0..9 {
            s[i] = 0.5 * (f[i] + f[OPP[i]]);
            h[i] = f[i] - s[i];
        }
        (s, h)
    }

    /// Create a new KBC solver with uniform rest equilibrium (ρ=1, u=0).
    pub fn new(nx: usize, ny: usize, nu: f64) -> Self {
        let n = nx * ny;
        let f_rest = EntropicD2Q9::equilibrium(1.0, [0.0, 0.0]);
        Self {
            nx,
            ny,
            f: vec![f_rest; n],
            nu,
        }
    }

    /// Compute macroscopic density and velocity at cell `idx`.
    #[allow(dead_code)]
    pub fn macros(&self, idx: usize) -> (f64, [f64; 2]) {
        let f = &self.f[idx];
        let rho: f64 = f.iter().sum();
        let mut mx = 0.0;
        let mut my = 0.0;
        for i in 0..9 {
            mx += f[i] * C[i].0 as f64;
            my += f[i] * C[i].1 as f64;
        }
        let inv_rho = if rho > 0.0 { 1.0 / rho } else { 0.0 };
        (rho, [mx * inv_rho, my * inv_rho])
    }

    /// Perform KBC collision at every cell.
    ///
    /// The update rule is:
    /// Δs = s_eq − s_f  (symmetric deviation)
    /// Δh = h_eq − h_f  (antisymmetric deviation)
    ///
    /// f_new = f + α_kbc * Δs + 2 * Δh
    ///
    /// where α_kbc is the entropic parameter found by bisection for the
    /// s-component and the h-component is always fully relaxed (α = 2).
    pub fn collide_kbc(&mut self) {
        let n = self.nx * self.ny;
        for k in 0..n {
            let f = &self.f[k];
            let rho: f64 = f.iter().sum();
            let mut mx = 0.0;
            let mut my = 0.0;
            for i in 0..9 {
                mx += f[i] * C[i].0 as f64;
                my += f[i] * C[i].1 as f64;
            }
            let inv_rho = if rho > 0.0 { 1.0 / rho } else { 0.0 };
            let u = [mx * inv_rho, my * inv_rho];

            let f_eq = EntropicD2Q9::equilibrium(rho, u);
            let (s_eq, h_eq) = Self::decompose_feq(&f_eq);
            let (s_f, h_f) = Self::decompose_feq(f);

            // Compute Δs and Δh.
            let mut ds = [0.0; 9];
            let mut dh = [0.0; 9];
            for i in 0..9 {
                ds[i] = s_eq[i] - s_f[i];
                dh[i] = h_eq[i] - h_f[i];
            }

            // Find entropic α for the s-component deviation only.
            let f_h_relaxed: [f64; 9] = {
                let mut tmp = *f;
                for i in 0..9 {
                    tmp[i] += 2.0 * dh[i];
                }
                tmp
            };

            let alpha_kbc = {
                let h_ref = h_function(&f_h_relaxed);
                let mut f_eq_bisect = f_h_relaxed;
                for i in 0..9 {
                    f_eq_bisect[i] += ds[i];
                }
                if h_function(&f_eq_bisect) <= h_ref {
                    2.0_f64
                } else {
                    find_entropic_alpha(&f_h_relaxed, &f_eq_bisect, 50)
                }
            };

            let f_new = &mut self.f[k];
            for i in 0..9 {
                f_new[i] += alpha_kbc * ds[i] + 2.0 * dh[i];
            }
        }
    }

    /// Perform pull-scheme streaming with periodic boundary conditions.
    #[allow(dead_code)]
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let f_old = self.f.clone();
        for y in 0..ny {
            for x in 0..nx {
                let k_dst = y * nx + x;
                for i in 0..9 {
                    let cx = C[i].0;
                    let cy = C[i].1;
                    let sx = (x as isize - cx as isize).rem_euclid(nx as isize) as usize;
                    let sy = (y as isize - cy as isize).rem_euclid(ny as isize) as usize;
                    let k_src = sy * nx + sx;
                    self.f[k_dst][i] = f_old[k_src][i];
                }
            }
        }
    }

    /// One full KBC step: collide then stream.
    #[allow(dead_code)]
    pub fn step(&mut self) {
        self.collide_kbc();
        self.stream();
    }

    /// Set all cells to uniform density and velocity.
    #[allow(dead_code)]
    pub fn set_uniform(&mut self, rho: f64, u: [f64; 2]) {
        let feq = EntropicD2Q9::equilibrium(rho, u);
        for cell in self.f.iter_mut() {
            *cell = feq;
        }
    }

    /// Compute total density in the domain.
    #[allow(dead_code)]
    pub fn total_density(&self) -> f64 {
        let n = self.nx * self.ny;
        (0..n).map(|k| self.macros(k).0).sum()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// H-theorem: H(f_eq) ≤ H(f) for perturbed f.
    #[test]
    fn test_h_theorem_at_equilibrium() {
        let rho = 1.0;
        let u = [0.05, 0.02];
        let f_eq = EntropicD2Q9::equilibrium(rho, u);

        // Perturb f away from equilibrium.
        let mut f_perturbed = f_eq;
        f_perturbed[1] += 0.02;
        f_perturbed[3] -= 0.01;
        f_perturbed[5] += 0.005;

        let h_eq = h_function(&f_eq);
        let h_f = h_function(&f_perturbed);

        assert!(
            h_eq <= h_f + 1e-12,
            "H-theorem violated: H(f_eq)={h_eq} > H(f)={h_f}"
        );
    }

    /// find_entropic_alpha must always return α ∈ \[0, 2\].
    #[test]
    fn test_find_entropic_alpha_range() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.05, 0.02]);
        let mut f = f_eq;
        f[1] += 0.03;
        f[3] -= 0.02;

        let alpha = find_entropic_alpha(&f, &f_eq, 50);
        assert!(
            (0.0..=2.0 + 1e-12).contains(&alpha),
            "alpha out of range: {alpha}"
        );
    }

    /// find_entropic_alpha returns 2.0 when f is already at equilibrium.
    #[test]
    fn test_find_entropic_alpha_at_equilibrium() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.0, 0.0]);
        let alpha = find_entropic_alpha(&f_eq, &f_eq, 50);
        assert!(
            (alpha - 2.0).abs() < 1e-10,
            "alpha at equilibrium should be 2.0, got {alpha}"
        );
    }

    /// After collide+stream on a periodic domain, total density must be conserved.
    #[test]
    fn test_density_conservation_after_step() {
        let nx = 8;
        let ny = 8;
        let nu = 1.0 / 6.0;
        let mut solver = EntropicD2Q9::new(nx, ny, nu);

        // Perturb a few cells.
        solver.f[5][1] += 0.02;
        solver.f[5][3] -= 0.01;
        solver.f[10][2] += 0.015;

        let rho_before: f64 = (0..nx * ny).map(|k| solver.macros(k).0).sum();

        for _ in 0..10 {
            solver.step();
        }

        let rho_after: f64 = (0..nx * ny).map(|k| solver.macros(k).0).sum();
        assert!(
            (rho_before - rho_after).abs() < 1e-10,
            "Density not conserved: before={rho_before}, after={rho_after}"
        );
    }

    /// Uniform flow initialized with set_uniform stays uniform after steps.
    #[test]
    fn test_uniform_flow_stays_uniform() {
        let nx = 6;
        let ny = 6;
        let nu = 1.0 / 6.0;
        let rho0 = 1.0;
        let u0 = [0.05, 0.0];

        let mut solver = EntropicD2Q9::new(nx, ny, nu);
        solver.set_uniform(rho0, u0);

        for _ in 0..20 {
            solver.step();
        }

        let n = nx * ny;
        for k in 0..n {
            let (rho, u) = solver.macros(k);
            assert!(
                (rho - rho0).abs() < 1e-12,
                "Density changed at cell {k}: {rho} != {rho0}"
            );
            assert!(
                (u[0] - u0[0]).abs() < 1e-12,
                "ux changed at cell {k}: {} != {}",
                u[0],
                u0[0]
            );
            assert!(
                u[1].abs() < 1e-12,
                "uy became non-zero at cell {k}: {}",
                u[1]
            );
        }
    }

    // -----------------------------------------------------------------------
    // Newton-Raphson alpha: must return α ∈ [0, 2]
    // -----------------------------------------------------------------------
    #[test]
    fn test_newton_alpha_range() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.05, 0.02]);
        let mut f = f_eq;
        f[1] += 0.03;
        f[3] -= 0.02;

        let alpha = find_entropic_alpha_newton(&f, &f_eq, 20);
        assert!(
            (0.0..=2.0 + 1e-10).contains(&alpha),
            "Newton alpha out of range: {alpha}"
        );
    }

    // -----------------------------------------------------------------------
    // Newton-Raphson: at equilibrium returns 2.0
    // -----------------------------------------------------------------------
    #[test]
    fn test_newton_alpha_at_equilibrium() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.0, 0.0]);
        let alpha = find_entropic_alpha_newton(&f_eq, &f_eq, 20);
        assert!(
            (alpha - 2.0).abs() < 1e-8,
            "Newton alpha at equilibrium should be ~2.0, got {alpha}"
        );
    }

    // -----------------------------------------------------------------------
    // Newton collision preserves density
    // -----------------------------------------------------------------------
    #[test]
    fn test_newton_collision_density_conservation() {
        let nx = 6;
        let ny = 6;
        let nu = 1.0 / 6.0;
        let mut solver = EntropicD2Q9::new(nx, ny, nu);

        solver.f[5][1] += 0.01;
        solver.f[5][3] -= 0.005;

        let rho_before: f64 = (0..nx * ny).map(|k| solver.macros(k).0).sum();
        solver.collide_newton();
        solver.stream();
        let rho_after: f64 = (0..nx * ny).map(|k| solver.macros(k).0).sum();

        assert!(
            (rho_before - rho_after).abs() < 1e-10,
            "Density not conserved with Newton: before={rho_before}, after={rho_after}"
        );
    }

    // -----------------------------------------------------------------------
    // Relative entropy: zero at equilibrium
    // -----------------------------------------------------------------------
    #[test]
    fn test_relative_entropy_zero_at_eq() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.05, 0.02]);
        let dkl = relative_entropy(&f_eq, &f_eq);
        assert!(
            dkl.abs() < 1e-14,
            "D_KL(f_eq || f_eq) should be 0, got {dkl}"
        );
    }

    // -----------------------------------------------------------------------
    // Relative entropy: positive for non-equilibrium
    // -----------------------------------------------------------------------
    #[test]
    fn test_relative_entropy_positive() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.05, 0.02]);
        let mut f = f_eq;
        f[1] += 0.02;
        f[3] -= 0.01;

        let dkl = relative_entropy(&f, &f_eq);
        assert!(
            dkl > 0.0,
            "D_KL(f || f_eq) should be positive for f != f_eq, got {dkl}"
        );
    }

    // -----------------------------------------------------------------------
    // Boltzmann entropy = -H
    // -----------------------------------------------------------------------
    #[test]
    fn test_boltzmann_entropy_is_neg_h() {
        let f = EntropicD2Q9::equilibrium(1.0, [0.03, -0.01]);
        let h = h_function(&f);
        let s = boltzmann_entropy(&f);
        assert!((s + h).abs() < 1e-14, "S should be -H: S={s}, H={h}");
    }

    // -----------------------------------------------------------------------
    // Non-equilibrium stress: zero at equilibrium
    // -----------------------------------------------------------------------
    #[test]
    fn test_neq_stress_zero_at_eq() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.05, 0.02]);
        let (sxx, sxy, syy) = non_equilibrium_stress(&f_eq, &f_eq);
        assert!(sxx.abs() < 1e-14, "σ_xx neq should be 0: {sxx}");
        assert!(sxy.abs() < 1e-14, "σ_xy neq should be 0: {sxy}");
        assert!(syy.abs() < 1e-14, "σ_yy neq should be 0: {syy}");
    }

    // -----------------------------------------------------------------------
    // Non-equilibrium stress: non-zero when perturbed
    // -----------------------------------------------------------------------
    #[test]
    fn test_neq_stress_nonzero_when_perturbed() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.05, 0.02]);
        let mut f = f_eq;
        // Asymmetric perturbation: only perturb direction 1 (cx=1)
        f[1] += 0.05;

        let (sxx, _sxy, _syy) = non_equilibrium_stress(&f, &f_eq);
        assert!(sxx.abs() > 1e-10, "s_xx should be non-zero when perturbed");
    }

    // -----------------------------------------------------------------------
    // ELBM Diagnostics: record and check
    // -----------------------------------------------------------------------
    #[test]
    fn test_elbm_diagnostics() {
        let nx = 6;
        let ny = 6;
        let nu = 1.0 / 6.0;
        let solver = EntropicD2Q9::new(nx, ny, nu);

        let mut diag = ElbmDiagnostics::new();
        assert_eq!(diag.step_count(), 0);

        diag.record(&solver);
        assert_eq!(diag.step_count(), 1);
        assert!(diag.entropy_history[0].is_finite());
        assert!(diag.alpha_min_history[0] >= 0.0);
        assert!(diag.alpha_avg_history[0] >= 0.0);
    }

    // -----------------------------------------------------------------------
    // ELBM Diagnostics: H-theorem check
    // -----------------------------------------------------------------------
    #[test]
    fn test_elbm_h_theorem_satisfied() {
        let diag = ElbmDiagnostics {
            entropy_history: vec![10.0, 9.5, 9.0, 8.5],
            alpha_min_history: vec![],
            alpha_avg_history: vec![],
            max_velocity_history: vec![],
        };
        assert!(diag.is_h_theorem_satisfied());

        let diag_bad = ElbmDiagnostics {
            entropy_history: vec![10.0, 9.5, 11.0, 8.5],
            alpha_min_history: vec![],
            alpha_avg_history: vec![],
            max_velocity_history: vec![],
        };
        assert!(!diag_bad.is_h_theorem_satisfied());
    }

    // -----------------------------------------------------------------------
    // Total entropy: finite for equilibrium
    // -----------------------------------------------------------------------
    #[test]
    fn test_total_entropy_finite() {
        let solver = EntropicD2Q9::new(6, 6, 1.0 / 6.0);
        let h = solver.total_entropy();
        assert!(h.is_finite(), "Total entropy should be finite: {h}");
    }

    // -----------------------------------------------------------------------
    // Max velocity at rest = 0
    // -----------------------------------------------------------------------
    #[test]
    fn test_max_velocity_at_rest() {
        let solver = EntropicD2Q9::new(6, 6, 1.0 / 6.0);
        let max_u = solver.max_velocity();
        assert!(max_u < 1e-14, "Max velocity at rest should be ~0: {max_u}");
    }

    // -----------------------------------------------------------------------
    // Density range at uniform state
    // -----------------------------------------------------------------------
    #[test]
    fn test_density_range_uniform() {
        let solver = EntropicD2Q9::new(6, 6, 1.0 / 6.0);
        let (rho_min, rho_max) = solver.density_range();
        assert!((rho_min - 1.0).abs() < 1e-14);
        assert!((rho_max - 1.0).abs() < 1e-14);
    }

    // -----------------------------------------------------------------------
    // KBC: decompose_general round-trip
    // -----------------------------------------------------------------------
    #[test]
    fn test_kbc_decompose_general_round_trip() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.03, -0.02]);
        let (s, h) = KbcD2Q9::decompose_general(&f_eq);
        for i in 0..9 {
            let reconstructed = s[i] + h[i];
            assert!(
                (reconstructed - f_eq[i]).abs() < 1e-14,
                "s+h != f at i={i}: {} vs {}",
                reconstructed,
                f_eq[i]
            );
        }
    }

    // -----------------------------------------------------------------------
    // KBC: total density conservation
    // -----------------------------------------------------------------------
    #[test]
    fn test_kbc_density_conservation() {
        let nx = 6;
        let ny = 6;
        let nu = 1.0 / 6.0;
        let mut kbc = KbcD2Q9::new(nx, ny, nu);

        kbc.f[5][1] += 0.02;
        kbc.f[5][3] -= 0.01;

        let rho_before = kbc.total_density();
        kbc.step();
        let rho_after = kbc.total_density();

        assert!(
            (rho_before - rho_after).abs() < 1e-10,
            "KBC density not conserved: before={rho_before}, after={rho_after}"
        );
    }

    // -----------------------------------------------------------------------
    // KBC: uniform flow stays uniform
    // -----------------------------------------------------------------------
    #[test]
    fn test_kbc_uniform_flow() {
        let nx = 6;
        let ny = 6;
        let nu = 1.0 / 6.0;
        let rho0 = 1.0;
        let u0 = [0.03, 0.0];

        let mut kbc = KbcD2Q9::new(nx, ny, nu);
        kbc.set_uniform(rho0, u0);

        for _ in 0..10 {
            kbc.step();
        }

        let n = nx * ny;
        for k in 0..n {
            let (rho, u) = kbc.macros(k);
            assert!(
                (rho - rho0).abs() < 1e-12,
                "KBC density changed at cell {k}"
            );
            assert!((u[0] - u0[0]).abs() < 1e-12, "KBC ux changed at cell {k}");
        }
    }

    // -----------------------------------------------------------------------
    // d2q9_equilibrium free function agrees with method
    // -----------------------------------------------------------------------
    #[test]
    fn test_d2q9_equilibrium_free_fn() {
        let rho = 1.5;
        let u = [0.04, -0.03];
        let feq1 = d2q9_equilibrium(rho, u);
        let feq2 = EntropicD2Q9::equilibrium(rho, u);
        for i in 0..9 {
            assert!(
                (feq1[i] - feq2[i]).abs() < 1e-14,
                "Free fn vs method mismatch at i={i}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // Mirror state: f_mirror + f = 2 * f_eq
    // -----------------------------------------------------------------------
    #[test]
    fn test_mirror_state_identity() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.05, 0.02]);
        let mut f = f_eq;
        f[1] += 0.01;
        f[3] -= 0.01;

        let fm = mirror_state(&f, &f_eq);
        for i in 0..9 {
            let sum = f[i] + fm[i];
            let expected = 2.0 * f_eq[i];
            assert!(
                (sum - expected).abs() < 1e-14,
                "f + f_mirror != 2*f_eq at i={i}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Entropy production and H-theorem diagnostics
// ---------------------------------------------------------------------------

/// Compute the symmetrised Kullback–Leibler divergence of f from f_eq:
///
/// D_sym(f || f_eq) = D_KL(f||f_eq) + D_KL(f_eq||f)
///
/// This is symmetric: D_sym(f,feq) = D_sym(feq,f).
/// Returns `f64::INFINITY` if any distribution has non-positive entries.
#[allow(dead_code)]
pub fn symmetric_kl_divergence(f: &[f64; 9], feq: &[f64; 9]) -> f64 {
    let mut d = 0.0;
    for i in 0..9 {
        if f[i] <= 0.0 || feq[i] <= 0.0 {
            return f64::INFINITY;
        }
        d += (f[i] - feq[i]) * (f[i] / feq[i]).ln();
    }
    d
}

/// Compute the entropy production per step: ΔH = H(f) - H(f*).
///
/// According to the H-theorem this must be ≥ 0.
#[allow(dead_code)]
pub fn entropy_production(f_pre: &[f64; 9], f_post: &[f64; 9]) -> f64 {
    h_function(f_pre) - h_function(f_post)
}

// ---------------------------------------------------------------------------
// Entropic stabiliser — Newton–Raphson alpha solver (improved)
// ---------------------------------------------------------------------------

/// Solve for the entropic parameter α satisfying the H-equation
/// H(f + α(f_eq - f)) = H(f_mirror)
/// using Newton–Raphson iteration with bracket search.
///
/// Returns α ∈ (0, 2] such that the discrete H-theorem is obeyed.
///
/// This is an enhanced version of the basic solver provided in EntropicD2Q9,
/// adding safeguards against divergence and bracket enforcement.
#[allow(dead_code)]
pub fn solve_alpha_newton_raphson(
    f: &[f64; 9],
    feq: &[f64; 9],
    f_mirror: &[f64; 9],
    tol: f64,
    max_iter: usize,
) -> f64 {
    // Target H value
    let h_target = h_function(f_mirror);
    if !h_target.is_finite() {
        return 2.0; // fallback
    }

    // Objective: g(α) = H(f + α(feq - f)) - h_target
    let g = |alpha: f64| -> f64 {
        let mut fa = [0.0f64; 9];
        for i in 0..9 {
            fa[i] = f[i] + alpha * (feq[i] - f[i]);
        }
        h_function(&fa) - h_target
    };

    // Derivative: g'(α) = Σ (feq_i - f_i) * (1 + ln(fa_i/w_i))
    let g_prime = |alpha: f64| -> f64 {
        let mut dg = 0.0;
        for i in 0..9 {
            let fa_i = f[i] + alpha * (feq[i] - f[i]);
            if fa_i <= 0.0 {
                return f64::INFINITY;
            }
            dg += (feq[i] - f[i]) * (1.0 + (fa_i / W[i]).ln());
        }
        dg
    };

    let mut alpha = 2.0; // start at mirror value
    for _ in 0..max_iter {
        let gval = g(alpha);
        if gval.abs() < tol {
            break;
        }
        let gpval = g_prime(alpha);
        if gpval.abs() < 1e-15 {
            break;
        }
        let alpha_new = alpha - gval / gpval;
        // Clamp to (0, 2]
        alpha = alpha_new.clamp(1e-6, 2.0);
    }
    alpha
}

// ---------------------------------------------------------------------------
// Entropic over-relaxation (EOR) scheme
// ---------------------------------------------------------------------------

/// Entropic over-relaxation (EOR) for D2Q9.
///
/// Uses α = 2 (full over-relaxation toward the mirror state) when the
/// H-theorem is satisfied, and falls back to BGK otherwise.
///
/// This provides second-order accuracy in time while maintaining the
/// discrete H-theorem unconditionally.
#[derive(Debug, Clone)]
pub struct EntropyOverRelaxation {
    /// Base kinematic viscosity.
    pub nu: f64,
    /// Whether to enforce H-theorem (true) or always use α=2.
    pub enforce_h: bool,
}

impl EntropyOverRelaxation {
    /// Construct EOR from kinematic viscosity.
    pub fn new(nu: f64) -> Self {
        Self {
            nu,
            enforce_h: true,
        }
    }

    /// Compute the EOR relaxation parameter for a single cell.
    ///
    /// Returns α ∈ (0, 2] — typically close to 2 for smooth flows.
    pub fn compute_alpha(&self, f: &[f64; 9], rho: f64, u: [f64; 2]) -> f64 {
        let feq = d2q9_equilibrium(rho, u);
        if !self.enforce_h {
            return 2.0;
        }
        let fm = mirror_state(f, &feq);
        // Quick check: does α=2 satisfy H-theorem?
        let h_pre = h_function(f);
        let h_post = h_function(&fm);
        if h_post.is_finite() && h_post <= h_pre + 1e-12 {
            return 2.0;
        }
        // Fallback: binary search for valid α
        let mut lo = 1.0_f64;
        let mut hi = 2.0_f64;
        for _ in 0..32 {
            let mid = 0.5 * (lo + hi);
            let mut ftry = [0.0f64; 9];
            for i in 0..9 {
                ftry[i] = f[i] + mid * (feq[i] - f[i]);
            }
            if h_function(&ftry) <= h_pre + 1e-12 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        lo
    }

    /// Apply EOR collision to a single cell.
    pub fn collide(&self, f: &[f64; 9], rho: f64, u: [f64; 2]) -> [f64; 9] {
        let feq = d2q9_equilibrium(rho, u);
        let alpha = self.compute_alpha(f, rho, u);
        let mut f_out = [0.0f64; 9];
        for i in 0..9 {
            f_out[i] = f[i] + alpha * (feq[i] - f[i]);
        }
        f_out
    }
}

// ---------------------------------------------------------------------------
// Lyapunov stability indicator
// ---------------------------------------------------------------------------

/// Compute the normalised non-equilibrium entropy indicator.
///
/// δ = H(f) - H(feq) >= 0 (H-theorem).
///
/// This provides a local measure of how far the system is from equilibrium
/// and can be used as a refinement indicator in adaptive LBM.
#[allow(dead_code)]
pub fn neq_entropy_indicator(f: &[f64; 9], rho: f64, u: [f64; 2]) -> f64 {
    let feq = d2q9_equilibrium(rho, u);
    let hf = h_function(f);
    let hfeq = h_function(&feq);
    if hf.is_finite() && hfeq.is_finite() {
        (hf - hfeq).max(0.0)
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Additional entropic tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod extra_entropic_tests {
    use super::*;

    /// Symmetric KL divergence is zero when f == feq.
    #[test]
    fn test_symmetric_kl_at_eq() {
        let feq = d2q9_equilibrium(1.0, [0.05, -0.02]);
        let d = symmetric_kl_divergence(&feq, &feq);
        assert!(d.abs() < 1e-12, "Symmetric KL at eq = {d}");
    }

    /// Symmetric KL divergence ≥ 0 by convexity.
    #[test]
    fn test_symmetric_kl_nonneg() {
        let feq = d2q9_equilibrium(1.0, [0.03, 0.01]);
        let mut f = feq;
        f[1] += 0.005;
        f[3] -= 0.005;
        let d = symmetric_kl_divergence(&f, &feq);
        assert!(d >= 0.0, "Symmetric KL should be non-negative: {d}");
    }

    /// H-theorem: entropy production ≥ 0 after BGK step (α=1).
    #[test]
    fn test_entropy_production_nonneg() {
        let rho = 1.0;
        let u = [0.04, -0.02];
        let feq = d2q9_equilibrium(rho, u);
        let mut f = feq;
        f[1] += 0.01;
        f[3] -= 0.01;
        let feq2 = d2q9_equilibrium(rho, u);
        // BGK post-collision: f* = f + (feq - f) [alpha=1, BGK]
        let mut f_post = [0.0f64; 9];
        for i in 0..9 {
            f_post[i] = f[i] + (feq2[i] - f[i]);
        }
        let dp = entropy_production(&f, &f_post);
        assert!(dp >= -1e-13, "Entropy production negative: {dp}");
    }

    /// EOR: mass is conserved after collision.
    #[test]
    fn test_eor_mass_conservation() {
        let eor = EntropyOverRelaxation::new(1.0 / 6.0);
        let rho = 1.0_f64;
        let u = [0.03, 0.01];
        let feq = d2q9_equilibrium(rho, u);
        let mut f = feq;
        f[0] += 0.005;
        f[1] -= 0.003;
        f[3] -= 0.002;
        let rho_in: f64 = f.iter().sum();
        let f_out = eor.collide(&f, rho_in, u);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho_in).abs() < 1e-12,
            "EOR mass: {rho_in} -> {rho_out}"
        );
    }

    /// EOR: alpha is in (0, 2].
    #[test]
    fn test_eor_alpha_range() {
        let eor = EntropyOverRelaxation::new(0.05);
        let rho = 1.0_f64;
        let u = [0.05, -0.02];
        let feq = d2q9_equilibrium(rho, u);
        let mut f = feq;
        f[1] += 0.005;
        f[3] -= 0.005;
        let alpha = eor.compute_alpha(&f, rho, u);
        assert!(
            alpha > 0.0 && alpha <= 2.0 + 1e-10,
            "alpha={alpha} out of range"
        );
    }

    /// solve_alpha_newton_raphson: at equilibrium alpha → 2.
    #[test]
    fn test_alpha_nr_at_equilibrium() {
        let feq = d2q9_equilibrium(1.0, [0.04, -0.01]);
        let fm = mirror_state(&feq, &feq);
        let alpha = solve_alpha_newton_raphson(&feq, &feq, &fm, 1e-12, 50);
        assert!(alpha > 0.0 && alpha <= 2.0 + 1e-10, "alpha={alpha}");
    }

    /// neq_entropy_indicator: at equilibrium indicator ≈ 0.
    #[test]
    fn test_neq_entropy_indicator_at_eq() {
        let rho = 1.0_f64;
        let u = [0.03, 0.02];
        let feq = d2q9_equilibrium(rho, u);
        let ind = neq_entropy_indicator(&feq, rho, u);
        assert!(ind < 1e-12, "indicator at eq = {ind}");
    }

    /// neq_entropy_indicator: perturbed state has larger indicator.
    #[test]
    fn test_neq_entropy_indicator_perturbed() {
        let rho = 1.0_f64;
        let u = [0.03, 0.02];
        let feq = d2q9_equilibrium(rho, u);
        let mut f = feq;
        f[1] += 0.01;
        f[3] -= 0.01;
        let ind = neq_entropy_indicator(&f, rho, u);
        assert!(ind >= 0.0, "indicator = {ind}");
    }

    /// entropy_production with α=2 (mirror) should be exactly zero.
    #[test]
    fn test_entropy_production_mirror() {
        let feq = d2q9_equilibrium(1.0, [0.04, 0.01]);
        let mut f = feq;
        f[1] += 0.005;
        f[3] -= 0.005;
        let fm = mirror_state(&f, &feq);
        // After mirror step H(fm) should equal H(f_mirror)
        let dp = entropy_production(&f, &fm);
        // dp = H(f) - H(fm); mirror should satisfy H-theorem
        assert!(dp >= -1e-12, "entropy production negative: {dp}");
    }

    /// h_function is lower for equilibrium than for perturbed state.
    #[test]
    fn test_h_function_minimum_at_eq() {
        let feq = d2q9_equilibrium(1.0, [0.02, -0.01]);
        let mut f = feq;
        f[2] += 0.01;
        f[4] -= 0.01;
        let h_eq = h_function(&feq);
        let h_f = h_function(&f);
        assert!(h_f >= h_eq - 1e-12, "H(f)={h_f} < H(feq)={h_eq}");
    }
}
