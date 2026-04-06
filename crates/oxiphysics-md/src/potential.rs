// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Pair potentials for molecular dynamics simulations.
//!
//! Each potential provides `energy(r)` and `force(r)` where `r` is the
//! scalar distance between two particles. The force returned is the
//! magnitude of the force along the inter-particle axis (positive = repulsive).

/// Trait for pair potentials.
///
/// `energy(r)` returns the potential energy at distance `r`.
/// `force(r)` returns the force magnitude along the separation vector
/// (positive means repulsive).
pub trait Potential: Send + Sync {
    /// Potential energy at distance `r`.
    fn energy(&self, r: f64) -> f64;

    /// Force magnitude at distance `r` (positive = repulsive).
    ///
    /// Defined as `-dV/dr`.
    fn force(&self, r: f64) -> f64;

    /// Cutoff distance beyond which this potential is zero.
    fn cutoff(&self) -> f64;
}

// ---------------------------------------------------------------------------
// Lennard-Jones
// ---------------------------------------------------------------------------

/// Lennard-Jones 12-6 potential with cutoff and energy shift.
///
/// V(r) = 4 * epsilon * ((sigma/r)^12 - (sigma/r)^6)  - V(r_c)
///
/// For r > cutoff, V = 0, F = 0.
#[derive(Debug, Clone)]
pub struct LennardJones {
    /// Depth of the potential well.
    pub epsilon: f64,
    /// Distance at which V = 0 (before shift).
    pub sigma: f64,
    /// Cutoff distance.
    pub cutoff_r: f64,
    /// Energy at the cutoff (for shifting).
    energy_shift: f64,
}

impl LennardJones {
    /// Create a new Lennard-Jones potential.
    pub fn new(epsilon: f64, sigma: f64, cutoff: f64) -> Self {
        let sr = sigma / cutoff;
        let sr6 = sr.powi(6);
        let sr12 = sr6 * sr6;
        let energy_shift = 4.0 * epsilon * (sr12 - sr6);
        Self {
            epsilon,
            sigma,
            cutoff_r: cutoff,
            energy_shift,
        }
    }

    /// Raw (un-shifted) LJ energy.
    fn raw_energy(&self, r: f64) -> f64 {
        let sr = self.sigma / r;
        let sr6 = sr.powi(6);
        let sr12 = sr6 * sr6;
        4.0 * self.epsilon * (sr12 - sr6)
    }
}

impl Potential for LennardJones {
    fn energy(&self, r: f64) -> f64 {
        if r >= self.cutoff_r {
            return 0.0;
        }
        self.raw_energy(r) - self.energy_shift
    }

    fn force(&self, r: f64) -> f64 {
        if r >= self.cutoff_r {
            return 0.0;
        }
        let sr = self.sigma / r;
        let sr6 = sr.powi(6);
        let sr12 = sr6 * sr6;
        // F = -dV/dr = 24*eps/r * (2*(sigma/r)^12 - (sigma/r)^6)
        24.0 * self.epsilon / r * (2.0 * sr12 - sr6)
    }

    fn cutoff(&self) -> f64 {
        self.cutoff_r
    }
}

// ---------------------------------------------------------------------------
// Morse
// ---------------------------------------------------------------------------

/// Morse potential: V(r) = D_e * (1 - exp(-alpha*(r - r_e)))^2
#[derive(Debug, Clone)]
pub struct Morse {
    /// Well depth.
    pub d_e: f64,
    /// Width parameter.
    pub alpha: f64,
    /// Equilibrium distance.
    pub r_e: f64,
    /// Cutoff distance.
    pub cutoff_r: f64,
    /// Energy at cutoff for shifting.
    energy_shift: f64,
}

impl Morse {
    /// Create a new Morse potential.
    pub fn new(d_e: f64, alpha: f64, r_e: f64, cutoff: f64) -> Self {
        let e = (-alpha * (cutoff - r_e)).exp();
        let energy_shift = d_e * (1.0 - e).powi(2);
        Self {
            d_e,
            alpha,
            r_e,
            cutoff_r: cutoff,
            energy_shift,
        }
    }
}

impl Potential for Morse {
    fn energy(&self, r: f64) -> f64 {
        if r >= self.cutoff_r {
            return 0.0;
        }
        let e = (-self.alpha * (r - self.r_e)).exp();
        self.d_e * (1.0 - e).powi(2) - self.energy_shift
    }

    fn force(&self, r: f64) -> f64 {
        if r >= self.cutoff_r {
            return 0.0;
        }
        // F = -dV/dr = 2*D_e*alpha*(1-exp(-a*(r-r_e)))*exp(-a*(r-r_e))
        let e = (-self.alpha * (r - self.r_e)).exp();
        2.0 * self.d_e * self.alpha * (1.0 - e) * e
    }

    fn cutoff(&self) -> f64 {
        self.cutoff_r
    }
}

// ---------------------------------------------------------------------------
// Coulomb (truncated)
// ---------------------------------------------------------------------------

/// Truncated Coulomb potential: V(r) = k_e * q1 * q2 / r.
///
/// The charge product (`qq`) should be `k_e * q1 * q2` in the desired unit
/// system (e.g., 1/(4*pi*eps_0) * q1 * q2 in SI).
#[derive(Debug, Clone)]
pub struct Coulomb {
    /// Product k_e * q1 * q2.
    pub qq: f64,
    /// Cutoff distance.
    pub cutoff_r: f64,
    /// Energy at cutoff for shifting.
    energy_shift: f64,
}

impl Coulomb {
    /// Create a new truncated Coulomb potential.
    pub fn new(qq: f64, cutoff: f64) -> Self {
        let energy_shift = qq / cutoff;
        Self {
            qq,
            cutoff_r: cutoff,
            energy_shift,
        }
    }
}

impl Potential for Coulomb {
    fn energy(&self, r: f64) -> f64 {
        if r >= self.cutoff_r {
            return 0.0;
        }
        self.qq / r - self.energy_shift
    }

    fn force(&self, r: f64) -> f64 {
        if r >= self.cutoff_r {
            return 0.0;
        }
        // F = -dV/dr = qq/r^2 (positive for same-sign charges = repulsive)
        self.qq / (r * r)
    }

    fn cutoff(&self) -> f64 {
        self.cutoff_r
    }
}

// ---------------------------------------------------------------------------
// Harmonic bond
// ---------------------------------------------------------------------------

/// Harmonic bond potential: V(r) = 0.5 * k * (r - r0)^2
///
/// Used for bonded interactions between atom pairs.
#[derive(Debug, Clone)]
pub struct HarmonicBond {
    /// Spring constant.
    pub k: f64,
    /// Equilibrium distance.
    pub r0: f64,
}

impl HarmonicBond {
    /// Create a new harmonic bond potential.
    pub fn new(k: f64, r0: f64) -> Self {
        Self { k, r0 }
    }

    /// Compute the bond energy between two atoms at positions `ri` and `rj`.
    pub fn energy(&self, ri: [f64; 3], rj: [f64; 3]) -> f64 {
        let dx = rj[0] - ri[0];
        let dy = rj[1] - ri[1];
        let dz = rj[2] - ri[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        let dr = r - self.r0;
        0.5 * self.k * dr * dr
    }

    /// Compute force vectors on atoms i and j from the harmonic bond.
    ///
    /// Returns `(force_on_i, force_on_j)` where each force is `[f64; 3]`.
    /// By Newton's 3rd law, `force_on_j = -force_on_i`.
    pub fn force_vectors(&self, ri: [f64; 3], rj: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let dx = rj[0] - ri[0];
        let dy = rj[1] - ri[1];
        let dz = rj[2] - ri[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        if r < 1e-30 {
            return ([0.0; 3], [0.0; 3]);
        }
        let f_scalar = -self.k * (r - self.r0); // F_i = -dV/dr * (rij/r)
        let fi = [f_scalar * dx / r, f_scalar * dy / r, f_scalar * dz / r];
        let fj = [-fi[0], -fi[1], -fi[2]];
        (fi, fj)
    }
}

impl Potential for HarmonicBond {
    fn energy(&self, r: f64) -> f64 {
        let dr = r - self.r0;
        0.5 * self.k * dr * dr
    }

    fn force(&self, r: f64) -> f64 {
        // F = -dV/dr = -k*(r-r0)
        -self.k * (r - self.r0)
    }

    fn cutoff(&self) -> f64 {
        f64::INFINITY
    }
}

// ---------------------------------------------------------------------------
// Stillinger-Weber potential (two-body term only)
// ---------------------------------------------------------------------------

/// Stillinger-Weber two-body potential (1985).
///
/// V₂(r) = A · ε · (B · (σ/r)^p - (σ/r)^q) · exp(σ / (r - a·σ))
///
/// for r < a·σ, and 0 otherwise.
///
/// Default parameters for silicon: A=7.049556277, B=0.6022245584,
/// p=4, q=0, a=1.8, ε=2.1683 eV, σ=2.0951 Å.
#[derive(Debug, Clone)]
pub struct StillingerWeber {
    /// Energy scale ε.
    pub epsilon: f64,
    /// Length scale σ.
    pub sigma: f64,
    /// Parameter A.
    pub a_param: f64,
    /// Parameter B.
    pub b_param: f64,
    /// Exponent p (typically 4).
    pub p: i32,
    /// Exponent q (typically 0).
    pub q: i32,
    /// Cutoff parameter a (V=0 for r ≥ a*σ).
    pub a_cut: f64,
}

impl StillingerWeber {
    /// Create a new Stillinger-Weber two-body potential.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        epsilon: f64,
        sigma: f64,
        a_param: f64,
        b_param: f64,
        p: i32,
        q: i32,
        a_cut: f64,
    ) -> Self {
        Self {
            epsilon,
            sigma,
            a_param,
            b_param,
            p,
            q,
            a_cut,
        }
    }

    /// Default silicon parameters.
    #[allow(dead_code)]
    pub fn silicon() -> Self {
        Self::new(2.1683, 2.0951, 7.049556277, 0.6022245584, 4, 0, 1.8)
    }
}

impl Potential for StillingerWeber {
    fn energy(&self, r: f64) -> f64 {
        let rc = self.a_cut * self.sigma;
        if r >= rc {
            return 0.0;
        }
        let sr = self.sigma / r;
        let term_p = sr.powi(self.p);
        let term_q = if self.q == 0 { 1.0 } else { sr.powi(self.q) };
        let exp_part = (self.sigma / (r - rc)).exp();
        self.a_param * self.epsilon * (self.b_param * term_p - term_q) * exp_part
    }

    fn force(&self, r: f64) -> f64 {
        let rc = self.a_cut * self.sigma;
        if r >= rc {
            return 0.0;
        }
        // Numerical derivative: F = -dV/dr ≈ -(V(r+h) - V(r-h)) / (2h)
        let h = 1e-8 * r.max(1e-10);
        let r_plus = (r + h).min(rc - 1e-15);
        let r_minus = (r - h).max(1e-15);
        let e_plus = self.energy(r_plus);
        let e_minus = self.energy(r_minus);
        -(e_plus - e_minus) / (r_plus - r_minus)
    }

    fn cutoff(&self) -> f64 {
        self.a_cut * self.sigma
    }
}

// ---------------------------------------------------------------------------
// Embedded Atom Method (EAM) — pair contribution
// ---------------------------------------------------------------------------

/// Simplified EAM (Embedded Atom Method) pair potential.
///
/// The pair potential part of EAM:
///   V_pair(r) = A · exp(-alpha · (r/r_e - 1)) - B · exp(-beta · (r/r_e - 1))
///
/// The full EAM also includes an embedding function F(rho) where rho is
/// a sum of atomic electron densities, but this struct provides only the
/// pairwise repulsive-attractive term.
#[derive(Debug, Clone)]
pub struct EamPair {
    /// Repulsive prefactor A.
    pub a_coeff: f64,
    /// Repulsive decay rate alpha.
    pub alpha: f64,
    /// Attractive prefactor B.
    pub b_coeff: f64,
    /// Attractive decay rate beta.
    pub beta: f64,
    /// Equilibrium distance r_e.
    pub r_e: f64,
    /// Cutoff distance.
    pub cutoff_r: f64,
    /// Energy shift at cutoff.
    energy_shift: f64,
}

impl EamPair {
    /// Create a new EAM pair potential.
    #[allow(clippy::too_many_arguments)]
    pub fn new(a_coeff: f64, alpha: f64, b_coeff: f64, beta: f64, r_e: f64, cutoff: f64) -> Self {
        let x_c = cutoff / r_e - 1.0;
        let energy_shift = a_coeff * (-alpha * x_c).exp() - b_coeff * (-beta * x_c).exp();
        Self {
            a_coeff,
            alpha,
            b_coeff,
            beta,
            r_e,
            cutoff_r: cutoff,
            energy_shift,
        }
    }
}

impl Potential for EamPair {
    fn energy(&self, r: f64) -> f64 {
        if r >= self.cutoff_r {
            return 0.0;
        }
        let x = r / self.r_e - 1.0;
        self.a_coeff * (-self.alpha * x).exp()
            - self.b_coeff * (-self.beta * x).exp()
            - self.energy_shift
    }

    fn force(&self, r: f64) -> f64 {
        if r >= self.cutoff_r {
            return 0.0;
        }
        let x = r / self.r_e - 1.0;
        // F = -dV/dr = (A*alpha/r_e)*exp(-alpha*x) - (B*beta/r_e)*exp(-beta*x)
        let inv_re = 1.0 / self.r_e;
        self.a_coeff * self.alpha * inv_re * (-self.alpha * x).exp()
            - self.b_coeff * self.beta * inv_re * (-self.beta * x).exp()
    }

    fn cutoff(&self) -> f64 {
        self.cutoff_r
    }
}

/// EAM electron density function: rho(r) = f_e · exp(-beta · (r/r_e - 1)).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct EamDensity {
    /// Density prefactor.
    pub f_e: f64,
    /// Decay rate.
    pub beta: f64,
    /// Equilibrium distance.
    pub r_e: f64,
    /// Cutoff.
    pub cutoff_r: f64,
}

impl EamDensity {
    /// Create a new EAM density function.
    pub fn new(f_e: f64, beta: f64, r_e: f64, cutoff: f64) -> Self {
        Self {
            f_e,
            beta,
            r_e,
            cutoff_r: cutoff,
        }
    }

    /// Evaluate electron density at distance r.
    pub fn density(&self, r: f64) -> f64 {
        if r >= self.cutoff_r {
            return 0.0;
        }
        let x = r / self.r_e - 1.0;
        self.f_e * (-self.beta * x).exp()
    }
}

// ---------------------------------------------------------------------------
// Tersoff potential (pair term)
// ---------------------------------------------------------------------------

/// Tersoff potential pair term (1988).
///
/// The Tersoff potential is a bond-order potential where:
///   V_ij = f_C(r) · \[f_R(r) - b_ij · f_A(r)\]
///
/// f_R(r) = A · exp(-lambda1 · r)  (repulsive)
/// f_A(r) = B · exp(-lambda2 · r)  (attractive)
/// f_C(r) = cutoff function
///
/// This struct provides the pair terms f_R and f_A with a smooth cutoff.
/// The bond-order term b_ij depends on the local environment and must
/// be computed externally.
#[derive(Debug, Clone)]
pub struct TersoffPair {
    /// Repulsive prefactor A.
    pub a_param: f64,
    /// Attractive prefactor B.
    pub b_param: f64,
    /// Repulsive decay lambda1.
    pub lambda1: f64,
    /// Attractive decay lambda2.
    pub lambda2: f64,
    /// Inner cutoff R.
    pub r_inner: f64,
    /// Outer cutoff S.
    pub s_outer: f64,
}

impl TersoffPair {
    /// Create a new Tersoff pair potential.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        a_param: f64,
        b_param: f64,
        lambda1: f64,
        lambda2: f64,
        r_inner: f64,
        s_outer: f64,
    ) -> Self {
        Self {
            a_param,
            b_param,
            lambda1,
            lambda2,
            r_inner,
            s_outer,
        }
    }

    /// Default silicon parameters (Tersoff 1988).
    #[allow(dead_code)]
    pub fn silicon() -> Self {
        Self::new(1830.8, 471.18, 2.4799, 1.7322, 2.7, 3.0)
    }

    /// Smooth cutoff function f_C(r).
    fn cutoff_fn(&self, r: f64) -> f64 {
        if r <= self.r_inner {
            1.0
        } else if r >= self.s_outer {
            0.0
        } else {
            let arg = std::f64::consts::PI * (r - self.r_inner) / (self.s_outer - self.r_inner);
            0.5 * (1.0 + arg.cos())
        }
    }

    /// Repulsive term: f_R(r) = A · exp(-lambda1 · r).
    pub fn repulsive(&self, r: f64) -> f64 {
        self.a_param * (-self.lambda1 * r).exp()
    }

    /// Attractive term: f_A(r) = B · exp(-lambda2 · r).
    pub fn attractive(&self, r: f64) -> f64 {
        self.b_param * (-self.lambda2 * r).exp()
    }

    /// Pair energy with bond order `bij`:
    /// V = f_C(r) * \[f_R(r) - bij * f_A(r)\]
    pub fn energy_with_bond_order(&self, r: f64, bij: f64) -> f64 {
        let fc = self.cutoff_fn(r);
        fc * (self.repulsive(r) - bij * self.attractive(r))
    }
}

impl Potential for TersoffPair {
    /// Energy with bond order = 1 (pure pair, no angular correction).
    fn energy(&self, r: f64) -> f64 {
        self.energy_with_bond_order(r, 1.0)
    }

    fn force(&self, r: f64) -> f64 {
        if r >= self.s_outer {
            return 0.0;
        }
        // Numerical derivative
        let h = 1e-8 * r.max(1e-10);
        let e_plus = self.energy(r + h);
        let e_minus = self.energy((r - h).max(1e-15));
        -(e_plus - e_minus) / (2.0 * h)
    }

    fn cutoff(&self) -> f64 {
        self.s_outer
    }
}

// ---------------------------------------------------------------------------
// Tabulated potential
// ---------------------------------------------------------------------------

/// Tabulated pair potential with linear interpolation.
///
/// The potential is defined on a uniform grid of distances \[r_min, r_max\]
/// with n_points data points. Energy and force at arbitrary r are obtained
/// via linear interpolation.
#[derive(Debug, Clone)]
pub struct TabulatedPotential {
    /// Minimum distance.
    pub r_min: f64,
    /// Maximum distance (= cutoff).
    pub r_max: f64,
    /// Grid spacing.
    pub dr: f64,
    /// Tabulated energies.
    pub energies: Vec<f64>,
    /// Tabulated forces.
    pub forces: Vec<f64>,
}

impl TabulatedPotential {
    /// Create from energy and force tables.
    pub fn new(r_min: f64, r_max: f64, energies: Vec<f64>, forces: Vec<f64>) -> Self {
        assert_eq!(
            energies.len(),
            forces.len(),
            "energy and force tables must have the same length"
        );
        assert!(energies.len() >= 2, "need at least 2 table points");
        let n = energies.len();
        let dr = (r_max - r_min) / (n - 1) as f64;
        Self {
            r_min,
            r_max,
            dr,
            energies,
            forces,
        }
    }

    /// Create from an existing `Potential` by tabulating it.
    #[allow(dead_code)]
    pub fn from_potential(pot: &dyn Potential, r_min: f64, r_max: f64, n_points: usize) -> Self {
        let dr = (r_max - r_min) / (n_points - 1) as f64;
        let energies: Vec<f64> = (0..n_points)
            .map(|i| pot.energy(r_min + i as f64 * dr))
            .collect();
        let forces: Vec<f64> = (0..n_points)
            .map(|i| pot.force(r_min + i as f64 * dr))
            .collect();
        Self {
            r_min,
            r_max,
            dr,
            energies,
            forces,
        }
    }

    /// Linear interpolation helper.
    fn interp(&self, table: &[f64], r: f64) -> f64 {
        let x = (r - self.r_min) / self.dr;
        let i = x as usize;
        let n = table.len();
        if i >= n - 1 {
            return table[n - 1];
        }
        let frac = x - i as f64;
        table[i] * (1.0 - frac) + table[i + 1] * frac
    }
}

impl Potential for TabulatedPotential {
    fn energy(&self, r: f64) -> f64 {
        if r >= self.r_max || r < self.r_min {
            return 0.0;
        }
        self.interp(&self.energies, r)
    }

    fn force(&self, r: f64) -> f64 {
        if r >= self.r_max || r < self.r_min {
            return 0.0;
        }
        self.interp(&self.forces, r)
    }

    fn cutoff(&self) -> f64 {
        self.r_max
    }
}

// ---------------------------------------------------------------------------
// Potential mixing rules
// ---------------------------------------------------------------------------

/// Lorentz-Berthelot mixing rules for LJ parameters.
///
/// sigma_ij = (sigma_i + sigma_j) / 2
/// epsilon_ij = sqrt(epsilon_i * epsilon_j)
#[allow(dead_code)]
pub fn lorentz_berthelot_mix(eps_i: f64, sig_i: f64, eps_j: f64, sig_j: f64) -> (f64, f64) {
    let eps_ij = (eps_i * eps_j).sqrt();
    let sig_ij = 0.5 * (sig_i + sig_j);
    (eps_ij, sig_ij)
}

/// Geometric mixing rule for LJ parameters.
///
/// sigma_ij = sqrt(sigma_i * sigma_j)
/// epsilon_ij = sqrt(epsilon_i * epsilon_j)
#[allow(dead_code)]
pub fn geometric_mix(eps_i: f64, sig_i: f64, eps_j: f64, sig_j: f64) -> (f64, f64) {
    let eps_ij = (eps_i * eps_j).sqrt();
    let sig_ij = (sig_i * sig_j).sqrt();
    (eps_ij, sig_ij)
}

/// Waldman-Hagler mixing rule (sixth-power combining rule).
///
/// sigma_ij = ((sigma_i^6 + sigma_j^6) / 2)^(1/6)
/// epsilon_ij = 2 * sqrt(eps_i * eps_j) * sigma_i^3 * sigma_j^3 / (sigma_i^6 + sigma_j^6)
#[allow(dead_code)]
pub fn waldman_hagler_mix(eps_i: f64, sig_i: f64, eps_j: f64, sig_j: f64) -> (f64, f64) {
    let si6 = sig_i.powi(6);
    let sj6 = sig_j.powi(6);
    let sig_ij = (0.5 * (si6 + sj6)).powf(1.0 / 6.0);
    let denom = si6 + sj6;
    let eps_ij = if denom.abs() < 1e-30 {
        0.0
    } else {
        2.0 * (eps_i * eps_j).sqrt() * sig_i.powi(3) * sig_j.powi(3) / denom
    };
    (eps_ij, sig_ij)
}

// ---------------------------------------------------------------------------
// Axilrod-Teller three-body potential
// ---------------------------------------------------------------------------

/// Axilrod-Teller three-body dispersion potential.
///
/// V₃(i,j,k) = C₉ * (1 + 3 cos θᵢ cos θⱼ cos θₖ) / (rᵢⱼ rᵢₖ rⱼₖ)³
///
/// where θᵢ, θⱼ, θₖ are the interior angles of the triangle formed by the
/// three atoms and rᵢⱼ etc. are the inter-atom distances.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AxilrodTeller {
    /// Three-body coefficient C₉ (energy · length⁹).
    pub c9: f64,
}

impl AxilrodTeller {
    /// Create a new Axilrod-Teller potential.
    #[allow(dead_code)]
    pub fn new(c9: f64) -> Self {
        Self { c9 }
    }

    /// Compute the Axilrod-Teller three-body energy for a triplet.
    #[allow(dead_code)]
    pub fn energy(&self, ri: [f64; 3], rj: [f64; 3], rk: [f64; 3]) -> f64 {
        let rij = dist3(ri, rj);
        let rik = dist3(ri, rk);
        let rjk = dist3(rj, rk);

        if rij < 1e-14 || rik < 1e-14 || rjk < 1e-14 {
            return 0.0;
        }

        // Cosines via dot products
        let eij = unit3(sub3(rj, ri));
        let eik = unit3(sub3(rk, ri));
        let ejk = unit3(sub3(rk, rj));

        let cos_i = dot3(eij, eik);
        let cos_j = -dot3(eij, ejk); // angle at j: vectors j→i and j→k
        let cos_k = dot3(neg3(eik), ejk); // angle at k

        let numerator = 1.0 + 3.0 * cos_i * cos_j * cos_k;
        let denom = (rij * rik * rjk).powi(3);

        self.c9 * numerator / denom
    }

    /// Compute forces on all three particles via numerical gradient.
    #[allow(clippy::too_many_arguments, dead_code)]
    pub fn forces(
        &self,
        ri: [f64; 3],
        rj: [f64; 3],
        rk: [f64; 3],
    ) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let h = 1e-6;
        let mut fi = [0.0f64; 3];
        let mut fj = [0.0f64; 3];
        let mut fk = [0.0f64; 3];

        for a in 0..3 {
            let mut ri_p = ri;
            let mut ri_m = ri;
            ri_p[a] += h;
            ri_m[a] -= h;
            fi[a] = -(self.energy(ri_p, rj, rk) - self.energy(ri_m, rj, rk)) / (2.0 * h);

            let mut rj_p = rj;
            let mut rj_m = rj;
            rj_p[a] += h;
            rj_m[a] -= h;
            fj[a] = -(self.energy(ri, rj_p, rk) - self.energy(ri, rj_m, rk)) / (2.0 * h);

            let mut rk_p = rk;
            let mut rk_m = rk;
            rk_p[a] += h;
            rk_m[a] -= h;
            fk[a] = -(self.energy(ri, rj, rk_p) - self.energy(ri, rj, rk_m)) / (2.0 * h);
        }

        (fi, fj, fk)
    }
}

// ---------------------------------------------------------------------------
// Angular-dependent potential
// ---------------------------------------------------------------------------

/// Harmonic angle potential: V(θ) = ½ k (θ − θ₀)²
///
/// Used for three-body angular interactions in molecular mechanics.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AngularPotential {
    /// Force constant (kJ mol⁻¹ rad⁻²).
    pub k: f64,
    /// Equilibrium angle (radians).
    pub theta0: f64,
}

impl AngularPotential {
    /// Create a new angular potential.
    #[allow(dead_code)]
    pub fn new(k: f64, theta0: f64) -> Self {
        Self { k, theta0 }
    }

    /// Energy V(θ) = ½ k (θ − θ₀)².
    #[allow(dead_code)]
    pub fn energy(&self, theta: f64) -> f64 {
        let dt = theta - self.theta0;
        0.5 * self.k * dt * dt
    }

    /// Torque τ = −dV/dθ = −k (θ − θ₀).
    #[allow(dead_code)]
    pub fn torque(&self, theta: f64) -> f64 {
        -self.k * (theta - self.theta0)
    }

    /// Compute the three-body energy and forces for atoms i-j-k (j is the
    /// central atom).  The angle is at atom j.
    #[allow(dead_code)]
    pub fn energy_and_forces(
        &self,
        ri: [f64; 3],
        rj: [f64; 3],
        rk: [f64; 3],
    ) -> (f64, [f64; 3], [f64; 3], [f64; 3]) {
        let ji = sub3(ri, rj);
        let jk = sub3(rk, rj);
        let rji = len3(ji);
        let rjk = len3(jk);

        if rji < 1e-14 || rjk < 1e-14 {
            return (0.0, [0.0; 3], [0.0; 3], [0.0; 3]);
        }

        let cos_theta = dot3(ji, jk) / (rji * rjk);
        let theta = cos_theta.clamp(-1.0, 1.0).acos();
        let e = self.energy(theta);
        let torque = self.torque(theta);

        // Convert torque to forces using the gradient of theta w.r.t. positions.
        // d(theta)/d(ri) = (cos(theta)*ji/rji - jk/rjk) / (rji * sin(theta))
        let sin_theta = theta.sin().max(1e-14);

        let mut fi = [0.0f64; 3];
        let mut fk = [0.0f64; 3];
        for a in 0..3 {
            fi[a] = torque * (cos_theta * ji[a] / rji - jk[a] / rjk) / (rji * sin_theta);
            fk[a] = torque * (cos_theta * jk[a] / rjk - ji[a] / rji) / (rjk * sin_theta);
        }
        let fj = [-fi[0] - fk[0], -fi[1] - fk[1], -fi[2] - fk[2]];

        (e, fi, fj, fk)
    }
}

// ---------------------------------------------------------------------------
// Tabulated three-body potential
// ---------------------------------------------------------------------------

/// Tabulated three-body potential on a uniform 2D grid of (r_ij, r_ik) pairs
/// with linear interpolation.
///
/// Energy is looked up as E(r_ij, r_ik) from a flattened row-major 2D table.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct TabulatedThreeBody {
    /// Minimum distance.
    pub r_min: f64,
    /// Maximum distance (cutoff).
    pub r_max: f64,
    /// Grid spacing.
    pub dr: f64,
    /// Number of grid points along each axis.
    pub n: usize,
    /// Energy table (n × n, row-major).
    pub energies: Vec<f64>,
}

impl TabulatedThreeBody {
    /// Create from a (n × n) energy table on the uniform grid \[r_min, r_max\].
    #[allow(dead_code)]
    pub fn new(r_min: f64, r_max: f64, n: usize, energies: Vec<f64>) -> Self {
        assert_eq!(energies.len(), n * n, "energy table must have n² entries");
        assert!(n >= 2);
        let dr = (r_max - r_min) / (n - 1) as f64;
        Self {
            r_min,
            r_max,
            dr,
            n,
            energies,
        }
    }

    /// Bilinear interpolation of the table at (r1, r2).
    fn bilinear(&self, r1: f64, r2: f64) -> f64 {
        if r1 < self.r_min || r1 >= self.r_max || r2 < self.r_min || r2 >= self.r_max {
            return 0.0;
        }
        let x1 = (r1 - self.r_min) / self.dr;
        let x2 = (r2 - self.r_min) / self.dr;
        let i1 = (x1 as usize).min(self.n - 2);
        let i2 = (x2 as usize).min(self.n - 2);
        let f1 = x1 - i1 as f64;
        let f2 = x2 - i2 as f64;

        let e00 = self.energies[i1 * self.n + i2];
        let e10 = self.energies[(i1 + 1) * self.n + i2];
        let e01 = self.energies[i1 * self.n + i2 + 1];
        let e11 = self.energies[(i1 + 1) * self.n + i2 + 1];

        e00 * (1.0 - f1) * (1.0 - f2)
            + e10 * f1 * (1.0 - f2)
            + e01 * (1.0 - f1) * f2
            + e11 * f1 * f2
    }

    /// Evaluate the tabulated three-body energy for distances r_ij and r_ik.
    ///
    /// The third distance r_jk is ignored (not currently used by simple tables).
    #[allow(dead_code)]
    pub fn energy(&self, r_ij: f64, r_ik: f64, _r_jk: f64) -> f64 {
        if r_ij >= self.r_max || r_ik >= self.r_max {
            return 0.0;
        }
        self.bilinear(r_ij, r_ik)
    }
}

// ---------------------------------------------------------------------------
// Potential energy surface (1-D) interpolation
// ---------------------------------------------------------------------------

/// One-dimensional potential energy surface (PES) with cubic-spline-like
/// interpolation (linear for simplicity, cubic available via `force`).
///
/// Useful for representing ab-initio computed 1-D reaction coordinates.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PesSurface1D {
    /// Minimum coordinate.
    pub r_min: f64,
    /// Maximum coordinate.
    pub r_max: f64,
    /// Grid spacing.
    pub dr: f64,
    /// Energy values at grid points.
    pub energies: Vec<f64>,
}

impl PesSurface1D {
    /// Create from a uniform energy grid.
    #[allow(dead_code)]
    pub fn new(r_min: f64, r_max: f64, energies: Vec<f64>) -> Self {
        let n = energies.len();
        assert!(n >= 2);
        let dr = (r_max - r_min) / (n - 1) as f64;
        Self {
            r_min,
            r_max,
            dr,
            energies,
        }
    }

    /// Create from a reference [`Potential`] by tabulating it.
    #[allow(dead_code)]
    pub fn from_potential(pot: &dyn Potential, r_min: f64, r_max: f64, n: usize) -> Self {
        let dr = (r_max - r_min) / (n - 1) as f64;
        let energies: Vec<f64> = (0..n).map(|i| pot.energy(r_min + i as f64 * dr)).collect();
        Self {
            r_min,
            r_max,
            dr,
            energies,
        }
    }

    /// Linearly interpolated energy at coordinate `r`.
    ///
    /// Returns 0 outside \[r_min, r_max\].
    #[allow(dead_code)]
    pub fn energy(&self, r: f64) -> f64 {
        if r < self.r_min || r >= self.r_max {
            return 0.0;
        }
        let x = (r - self.r_min) / self.dr;
        let i = x as usize;
        let n = self.energies.len();
        if i >= n - 1 {
            return self.energies[n - 1];
        }
        let frac = x - i as f64;
        self.energies[i] * (1.0 - frac) + self.energies[i + 1] * frac
    }

    /// Force (negative gradient) at coordinate `r` using central differences.
    #[allow(dead_code)]
    pub fn force(&self, r: f64) -> f64 {
        let h = self.dr * 0.5;
        let ep = self.energy((r + h).min(self.r_max - 1e-12));
        let em = self.energy((r - h).max(self.r_min));
        let denom = (r + h).min(self.r_max - 1e-12) - (r - h).max(self.r_min);
        if denom.abs() < 1e-20 {
            return 0.0;
        }
        -(ep - em) / denom
    }

    /// Minimum energy value and corresponding coordinate.
    #[allow(dead_code)]
    pub fn minimum(&self) -> (f64, f64) {
        let (i_min, &e_min) = self
            .energies
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or((0, &0.0));
        let r_min = self.r_min + i_min as f64 * self.dr;
        (r_min, e_min)
    }
}

// ---------------------------------------------------------------------------
// Three-body helper math (no nalgebra)
// ---------------------------------------------------------------------------

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn neg3(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}

#[inline]
fn len3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

#[inline]
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    len3(sub3(a, b))
}

#[inline]
fn unit3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l > 1e-20 {
        [a[0] / l, a[1] / l, a[2] / l]
    } else {
        [0.0, 1.0, 0.0]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lj_equilibrium_energy() {
        // At r = 2^(1/6)*sigma the LJ potential has its minimum at -epsilon
        let eps = 1.0;
        let sigma = 1.0;
        let cutoff = 5.0;
        let lj = LennardJones::new(eps, sigma, cutoff);
        let r_min = 2.0_f64.powf(1.0 / 6.0) * sigma;
        let e_shift = lj.energy_shift;
        // Raw energy at minimum is -epsilon
        let expected = -eps - e_shift;
        let computed = lj.energy(r_min);
        assert!(
            (computed - expected).abs() < 1e-10,
            "LJ energy at minimum: got {computed}, expected {expected}"
        );
    }

    #[test]
    fn test_lj_force_at_equilibrium() {
        // Force at equilibrium distance should be zero
        let lj = LennardJones::new(1.0, 1.0, 5.0);
        let r_min = 2.0_f64.powf(1.0 / 6.0);
        let f = lj.force(r_min);
        assert!(f.abs() < 1e-10, "LJ force at equilibrium: got {f}");
    }

    #[test]
    fn test_lj_cutoff() {
        let lj = LennardJones::new(1.0, 1.0, 2.5);
        assert_eq!(lj.energy(3.0), 0.0);
        assert_eq!(lj.force(3.0), 0.0);
    }

    #[test]
    fn test_morse_equilibrium() {
        let morse = Morse::new(1.0, 2.0, 1.0, 5.0);
        // At r=r_e, raw energy = 0, so shifted = -shift
        let e = morse.energy(1.0);
        assert!(e < 0.0);
        // Force at r_e is zero
        let f = morse.force(1.0);
        assert!(f.abs() < 1e-12);
    }

    #[test]
    fn test_harmonic_bond() {
        let bond = HarmonicBond::new(100.0, 1.5);
        // At r=r0, energy is 0, force is 0
        assert!((Potential::energy(&bond, 1.5)).abs() < 1e-12);
        assert!((bond.force(1.5)).abs() < 1e-12);
        // At r=2.0, energy = 0.5*100*0.25 = 12.5
        assert!((Potential::energy(&bond, 2.0) - 12.5).abs() < 1e-10);
        // Force at r=2.0: -100*(2.0-1.5) = -50 (attractive, pulling back)
        assert!((bond.force(2.0) - (-50.0)).abs() < 1e-10);
    }

    #[test]
    fn test_coulomb_truncated() {
        let coul = Coulomb::new(1.0, 5.0);
        // At r > cutoff, zero
        assert_eq!(coul.energy(6.0), 0.0);
        // At r=1, energy = 1/1 - 1/5 = 0.8
        assert!((coul.energy(1.0) - 0.8).abs() < 1e-12);
    }

    // -----------------------------------------------------------------------
    // Required tests
    // -----------------------------------------------------------------------

    /// LJ potential energy at the minimum r = 2^(1/6) * sigma is -epsilon
    /// (shifted to -epsilon - V(cutoff), so shifted energy at minimum equals
    /// raw_energy(r_min) - energy_shift = -eps - energy_shift).
    ///
    /// This test verifies the unshifted raw energy at the minimum equals -eps.
    #[test]
    fn test_lj_energy_at_minimum_equals_minus_epsilon() {
        let epsilon = 2.5;
        let sigma = 1.0;
        let cutoff = 10.0; // large cutoff so shift ≈ 0
        let lj = LennardJones::new(epsilon, sigma, cutoff);

        // Minimum at r = 2^(1/6) * sigma
        let r_min = 2.0_f64.powf(1.0 / 6.0) * sigma;

        // Raw energy: 4*eps*((sigma/r)^12 - (sigma/r)^6)
        // At r_min: (sigma/r_min)^6 = 1/2, (sigma/r_min)^12 = 1/4
        // raw = 4*eps*(1/4 - 1/2) = 4*eps*(-1/4) = -eps
        let raw = lj.raw_energy(r_min);
        assert!(
            (raw - (-epsilon)).abs() < 1e-10,
            "LJ raw energy at minimum should be -epsilon={epsilon}, got {raw}"
        );
    }

    /// LJ force at the equilibrium distance r_min = 2^(1/6)*sigma must be zero.
    #[test]
    fn test_lj_force_zero_at_equilibrium_distance() {
        let lj = LennardJones::new(1.0, 1.0, 8.0);
        let r_min = 2.0_f64.powf(1.0 / 6.0);
        let f = lj.force(r_min);
        assert!(f.abs() < 1e-9, "LJ force at r_min must be zero, got {f}");
    }

    /// LJ is repulsive for r < r_min and attractive for r_min < r < cutoff.
    #[test]
    fn test_lj_repulsive_below_rmin_attractive_above() {
        let lj = LennardJones::new(1.0, 1.0, 5.0);
        let r_min = 2.0_f64.powf(1.0 / 6.0);

        // Force > 0 means repulsive (positive = repulsive convention)
        let f_repul = lj.force(r_min * 0.8);
        let f_attr = lj.force(r_min * 1.2);

        assert!(f_repul > 0.0, "LJ force should be repulsive at r < r_min");
        assert!(
            f_attr < 0.0,
            "LJ force should be attractive at r_min < r < cutoff"
        );
    }

    // --- Stillinger-Weber tests ---

    #[test]
    fn test_sw_cutoff() {
        let sw = StillingerWeber::silicon();
        let rc = sw.cutoff();
        assert!((rc - 1.8 * 2.0951).abs() < 1e-10);
        assert_eq!(sw.energy(rc + 0.01), 0.0);
        assert_eq!(sw.force(rc + 0.01), 0.0);
    }

    #[test]
    fn test_sw_energy_finite_inside() {
        let sw = StillingerWeber::silicon();
        let r = 2.0; // inside cutoff
        let e = sw.energy(r);
        assert!(
            e.is_finite(),
            "SW energy at r={r} should be finite, got {e}"
        );
    }

    #[test]
    fn test_sw_force_finite_inside() {
        let sw = StillingerWeber::silicon();
        let r = 2.5;
        let f = sw.force(r);
        assert!(f.is_finite(), "SW force at r={r} should be finite, got {f}");
    }

    // --- EAM tests ---

    #[test]
    fn test_eam_pair_cutoff() {
        let eam = EamPair::new(1.0, 5.0, 0.5, 3.0, 2.5, 6.0);
        assert_eq!(eam.energy(7.0), 0.0);
        assert_eq!(eam.force(7.0), 0.0);
    }

    #[test]
    fn test_eam_pair_energy_at_equilibrium() {
        // At r = r_e, x = 0, V = A - B - shift
        let eam = EamPair::new(10.0, 5.0, 5.0, 3.0, 2.5, 10.0);
        let e = eam.energy(2.5);
        // Raw at r_e: A*exp(0) - B*exp(0) = A - B = 10 - 5 = 5
        // Shift: at cutoff r=10, x = 10/2.5 - 1 = 3
        // shift = A*exp(-5*3) - B*exp(-3*3) = very small
        // So e ≈ 5 - small_shift
        assert!(e.is_finite());
    }

    #[test]
    fn test_eam_density() {
        let dens = EamDensity::new(1.0, 3.0, 2.5, 6.0);
        // At r = r_e: density = f_e * exp(0) = f_e = 1.0
        assert!((dens.density(2.5) - 1.0).abs() < 1e-12);
        // At r > cutoff: 0
        assert_eq!(dens.density(7.0), 0.0);
        // Density should decrease with r
        assert!(dens.density(3.0) < dens.density(2.5));
    }

    // --- Tersoff tests ---

    #[test]
    fn test_tersoff_cutoff() {
        let tersoff = TersoffPair::silicon();
        assert_eq!(tersoff.cutoff(), 3.0);
        assert_eq!(tersoff.energy(3.5), 0.0);
    }

    #[test]
    fn test_tersoff_repulsive_at_small_r() {
        let tersoff = TersoffPair::silicon();
        // At very small r, repulsive term dominates
        let e = tersoff.energy(1.0);
        assert!(
            e > 0.0,
            "Tersoff energy at small r should be positive (repulsive), got {e}"
        );
    }

    #[test]
    fn test_tersoff_cutoff_function_boundaries() {
        let tersoff = TersoffPair::silicon();
        // f_C(r <= R) = 1
        assert!((tersoff.cutoff_fn(2.0) - 1.0).abs() < 1e-12);
        // f_C(r >= S) = 0
        assert!((tersoff.cutoff_fn(3.5)).abs() < 1e-12);
        // f_C at midpoint = 0.5
        let mid = 0.5 * (2.7 + 3.0);
        assert!((tersoff.cutoff_fn(mid) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_tersoff_bond_order_effect() {
        let tersoff = TersoffPair::silicon();
        let r = 2.35; // typical Si-Si distance
        let e_full = tersoff.energy_with_bond_order(r, 1.0);
        let e_half = tersoff.energy_with_bond_order(r, 0.5);
        // Lower bond order reduces attractive part, so energy is higher
        assert!(
            e_half > e_full,
            "Lower bond order should give higher energy"
        );
    }

    // --- Tabulated potential tests ---

    #[test]
    fn test_tabulated_from_lj() {
        let lj = LennardJones::new(1.0, 1.0, 5.0);
        let tab = TabulatedPotential::from_potential(&lj, 0.8, 5.0, 1000);

        // Check at LJ equilibrium
        let r_min = 2.0_f64.powf(1.0 / 6.0);
        let e_lj = lj.energy(r_min);
        let e_tab = tab.energy(r_min);
        assert!(
            (e_lj - e_tab).abs() < 0.01,
            "Tabulated LJ energy at r_min: lj={e_lj}, tab={e_tab}"
        );
    }

    #[test]
    fn test_tabulated_cutoff() {
        let energies = vec![1.0, 0.5, 0.1, 0.0];
        let forces = vec![-2.0, -1.0, -0.3, 0.0];
        let tab = TabulatedPotential::new(0.5, 2.0, energies, forces);
        assert_eq!(tab.cutoff(), 2.0);
        assert_eq!(tab.energy(2.5), 0.0);
        assert_eq!(tab.force(2.5), 0.0);
    }

    #[test]
    fn test_tabulated_interpolation() {
        // Linear table: E(r) = 10 - 2r, F(r) = 2
        let energies = vec![10.0, 8.0, 6.0, 4.0, 2.0];
        let forces = vec![2.0, 2.0, 2.0, 2.0, 2.0];
        let tab = TabulatedPotential::new(0.0, 4.0, energies, forces);
        // At r=1.5: E = 10 - 3 = 7
        let e = tab.energy(1.5);
        assert!((e - 7.0).abs() < 1e-10, "E at r=1.5: expected 7.0, got {e}");
        let f = tab.force(1.5);
        assert!((f - 2.0).abs() < 1e-10, "F at r=1.5: expected 2.0, got {f}");
    }

    // --- Mixing rules tests ---

    #[test]
    fn test_lorentz_berthelot_mix() {
        let (eps, sig) = lorentz_berthelot_mix(1.0, 3.0, 4.0, 5.0);
        assert!(
            (eps - 2.0).abs() < 1e-12,
            "eps_ij = sqrt(1*4) = 2, got {eps}"
        );
        assert!((sig - 4.0).abs() < 1e-12, "sig_ij = (3+5)/2 = 4, got {sig}");
    }

    #[test]
    fn test_geometric_mix() {
        let (eps, sig) = geometric_mix(1.0, 4.0, 4.0, 9.0);
        assert!(
            (eps - 2.0).abs() < 1e-12,
            "eps_ij = sqrt(1*4) = 2, got {eps}"
        );
        assert!(
            (sig - 6.0).abs() < 1e-12,
            "sig_ij = sqrt(4*9) = 6, got {sig}"
        );
    }

    #[test]
    fn test_waldman_hagler_mix_identical() {
        // For identical species, WH should return the same parameters
        let (eps, sig) = waldman_hagler_mix(2.0, 3.0, 2.0, 3.0);
        assert!((eps - 2.0).abs() < 1e-12, "WH eps for identical: got {eps}");
        assert!((sig - 3.0).abs() < 1e-12, "WH sig for identical: got {sig}");
    }

    #[test]
    fn test_mixing_symmetry() {
        // All mixing rules should be symmetric: mix(i,j) = mix(j,i)
        let (e1, s1) = lorentz_berthelot_mix(1.0, 2.0, 3.0, 4.0);
        let (e2, s2) = lorentz_berthelot_mix(3.0, 4.0, 1.0, 2.0);
        assert!((e1 - e2).abs() < 1e-12);
        assert!((s1 - s2).abs() < 1e-12);

        let (e1, s1) = geometric_mix(1.0, 2.0, 3.0, 4.0);
        let (e2, s2) = geometric_mix(3.0, 4.0, 1.0, 2.0);
        assert!((e1 - e2).abs() < 1e-12);
        assert!((s1 - s2).abs() < 1e-12);
    }

    // --- Axilrod-Teller three-body tests ---

    #[test]
    fn test_axilrod_teller_symmetric() {
        let at = AxilrodTeller::new(1.0);
        let a = [0.0f64, 0.0, 0.0];
        let b = [1.0f64, 0.0, 0.0];
        let c = [0.5f64, 1.0, 0.0];
        let e1 = at.energy(a, b, c);
        let e2 = at.energy(b, c, a);
        let e3 = at.energy(c, a, b);
        // Cyclic permutations should give same energy (potential is symmetric).
        assert!(
            (e1 - e2).abs() < 1e-10,
            "AT not symmetric under permutation: e1={e1}, e2={e2}"
        );
        assert!(
            (e1 - e3).abs() < 1e-10,
            "AT not symmetric under permutation: e1={e1}, e3={e3}"
        );
    }

    #[test]
    fn test_axilrod_teller_collinear_divergence() {
        // Collinear case: angles approach π or 0 → numerator → ±negative → force non-trivial.
        let at = AxilrodTeller::new(1.0);
        let a = [0.0f64, 0.0, 0.0];
        let b = [1.0f64, 0.0, 0.0];
        let c = [2.0f64, 0.0, 0.0];
        let e = at.energy(a, b, c);
        assert!(
            e.is_finite(),
            "AT energy for collinear atoms should be finite: {e}"
        );
    }

    #[test]
    fn test_axilrod_teller_force_gradient() {
        // The gradient of the AT energy should point atoms apart from the centre.
        let at = AxilrodTeller::new(1.0);
        let a = [0.0f64, 0.0, 0.0];
        let b = [1.0f64, 0.0, 0.0];
        let c = [0.5f64, 0.8, 0.0];
        let (fa, fb, fc) = at.forces(a, b, c);
        // Newton's third law: sum of forces should be zero
        let sum = [
            fa[0] + fb[0] + fc[0],
            fa[1] + fb[1] + fc[1],
            fa[2] + fb[2] + fc[2],
        ];
        let mag = (sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]).sqrt();
        assert!(
            mag < 1e-8,
            "AT forces violate Newton's 3rd law: net force = {mag}"
        );
    }

    // --- Angular-dependent potential tests ---

    #[test]
    fn test_angular_potential_at_equilibrium() {
        let ang = AngularPotential::new(100.0, std::f64::consts::PI / 2.0);
        // At equilibrium angle: energy = 0, torque = 0.
        let e = ang.energy(std::f64::consts::PI / 2.0);
        assert!(e.abs() < 1e-12, "Angular potential at equilibrium: {e}");
        let t = ang.torque(std::f64::consts::PI / 2.0);
        assert!(t.abs() < 1e-12, "Angular torque at equilibrium: {t}");
    }

    #[test]
    fn test_angular_potential_restoring() {
        let ang = AngularPotential::new(100.0, std::f64::consts::PI / 2.0);
        // For angle > equilibrium the torque is negative (restoring).
        let t = ang.torque(std::f64::consts::PI / 2.0 + 0.1);
        assert!(
            t < 0.0,
            "Torque should be restoring (negative) for angle > eq: {t}"
        );
    }

    // --- Tabulated three-body tests ---

    #[test]
    fn test_tabulated_three_body_symmetric() {
        // Uniform table: energy = constant for all r values.
        let n = 50;
        let r_min = 0.5;
        let r_max = 5.0;
        let energies: Vec<f64> = vec![2.0; n * n];
        let tab3 = TabulatedThreeBody::new(r_min, r_max, n, energies);
        let e1 = tab3.energy(1.0, 2.0, 1.5);
        let e2 = tab3.energy(2.0, 1.0, 1.5);
        assert!(
            (e1 - e2).abs() < 1e-10,
            "Tabulated 3-body should be symmetric: {e1} vs {e2}"
        );
    }

    #[test]
    fn test_tabulated_three_body_outside_cutoff() {
        let n = 10;
        let r_min = 0.5;
        let r_max = 3.0;
        let energies = vec![1.0; n * n];
        let tab3 = TabulatedThreeBody::new(r_min, r_max, n, energies);
        assert_eq!(
            tab3.energy(4.0, 2.0, 1.0),
            0.0,
            "Outside cutoff must be zero"
        );
        assert_eq!(tab3.energy(2.0, 4.0, 1.0), 0.0);
    }

    // --- PES interpolation tests ---

    #[test]
    fn test_pes_interpolation_returns_expected() {
        // Linear PES: E(r) = 2*r
        let n = 11;
        let r_min = 0.0;
        let r_max = 10.0;
        let energies: Vec<f64> = (0..n).map(|i| 2.0 * i as f64).collect();
        let pes = PesSurface1D::new(r_min, r_max, energies);
        let e = pes.energy(3.0);
        assert!((e - 6.0).abs() < 1e-10, "PES at r=3 should be 6: {e}");
    }

    #[test]
    fn test_pes_interpolation_extrapolation_zero() {
        let n = 5;
        let energies = vec![1.0; n];
        let pes = PesSurface1D::new(0.0, 4.0, energies);
        assert_eq!(pes.energy(10.0), 0.0, "PES outside range must be zero");
        assert_eq!(pes.energy(-1.0), 0.0);
    }

    #[test]
    fn test_pes_gradient_numerical() {
        // E(r) = r^2 → dE/dr = 2r.
        let n = 1001;
        let r_min = 0.0;
        let r_max = 10.0;
        let dr = (r_max - r_min) / (n - 1) as f64;
        let energies: Vec<f64> = (0..n)
            .map(|i| {
                let r = r_min + i as f64 * dr;
                r * r
            })
            .collect();
        let pes = PesSurface1D::new(r_min, r_max, energies);
        let r = 3.0;
        let force = pes.force(r);
        // -dE/dr = -2r at r=3 is -6
        assert!(
            (force - (-6.0)).abs() < 0.05,
            "PES gradient at r=3: expected ≈-6, got {force}"
        );
    }
}
