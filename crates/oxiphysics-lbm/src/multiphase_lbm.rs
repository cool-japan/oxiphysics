// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Multiphase LBM: Shan-Chen and free-energy models.
//!
//! This module provides:
//! - [`ShanChenModel`]: Shan-Chen pseudo-potential multiphase LBM
//! - [`FreeEnergyModel`]: Thermodynamically consistent free-energy multiphase LBM
//! - Phase diagram construction via Maxwell equal-area construction
//! - Interface width estimation from tanh profiles

// ---------------------------------------------------------------------------
// Shan-Chen model
// ---------------------------------------------------------------------------

/// Shan-Chen multiphase LBM model parameters.
///
/// The interaction strength `g_coupling` drives phase separation when
/// it exceeds a critical value dependent on the pseudo-potential.
#[derive(Debug, Clone)]
pub struct ShanChenModel {
    /// Interaction coupling constant G (negative = attractive).
    pub g_coupling: f64,
    /// Reference density for the pseudo-potential.
    pub rho_ref: f64,
    /// BGK relaxation frequency for phase 1.
    pub omega1: f64,
    /// BGK relaxation frequency for phase 2.
    pub omega2: f64,
}

impl ShanChenModel {
    /// Create a new Shan-Chen model.
    ///
    /// # Arguments
    /// - `g_coupling`: interaction strength (typically negative)
    /// - `rho_ref`: reference density (e.g. 1.0)
    /// - `omega1`, `omega2`: relaxation frequencies for the two phases
    pub fn new(g_coupling: f64, rho_ref: f64, omega1: f64, omega2: f64) -> Self {
        Self {
            g_coupling,
            rho_ref,
            omega1,
            omega2,
        }
    }

    /// Compute the pseudo-potential psi(rho) = rho_ref * (1 - exp(-rho / rho_ref)).
    ///
    /// This is the original Shan-Chen exponential form.
    pub fn psi_potential(&self, rho: f64) -> f64 {
        psi_exponential(rho, self.rho_ref)
    }

    /// Compute the Shan-Chen interaction force at a node.
    ///
    /// Uses D2Q9 nearest-neighbor stencil.  The force is:
    ///
    /// `F_alpha = -G * psi(rho) * sum_i w_i * psi(rho(x + e_i)) * e_i_alpha`
    ///
    /// # Arguments
    /// - `psi_center`: psi at the current node
    /// - `psi_neighbors`: psi values at the 8 D2Q9 neighbors
    ///   (ordered: E, N, W, S, NE, NW, SW, SE)
    ///
    /// Returns `(fx, fy)`.
    pub fn shan_chen_force(&self, psi_center: f64, psi_neighbors: &[f64; 8]) -> (f64, f64) {
        shan_chen_force_d2q9(self.g_coupling, psi_center, psi_neighbors)
    }
}

/// Shan-Chen exponential pseudo-potential.
///
/// `psi(rho) = rho_0 * (1 - exp(-rho / rho_0))`
pub fn psi_exponential(rho: f64, rho_0: f64) -> f64 {
    rho_0 * (1.0 - (-rho / rho_0.max(1e-30)).exp())
}

/// Shan-Chen linear pseudo-potential: psi(rho) = rho.
pub fn psi_linear(rho: f64) -> f64 {
    rho
}

/// Sukop-Thorne pseudo-potential: psi(rho) = exp(-rho_0 / rho).
pub fn psi_sukop_thorne(rho: f64, rho_0: f64) -> f64 {
    (-rho_0 / rho.max(1e-30)).exp()
}

/// Compute the D2Q9 Shan-Chen interaction force.
///
/// D2Q9 velocity set (excluding rest):
/// - i=0: E  (1,0), w=1/9
/// - i=1: N  (0,1), w=1/9
/// - i=2: W  (-1,0), w=1/9
/// - i=3: S  (0,-1), w=1/9
/// - i=4: NE (1,1), w=1/36
/// - i=5: NW (-1,1), w=1/36
/// - i=6: SW (-1,-1), w=1/36
/// - i=7: SE (1,-1), w=1/36
///
/// # Arguments
/// - `g`: coupling constant
/// - `psi_c`: psi at center node
/// - `psi_nb`: psi at 8 neighbors in the order above
///
/// Returns `(fx, fy)`.
pub fn shan_chen_force_d2q9(g: f64, psi_c: f64, psi_nb: &[f64; 8]) -> (f64, f64) {
    // Velocity directions: (ex, ey, weight)
    let dirs: [(f64, f64, f64); 8] = [
        (1.0, 0.0, 1.0 / 9.0),
        (0.0, 1.0, 1.0 / 9.0),
        (-1.0, 0.0, 1.0 / 9.0),
        (0.0, -1.0, 1.0 / 9.0),
        (1.0, 1.0, 1.0 / 36.0),
        (-1.0, 1.0, 1.0 / 36.0),
        (-1.0, -1.0, 1.0 / 36.0),
        (1.0, -1.0, 1.0 / 36.0),
    ];
    let mut fx = 0.0_f64;
    let mut fy = 0.0_f64;
    for (i, &(ex, ey, w)) in dirs.iter().enumerate() {
        let psi_i = psi_nb[i];
        fx += w * psi_i * ex;
        fy += w * psi_i * ey;
    }
    fx *= -g * psi_c;
    fy *= -g * psi_c;
    (fx, fy)
}

// ---------------------------------------------------------------------------
// Free-energy model
// ---------------------------------------------------------------------------

/// Free-energy multiphase LBM model parameters.
///
/// Based on the van der Waals / Swift *et al.* free-energy model.
#[derive(Debug, Clone)]
pub struct FreeEnergyModel {
    /// Surface tension coefficient kappa (controls interface energy).
    pub kappa: f64,
    /// Bulk free-energy parameter a (controls phase coexistence).
    pub a_coeff: f64,
    /// Bulk free-energy parameter b (excluded volume).
    pub b_coeff: f64,
    /// Temperature (in units of critical temperature).
    pub temperature: f64,
}

impl FreeEnergyModel {
    /// Create a new free-energy model.
    pub fn new(kappa: f64, a_coeff: f64, b_coeff: f64, temperature: f64) -> Self {
        Self {
            kappa,
            a_coeff,
            b_coeff,
            temperature,
        }
    }

    /// Compute the chemical potential mu(rho) from the van der Waals EOS.
    ///
    /// mu = d(f_bulk)/d(rho) = T * b / (b - rho*b_mol) - 2*a*rho
    ///
    /// For the simple Landau form: f = a*(rho - rho_c)^2 + b*(rho - rho_c)^4
    /// mu = 2*a*(rho - rho_c) + 4*b*(rho - rho_c)^3
    pub fn chemical_potential(&self, rho: f64) -> f64 {
        chemical_potential_landau(rho, self.a_coeff, self.b_coeff)
    }

    /// Estimate the interface width from the tanh profile.
    ///
    /// For a planar interface: rho(x) = rho_avg + delta_rho/2 * tanh(x / (2*xi))
    /// where xi = sqrt(kappa / |a|)
    pub fn interface_width(&self) -> f64 {
        interface_width(self.kappa, self.a_coeff)
    }
}

/// Chemical potential from a Landau double-well free energy.
///
/// `f(rho) = a*(rho - rho_c)^2 + b*(rho - rho_c)^4` with rho_c = 0
/// `mu = df/d(rho) = 2*a*rho + 4*b*rho^3`
pub fn chemical_potential_landau(rho: f64, a: f64, b: f64) -> f64 {
    2.0 * a * rho + 4.0 * b * rho * rho * rho
}

/// Chemical potential from van der Waals equation of state.
///
/// `p_vdw = rho*T / (1 - rho*b) - a*rho^2`
///
/// Returns `mu = integral(dp/rho)` simplified:
/// `mu = T * ln(rho/(1-rho*b)) + T*b*rho/(1-rho*b) - 2*a*rho`
pub fn chemical_potential_vdw(rho: f64, temperature: f64, a: f64, b: f64) -> f64 {
    let denom = 1.0 - rho * b;
    if denom <= 1e-12 || rho <= 1e-30 {
        return 0.0;
    }
    vdw_chemical_potential_raw(rho, temperature, a, b)
}

/// Van der Waals equation of state pressure.
///
/// `p = rho * T / (1 - rho*b) - a * rho^2`
pub fn vdw_pressure(rho: f64, temperature: f64, a: f64, b: f64) -> f64 {
    let denom = 1.0 - rho * b;
    if denom <= 1e-12 {
        return 0.0;
    }
    vdw_pressure_raw(rho, temperature, a, b)
}

/// Van der Waals pressure without the near-singularity clamp of [`vdw_pressure`].
///
/// Valid for `0 <= rho < 1/b`; diverges to `+inf` at `rho = 1/b`.
fn vdw_pressure_raw(rho: f64, temperature: f64, a: f64, b: f64) -> f64 {
    rho * temperature / (1.0 - rho * b) - a * rho * rho
}

/// Van der Waals chemical potential without the clamps of [`chemical_potential_vdw`].
///
/// Valid for `0 < rho < 1/b`.
fn vdw_chemical_potential_raw(rho: f64, temperature: f64, a: f64, b: f64) -> f64 {
    let denom = 1.0 - rho * b;
    temperature * (rho / denom).ln() + temperature * b * rho / denom - 2.0 * a * rho
}

/// Slope `dp/drho = T / (1 - rho*b)^2 - 2*a*rho` of the van der Waals isotherm.
fn vdw_pressure_derivative(rho: f64, temperature: f64, a: f64, b: f64) -> f64 {
    let denom = 1.0 - rho * b;
    temperature / (denom * denom) - 2.0 * a * rho
}

/// Interface width parameter xi = sqrt(kappa / |a|).
///
/// For a tanh interface profile: the width is proportional to xi.
pub fn interface_width(kappa: f64, a: f64) -> f64 {
    if a.abs() < 1e-30 {
        return f64::INFINITY;
    }
    (kappa / a.abs()).sqrt()
}

/// Result of the Maxwell construction on a van der Waals isotherm.
///
/// The two coexisting phases have equal pressure and equal chemical potential
/// (equivalently: the areas enclosed by the isotherm above and below the tie
/// line are equal in the `p`-`v` plane).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaxwellCoexistence {
    /// Density of the vapour phase (below the lower spinodal).
    pub rho_vapor: f64,
    /// Density of the liquid phase (above the upper spinodal).
    pub rho_liquid: f64,
    /// Common coexistence (saturation) pressure.
    pub pressure: f64,
    /// Common coexistence chemical potential, in the convention of
    /// [`chemical_potential_vdw`].
    pub chemical_potential: f64,
}

/// Upper bound on the number of bisection halvings used per root.
///
/// Fifty-five halvings already exhaust double precision; anything beyond a few
/// hundred only burns time.
const MAX_BISECTION_STEPS: usize = 200;

/// Bisection for a continuous function with a sign change between `lo` and `hi`.
///
/// The sign of `f(lo)` is taken as the "low side" sign.  Iteration stops after
/// `max_steps` halvings or as soon as the bracket cannot be halved in floating
/// point; the midpoint of the final bracket is returned.
fn bisect_sign_change<F: Fn(f64) -> f64>(f: F, mut lo: f64, mut hi: f64, max_steps: usize) -> f64 {
    let low_side_positive = f(lo) > 0.0;
    for _ in 0..max_steps {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            break;
        }
        if (f(mid) > 0.0) == low_side_positive {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Maxwell construction for the van der Waals isotherm at `temperature`.
///
/// Solves for the vapour and liquid densities `rho_v < rho_l` with
///
/// - equal pressure  `p(rho_v) = p(rho_l) = p*` and
/// - equal chemical potential `mu(rho_v) = mu(rho_l)`
///
/// which is equivalent to the equal-area rule `integral (p - p*) dv = 0`
/// over the specific volume `v = 1/rho`.
///
/// Algorithm: the two spinodal densities (where `dp/drho = 0`) are located by
/// bisection on the analytic slope; the saturation pressure `p*` is then
/// bisected (on a logarithmic scale, so that the exponentially small vapour
/// pressures at low temperature are resolved) between `max(p(rho_spinodal_hi), 0)`
/// and `p(rho_spinodal_lo)`.  For every trial `p*` the outer roots of
/// `p(rho) = p*` on the vapour and liquid branches are found by bisection, and
/// the residual is `mu(rho_l) - mu(rho_v)`, which decreases monotonically
/// with `p*`.
///
/// # Arguments
/// - `temperature`: temperature in the units of the equation of state
///   (must satisfy `0 < T < T_c = 8a/(27b)`)
/// - `a`, `b`: van der Waals attraction and excluded-volume parameters (> 0)
/// - `max_steps`: number of bisection halvings applied to every bracket
///   (at least 4; 55 or more reach double precision, values above 200 are
///   clamped)
///
/// Returns `None` when `T` is not strictly between 0 and `T_c`, when a
/// parameter is not positive and finite, or when the construction fails.
pub fn maxwell_construction(
    temperature: f64,
    a: f64,
    b: f64,
    max_steps: usize,
) -> Option<MaxwellCoexistence> {
    let parameters_valid = temperature.is_finite()
        && a.is_finite()
        && b.is_finite()
        && temperature > 0.0
        && a > 0.0
        && b > 0.0;
    if !parameters_valid || max_steps < 4 {
        return None;
    }
    let steps = max_steps.min(MAX_BISECTION_STEPS);
    // Critical point: rho_c = 1/(3b), T_c = 8a/(27b).
    let t_c = 8.0 * a / (27.0 * b);
    if temperature >= t_c {
        return None;
    }
    let rho_c = 1.0 / (3.0 * b);
    let rho_close = 1.0 / b;
    let slope = |rho: f64| vdw_pressure_derivative(rho, temperature, a, b);
    let pressure = |rho: f64| vdw_pressure_raw(rho, temperature, a, b);
    let potential = |rho: f64| vdw_chemical_potential_raw(rho, temperature, a, b);

    // Spinodals: dp/drho > 0 near rho = 0, < 0 at rho_c (T < T_c), > 0 near 1/b.
    if slope(rho_c) >= 0.0 {
        return None;
    }
    let rho_spinodal_lo = bisect_sign_change(slope, 0.0, rho_c, steps);
    let mut bracket_lo = rho_c;
    let mut gap = 0.5 * (rho_close - rho_c);
    let mut bracket_hi = rho_close - gap;
    let mut attempts = 0usize;
    while slope(bracket_hi) <= 0.0 {
        bracket_lo = bracket_hi;
        gap *= 0.5;
        bracket_hi = rho_close - gap;
        attempts += 1;
        if attempts > 1000 || bracket_hi >= rho_close {
            return None;
        }
    }
    let rho_spinodal_hi = bisect_sign_change(slope, bracket_lo, bracket_hi, steps);
    let p_max = pressure(rho_spinodal_lo);
    let p_min = pressure(rho_spinodal_hi);
    if !(p_max.is_finite() && p_max > 0.0 && p_max > p_min) {
        return None;
    }

    // Outer root of p(rho) = p_star on the vapour branch (0, rho_spinodal_lo):
    // the pressure rises monotonically from 0 to p_max.  The bracket is first
    // tightened to a factor two so that tiny densities keep full relative accuracy.
    let vapor_root = |p_star: f64| -> f64 {
        let mut hi = rho_spinodal_lo;
        let mut lo = rho_spinodal_lo;
        let mut halvings = 0usize;
        while pressure(lo) > p_star && halvings < 1100 {
            hi = lo;
            lo *= 0.5;
            halvings += 1;
        }
        bisect_sign_change(|rho| pressure(rho) - p_star, lo, hi, steps)
    };
    // Outer root on the liquid branch (rho_spinodal_hi, 1/b): the pressure rises
    // monotonically from p_min to +infinity.
    let liquid_root = |p_star: f64| -> f64 {
        let mut lo = rho_spinodal_hi;
        let mut gap = 0.5 * (rho_close - rho_spinodal_hi);
        let mut hi = rho_close - gap;
        let mut halvings = 0usize;
        while pressure(hi) < p_star && halvings < 1100 && hi < rho_close {
            lo = hi;
            gap *= 0.5;
            hi = rho_close - gap;
            halvings += 1;
        }
        bisect_sign_change(|rho| pressure(rho) - p_star, lo, hi, steps)
    };

    // mu(rho_l) - mu(rho_v) as a function of ln(p_star); it decreases with p_star
    // because d(mu_l - mu_v)/dp = 1/rho_l - 1/rho_v < 0.
    let residual = |log_p: f64| -> f64 {
        let p_star = log_p.exp();
        potential(liquid_root(p_star)) - potential(vapor_root(p_star))
    };
    let log_p_hi = p_max.ln();
    let p_floor = p_min.max(p_max * 1e-300).max(f64::MIN_POSITIVE);
    let log_p_lo = p_floor.ln();
    let residual_at_floor = residual(log_p_lo);
    if residual_at_floor.is_nan() || residual_at_floor <= 0.0 {
        return None;
    }
    let p_star = bisect_sign_change(residual, log_p_lo, log_p_hi, steps).exp();

    let rho_vapor = vapor_root(p_star);
    let rho_liquid = liquid_root(p_star);
    let ordered = rho_vapor.is_finite()
        && rho_liquid.is_finite()
        && rho_vapor > 0.0
        && rho_vapor <= rho_spinodal_lo
        && rho_spinodal_lo < rho_c
        && rho_c < rho_spinodal_hi
        && rho_spinodal_hi <= rho_liquid
        && rho_liquid < rho_close;
    if !ordered {
        return None;
    }
    Some(MaxwellCoexistence {
        rho_vapor,
        rho_liquid,
        pressure: p_star,
        chemical_potential: 0.5 * (potential(rho_vapor) + potential(rho_liquid)),
    })
}

/// Coexistence densities of a van der Waals fluid from the Maxwell construction.
///
/// Finds `rho_v` and `rho_l` with equal pressure and equal chemical potential
/// at the given temperature (see [`maxwell_construction`], which also returns
/// the saturation pressure and chemical potential).
///
/// # Arguments
/// - `temperature`: temperature (must satisfy `0 < T < T_c = 8a/(27b)`)
/// - `a`, `b`: vdW parameters (> 0)
/// - `n_points`: resolution of the numerical solve, i.e. the number of
///   bisection halvings applied to every bracket (at least 4; 55 or more
///   reach double precision)
///
/// Returns `(rho_vapor, rho_liquid)` or `(f64::NAN, f64::NAN)` if there is no
/// coexistence (`T >= T_c`, non-positive or non-finite inputs, `n_points < 4`)
/// or the construction fails.
pub fn phase_diagram(temperature: f64, a: f64, b: f64, n_points: usize) -> (f64, f64) {
    match maxwell_construction(temperature, a, b, n_points) {
        Some(coexistence) => (coexistence.rho_vapor, coexistence.rho_liquid),
        None => (f64::NAN, f64::NAN),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- psi potentials ---

    #[test]
    fn test_psi_exponential_zero_density() {
        // psi(0) = rho_0 * (1 - exp(0)) = 0
        let val = psi_exponential(0.0, 1.0);
        assert!(val.abs() < 1e-14);
    }

    #[test]
    fn test_psi_exponential_large_density() {
        // psi -> rho_0 as rho -> infinity
        let val = psi_exponential(1000.0, 1.0);
        assert!((val - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_psi_exponential_positive() {
        for &rho in &[0.1, 0.5, 1.0, 2.0] {
            assert!(psi_exponential(rho, 1.0) >= 0.0);
        }
    }

    #[test]
    fn test_psi_linear() {
        assert_eq!(psi_linear(0.5), 0.5);
        assert_eq!(psi_linear(0.0), 0.0);
        assert_eq!(psi_linear(2.0), 2.0);
    }

    #[test]
    fn test_psi_sukop_thorne_positive() {
        let val = psi_sukop_thorne(1.0, 1.0);
        assert!(val > 0.0 && val <= 1.0);
    }

    #[test]
    fn test_psi_sukop_thorne_increases_with_density() {
        let v1 = psi_sukop_thorne(0.5, 1.0);
        let v2 = psi_sukop_thorne(2.0, 1.0);
        assert!(v2 > v1);
    }

    // --- ShanChenModel ---

    #[test]
    fn test_shan_chen_model_new() {
        let sc = ShanChenModel::new(-1.0, 1.0, 1.0, 1.0);
        assert_eq!(sc.g_coupling, -1.0);
        assert_eq!(sc.rho_ref, 1.0);
    }

    #[test]
    fn test_shan_chen_psi() {
        let sc = ShanChenModel::new(-1.0, 1.0, 1.0, 1.0);
        let p = sc.psi_potential(0.0);
        assert!(p.abs() < 1e-14);
    }

    #[test]
    fn test_shan_chen_force_symmetric_zero() {
        // With all neighbor psi equal to psi_center, net force should be zero
        let sc = ShanChenModel::new(-1.0, 1.0, 1.0, 1.0);
        let psi_c = 0.5;
        let psi_nb = [psi_c; 8];
        let (fx, fy) = sc.shan_chen_force(psi_c, &psi_nb);
        assert!(fx.abs() < 1e-12);
        assert!(fy.abs() < 1e-12);
    }

    #[test]
    fn test_shan_chen_force_x_gradient() {
        // Higher psi on E side => net force in +x for negative G
        let g = -1.0;
        let psi_c = 0.5;
        let mut psi_nb = [0.5_f64; 8];
        psi_nb[0] = 0.8; // E neighbor has higher psi
        let (fx, _fy) = shan_chen_force_d2q9(g, psi_c, &psi_nb);
        // force = -g * psi_c * sum(w * psi * e_x)
        // sum contribution from E: (1/9)*0.8*1 + ... others at 0.5 cancel in pairs except E
        // net sum_x = (1/9)*(0.8 - 0.5) > 0, force_x = -(-1)*0.5*positive > 0
        assert!(fx > 0.0);
    }

    #[test]
    fn test_shan_chen_force_units() {
        // With all zeros, force is zero
        let psi_nb = [0.0_f64; 8];
        let (fx, fy) = shan_chen_force_d2q9(-1.0, 0.5, &psi_nb);
        assert_eq!(fx, 0.0);
        assert_eq!(fy, 0.0);
    }

    // --- FreeEnergyModel ---

    #[test]
    fn test_free_energy_model_new() {
        let fe = FreeEnergyModel::new(0.01, -0.1, 0.1, 0.7);
        assert_eq!(fe.kappa, 0.01);
        assert_eq!(fe.temperature, 0.7);
    }

    #[test]
    fn test_chemical_potential_landau_zero() {
        // mu(0) = 0 for Landau
        let mu = chemical_potential_landau(0.0, -0.1, 0.1);
        assert_eq!(mu, 0.0);
    }

    #[test]
    fn test_chemical_potential_landau_antisymmetric() {
        // mu(-rho) = -mu(rho) because mu = 2a*rho + 4b*rho^3 is odd
        let mu_pos = chemical_potential_landau(0.5, -1.0, 0.5);
        let mu_neg = chemical_potential_landau(-0.5, -1.0, 0.5);
        assert!((mu_pos + mu_neg).abs() < 1e-14);
    }

    #[test]
    fn test_chemical_potential_landau_double_well() {
        // For a < 0, b > 0: two minima of f at rho = +-sqrt(-a/(2b))
        // At those points mu = 0
        let a = -1.0;
        let b = 1.0;
        let rho_min = (-a / (2.0_f64 * b)).sqrt(); // = sqrt(0.5)
        let mu = chemical_potential_landau(rho_min, a, b);
        assert!(mu.abs() < 1e-12);
    }

    #[test]
    fn test_vdw_pressure_basic() {
        // At rho=0, pressure=0
        let p = vdw_pressure(0.0, 1.0, 1.0, 0.3);
        assert_eq!(p, 0.0);
    }

    #[test]
    fn test_vdw_pressure_positive_at_low_density() {
        let p = vdw_pressure(0.1, 1.0, 1.0, 0.3);
        // rho*T/(1-rho*b) - a*rho^2 = 0.1*1.0/0.97 - 1.0*0.01 ≈ 0.103 - 0.01 > 0
        assert!(p > 0.0);
    }

    #[test]
    fn test_interface_width_basic() {
        let xi = interface_width(0.01, 0.1);
        let expected = (0.01 / 0.1_f64).sqrt();
        assert!((xi - expected).abs() < 1e-14);
    }

    #[test]
    fn test_interface_width_zero_a() {
        let xi = interface_width(0.01, 0.0);
        assert_eq!(xi, f64::INFINITY);
    }

    #[test]
    fn test_interface_width_larger_kappa() {
        let xi1 = interface_width(0.01, 0.1);
        let xi2 = interface_width(0.04, 0.1);
        assert!(xi2 > xi1);
    }

    #[test]
    fn test_free_energy_model_interface_width() {
        let fe = FreeEnergyModel::new(0.01, -0.1, 0.1, 0.7);
        let xi = fe.interface_width();
        let expected = interface_width(0.01, 0.1);
        assert!((xi - expected).abs() < 1e-14);
    }

    #[test]
    fn test_phase_diagram_above_tc() {
        // Above critical temperature: no phase separation
        let (rv, rl) = phase_diagram(1.5, 1.0, 0.3, 100);
        assert!(rv.is_nan() || rl.is_nan());
    }

    #[test]
    fn test_phase_diagram_below_tc() {
        // T < T_c = 8a/(27b) = 8*1.0/(27*0.3) ≈ 0.988
        let t_c = 8.0 / (27.0 * 0.3);
        let t = 0.8 * t_c;
        let (rv, rl) = phase_diagram(t, 1.0, 0.3, 200);
        assert!(rv.is_finite() && rl.is_finite(), "rv = {rv}, rl = {rl}");
        assert!(rv > 0.0 && rv < rl);
        // High-precision (mpmath, 50 digits) Maxwell construction at T = 0.8 T_c.
        assert!(((rv - 0.26629657982343397) / rv).abs() < 1e-9);
        assert!(((rl - 2.147450920666331) / rl).abs() < 1e-9);
    }

    #[test]
    fn test_phase_diagram_rho_ordering() {
        // Below T_c the vapour is always strictly lighter than the liquid and
        // both lie on either side of the critical density 1/(3b).
        let t_c = 8.0 * 2.0 / (27.0 * 0.5);
        let t = 0.7 * t_c;
        let (rv, rl) = phase_diagram(t, 2.0, 0.5, 500);
        assert!(rv.is_finite() && rl.is_finite(), "rv = {rv}, rl = {rl}");
        let rho_c = 1.0 / (3.0 * 0.5);
        assert!(0.0 < rv && rv < rho_c && rho_c < rl);
        // High-precision (mpmath, 50 digits) Maxwell construction at T = 0.7 T_c.
        assert!(((rv - 0.08534820111052453) / rv).abs() < 1e-9);
        assert!(((rl - 1.4269616990038085) / rl).abs() < 1e-9);
    }

    #[test]
    fn test_chemical_potential_vdw_zero_density() {
        let mu = chemical_potential_vdw(0.0, 1.0, 1.0, 0.3);
        assert_eq!(mu, 0.0);
    }

    // --- Maxwell construction (phase coexistence) ---

    /// Composite Simpson integral of `f` over `[lo, hi]` with `panels` (even) panels.
    fn simpson<F: Fn(f64) -> f64>(f: F, lo: f64, hi: f64, panels: usize) -> f64 {
        let h = (hi - lo) / panels as f64;
        let mut sum = f(lo) + f(hi);
        for k in 1..panels {
            let weight = if k % 2 == 1 { 4.0 } else { 2.0 };
            sum += weight * f(lo + k as f64 * h);
        }
        sum * h / 3.0
    }

    #[test]
    fn test_phase_diagram_matches_reference_coexistence_values() {
        // High-precision (mpmath, 50 digits) solutions of p(rho_v) = p(rho_l),
        // mu(rho_v) = mu(rho_l) for a = 1, b = 0.3 (T_c = 8a/(27b)).
        let (a, b) = (1.0, 0.3);
        let t_c = 8.0 * a / (27.0 * b);
        let cases = [
            // (T/T_c, rho_v, rho_l, p*)
            (
                0.9,
                0.4730462641378402,
                1.8414113466648023,
                0.26625446579104983,
            ),
            (
                0.5,
                0.02416311905317121,
                2.7316577781668205,
                0.011435676972514519,
            ),
            (
                0.3,
                0.00044340585166953054,
                3.00462698376084,
                0.00013120038151480452,
            ),
            (
                0.999,
                1.0413011549174354,
                1.181810281369748,
                0.40987851813447684,
            ),
        ];
        for &(t_red, rho_v_ref, rho_l_ref, p_ref) in &cases {
            let t = t_red * t_c;
            let (rho_v, rho_l) = phase_diagram(t, a, b, 200);
            assert!(
                ((rho_v - rho_v_ref) / rho_v_ref).abs() < 1e-9,
                "T/Tc={t_red}: rho_v = {rho_v}, expected {rho_v_ref}"
            );
            assert!(
                ((rho_l - rho_l_ref) / rho_l_ref).abs() < 1e-9,
                "T/Tc={t_red}: rho_l = {rho_l}, expected {rho_l_ref}"
            );
            let coexistence =
                maxwell_construction(t, a, b, 200).expect("coexistence must exist below T_c");
            assert!(
                ((coexistence.pressure - p_ref) / p_ref).abs() < 1e-9,
                "T/Tc={t_red}: p* = {}, expected {p_ref}",
                coexistence.pressure
            );
        }
    }

    #[test]
    fn test_phase_diagram_equal_pressure_and_chemical_potential() {
        let (a, b) = (1.0, 0.3);
        let t_c = 8.0 * a / (27.0 * b);
        let rho_c = 1.0 / (3.0 * b);
        for &t_red in &[0.5, 0.7, 0.8, 0.9, 0.95] {
            let t = t_red * t_c;
            let (rho_v, rho_l) = phase_diagram(t, a, b, 200);
            assert!(
                rho_v < rho_c && rho_c < rho_l,
                "T/Tc={t_red}: rho_v={rho_v}, rho_c={rho_c}, rho_l={rho_l}"
            );
            let p_v = vdw_pressure(rho_v, t, a, b);
            let p_l = vdw_pressure(rho_l, t, a, b);
            assert!(
                ((p_v - p_l) / p_l).abs() < 1e-8,
                "T/Tc={t_red}: p_v = {p_v}, p_l = {p_l}"
            );
            let mu_v = chemical_potential_vdw(rho_v, t, a, b);
            let mu_l = chemical_potential_vdw(rho_l, t, a, b);
            assert!(
                (mu_v - mu_l).abs() < 1e-8,
                "T/Tc={t_red}: mu_v = {mu_v}, mu_l = {mu_l}"
            );
        }
    }

    #[test]
    fn test_phase_diagram_satisfies_equal_area_rule() {
        // Independent check of the classical Maxwell rule in the (p, v) plane:
        // integral_{rho_v}^{rho_l} (p(rho) - p*) / rho^2 d(rho) = 0  with v = 1/rho.
        let (a, b) = (1.0, 0.3);
        let t = 0.9 * 8.0 * a / (27.0 * b);
        let coexistence = maxwell_construction(t, a, b, 200).expect("T < T_c");
        let area = simpson(
            |rho| (vdw_pressure(rho, t, a, b) - coexistence.pressure) / (rho * rho),
            coexistence.rho_vapor,
            coexistence.rho_liquid,
            20_000,
        );
        let scale = simpson(
            |rho| (vdw_pressure(rho, t, a, b) - coexistence.pressure).abs() / (rho * rho),
            coexistence.rho_vapor,
            coexistence.rho_liquid,
            20_000,
        );
        assert!(scale > 1e-3, "the isotherm must oscillate around p*");
        assert!(
            area.abs() < 1e-9 * scale,
            "unbalanced area {area} (total area {scale})"
        );
        // The coexistence densities lie strictly outside the spinodal region
        // (where dp/drho < 0); the old heuristic returned the spinodal itself.
        let slope = |rho: f64| t / ((1.0 - b * rho) * (1.0 - b * rho)) - 2.0 * a * rho;
        assert!(slope(coexistence.rho_vapor) > 0.0);
        assert!(slope(coexistence.rho_liquid) > 0.0);
        assert!(coexistence.rho_vapor < 1.0 / (3.0 * b));
        assert!(coexistence.rho_liquid > 1.0 / (3.0 * b));
    }

    #[test]
    fn test_phase_diagram_scales_with_van_der_waals_parameters() {
        // rho * b and p * b^2 / a are universal functions of T / T_c.
        let (rho_v1, rho_l1) = {
            let (a, b) = (1.0, 0.3);
            phase_diagram(0.9 * 8.0 * a / (27.0 * b), a, b, 200)
        };
        let (a2, b2) = (2.0, 0.5);
        let t2 = 0.9 * 8.0 * a2 / (27.0 * b2);
        let (rho_v2, rho_l2) = phase_diagram(t2, a2, b2, 200);
        assert!((rho_v1 * 0.3 - rho_v2 * b2).abs() < 1e-12);
        assert!((rho_l1 * 0.3 - rho_l2 * b2).abs() < 1e-12);
        assert!((rho_v2 - 0.2838277584827041).abs() < 1e-9);
        assert!((rho_l2 - 1.1048468079988814).abs() < 1e-9);
    }

    #[test]
    fn test_phase_diagram_follows_mean_field_critical_scaling() {
        // Close to T_c the van der Waals binodal obeys
        // (rho_l - rho_v) / rho_c = 4 sqrt(1 - T / T_c).
        let (a, b) = (1.0, 0.3);
        let t_c = 8.0 * a / (27.0 * b);
        let rho_c = 1.0 / (3.0 * b);
        for &t_red in &[0.99, 0.999] {
            let (rho_v, rho_l) = phase_diagram(t_red * t_c, a, b, 200);
            let width = (rho_l - rho_v) / rho_c;
            let expected = 4.0 * (1.0 - t_red).sqrt();
            assert!(
                (width / expected - 1.0).abs() < 0.03,
                "T/Tc={t_red}: width {width}, mean-field {expected}"
            );
        }
    }

    #[test]
    fn test_phase_diagram_widens_on_cooling() {
        let (a, b) = (1.0, 0.3);
        let t_c = 8.0 * a / (27.0 * b);
        let mut previous = phase_diagram(0.95 * t_c, a, b, 200);
        for &t_red in &[0.9, 0.8, 0.7, 0.6, 0.5] {
            let current = phase_diagram(t_red * t_c, a, b, 200);
            assert!(
                current.0 < previous.0,
                "vapour density must fall on cooling"
            );
            assert!(
                current.1 > previous.1,
                "liquid density must rise on cooling"
            );
            previous = current;
        }
    }

    #[test]
    fn test_phase_diagram_rejects_invalid_input() {
        let (a, b) = (1.0, 0.3);
        let t_c = 8.0 * a / (27.0 * b);
        let invalid = [
            phase_diagram(t_c, a, b, 200),
            phase_diagram(1.1 * t_c, a, b, 200),
            phase_diagram(0.0, a, b, 200),
            phase_diagram(-0.1, a, b, 200),
            phase_diagram(f64::NAN, a, b, 200),
            phase_diagram(0.8 * t_c, 0.0, b, 200),
            phase_diagram(0.8 * t_c, a, 0.0, 200),
            phase_diagram(0.8 * t_c, a, -b, 200),
            phase_diagram(0.8 * t_c, a, b, 3),
        ];
        for (rv, rl) in invalid {
            assert!(rv.is_nan() && rl.is_nan());
        }
        assert!(maxwell_construction(0.8 * t_c, a, b, 3).is_none());
    }

    #[test]
    fn test_phase_diagram_resolution_parameter_converges() {
        // 60 bisection halvings already reach full double precision.
        let (a, b) = (1.0, 0.3);
        let t = 0.85 * 8.0 * a / (27.0 * b);
        let (rv_fine, rl_fine) = phase_diagram(t, a, b, 200);
        let (rv_60, rl_60) = phase_diagram(t, a, b, 60);
        assert!(((rv_60 - rv_fine) / rv_fine).abs() < 1e-12);
        assert!(((rl_60 - rl_fine) / rl_fine).abs() < 1e-12);
    }

    #[test]
    fn test_maxwell_construction_reports_common_chemical_potential() {
        let (a, b) = (1.0, 0.3);
        let t = 0.9 * 8.0 * a / (27.0 * b);
        let c = maxwell_construction(t, a, b, 200).expect("T < T_c");
        assert!((c.chemical_potential - (-1.3284277133072706)).abs() < 1e-9);
        assert!((chemical_potential_vdw(c.rho_vapor, t, a, b) - c.chemical_potential).abs() < 1e-9);
        assert!(
            (chemical_potential_vdw(c.rho_liquid, t, a, b) - c.chemical_potential).abs() < 1e-9
        );
    }
}
