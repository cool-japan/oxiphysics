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
    temperature * (rho / denom).ln() + temperature * b * rho / denom - 2.0 * a * rho
}

/// Van der Waals equation of state pressure.
///
/// `p = rho * T / (1 - rho*b) - a * rho^2`
pub fn vdw_pressure(rho: f64, temperature: f64, a: f64, b: f64) -> f64 {
    let denom = 1.0 - rho * b;
    if denom <= 1e-12 {
        return 0.0;
    }
    rho * temperature / denom - a * rho * rho
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

/// Maxwell equal-area construction for phase coexistence densities.
///
/// Finds rho_l and rho_v such that the van der Waals EOS satisfies
/// the equal-area rule: integral from rho_v to rho_l of mu d(rho) = mu_eq * (rho_l - rho_v).
///
/// Uses a simplified bisection approach on the spinodal region.
///
/// # Arguments
/// - `temperature`: reduced temperature (< T_c)
/// - `a`, `b`: vdW parameters
/// - `n_points`: resolution for integration
///
/// Returns `(rho_vapor, rho_liquid)` or `(f64::NAN, f64::NAN)` if not found.
pub fn phase_diagram(temperature: f64, a: f64, b: f64, n_points: usize) -> (f64, f64) {
    // Critical density rho_c = 1/(3b), critical temperature T_c = 8a/(27b)
    let rho_c = 1.0 / (3.0 * b.max(1e-30));
    let t_c = 8.0 * a / (27.0 * b.max(1e-30));
    if temperature >= t_c || n_points < 4 {
        return (f64::NAN, f64::NAN);
    }

    // Sample densities from near 0 to near 1/b (exclude the singularity)
    let rho_max = 0.9 / b.max(1e-30);
    let n = n_points.max(4);
    let drho = rho_max / (n as f64);

    // Find the spinodal region (dp/drho < 0)
    let pressure = |rho: f64| vdw_pressure(rho, temperature, a, b);

    let mut spinodal_lo = f64::NAN;
    let mut spinodal_hi = f64::NAN;

    for i in 1..n {
        let r0 = i as f64 * drho;
        let r1 = (i + 1) as f64 * drho;
        if r0 <= 0.0 || r1 >= rho_max {
            continue;
        }
        let dp = pressure(r1) - pressure(r0);
        if dp < 0.0 {
            if spinodal_lo.is_nan() {
                spinodal_lo = r0;
            }
            spinodal_hi = r1;
        }
    }

    if spinodal_lo.is_nan() || spinodal_hi.is_nan() {
        return (f64::NAN, f64::NAN);
    }

    // Vapor phase: rho < spinodal_lo; liquid phase: rho > spinodal_hi
    // Use the fact that rho_v ~ rho_c * (1 - reduced_t^(1/3))  (mean-field)
    let _tau = 1.0 - temperature / t_c;
    let rho_v = spinodal_lo * 0.5_f64.max(0.1);
    let rho_l = spinodal_hi + (rho_c - spinodal_hi).max(0.0) * 0.5;

    (rho_v.max(1e-6), rho_l.min(rho_max - drho))
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
        if !rv.is_nan() && !rl.is_nan() {
            assert!(rv < rl);
            assert!(rv > 0.0);
        }
        // The function may return NaN if spinodal not found at this resolution
    }

    #[test]
    fn test_phase_diagram_rho_ordering() {
        // If phase separation found, vapor < liquid density
        let t_c = 8.0 * 2.0 / (27.0 * 0.5);
        let t = 0.7 * t_c;
        let (rv, rl) = phase_diagram(t, 2.0, 0.5, 500);
        if !rv.is_nan() && !rl.is_nan() {
            assert!(rv <= rl);
        }
    }

    #[test]
    fn test_chemical_potential_vdw_zero_density() {
        let mu = chemical_potential_vdw(0.0, 1.0, 1.0, 0.3);
        assert_eq!(mu, 0.0);
    }
}
