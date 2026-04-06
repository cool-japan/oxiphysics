//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::R_UNIVERSAL;
#[allow(unused_imports)]
use super::functions::*;

/// Rankine-Hugoniot shock jump conditions.
///
/// Relates shock velocity U_s to particle velocity U_p via the linear fit:
/// U_s = c₀ + s * U_p
#[derive(Debug, Clone, Copy)]
pub struct ShockHugoniot {
    /// Reference density ρ₀ (kg/m³).
    pub rho0: f64,
    /// Reference bulk sound speed c₀ (m/s).
    pub c0: f64,
    /// Linear Hugoniot slope s.
    pub s: f64,
    /// Reference pressure p₀ (Pa).
    pub p0: f64,
}
impl ShockHugoniot {
    /// Create a Hugoniot.
    pub fn new(rho0: f64, c0: f64, s: f64) -> Self {
        Self {
            rho0,
            c0,
            s,
            p0: 0.0,
        }
    }
    /// Aluminium: ρ₀=2700, c₀=5240 m/s, s=1.4.
    pub fn aluminum() -> Self {
        Self::new(2700.0, 5240.0, 1.4)
    }
    /// Iron: ρ₀=7874, c₀=3574 m/s, s=1.92.
    pub fn iron() -> Self {
        Self::new(7874.0, 3574.0, 1.92)
    }
    /// Shock velocity from particle velocity.
    pub fn shock_velocity(&self, u_p: f64) -> f64 {
        self.c0 + self.s * u_p
    }
    /// Hugoniot pressure from particle velocity.
    ///
    /// p_H = ρ₀ * U_s * U_p
    pub fn hugoniot_pressure(&self, u_p: f64) -> f64 {
        let u_s = self.shock_velocity(u_p);
        self.rho0 * u_s * u_p + self.p0
    }
    /// Post-shock density from particle velocity.
    ///
    /// ρ = ρ₀ * U_s / (U_s - U_p)
    pub fn post_shock_density(&self, u_p: f64) -> f64 {
        let u_s = self.shock_velocity(u_p);
        if (u_s - u_p).abs() < 1e-10 {
            return f64::INFINITY;
        }
        self.rho0 * u_s / (u_s - u_p)
    }
    /// Hugoniot specific internal energy from particle velocity.
    ///
    /// e_H = 0.5 * (p₀ + p_H) * (V₀ - V)
    pub fn hugoniot_energy(&self, u_p: f64) -> f64 {
        let p_h = self.hugoniot_pressure(u_p);
        let rho_h = self.post_shock_density(u_p);
        0.5 * (self.p0 + p_h) * (1.0 / self.rho0 - 1.0 / rho_h)
    }
}
/// Vinet (Universal) equation of state for solids.
///
/// p(V) = 3 K₀ * x^(-2) * (1 - x) * exp(η(1-x))
/// where x = (V/V₀)^(1/3), η = (3/2)(K₀' - 1).
#[derive(Debug, Clone, Copy)]
pub struct VinetEos {
    /// Reference specific volume V₀ (m³/kg).
    pub v0: f64,
    /// Bulk modulus at zero pressure K₀ (Pa).
    pub k0: f64,
    /// Pressure derivative K₀'.
    pub k0_prime: f64,
}
impl VinetEos {
    /// Create a Vinet EOS.
    pub fn new(v0: f64, k0: f64, k0_prime: f64) -> Self {
        Self { v0, k0, k0_prime }
    }
    /// Diamond (cubic): V₀ = 1/3515, K₀ = 443 GPa, K₀' = 3.6.
    pub fn diamond() -> Self {
        Self::new(1.0 / 3515.0, 443.0e9, 3.6)
    }
    /// Copper: V₀ = 1/8960, K₀ = 136.7 GPa, K₀' = 5.3.
    pub fn copper() -> Self {
        Self::new(1.0 / 8960.0, 136.7e9, 5.3)
    }
    /// Pressure from specific volume.
    pub fn pressure_from_volume(&self, v: f64) -> f64 {
        let x = (v / self.v0).powf(1.0 / 3.0);
        if (x - 1.0).abs() < 1e-15 {
            return 0.0;
        }
        let eta = 1.5 * (self.k0_prime - 1.0);
        3.0 * self.k0 * (1.0 - x) / (x * x) * (eta * (1.0 - x)).exp()
    }
    /// Volume from pressure.
    pub fn volume(&self, pressure: f64) -> f64 {
        let mut lo = self.v0 * 0.3;
        let mut hi = self.v0 * 2.0;
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if self.pressure_from_volume(mid) > pressure {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }
}
/// Jones-Wilkins-Lee equation of state for detonation products.
///
/// p = A*(1 - ω/(R1*v))*exp(-R1*v) + B*(1 - ω/(R2*v))*exp(-R2*v) + ω*e/v
///
/// where v = ρ₀/ρ (relative specific volume = 1/η).
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct JwlEos {
    /// JWL coefficient A (Pa).
    pub big_a: f64,
    /// JWL coefficient B (Pa).
    pub big_b: f64,
    /// JWL constant R1.
    pub r1: f64,
    /// JWL constant R2.
    pub r2: f64,
    /// Grüneisen-like constant ω.
    pub omega: f64,
    /// Reference detonation energy E0 (J/m³).
    pub e0: f64,
    /// Reference density ρ₀ (kg/m³).
    pub rho0: f64,
}
#[allow(dead_code)]
impl JwlEos {
    /// Create a JWL EOS.
    #[allow(clippy::too_many_arguments)]
    pub fn new(big_a: f64, big_b: f64, r1: f64, r2: f64, omega: f64, e0: f64, rho0: f64) -> Self {
        Self {
            big_a,
            big_b,
            r1,
            r2,
            omega,
            e0,
            rho0,
        }
    }
    /// TNT parameters (approximate literature values).
    pub fn tnt() -> Self {
        Self::new(3.712e11, 3.231e9, 4.15, 0.95, 0.30, 7.0e9, 1630.0)
    }
    /// Pressure from density `rho` (kg/m³) and specific internal energy `e` (J/m³).
    ///
    /// v = ρ₀/ρ (relative volume)
    /// p = A*(1 - ω/(R1*v))*exp(-R1*v) + B*(1 - ω/(R2*v))*exp(-R2*v) + ω*e/v
    pub fn pressure(&self, rho: f64, e: f64) -> f64 {
        let v = self.rho0 / rho.max(1e-30);
        let term1 = self.big_a * (1.0 - self.omega / (self.r1 * v)) * (-self.r1 * v).exp();
        let term2 = self.big_b * (1.0 - self.omega / (self.r2 * v)) * (-self.r2 * v).exp();
        let term3 = self.omega * e / v;
        term1 + term2 + term3
    }
}
/// Polynomial EOS commonly used in hydrocode simulations.
///
/// p = a₁μ + a₂μ² + a₃μ³ + (b₀ + b₁μ + b₂μ²)ρe  (compressed)
/// p = t₁μ + t₂μ²                                  (expanded, e<0)
/// where μ = ρ/ρ₀ - 1.
#[derive(Debug, Clone, Copy)]
pub struct PolynomialEos {
    /// Reference density ρ₀ (kg/m³).
    pub rho0: f64,
    /// Compression coefficients a₁, a₂, a₃ (Pa).
    pub a: [f64; 3],
    /// Energy coupling b₀, b₁, b₂.
    pub b: [f64; 3],
    /// Tension coefficients t₁, t₂ (Pa).
    pub t: [f64; 2],
}
impl PolynomialEos {
    /// Create a polynomial EOS.
    pub fn new(rho0: f64, a: [f64; 3], b: [f64; 3], t: [f64; 2]) -> Self {
        Self { rho0, a, b, t }
    }
    /// Compression μ = ρ/ρ₀ - 1.
    pub fn mu(&self, density: f64) -> f64 {
        density / self.rho0 - 1.0
    }
    /// Pressure from density and energy.
    pub fn pressure_energy(&self, density: f64, energy: f64) -> f64 {
        let mu = self.mu(density);
        if mu >= 0.0 {
            let poly = self.a[0] * mu + self.a[1] * mu * mu + self.a[2] * mu * mu * mu;
            let energy_term = (self.b[0] + self.b[1] * mu + self.b[2] * mu * mu) * density * energy;
            poly + energy_term
        } else {
            self.t[0] * mu + self.t[1] * mu * mu
        }
    }
}
/// Volume-dependent Grüneisen parameter model.
///
/// Γ(V) = Γ₀ * (V/V₀)^q  (power-law, common approximation).
#[derive(Debug, Clone, Copy)]
pub struct GruneisenParameter {
    /// Grüneisen parameter at reference volume.
    pub gamma_0: f64,
    /// Reference specific volume V₀.
    pub v0: f64,
    /// Exponent q (often taken as 1 or 2/3).
    pub q: f64,
}
impl GruneisenParameter {
    /// Create a Grüneisen parameter model.
    pub fn new(gamma_0: f64, v0: f64, q: f64) -> Self {
        Self { gamma_0, v0, q }
    }
    /// Evaluate Γ(V).
    pub fn evaluate(&self, v: f64) -> f64 {
        self.gamma_0 * (v / self.v0).powf(self.q)
    }
    /// Thermal pressure contribution: p_th = Γ(V) * ρ * cv * (T - T₀).
    pub fn thermal_pressure(&self, density: f64, cv: f64, temperature: f64, t0: f64) -> f64 {
        let v = 1.0 / density;
        self.evaluate(v) * density * cv * (temperature - t0)
    }
}
/// Peng-Robinson equation of state.
///
/// P = R*T/(V_m - b) - a(T) / \[V_m*(V_m + b) + b*(V_m - b)\]
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct PengRobinson {
    /// Critical temperature Tc (K).
    pub tc: f64,
    /// Critical pressure Pc (Pa).
    pub pc: f64,
    /// Acentric factor ω.
    pub omega: f64,
}
#[allow(dead_code)]
impl PengRobinson {
    /// Create a Peng-Robinson EOS.
    pub fn new(tc: f64, pc: f64, omega: f64) -> Self {
        Self { tc, pc, omega }
    }
    /// PR `b` parameter: b = 0.07780 * R * Tc / Pc
    pub fn b_param(&self) -> f64 {
        0.07780 * R_UNIVERSAL * self.tc / self.pc
    }
    /// PR `a` parameter at critical: a_c = 0.45724 * R² * Tc² / Pc
    pub fn a_critical(&self) -> f64 {
        0.45724 * R_UNIVERSAL * R_UNIVERSAL * self.tc * self.tc / self.pc
    }
    /// α(T) correction function for Peng-Robinson.
    ///
    /// α = \[1 + κ*(1 - sqrt(T/Tc))\]²
    /// κ = 0.37464 + 1.54226*ω - 0.26992*ω²
    pub fn acentric_factor_correction(&self, t: f64) -> f64 {
        let kappa = 0.37464 + 1.54226 * self.omega - 0.26992 * self.omega * self.omega;

        (1.0 + kappa * (1.0 - (t / self.tc).sqrt())).powi(2)
    }
    /// Temperature-dependent `a(T) = a_c * α(T)`.
    pub fn a_param(&self, t: f64) -> f64 {
        self.a_critical() * self.acentric_factor_correction(t)
    }
    /// Pressure at molar volume `v_mol` and temperature `t`.
    pub fn pressure(&self, t: f64, v_mol: f64) -> f64 {
        let b = self.b_param();
        let a_t = self.a_param(t);
        if (v_mol - b).abs() < 1e-20 {
            return f64::INFINITY;
        }
        R_UNIVERSAL * t / (v_mol - b) - a_t / (v_mol * (v_mol + b) + b * (v_mol - b))
    }
}
/// Murnaghan equation of state (first-order pressure derivative of bulk modulus).
///
/// P(V) = (K₀/n) * \[ (V₀/V)^n − 1 \]
///
/// This is the original Murnaghan 1944 form, simpler than Birch-Murnaghan.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct MurnaghanEos {
    /// Reference specific volume V₀ (m³/kg).
    pub v0: f64,
    /// Bulk modulus at zero pressure K₀ (Pa).
    pub k0: f64,
    /// Pressure derivative n = K₀' (dimensionless, typically 3–5).
    pub n: f64,
}
#[allow(dead_code)]
impl MurnaghanEos {
    /// Create a Murnaghan EOS.
    pub fn new(v0: f64, k0: f64, n: f64) -> Self {
        Self { v0, k0, n }
    }
    /// Aluminium: V₀ = 1/2700, K₀ = 76 GPa, n = 4.5.
    pub fn aluminum() -> Self {
        Self::new(1.0 / 2700.0, 76.0e9, 4.5)
    }
    /// Iron: V₀ = 1/7874, K₀ = 170 GPa, n = 4.9.
    pub fn iron() -> Self {
        Self::new(1.0 / 7874.0, 170.0e9, 4.9)
    }
    /// Pressure from specific volume V (m³/kg).
    ///
    /// P = (K₀/n) * \[(V₀/V)^n − 1\]
    pub fn pressure_from_volume(&self, v: f64) -> f64 {
        if self.n.abs() < f64::EPSILON {
            return 0.0;
        }
        (self.k0 / self.n) * ((self.v0 / v).powf(self.n) - 1.0)
    }
    /// Volume from pressure (analytical inversion).
    ///
    /// V = V₀ / (1 + n*P/K₀)^(1/n)
    pub fn volume_from_pressure(&self, pressure: f64) -> f64 {
        let x = 1.0 + self.n * pressure / self.k0;
        if x <= 0.0 {
            return self.v0 * 0.01;
        }
        self.v0 / x.powf(1.0 / self.n)
    }
    /// Bulk modulus at volume V: K(V) = K₀ * (V₀/V)^n.
    pub fn bulk_modulus(&self, v: f64) -> f64 {
        self.k0 * (self.v0 / v).powf(self.n)
    }
}
/// Tait equation of state in bulk-modulus / isothermal form.
///
/// Often used for water and other liquids:
/// P = K0/K0' * \[ (ρ/ρ0)^K0' - 1 \]
///
/// Different from `TaitEos` which uses the c₀-based SPH form.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct TaitBulkEos {
    /// Reference bulk modulus K₀ (Pa).
    pub k0: f64,
    /// Pressure derivative of bulk modulus K₀' (dimensionless).
    pub k0_prime: f64,
    /// Reference density ρ₀ (kg/m³).
    pub rho0: f64,
}
#[allow(dead_code)]
impl TaitBulkEos {
    /// Create a Tait bulk-modulus EOS.
    pub fn new(k0: f64, k0_prime: f64, rho0: f64) -> Self {
        Self { k0, k0_prime, rho0 }
    }
    /// Water at 300 K: K₀ = 2.2 GPa, K₀' = 6.0, ρ₀ = 997 kg/m³.
    pub fn water() -> Self {
        Self::new(2.2e9, 6.0, 997.0)
    }
    /// Pressure from density ρ (kg/m³).
    ///
    /// P = K0/K0' * \[ (ρ/ρ0)^K0' - 1 \]
    ///
    /// Returns 0 at reference density.
    pub fn pressure(&self, rho: f64) -> f64 {
        let ratio = rho / self.rho0;
        self.k0 / self.k0_prime * (ratio.powf(self.k0_prime) - 1.0)
    }
    /// Sound speed: c = sqrt(K₀ * (ρ/ρ₀)^(K₀'-1) / ρ₀) \[approx\].
    #[allow(dead_code)]
    pub fn sound_speed(&self, rho: f64) -> f64 {
        let ratio = rho / self.rho0;
        (self.k0 * ratio.powf(self.k0_prime - 1.0) / self.rho0).sqrt()
    }
}
/// Ideal gas equation of state: p = ρ R_specific T.
#[derive(Debug, Clone, Copy)]
pub struct IdealGasEos {
    /// Adiabatic index γ = Cp/Cv.
    pub gamma: f64,
    /// Specific gas constant R_specific (J/(kg·K)).
    pub specific_gas_constant: f64,
}
impl IdealGasEos {
    /// Create a new ideal gas EOS.
    pub fn new(gamma: f64, specific_gas_constant: f64) -> Self {
        Self {
            gamma,
            specific_gas_constant,
        }
    }
    /// Dry air (γ=1.4, R=287 J/(kg·K)).
    pub fn air() -> Self {
        Self::new(1.4, 287.0)
    }
    /// Diatomic hydrogen (γ≈1.4, R=4157 J/(kg·K)).
    pub fn hydrogen() -> Self {
        Self::new(1.4, 4157.0)
    }
    /// Compute pressure at a given density and temperature.
    pub fn pressure_at_temperature(&self, density: f64, temperature: f64) -> f64 {
        density * self.specific_gas_constant * temperature
    }
    /// Temperature from density and pressure.
    pub fn temperature(&self, density: f64, pressure: f64) -> f64 {
        if density.abs() < f64::EPSILON {
            return 0.0;
        }
        pressure / (density * self.specific_gas_constant)
    }
}
/// 3rd-order Birch-Murnaghan equation of state for solids.
///
/// p(V) = (3/2) * K₀ * \[ (V₀/V)^(7/3) - (V₀/V)^(5/3) \]
///          * { 1 + (3/4)*(K₀' - 4)*\[ (V₀/V)^(2/3) - 1 \] }
#[derive(Debug, Clone, Copy)]
pub struct BirchMurnaghan3Eos {
    /// Reference volume V₀ (m³/kg, i.e. specific volume).
    pub v0: f64,
    /// Isothermal bulk modulus K₀ at zero pressure (Pa).
    pub k0: f64,
    /// Pressure derivative of K₀: dK/dP at zero pressure.
    pub k0_prime: f64,
}
impl BirchMurnaghan3Eos {
    /// Create a 3rd-order BM EOS.
    pub fn new(v0: f64, k0: f64, k0_prime: f64) -> Self {
        Self { v0, k0, k0_prime }
    }
    /// Iron at ambient conditions (approximate).
    pub fn iron() -> Self {
        Self::new(1.0 / 7874.0, 170.0e9, 4.9)
    }
    /// MgO (periclase): V₀ = 1/3580, K₀ = 160.2 GPa, K₀' = 3.99.
    pub fn mgo() -> Self {
        Self::new(1.0 / 3580.0, 160.2e9, 3.99)
    }
    /// Pressure from specific volume V (m³/kg).
    pub fn pressure_from_volume(&self, v: f64) -> f64 {
        let x = self.v0 / v;
        let x23 = x.powf(2.0 / 3.0);
        let x73 = x.powf(7.0 / 3.0);
        let x53 = x.powf(5.0 / 3.0);
        let correction = 1.0 + 0.75 * (self.k0_prime - 4.0) * (x23 - 1.0);
        1.5 * self.k0 * (x73 - x53) * correction
    }
    /// Volume from pressure (numerical bisection).
    pub fn volume(&self, pressure: f64) -> f64 {
        let mut lo = self.v0 * 0.1;
        let mut hi = self.v0 * 10.0;
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if self.pressure_from_volume(mid) > pressure {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }
    /// Isothermal bulk modulus K_T(V) at volume V.
    pub fn bulk_modulus(&self, v: f64) -> f64 {
        let h = v * 1e-6;
        let dp = self.pressure_from_volume(v - h) - self.pressure_from_volume(v + h);
        v * dp / (2.0 * h)
    }
}
/// Van der Waals equation of state in molar form.
///
/// (P + a/V_m²)(V_m - b) = R * T
/// P = R*T/(V_m - b) - a/V_m²
///
/// Different from `VanDerWaalsEos` (which works with mass density).
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct VanDerWaals {
    /// Attractive interaction constant a (Pa·m⁶/mol²).
    pub a: f64,
    /// Excluded volume constant b (m³/mol).
    pub b: f64,
    /// Universal gas constant R (J/(mol·K)).
    pub r_gas: f64,
}
#[allow(dead_code)]
impl VanDerWaals {
    /// Create a new VdW EOS.  `R` is the universal gas constant (8.314 J/mol/K).
    pub fn new(a: f64, b: f64, r_gas: f64) -> Self {
        Self { a, b, r_gas }
    }
    /// Standard instance with R = 8.314 J/(mol·K).
    pub fn with_r(a: f64, b: f64) -> Self {
        Self::new(a, b, 8.314_462_618)
    }
    /// Pressure: P = R*T/(V_m - b) - a/V_m²
    ///
    /// # Arguments
    /// * `t`   - Temperature (K).
    /// * `v`   - Molar volume (m³/mol).
    /// * `_n`  - Amount of substance (mol) — unused in molar form, kept for API.
    pub fn pressure(&self, t: f64, v: f64, _n: f64) -> f64 {
        if (v - self.b).abs() < 1e-20 {
            return f64::INFINITY;
        }
        self.r_gas * t / (v - self.b) - self.a / (v * v)
    }
    /// Compressibility factor Z = P V_m / (R T).
    ///
    /// # Arguments
    /// * `t` - Temperature (K).
    /// * `p` - Pressure (Pa).
    /// * `n` - Amount (mol), used as hint to select molar volume root via bisection.
    pub fn compressibility(&self, t: f64, p: f64, _n: f64) -> f64 {
        let rt = self.r_gas * t;
        if rt.abs() < 1e-30 {
            return 0.0;
        }
        let mut lo = self.b + 1e-10;
        let mut hi = rt / p.max(1.0) * 100.0 + 1.0;
        hi = hi.max(lo * 1000.0);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if self.pressure(t, mid, 1.0) > p {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let v_m = 0.5 * (lo + hi);
        p * v_m / rt
    }
}
/// Mie-Grüneisen equation of state for shock-compressed solids.
///
/// p = p_H + Γ*ρ*(e - e_H)
/// where p_H and e_H are the Hugoniot pressure and energy.
#[derive(Debug, Clone, Copy)]
pub struct MieGruneisenEos {
    /// Reference density ρ₀ (kg/m³).
    pub rho0: f64,
    /// Grüneisen parameter Γ (dimensionless).
    pub gamma_0: f64,
    /// Linear shock velocity coefficient c₀ (m/s).
    pub c0: f64,
    /// Slope of shock velocity - particle velocity (S parameter).
    pub s: f64,
}
impl MieGruneisenEos {
    /// Create a Mie-Grüneisen EOS.
    pub fn new(rho0: f64, gamma_0: f64, c0: f64, s: f64) -> Self {
        Self {
            rho0,
            gamma_0,
            c0,
            s,
        }
    }
    /// Copper (ρ₀=8930, Γ=2.0, c₀=3933, S=1.5).
    pub fn copper() -> Self {
        Self::new(8930.0, 2.0, 3933.0, 1.5)
    }
    /// Aluminium (ρ₀=2785, Γ=2.0, c₀=5386, S=1.339).
    pub fn aluminum() -> Self {
        Self::new(2785.0, 2.0, 5386.0, 1.339)
    }
    /// Iron (ρ₀=7896, Γ=1.81, c₀=3574, S=1.92).
    pub fn iron() -> Self {
        Self::new(7896.0, 1.81, 3574.0, 1.92)
    }
    /// Volumetric compression η = 1 - ρ₀/ρ.
    pub fn compression(&self, density: f64) -> f64 {
        1.0 - self.rho0 / density
    }
    /// Hugoniot pressure (reference Rankine-Hugoniot curve).
    pub fn hugoniot_pressure(&self, density: f64) -> f64 {
        let mu = density / self.rho0 - 1.0;
        if mu <= 0.0 {
            return 0.0;
        }
        let denom = (1.0 - self.s * mu).powi(2);
        if denom.abs() < f64::EPSILON {
            return 0.0;
        }
        self.rho0 * self.c0 * self.c0 * mu / denom
    }
    /// Pressure from density and specific internal energy (J/kg).
    pub fn pressure_from_energy(&self, density: f64, energy: f64) -> f64 {
        let ph = self.hugoniot_pressure(density);
        let eh = if density > self.rho0 {
            ph * 0.5 * (1.0 / self.rho0 - 1.0 / density)
        } else {
            0.0
        };
        ph + self.gamma_0 * density * (energy - eh)
    }
}
/// Tait equation of state, commonly used for weakly-compressible fluids (e.g. water in SPH).
///
/// pressure = B * ((ρ/ρ₀)^γ − 1) + p_ref
/// where B = ρ₀ c₀² / γ
#[derive(Debug, Clone, Copy)]
pub struct TaitEos {
    /// Reference density ρ₀ (kg/m³).
    pub reference_density: f64,
    /// Reference pressure p_ref (Pa).
    pub reference_pressure: f64,
    /// Reference speed of sound c₀ (m/s).
    pub speed_of_sound: f64,
    /// Polytropic exponent γ (typically 7 for water).
    pub gamma: f64,
}
impl TaitEos {
    /// Create a new Tait equation of state.
    pub fn new(
        reference_density: f64,
        reference_pressure: f64,
        speed_of_sound: f64,
        gamma: f64,
    ) -> Self {
        Self {
            reference_density,
            reference_pressure,
            speed_of_sound,
            gamma,
        }
    }
    /// Water at standard conditions: ρ₀=1000, p_ref=0, c₀=1500, γ=7.
    pub fn water() -> Self {
        Self::new(1000.0, 0.0, 1500.0, 7.0)
    }
    /// B coefficient = ρ₀ c₀² / γ.
    pub fn b_coefficient(&self) -> f64 {
        self.reference_density * self.speed_of_sound.powi(2) / self.gamma
    }
}
/// Redlich-Kwong equation of state.
///
/// P = R*T/(V_m - b) - a / (T^0.5 * V_m * (V_m + b))
///
/// Uses R = 8.314 J/(mol·K).
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct RedlichKwong {
    /// RK attraction parameter a (Pa·m⁶·K^0.5/mol²).
    pub a: f64,
    /// RK co-volume parameter b (m³/mol).
    pub b: f64,
}
#[allow(dead_code)]
impl RedlichKwong {
    /// Create a Redlich-Kwong EOS.
    pub fn new(a: f64, b: f64) -> Self {
        Self { a, b }
    }
    /// Compute a and b from critical properties.
    ///
    /// a = 0.42748 * R² * Tc^2.5 / Pc
    /// b = 0.08664 * R * Tc / Pc
    pub fn from_critical(tc: f64, pc: f64) -> Self {
        let a = 0.42748 * R_UNIVERSAL * R_UNIVERSAL * tc.powf(2.5) / pc;
        let b = 0.08664 * R_UNIVERSAL * tc / pc;
        Self { a, b }
    }
    /// Pressure at molar volume `v_mol` and temperature `t`.
    pub fn pressure(&self, t: f64, v_mol: f64) -> f64 {
        if (v_mol - self.b).abs() < 1e-20 {
            return f64::INFINITY;
        }
        R_UNIVERSAL * t / (v_mol - self.b) - self.a / (t.sqrt() * v_mol * (v_mol + self.b))
    }
}
/// Van der Waals equation of state for real gases.
///
/// (p + a*ρ²/M²)(1/ρ - b/M) = R*T/M
/// where a and b are van der Waals constants.
#[derive(Debug, Clone, Copy)]
pub struct VanDerWaalsEos {
    /// Attractive interaction constant a (Pa·m⁶/mol²).
    pub a: f64,
    /// Excluded volume constant b (m³/mol).
    pub b: f64,
    /// Molar mass M (kg/mol).
    pub molar_mass: f64,
    /// Universal gas constant R (8.314 J/(mol·K)).
    pub r_gas: f64,
}
impl VanDerWaalsEos {
    /// Create a van der Waals EOS.
    pub fn new(a: f64, b: f64, molar_mass: f64) -> Self {
        Self {
            a,
            b,
            molar_mass,
            r_gas: 8.314,
        }
    }
    /// Carbon dioxide (a=0.3658, b=42.9e-6, M=0.044).
    pub fn co2() -> Self {
        Self::new(0.3658, 42.9e-6, 0.044)
    }
    /// Water vapour (a=0.5537, b=30.5e-6, M=0.018).
    pub fn water_vapor() -> Self {
        Self::new(0.5537, 30.5e-6, 0.018)
    }
    /// Nitrogen (a=0.1370, b=38.7e-6, M=0.028).
    pub fn nitrogen() -> Self {
        Self::new(0.1370, 38.7e-6, 0.028)
    }
    /// Pressure at given density and temperature.
    pub fn pressure_at_temperature(&self, density: f64, temperature: f64) -> f64 {
        let specific_volume = 1.0 / density;
        let molar_volume = specific_volume * self.molar_mass;
        if (molar_volume - self.b).abs() < 1e-15 {
            return f64::INFINITY;
        }
        self.r_gas * temperature / (molar_volume - self.b) - self.a / (molar_volume * molar_volume)
    }
}
/// Tillotson equation of state for hypervelocity impact and planetary science.
///
/// Used for rocks, minerals, and ices under extreme shock conditions.
#[derive(Debug, Clone, Copy)]
pub struct TillotsonEos {
    /// Reference density ρ₀ (kg/m³).
    pub rho0: f64,
    /// Reference internal energy E₀ (J/kg).
    pub e0: f64,
    /// Cold pressure coefficient a.
    pub a: f64,
    /// Cold pressure coefficient b.
    pub b: f64,
    /// Bulk modulus-like coefficient A (Pa).
    pub big_a: f64,
    /// Bulk modulus-like coefficient B (Pa).
    pub big_b: f64,
    /// Incipient vaporization energy E_iv (J/kg).
    pub e_iv: f64,
    /// Complete vaporization energy E_cv (J/kg).
    pub e_cv: f64,
    /// Expansion parameters α and β.
    pub alpha: f64,
    /// Expansion parameter β.
    pub beta: f64,
}
impl TillotsonEos {
    /// Create a Tillotson EOS.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        rho0: f64,
        e0: f64,
        a: f64,
        b: f64,
        big_a: f64,
        big_b: f64,
        e_iv: f64,
        e_cv: f64,
        alpha: f64,
        beta: f64,
    ) -> Self {
        Self {
            rho0,
            e0,
            a,
            b,
            big_a,
            big_b,
            e_iv,
            e_cv,
            alpha,
            beta,
        }
    }
    /// Granite (approximate parameters).
    pub fn granite() -> Self {
        Self::new(
            2680.0, 16.0e6, 0.5, 1.3, 18.0e9, 18.0e9, 3.5e6, 18.0e6, 5.0, 5.0,
        )
    }
    /// Basalt (approximate parameters).
    pub fn basalt() -> Self {
        Self::new(
            2860.0, 49.0e6, 0.49, 1.39, 26.7e9, 26.7e9, 4.72e6, 14.0e6, 5.0, 5.0,
        )
    }
    /// Compression ratio η = ρ/ρ₀.
    pub fn eta(&self, density: f64) -> f64 {
        density / self.rho0
    }
    /// Pressure in compressed/hot region.
    pub fn pressure_compressed(&self, density: f64, energy: f64) -> f64 {
        let eta = self.eta(density);
        let mu = eta - 1.0;
        let z = energy / (self.e0 * eta * eta);
        (self.a + self.b / (z + 1.0)) * density * energy + self.big_a * mu + self.big_b * mu * mu
    }
    /// Pressure in expanded region (gas-like).
    pub fn pressure_expanded(&self, density: f64, energy: f64) -> f64 {
        let eta = self.eta(density);
        let mu = eta - 1.0;
        let z = energy / (self.e0 * eta * eta);
        let term1 = self.a * density * energy;
        let term2 = (self.b * density * energy / (z + 1.0)
            + self.big_a * mu * (-self.beta * (1.0 / eta - 1.0)).exp())
            * (-self.alpha * (1.0 / eta - 1.0).powi(2)).exp();
        term1 + term2
    }
    /// Full pressure from density and specific internal energy.
    pub fn pressure_from_energy(&self, density: f64, energy: f64) -> f64 {
        let eta = self.eta(density);
        if eta >= 1.0 || energy <= self.e_iv {
            self.pressure_compressed(density, energy)
        } else if energy >= self.e_cv {
            self.pressure_expanded(density, energy)
        } else {
            let t = (energy - self.e_iv) / (self.e_cv - self.e_iv);
            let pc = self.pressure_compressed(density, energy);
            let pe = self.pressure_expanded(density, energy);
            (1.0 - t) * pc + t * pe
        }
    }
}
/// Soave-Redlich-Kwong equation of state.
///
/// P = R*T/(V_m - b) - a*α(T) / \[V_m*(V_m + b)\]
/// α(T) = \[1 + m*(1 - sqrt(T/Tc))\]²
/// m = 0.480 + 1.574*ω - 0.176*ω²
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct SoaveRedlichKwong {
    /// SRK attraction parameter a (Pa·m⁶/mol²).
    pub a: f64,
    /// SRK co-volume b (m³/mol).
    pub b: f64,
    /// Acentric factor ω.
    pub omega: f64,
    /// Critical temperature Tc (K).
    pub tc: f64,
}
#[allow(dead_code)]
impl SoaveRedlichKwong {
    /// Create a Soave-Redlich-Kwong EOS.
    pub fn new(a: f64, b: f64, omega: f64, tc: f64) -> Self {
        Self { a, b, omega, tc }
    }
    /// Construct from critical properties.
    ///
    /// a = 0.42748 * R² * Tc² / Pc
    /// b = 0.08664 * R * Tc / Pc
    pub fn from_critical(tc: f64, pc: f64, omega: f64) -> Self {
        let a = 0.42748 * R_UNIVERSAL * R_UNIVERSAL * tc * tc / pc;
        let b = 0.08664 * R_UNIVERSAL * tc / pc;
        Self { a, b, omega, tc }
    }
    /// α(T) correction: α = \[1 + m*(1 - sqrt(T/Tc))\]²
    fn alpha(&self, t: f64) -> f64 {
        let m = 0.480 + 1.574 * self.omega - 0.176 * self.omega * self.omega;
        (1.0 + m * (1.0 - (t / self.tc).sqrt())).powi(2)
    }
    /// Pressure at molar volume `v_mol` and temperature `t`.
    pub fn pressure(&self, t: f64, v_mol: f64) -> f64 {
        let a_t = self.a * self.alpha(t);
        if (v_mol - self.b).abs() < 1e-20 {
            return f64::INFINITY;
        }
        R_UNIVERSAL * t / (v_mol - self.b) - a_t / (v_mol * (v_mol + self.b))
    }
}
/// Stiffened gas EOS: p + p_inf = (γ-1) ρ e.
///
/// Generalises ideal gas to account for liquid-like stiffness.
/// Used for water shock and cavitation in underwater explosions.
#[derive(Debug, Clone, Copy)]
pub struct StiffenedGasEos {
    /// Adiabatic index γ.
    pub gamma: f64,
    /// Stiffness pressure p_inf (Pa). For ideal gas, p_inf = 0.
    pub p_inf: f64,
    /// Reference density ρ₀ (kg/m³).
    pub reference_density: f64,
    /// Reference energy e₀ (J/kg).
    pub reference_energy: f64,
}
impl StiffenedGasEos {
    /// Create a new stiffened gas EOS.
    pub fn new(gamma: f64, p_inf: f64, reference_density: f64, reference_energy: f64) -> Self {
        Self {
            gamma,
            p_inf,
            reference_density,
            reference_energy,
        }
    }
    /// Water (Noble-Abel stiffened gas): γ=6.59, p_inf=4.02e8 Pa.
    pub fn water() -> Self {
        Self::new(6.59, 4.02e8, 1000.0, 0.0)
    }
    /// Air treated as stiffened gas (reduces to ideal gas with p_inf=0).
    pub fn air() -> Self {
        Self::new(1.4, 0.0, 1.225, 0.0)
    }
    /// Pressure from density and internal energy.
    pub fn pressure_from_energy(&self, density: f64, energy: f64) -> f64 {
        (self.gamma - 1.0) * density * (energy - self.reference_energy) - self.p_inf
    }
}
/// Noble-Abel equation of state for dense propellant gases.
///
/// P = ρ R T / (1 - b ρ)
///
/// Accounts for co-volume b (repulsive hard-core interactions).
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct NobleAbelEos {
    /// Specific gas constant R (J/(kg·K)).
    pub r_specific: f64,
    /// Co-volume b (m³/kg).
    pub b_covolume: f64,
}
#[allow(dead_code)]
impl NobleAbelEos {
    /// Create a Noble-Abel EOS.
    pub fn new(r_specific: f64, b_covolume: f64) -> Self {
        Self {
            r_specific,
            b_covolume,
        }
    }
    /// Gunpowder propellant gases (approximate).
    pub fn propellant() -> Self {
        Self::new(290.0, 0.001)
    }
    /// Pressure at density and temperature.
    pub fn pressure_at_temperature(&self, density: f64, temperature: f64) -> f64 {
        let denom = 1.0 - self.b_covolume * density;
        if denom <= 0.0 {
            return f64::INFINITY;
        }
        density * self.r_specific * temperature / denom
    }
    /// Density from pressure and temperature.
    pub fn density_at_temperature(&self, pressure: f64, temperature: f64) -> f64 {
        let rt = self.r_specific * temperature;
        if rt.abs() < 1e-30 {
            return 0.0;
        }
        pressure / (rt + pressure * self.b_covolume)
    }
}
/// Simple P-V-T EOS combining Murnaghan cold compression with Grüneisen thermal.
///
/// p(V,T) = p_cold(V) + Γ(V)/V * cv * (T - T₀)
pub struct PvtEos {
    /// Cold compression EOS (Birch-Murnaghan 3rd order).
    pub cold_eos: BirchMurnaghan3Eos,
    /// Grüneisen parameter model.
    pub gruneisen: GruneisenParameter,
    /// Specific heat at constant volume cv (J/(kg·K)).
    pub cv: f64,
    /// Reference temperature T₀ (K).
    pub t0: f64,
}
impl PvtEos {
    /// Create a P-V-T EOS.
    pub fn new(
        cold_eos: BirchMurnaghan3Eos,
        gruneisen: GruneisenParameter,
        cv: f64,
        t0: f64,
    ) -> Self {
        Self {
            cold_eos,
            gruneisen,
            cv,
            t0,
        }
    }
    /// Iron P-V-T EOS (simplified).
    pub fn iron() -> Self {
        let cold = BirchMurnaghan3Eos::iron();
        let v0 = cold.v0;
        let gruneisen = GruneisenParameter::new(1.81, v0, 1.0);
        Self::new(cold, gruneisen, 450.0, 300.0)
    }
    /// Pressure at specific volume V and temperature T.
    pub fn pressure(&self, v: f64, temperature: f64) -> f64 {
        let density = 1.0 / v;
        let p_cold = self.cold_eos.pressure_from_volume(v);
        let p_th = self
            .gruneisen
            .thermal_pressure(density, self.cv, temperature, self.t0);
        p_cold + p_th
    }
    /// Density from pressure and temperature (numerical bisection).
    pub fn density(&self, pressure: f64, temperature: f64) -> f64 {
        let mut lo = 1.0 / (self.cold_eos.v0 * 3.0);
        let mut hi = 1.0 / (self.cold_eos.v0 * 0.3);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if self.pressure(1.0 / mid, temperature) < pressure {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }
}
