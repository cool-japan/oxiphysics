// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Dielectric and electromagnetic material models.
//!
//! Implements:
//!
//! - Debye relaxation model for complex permittivity
//! - Skin depth, wave impedance, Fresnel reflection
//! - Dielectric loss tangent and Snell refraction

// ---------------------------------------------------------------------------
// DielectricMaterial
// ---------------------------------------------------------------------------

/// A dielectric / electromagnetic material defined by its permittivity,
/// permeability, and DC conductivity.
#[derive(Debug, Clone)]
pub struct DielectricMaterial {
    /// High-frequency (optical) relative permittivity ε_∞.
    pub epsilon_inf: f64,
    /// Static (DC) relative permittivity ε_s.
    pub epsilon_static: f64,
    /// Relative permeability μ_r.
    pub mu_r: f64,
    /// DC electrical conductivity σ \[S/m\].
    pub conductivity: f64,
    /// Debye relaxation time τ \[s\].
    pub relaxation_time: f64,
}

/// Vacuum permittivity ε₀ \[F/m\].
const EPSILON_0: f64 = 8.854_187_817e-12;
/// Vacuum permeability μ₀ \[H/m\].
const MU_0: f64 = 1.256_637_061_4e-6;

impl DielectricMaterial {
    /// Create a new dielectric material.
    ///
    /// # Parameters
    /// - `epsilon_inf`      : High-frequency relative permittivity.
    /// - `epsilon_static`   : Static (DC) relative permittivity.
    /// - `mu_r`             : Relative permeability.
    /// - `conductivity`     : DC conductivity \[S/m\].
    /// - `relaxation_time`  : Debye relaxation time τ \[s\].
    pub fn new(
        epsilon_inf: f64,
        epsilon_static: f64,
        mu_r: f64,
        conductivity: f64,
        relaxation_time: f64,
    ) -> Self {
        Self {
            epsilon_inf,
            epsilon_static,
            mu_r,
            conductivity,
            relaxation_time,
        }
    }

    /// Compute the **complex permittivity** using the Debye relaxation model.
    ///
    /// ε(ω) = ε_∞ + (ε_s − ε_∞) / (1 + iωτ)
    ///
    /// Returns `(ε', ε'')` — real and imaginary parts of the relative
    /// permittivity (the imaginary part is negative for passive materials).
    pub fn complex_permittivity(&self, omega: f64) -> (f64, f64) {
        let delta_eps = self.epsilon_static - self.epsilon_inf;
        let denom = 1.0 + (omega * self.relaxation_time).powi(2);
        let eps_real = self.epsilon_inf + delta_eps / denom;
        let eps_imag = -delta_eps * omega * self.relaxation_time / denom;
        (eps_real, eps_imag)
    }

    /// Compute the **skin depth** δ \[m\] at angular frequency ω.
    ///
    /// δ = √(2 / (ω μ σ))
    ///
    /// Returns `None` if ω or σ is zero (skin depth is infinite).
    pub fn skin_depth(&self, omega: f64) -> Option<f64> {
        if omega <= 0.0 || self.conductivity <= 0.0 {
            return None;
        }
        let mu = self.mu_r * MU_0;
        Some((2.0 / (omega * mu * self.conductivity)).sqrt())
    }

    /// Compute the **wave impedance** Z \[Ω\] at angular frequency ω.
    ///
    /// Z = √(μ / ε)  where both μ and ε are absolute (not relative).
    pub fn wave_impedance(&self, omega: f64) -> f64 {
        let (eps_real, eps_imag) = self.complex_permittivity(omega);
        let mu = self.mu_r * MU_0;
        let eps_abs_real = eps_real * EPSILON_0;
        let eps_abs_imag = eps_imag * EPSILON_0;
        // |ε| = sqrt(ε'^2 + ε''^2)
        let eps_mag = (eps_abs_real.powi(2) + eps_abs_imag.powi(2)).sqrt();
        (mu / eps_mag).sqrt()
    }

    /// Compute the **Fresnel reflection coefficient** R at normal incidence
    /// between this material and a second material.
    ///
    /// R = ((Z2 - Z1) / (Z2 + Z1))²
    pub fn reflection_coefficient(&self, other: &DielectricMaterial, omega: f64) -> f64 {
        let z1 = self.wave_impedance(omega);
        let z2 = other.wave_impedance(omega);
        let r = (z2 - z1) / (z2 + z1);
        r * r
    }

    /// Compute the **dielectric loss tangent** tan(δ) = ε'' / ε' at ω.
    pub fn dielectric_loss_tangent(&self, omega: f64) -> f64 {
        let (eps_real, eps_imag) = self.complex_permittivity(omega);
        if eps_real.abs() < f64::EPSILON {
            return f64::INFINITY;
        }
        (-eps_imag / eps_real).abs()
    }

    /// Compute the **refraction angle** θ_t \[radians\] using Snell's law.
    ///
    /// n1 sin(θ_i) = n2 sin(θ_t)
    ///
    /// where the refractive index n = √(ε_r · μ_r).
    ///
    /// Returns `None` for total internal reflection.
    pub fn snell_refraction(&self, other: &DielectricMaterial, theta_i: f64) -> Option<f64> {
        let n1 = (self.epsilon_static * self.mu_r).sqrt();
        let n2 = (other.epsilon_static * other.mu_r).sqrt();
        let sin_t = n1 * theta_i.sin() / n2;
        if sin_t.abs() > 1.0 {
            None // total internal reflection
        } else {
            Some(sin_t.asin())
        }
    }
}

/// Compute the **attenuation coefficient** α \[Np/m\] for a plane wave.
///
/// α = ω √(με/2) · √(√(1 + (σ/(ωε))²) − 1)
pub fn attenuation_coefficient(material: &DielectricMaterial, omega: f64) -> f64 {
    if omega <= 0.0 {
        return 0.0;
    }
    let (eps_real, _) = material.complex_permittivity(omega);
    let mu = material.mu_r * MU_0;
    let eps = eps_real * EPSILON_0;
    let sigma_over_oe = material.conductivity / (omega * eps.abs().max(f64::EPSILON));
    let inner = (1.0 + sigma_over_oe.powi(2)).sqrt() - 1.0;
    omega * (mu * eps.abs() / 2.0).sqrt() * inner.sqrt()
}

/// Compute the **phase velocity** v_p \[m/s\] in a non-conducting dielectric.
///
/// v_p = c / √(ε_r · μ_r)
pub fn phase_velocity(material: &DielectricMaterial) -> f64 {
    let c = 1.0 / (EPSILON_0 * MU_0).sqrt();
    c / (material.epsilon_static * material.mu_r).sqrt()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Water-like dielectric (ε_s ≈ 80, ε_∞ ≈ 5.5, τ ≈ 9.4 ps).
    fn water_like() -> DielectricMaterial {
        DielectricMaterial::new(5.5, 80.0, 1.0, 0.0, 9.4e-12)
    }

    /// Air/vacuum material.
    fn vacuum() -> DielectricMaterial {
        DielectricMaterial::new(1.0, 1.0, 1.0, 0.0, 0.0)
    }

    // 1. At ω→0 the real part of permittivity approaches ε_s.
    #[test]
    fn test_debye_dc_limit() {
        let mat = water_like();
        let (eps_r, eps_i) = mat.complex_permittivity(1.0); // very low frequency
        assert!(
            (eps_r - mat.epsilon_static).abs() < 1.0,
            "At low ω, ε' ≈ ε_s, got {eps_r}"
        );
        assert!(eps_i.abs() < 1.0, "At low ω, ε'' ≈ 0, got {eps_i}");
    }

    // 2. At ω→∞ the real part approaches ε_∞.
    #[test]
    fn test_debye_high_freq_limit() {
        let mat = water_like();
        let (eps_r, eps_i) = mat.complex_permittivity(1e20);
        assert!(
            (eps_r - mat.epsilon_inf).abs() < 0.1,
            "At high ω, ε' ≈ ε_∞ = {}, got {eps_r}",
            mat.epsilon_inf
        );
        assert!(eps_i.abs() < 0.1, "At high ω, ε'' ≈ 0, got {eps_i}");
    }

    // 3. At ωτ = 1 (relaxation frequency), loss is maximum.
    #[test]
    fn test_debye_relaxation_peak() {
        let mat = water_like();
        let omega_r = 1.0 / mat.relaxation_time;
        let (_, eps_i_r) = mat.complex_permittivity(omega_r);
        let (_, eps_i_low) = mat.complex_permittivity(omega_r * 0.01);
        let (_, eps_i_high) = mat.complex_permittivity(omega_r * 100.0);
        assert!(
            eps_i_r.abs() > eps_i_low.abs(),
            "Loss at relaxation should exceed loss at low ω"
        );
        assert!(
            eps_i_r.abs() > eps_i_high.abs(),
            "Loss at relaxation should exceed loss at high ω"
        );
    }

    // 4. Imaginary part is always ≤ 0 (passive material, energy dissipation).
    #[test]
    fn test_debye_imaginary_nonpositive() {
        let mat = water_like();
        for exp in 0..20 {
            let omega = 10f64.powi(exp);
            let (_, eps_i) = mat.complex_permittivity(omega);
            assert!(
                eps_i <= 1e-12,
                "Imaginary part must be ≤ 0 for passive Debye, got {eps_i} at ω={omega}"
            );
        }
    }

    // 5. Skin depth decreases with increasing frequency.
    #[test]
    fn test_skin_depth_decreases_with_freq() {
        let copper = DielectricMaterial::new(1.0, 1.0, 1.0, 5.8e7, 0.0);
        let d1 = copper.skin_depth(1e3).unwrap();
        let d2 = copper.skin_depth(1e6).unwrap();
        let d3 = copper.skin_depth(1e9).unwrap();
        assert!(d1 > d2, "Skin depth should decrease with frequency");
        assert!(d2 > d3, "Skin depth should decrease with frequency");
    }

    // 6. Skin depth scales as 1/√ω.
    #[test]
    fn test_skin_depth_scaling() {
        let mat = DielectricMaterial::new(1.0, 1.0, 1.0, 1e6, 0.0);
        let d1 = mat.skin_depth(1e6).unwrap();
        let d2 = mat.skin_depth(4e6).unwrap(); // 4× frequency → depth halves
        let ratio = d1 / d2;
        assert!(
            (ratio - 2.0).abs() < 1e-6,
            "Skin depth should scale as 1/√ω, ratio={ratio}"
        );
    }

    // 7. Skin depth returns None for zero conductivity.
    #[test]
    fn test_skin_depth_none_for_zero_conductivity() {
        let mat = water_like();
        assert!(
            mat.skin_depth(1e9).is_none(),
            "Skin depth should be None for zero conductivity"
        );
    }

    // 8. Wave impedance of vacuum ≈ 377 Ω.
    #[test]
    fn test_wave_impedance_vacuum() {
        let vac = vacuum();
        let z = vac.wave_impedance(1e9);
        assert!(
            (z - 377.0).abs() < 1.0,
            "Vacuum impedance should be ≈ 377 Ω, got {z}"
        );
    }

    // 9. Wave impedance scales inversely with √ε_r.
    #[test]
    fn test_wave_impedance_scales_with_eps() {
        let mat1 = DielectricMaterial::new(1.0, 1.0, 1.0, 0.0, 0.0);
        let mat4 = DielectricMaterial::new(1.0, 4.0, 1.0, 0.0, 0.0);
        let z1 = mat1.wave_impedance(1e9);
        let z4 = mat4.wave_impedance(1e9);
        assert!(
            (z1 / z4 - 2.0).abs() < 1e-4,
            "Impedance ratio should be 2 for ε_r ratio of 4, got {}",
            z1 / z4
        );
    }

    // 10. Reflection at vacuum-vacuum interface is zero.
    #[test]
    fn test_reflection_vacuum_vacuum() {
        let vac1 = vacuum();
        let vac2 = vacuum();
        let r = vac1.reflection_coefficient(&vac2, 1e9);
        assert!(
            r < 1e-20,
            "Reflection at same-material boundary should be 0, got {r}"
        );
    }

    // 11. Reflection coefficient is between 0 and 1.
    #[test]
    fn test_reflection_in_range() {
        let vac = vacuum();
        let water = water_like();
        let r = vac.reflection_coefficient(&water, 1e9);
        assert!(r >= 0.0, "Reflection coefficient must be non-negative");
        assert!(r <= 1.0, "Reflection coefficient must be ≤ 1, got {r}");
    }

    // 12. Fresnel reflection is symmetric (R(1→2) = R(2→1)).
    #[test]
    fn test_reflection_symmetry() {
        let vac = vacuum();
        let water = water_like();
        let r12 = vac.reflection_coefficient(&water, 1e9);
        let r21 = water.reflection_coefficient(&vac, 1e9);
        assert!(
            (r12 - r21).abs() < 1e-10,
            "Reflection should be symmetric: r12={r12}, r21={r21}"
        );
    }

    // 13. Loss tangent is zero for lossless Debye at very high frequency.
    #[test]
    fn test_loss_tangent_zero_at_high_freq() {
        let mat = water_like();
        let tan_d = mat.dielectric_loss_tangent(1e20);
        assert!(
            tan_d < 1e-3,
            "Loss tangent should be ~0 at very high ω, got {tan_d}"
        );
    }

    // 14. Loss tangent is positive.
    #[test]
    fn test_loss_tangent_positive() {
        let mat = water_like();
        let tan_d = mat.dielectric_loss_tangent(1.0 / mat.relaxation_time);
        assert!(tan_d > 0.0, "Loss tangent must be positive, got {tan_d}");
    }

    // 15. Snell's law: normal incidence gives zero refraction angle.
    #[test]
    fn test_snell_normal_incidence() {
        let vac = vacuum();
        let glass = DielectricMaterial::new(1.0, 2.25, 1.0, 0.0, 0.0);
        let theta_t = vac.snell_refraction(&glass, 0.0).unwrap();
        assert!(
            theta_t.abs() < 1e-10,
            "Normal incidence: θ_t should be 0, got {theta_t}"
        );
    }

    // 16. Snell's law from denser to less dense medium: correct bending.
    #[test]
    fn test_snell_refraction_bends_away() {
        let glass = DielectricMaterial::new(1.0, 2.25, 1.0, 0.0, 0.0);
        let vac = vacuum();
        let theta_i = 30.0_f64.to_radians();
        let theta_t = glass.snell_refraction(&vac, theta_i).unwrap();
        assert!(
            theta_t > theta_i,
            "Going denser→rarer, angle must increase: θ_i={theta_i:.4}, θ_t={theta_t:.4}"
        );
    }

    // 17. Snell's law: total internal reflection.
    #[test]
    fn test_snell_total_internal_reflection() {
        let glass = DielectricMaterial::new(1.0, 2.25, 1.0, 0.0, 0.0);
        let vac = vacuum();
        // Critical angle for n=1.5 glass → vacuum: θ_c = arcsin(1/1.5) ≈ 41.8°.
        let theta_i = 80.0_f64.to_radians(); // well above critical angle
        let result = glass.snell_refraction(&vac, theta_i);
        assert!(
            result.is_none(),
            "Should be total internal reflection at {theta_i:.2} rad"
        );
    }

    // 18. phase_velocity in vacuum ≈ speed of light.
    #[test]
    fn test_phase_velocity_vacuum() {
        let vac = vacuum();
        let v = phase_velocity(&vac);
        let c = 3e8;
        assert!(
            (v - c).abs() / c < 0.01,
            "Phase velocity in vacuum should be ≈ c, got {v:.3e}"
        );
    }

    // 19. phase_velocity in glass is slower.
    #[test]
    fn test_phase_velocity_glass_slower() {
        let vac = vacuum();
        let glass = DielectricMaterial::new(1.0, 2.25, 1.0, 0.0, 0.0);
        let v_vac = phase_velocity(&vac);
        let v_glass = phase_velocity(&glass);
        assert!(
            v_glass < v_vac,
            "Phase velocity in glass should be slower than vacuum"
        );
    }

    // 20. Debye ε' is monotone decreasing with ω.
    #[test]
    fn test_debye_eps_real_decreasing() {
        let mat = water_like();
        let freqs = [1e6, 1e9, 1e11, 1e13];
        let mut prev_eps = f64::INFINITY;
        for &omega in &freqs {
            let (eps_r, _) = mat.complex_permittivity(omega);
            assert!(
                eps_r < prev_eps + 1e-10,
                "ε' should decrease monotonically with ω"
            );
            prev_eps = eps_r;
        }
    }

    // 21. Skin depth is always positive when defined.
    #[test]
    fn test_skin_depth_positive() {
        let mat = DielectricMaterial::new(1.0, 1.0, 1.0, 1e4, 0.0);
        for exp in 1..10 {
            let omega = 10f64.powi(exp);
            let d = mat.skin_depth(omega).unwrap();
            assert!(d > 0.0, "Skin depth must be positive, got {d}");
        }
    }

    // 22. Skin depth formula: δ = sqrt(2/(ωμσ)) matches implementation.
    #[test]
    fn test_skin_depth_formula() {
        let sigma = 1e5_f64;
        let mat = DielectricMaterial::new(1.0, 1.0, 1.0, sigma, 0.0);
        let omega = 2.0 * PI * 1e6;
        let expected = (2.0 / (omega * MU_0 * sigma)).sqrt();
        let got = mat.skin_depth(omega).unwrap();
        assert!(
            (got - expected).abs() / expected < 1e-6,
            "Skin depth mismatch: expected {expected:.6e}, got {got:.6e}"
        );
    }

    // 23. wave_impedance grows with μ_r.
    #[test]
    fn test_wave_impedance_mu_dependence() {
        let mat1 = DielectricMaterial::new(1.0, 1.0, 1.0, 0.0, 0.0);
        let mat4 = DielectricMaterial::new(1.0, 1.0, 4.0, 0.0, 0.0);
        let z1 = mat1.wave_impedance(1e9);
        let z4 = mat4.wave_impedance(1e9);
        assert!(z4 > z1, "Higher μ_r should give higher wave impedance");
    }

    // 24. attenuation_coefficient is positive in a lossy medium.
    #[test]
    fn test_attenuation_positive() {
        let mat = DielectricMaterial::new(1.0, 1.0, 1.0, 1e6, 0.0);
        let alpha = attenuation_coefficient(&mat, 2.0 * PI * 1e9);
        assert!(
            alpha > 0.0,
            "Attenuation should be positive in lossy medium, got {alpha}"
        );
    }

    // 25. attenuation_coefficient is zero at zero frequency.
    #[test]
    fn test_attenuation_zero_at_zero_freq() {
        let mat = DielectricMaterial::new(1.0, 1.0, 1.0, 1e6, 0.0);
        let alpha = attenuation_coefficient(&mat, 0.0);
        assert!(
            alpha.abs() < 1e-20,
            "Attenuation at ω=0 should be 0, got {alpha}"
        );
    }
}
