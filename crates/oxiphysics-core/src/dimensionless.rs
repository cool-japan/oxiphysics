// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Dimensionless number calculations and similarity analysis.
//!
//! Provides standard dimensionless groups used in fluid mechanics, heat
//! transfer, and related disciplines.  All functions operate on SI-unit
//! inputs and return the corresponding dimensionless ratio.

// ---------------------------------------------------------------------------
// Fluid mechanics
// ---------------------------------------------------------------------------

/// Reynolds number: Re = ρ·v·L / μ.
///
/// Characterises the ratio of inertial to viscous forces.
///
/// - `rho` — fluid density \[kg/m³\].
/// - `v` — characteristic velocity \[m/s\].
/// - `l` — characteristic length \[m\].
/// - `mu` — dynamic viscosity \[Pa·s\].
#[allow(dead_code)]
pub fn reynolds_number(rho: f64, v: f64, l: f64, mu: f64) -> f64 {
    rho * v * l / mu
}

/// Mach number: Ma = v / c_sound.
///
/// Ratio of flow velocity to the local speed of sound.
///
/// - `v` — flow speed \[m/s\].
/// - `c_sound` — speed of sound in the medium \[m/s\].
#[allow(dead_code)]
pub fn mach_number(v: f64, c_sound: f64) -> f64 {
    v / c_sound
}

/// Froude number: Fr = v / √(g·L).
///
/// Ratio of inertial to gravitational forces in free-surface flows.
///
/// - `v` — characteristic velocity \[m/s\].
/// - `g` — gravitational acceleration \[m/s²\].
/// - `l` — characteristic length \[m\].
#[allow(dead_code)]
pub fn froude_number(v: f64, g: f64, l: f64) -> f64 {
    v / (g * l).sqrt()
}

/// Strouhal number: St = f·L / v.
///
/// Describes oscillating flow mechanisms.
///
/// - `f` — frequency of oscillation \[Hz\].
/// - `l` — characteristic length \[m\].
/// - `v` — characteristic velocity \[m/s\].
#[allow(dead_code)]
pub fn strouhal_number(f: f64, l: f64, v: f64) -> f64 {
    f * l / v
}

/// Weber number: We = ρ·v²·L / σ.
///
/// Ratio of inertial force to surface tension.
///
/// - `rho` — fluid density \[kg/m³\].
/// - `v` — characteristic velocity \[m/s\].
/// - `l` — characteristic length \[m\].
/// - `sigma` — surface tension \[N/m\].
#[allow(dead_code)]
pub fn weber_number(rho: f64, v: f64, l: f64, sigma: f64) -> f64 {
    rho * v * v * l / sigma
}

/// Knudsen number: Kn = λ / L.
///
/// Ratio of the molecular mean free path to the characteristic length.
/// Determines when continuum mechanics applies (Kn ≪ 1).
///
/// - `lambda_mfp` — mean free path \[m\].
/// - `l` — characteristic length \[m\].
#[allow(dead_code)]
pub fn knudsen_number(lambda_mfp: f64, l: f64) -> f64 {
    lambda_mfp / l
}

// ---------------------------------------------------------------------------
// Heat and mass transfer
// ---------------------------------------------------------------------------

/// Nusselt number: Nu = h·L / k_thermal.
///
/// Ratio of convective to conductive heat transfer.
///
/// - `h` — convective heat transfer coefficient \[W/(m²·K)\].
/// - `l` — characteristic length \[m\].
/// - `k_thermal` — thermal conductivity \[W/(m·K)\].
#[allow(dead_code)]
pub fn nusselt_number(h: f64, l: f64, k_thermal: f64) -> f64 {
    h * l / k_thermal
}

/// Prandtl number: Pr = μ·cp / k.
///
/// Ratio of momentum diffusivity to thermal diffusivity.
///
/// - `mu` — dynamic viscosity \[Pa·s\].
/// - `cp` — specific heat capacity at constant pressure \[J/(kg·K)\].
/// - `k` — thermal conductivity \[W/(m·K)\].
#[allow(dead_code)]
pub fn prandtl_number(mu: f64, cp: f64, k: f64) -> f64 {
    mu * cp / k
}

/// Grashof number: Gr = g·β·ΔT·L³ / ν².
///
/// Ratio of buoyancy to viscous forces in natural convection.
///
/// - `g` — gravitational acceleration \[m/s²\].
/// - `beta` — volumetric thermal expansion coefficient \[1/K\].
/// - `delta_t` — temperature difference \[K\].
/// - `l` — characteristic length \[m\].
/// - `nu` — kinematic viscosity \[m²/s\].
#[allow(dead_code)]
pub fn grashof_number(g: f64, beta: f64, delta_t: f64, l: f64, nu: f64) -> f64 {
    g * beta * delta_t * l * l * l / (nu * nu)
}

/// Péclet number: Pe = v·L / D.
///
/// Ratio of advective transport rate to diffusive transport rate.
///
/// - `v` — characteristic velocity \[m/s\].
/// - `l` — characteristic length \[m\].
/// - `d` — diffusion coefficient \[m²/s\].
#[allow(dead_code)]
pub fn peclet_number(v: f64, l: f64, d: f64) -> f64 {
    v * l / d
}

// ---------------------------------------------------------------------------
// Derived / combined numbers
// ---------------------------------------------------------------------------

/// Rayleigh number: Ra = Gr·Pr.
///
/// Used in natural convection; equals the product of the Grashof and Prandtl
/// numbers.
#[allow(dead_code)]
pub fn rayleigh_number(g: f64, beta: f64, delta_t: f64, l: f64, nu: f64, alpha: f64) -> f64 {
    g * beta * delta_t * l * l * l / (nu * alpha)
}

/// Biot number: Bi = h·L / k_s.
///
/// Ratio of internal thermal resistance to surface thermal resistance.
///
/// - `h` — convective heat transfer coefficient \[W/(m²·K)\].
/// - `l` — characteristic length \[m\].
/// - `k_s` — solid thermal conductivity \[W/(m·K)\].
#[allow(dead_code)]
pub fn biot_number(h: f64, l: f64, k_s: f64) -> f64 {
    h * l / k_s
}

/// Euler number: Eu = ΔP / (ρ·v²).
///
/// Ratio of pressure forces to inertial forces.
///
/// - `delta_p` — pressure drop \[Pa\].
/// - `rho` — fluid density \[kg/m³\].
/// - `v` — characteristic velocity \[m/s\].
#[allow(dead_code)]
pub fn euler_number(delta_p: f64, rho: f64, v: f64) -> f64 {
    delta_p / (rho * v * v)
}

/// Stokes number: Stk = t_p·v / l.
///
/// Characterises particle-flow coupling; high Stk means particles do not
/// follow streamlines.
///
/// - `t_p` — particle relaxation time \[s\].
/// - `v` — characteristic velocity \[m/s\].
/// - `l` — characteristic length \[m\].
#[allow(dead_code)]
pub fn stokes_number(t_p: f64, v: f64, l: f64) -> f64 {
    t_p * v / l
}

/// Archimedes number: Ar = g·L³·ρ_f·(ρ_p - ρ_f) / μ².
///
/// Determines the motion of fluids due to density differences.
///
/// - `g` — gravitational acceleration \[m/s²\].
/// - `l` — characteristic length \[m\].
/// - `rho_f` — fluid density \[kg/m³\].
/// - `rho_p` — particle/second-phase density \[kg/m³\].
/// - `mu` — dynamic viscosity \[Pa·s\].
#[allow(dead_code)]
pub fn archimedes_number(g: f64, l: f64, rho_f: f64, rho_p: f64, mu: f64) -> f64 {
    g * l * l * l * rho_f * (rho_p - rho_f) / (mu * mu)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Water at ~20 °C
    const RHO_WATER: f64 = 1000.0; // kg/m³
    const MU_WATER: f64 = 1e-3; // Pa·s
    const NU_WATER: f64 = 1e-6; // m²/s

    // Air at ~20 °C
    const C_SOUND_AIR: f64 = 343.0; // m/s
    const G: f64 = 9.81; // m/s²

    #[test]
    fn test_reynolds_pipe_flow() {
        // Pipe: D=0.01 m, v=1 m/s, water → Re = 1000*1*0.01/1e-3 = 10 000
        let re = reynolds_number(RHO_WATER, 1.0, 0.01, MU_WATER);
        assert!((re - 10_000.0).abs() < 1.0);
    }

    #[test]
    fn test_reynolds_proportional_to_velocity() {
        let re1 = reynolds_number(1.0, 1.0, 1.0, 1.0);
        let re2 = reynolds_number(1.0, 2.0, 1.0, 1.0);
        assert!((re2 - 2.0 * re1).abs() < 1e-12);
    }

    #[test]
    fn test_mach_subsonic() {
        let ma = mach_number(100.0, C_SOUND_AIR);
        assert!(ma < 1.0, "100 m/s in air is subsonic");
        assert!((ma - 100.0 / 343.0).abs() < 1e-10);
    }

    #[test]
    fn test_mach_sonic() {
        let ma = mach_number(C_SOUND_AIR, C_SOUND_AIR);
        assert!((ma - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_mach_supersonic() {
        let ma = mach_number(2.0 * C_SOUND_AIR, C_SOUND_AIR);
        assert!((ma - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_froude_number_open_channel() {
        // v=3 m/s, g=9.81, L=1 m → Fr = 3/√9.81 ≈ 0.9579
        let fr = froude_number(3.0, G, 1.0);
        let expected = 3.0 / G.sqrt();
        assert!((fr - expected).abs() < 1e-10);
    }

    #[test]
    fn test_froude_critical_flow() {
        // Critical flow: Fr = 1 when v = √(g·L)
        let l = 2.0;
        let v = (G * l).sqrt();
        let fr = froude_number(v, G, l);
        assert!((fr - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_strouhal_cylinder() {
        // For a cylinder: St ≈ 0.2 at moderate Re
        // f=2 Hz, L=0.1 m, v=1 m/s → St = 0.2
        let st = strouhal_number(2.0, 0.1, 1.0);
        assert!((st - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_nusselt_number_basic() {
        // h=100 W/(m²K), L=0.5 m, k=0.6 W/(mK) → Nu = 100*0.5/0.6 ≈ 83.33
        let nu = nusselt_number(100.0, 0.5, 0.6);
        assert!((nu - 100.0 * 0.5 / 0.6).abs() < 1e-10);
    }

    #[test]
    fn test_prandtl_water() {
        // Water ~20°C: μ=1e-3, cp=4182, k=0.598 → Pr ≈ 7.0
        let pr = prandtl_number(MU_WATER, 4182.0, 0.598);
        assert!(
            pr > 6.0 && pr < 8.0,
            "Prandtl for water should be ~7, got {pr}"
        );
    }

    #[test]
    fn test_prandtl_air() {
        // Air ~20°C: μ=1.81e-5, cp=1005, k=0.0257 → Pr ≈ 0.71
        let pr = prandtl_number(1.81e-5, 1005.0, 0.0257);
        assert!(
            pr > 0.6 && pr < 0.8,
            "Prandtl for air should be ~0.71, got {pr}"
        );
    }

    #[test]
    fn test_grashof_number_positive_delta_t() {
        // g=9.81, β=3.4e-3 K⁻¹ (water), ΔT=10 K, L=0.1 m, ν=1e-6 m²/s
        let gr = grashof_number(G, 3.4e-3, 10.0, 0.1, NU_WATER);
        assert!(gr > 0.0);
    }

    #[test]
    fn test_grashof_scales_with_l_cubed() {
        let gr1 = grashof_number(G, 1e-3, 10.0, 1.0, 1e-6);
        let gr2 = grashof_number(G, 1e-3, 10.0, 2.0, 1e-6);
        assert!((gr2 - 8.0 * gr1).abs() < 1e-6 * gr1.abs());
    }

    #[test]
    fn test_weber_number_basic() {
        // ρ=1000, v=1, L=0.01, σ=0.072 (water-air) → We = 1000*1*0.01/0.072 ≈ 138.9
        let we = weber_number(1000.0, 1.0, 0.01, 0.072);
        assert!((we - 1000.0 * 0.01 / 0.072).abs() < 0.1);
    }

    #[test]
    fn test_knudsen_continuum_limit() {
        // λ_mfp = 70 nm, L = 1 m → Kn ~ 7e-8 ≪ 1 (continuum)
        let kn = knudsen_number(70e-9, 1.0);
        assert!(kn < 1e-6, "Kn should be tiny in continuum regime");
    }

    #[test]
    fn test_knudsen_free_molecular() {
        // λ_mfp ~ L → Kn ~ 1 (free molecular regime)
        let kn = knudsen_number(1.0, 1.0);
        assert!((kn - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_peclet_number_basic() {
        // v=1, L=1, D=1e-3 → Pe = 1000
        let pe = peclet_number(1.0, 1.0, 1e-3);
        assert!((pe - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn test_peclet_proportional_to_velocity() {
        let pe1 = peclet_number(1.0, 1.0, 1.0);
        let pe2 = peclet_number(3.0, 1.0, 1.0);
        assert!((pe2 - 3.0 * pe1).abs() < 1e-12);
    }

    #[test]
    fn test_rayleigh_number_positive() {
        // α ≈ 1.4e-7 m²/s (water thermal diffusivity)
        let ra = rayleigh_number(G, 3.4e-3, 10.0, 0.1, NU_WATER, 1.4e-7);
        assert!(ra > 0.0);
    }

    #[test]
    fn test_biot_number_small_lump() {
        // Bi < 0.1 → lumped capacitance valid
        let bi = biot_number(10.0, 0.005, 200.0); // aluminium
        assert!(bi < 0.1, "Bi = {bi} should be < 0.1");
    }

    #[test]
    fn test_biot_number_large() {
        let bi = biot_number(1000.0, 0.1, 1.0);
        assert!(bi > 1.0);
    }

    #[test]
    fn test_euler_number_definition() {
        // ΔP = 100 Pa, ρ=1, v=10 → Eu = 100/(1*100) = 1
        let eu = euler_number(100.0, 1.0, 10.0);
        assert!((eu - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_stokes_number_unity() {
        // t_p=1, v=1, L=1 → Stk=1
        let stk = stokes_number(1.0, 1.0, 1.0);
        assert!((stk - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_stokes_number_small() {
        // Small particle: t_p very small → Stk ≪ 1
        let stk = stokes_number(1e-6, 1.0, 0.1);
        assert!(stk < 1e-4);
    }

    #[test]
    fn test_archimedes_number_positive_buoyancy() {
        // ρ_p > ρ_f → Ar > 0
        let ar = archimedes_number(G, 0.01, 1000.0, 2500.0, MU_WATER);
        assert!(ar > 0.0);
    }

    #[test]
    fn test_archimedes_number_negative_buoyancy() {
        // ρ_p < ρ_f → Ar < 0 (particle floats)
        let ar = archimedes_number(G, 0.01, 1000.0, 500.0, MU_WATER);
        assert!(ar < 0.0);
    }

    #[test]
    fn test_all_dimensionless_finite() {
        // Quick sanity: no NaN or infinity for positive inputs
        assert!(reynolds_number(1.0, 1.0, 1.0, 1.0).is_finite());
        assert!(mach_number(1.0, 1.0).is_finite());
        assert!(froude_number(1.0, 1.0, 1.0).is_finite());
        assert!(strouhal_number(1.0, 1.0, 1.0).is_finite());
        assert!(nusselt_number(1.0, 1.0, 1.0).is_finite());
        assert!(prandtl_number(1.0, 1.0, 1.0).is_finite());
        assert!(grashof_number(1.0, 1.0, 1.0, 1.0, 1.0).is_finite());
        assert!(weber_number(1.0, 1.0, 1.0, 1.0).is_finite());
        assert!(knudsen_number(1.0, 1.0).is_finite());
        assert!(peclet_number(1.0, 1.0, 1.0).is_finite());
    }

    #[test]
    fn test_reynolds_turbulent_threshold() {
        // Pipe flow transitions to turbulence at Re ≈ 4000
        let re_laminar = reynolds_number(1000.0, 0.1, 0.01, 1e-3); // Re=1000
        let re_turbulent = reynolds_number(1000.0, 1.0, 0.01, 1e-3); // Re=10000
        assert!(re_laminar < 4000.0);
        assert!(re_turbulent > 4000.0);
    }
}
