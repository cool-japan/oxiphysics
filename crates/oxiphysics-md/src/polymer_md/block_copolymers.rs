// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Block copolymer and crystallization models.

use super::block_copolymer_chi;

/// A-B block copolymer microphase separation.
///
/// Tracks chi*N parameter and computes lamellar spacing.
#[derive(Debug, Clone)]
pub struct BlockCopolymer {
    /// Number of A monomers.
    pub n_a: usize,
    /// Number of B monomers.
    pub n_b: usize,
    /// Flory-Huggins chi parameter.
    pub chi: f64,
    /// Statistical segment length b.
    pub b: f64,
    /// Temperature K.
    pub temperature: f64,
}

impl BlockCopolymer {
    /// Create a new block copolymer.
    pub fn new(n_a: usize, n_b: usize, chi: f64, b: f64, temperature: f64) -> Self {
        Self {
            n_a,
            n_b,
            chi,
            b,
            temperature,
        }
    }

    /// Total degree of polymerization.
    pub fn n_total(&self) -> usize {
        self.n_a + self.n_b
    }

    /// Volume fraction of A block.
    pub fn f_a(&self) -> f64 {
        self.n_a as f64 / self.n_total() as f64
    }

    /// Chi*N interaction parameter.
    pub fn chi_n(&self) -> f64 {
        block_copolymer_chi(self.chi, self.n_total())
    }

    /// Critical chi*N for microphase separation (symmetric: 10.5).
    pub fn chi_n_critical(&self) -> f64 {
        // For symmetric f=0.5: (chi*N)_c = 10.495
        // For asymmetric: varies
        let f = self.f_a();
        10.495 / (4.0 * f * (1.0 - f))
    }

    /// Is the system phase-separated?
    pub fn is_phase_separated(&self) -> bool {
        self.chi_n() > self.chi_n_critical()
    }

    /// Lamellar periodicity d = 1.95 * b * N^(2/3) (strong segregation).
    pub fn lamellar_spacing(&self) -> f64 {
        1.95 * self.b * (self.n_total() as f64).powf(2.0 / 3.0)
    }

    /// Interfacial width w = b / sqrt(6*chi).
    pub fn interfacial_width(&self) -> f64 {
        if self.chi < 1e-15 {
            return f64::INFINITY;
        }
        self.b / (6.0 * self.chi).sqrt()
    }

    /// Free energy of microphase separation (Landau).
    pub fn free_energy_density(&self) -> f64 {
        let f = self.f_a();
        let chi_n = self.chi_n();
        let chi_n_c = self.chi_n_critical();
        // Landau expansion: f ~ a*(chi_n - chi_nc)/chi_nc * m^2 + b*m^4
        let m = (f - 0.5).abs();
        let a = (chi_n - chi_n_c) / chi_n_c;
        a * m * m + m * m * m * m
    }
}

/// Polymer crystallization model.
///
/// Implements Lauritzen-Hoffman nucleation theory, crystal growth, and lamellar thickness.
#[derive(Debug, Clone)]
pub struct PolymerCrystallization {
    /// Temperature K.
    pub temperature: f64,
    /// Equilibrium melting temperature Tm0 in K.
    pub tm0: f64,
    /// Surface free energy sigma_e in J/m².
    pub sigma_e: f64,
    /// Lateral surface energy sigma in J/m².
    pub sigma: f64,
    /// Enthalpy of fusion delta_hf in J/mol.
    pub delta_hf: f64,
    /// Pre-exponential factor G0.
    pub g0: f64,
    /// Transport activation energy U* (J/mol).
    pub u_star: f64,
    /// Vogel temperature T_inf (K).
    pub t_inf: f64,
    /// Gas constant.
    pub r_gas: f64,
}

impl PolymerCrystallization {
    /// Create a new polymer crystallization model.
    pub fn new(
        temperature: f64,
        tm0: f64,
        sigma_e: f64,
        sigma: f64,
        delta_hf: f64,
        g0: f64,
        u_star: f64,
        t_inf: f64,
    ) -> Self {
        Self {
            temperature,
            tm0,
            sigma_e,
            sigma,
            delta_hf,
            g0,
            u_star,
            t_inf,
            r_gas: 8.314,
        }
    }

    /// Undercooling delta_T = Tm0 - T.
    pub fn undercooling(&self) -> f64 {
        self.tm0 - self.temperature
    }

    /// Critical lamellar thickness l* = 2*sigma_e*Tm0 / (delta_hf * delta_T).
    pub fn critical_lamellar_thickness(&self) -> f64 {
        let dt = self.undercooling().max(1.0);
        2.0 * self.sigma_e * self.tm0 / (self.delta_hf * dt)
    }

    /// Observed lamellar thickness l = l* + delta_l.
    pub fn lamellar_thickness(&self) -> f64 {
        let l_star = self.critical_lamellar_thickness();
        l_star * 1.1 // empirical constant
    }

    /// Lauritzen-Hoffman nucleation rate G = G0 * exp(-U*/(R*(T-T_inf))) * exp(-Kg/(T*delta_T)).
    pub fn nucleation_rate(&self) -> f64 {
        let t = self.temperature;
        let dt = self.undercooling().max(1.0);
        let transport = (-self.u_star / (self.r_gas * (t - self.t_inf))).exp();
        let kg = 4.0 * self.sigma * self.sigma_e * self.tm0 / (self.delta_hf * self.r_gas);
        let nucleation = (-kg / (t * dt)).exp();
        self.g0 * transport * nucleation
    }

    /// Observed crystallization exponent n (Avrami).
    pub fn avrami_exponent(&self) -> f64 {
        3.0 // typical for 3D spherulitic crystallization
    }

    /// Avrami crystallinity at time t.
    pub fn crystallinity(&self, t: f64, k_avrami: f64) -> f64 {
        let n = self.avrami_exponent();
        1.0 - (-k_avrami * t.powf(n)).exp()
    }
}
