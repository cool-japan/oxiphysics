// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Melt and network models: PolymerMelt, PolymerNetwork, FloryHuggins,
//! PolymerBrush, Polyelectrolyte, DebyeHuckel.

use std::f64::consts::PI;

/// Entangled polymer melt (tube model).
///
/// Implements reptation dynamics: relaxation time tau ~ N³, viscosity η ~ N³.
#[derive(Debug, Clone)]
pub struct PolymerMelt {
    /// Degree of polymerization.
    pub n: usize,
    /// Entanglement length Ne.
    pub n_entangle: usize,
    /// Segment friction coefficient zeta.
    pub zeta: f64,
    /// Bond length b.
    pub b: f64,
    /// Thermal energy kT.
    pub kt: f64,
    /// Monomer mass.
    pub mass: f64,
}

impl PolymerMelt {
    /// Create a new polymer melt model.
    pub fn new(n: usize, n_entangle: usize, zeta: f64, b: f64, kt: f64, mass: f64) -> Self {
        Self {
            n,
            n_entangle,
            zeta,
            b,
            kt,
            mass,
        }
    }

    /// Rouse time tau_R = zeta * N^2 * b^2 / (3*pi^2*kT).
    pub fn rouse_time(&self) -> f64 {
        self.zeta * (self.n as f64).powi(2) * self.b * self.b
            / (3.0 * std::f64::consts::PI.powi(2) * self.kt)
    }

    /// Reptation (disengagement) time tau_d ~ N^3.
    pub fn reptation_time(&self) -> f64 {
        let n = self.n as f64;
        let ne = self.n_entangle as f64;
        3.0 * self.rouse_time() * n / ne
    }

    /// Zero-shear viscosity eta_0 ~ N^3.
    pub fn viscosity(&self) -> f64 {
        let tau_d = self.reptation_time();
        let n = self.n as f64;
        let ne = self.n_entangle as f64;
        let rho = self.mass * n / (self.b * self.b * self.b * n);
        rho * self.kt * tau_d * n / ne
    }

    /// Diffusion coefficient D ~ 1/N^2.
    pub fn diffusion_coefficient(&self) -> f64 {
        let n = self.n as f64;
        let ne = self.n_entangle as f64;
        self.kt * ne / (self.zeta * n * n)
    }

    /// Number of entanglements per chain.
    pub fn z_entanglements(&self) -> f64 {
        self.n as f64 / self.n_entangle as f64
    }
}

/// Cross-linked polymer network.
///
/// Models network modulus G_N^0 and junction behavior.
#[derive(Debug, Clone)]
pub struct PolymerNetwork {
    /// Number of network strands.
    pub n_strands: usize,
    /// Number of junctions.
    pub n_junctions: usize,
    /// Strand molecular weight Mc.
    pub mc: f64,
    /// Density kg/m³.
    pub density: f64,
    /// Temperature K.
    pub temperature: f64,
    /// Gas constant J/(mol·K).
    pub r_gas: f64,
    /// Network modulus.
    pub g_n0: f64,
}

impl PolymerNetwork {
    /// Create a new polymer network.
    pub fn new(
        n_strands: usize,
        n_junctions: usize,
        mc: f64,
        density: f64,
        temperature: f64,
    ) -> Self {
        let r_gas = 8.314;
        let g_n0 = density * r_gas * temperature / mc;
        Self {
            n_strands,
            n_junctions,
            mc,
            density,
            temperature,
            r_gas,
            g_n0,
        }
    }

    /// Network modulus G_N^0 = rho * R * T / Mc.
    pub fn modulus(&self) -> f64 {
        self.g_n0
    }

    /// Crosslink density nu = rho / Mc.
    pub fn crosslink_density(&self) -> f64 {
        self.density / self.mc
    }

    /// Swelling ratio Q for crosslinked network.
    pub fn swelling_ratio(&self, chi: f64) -> f64 {
        // Flory-Rehner: estimate equilibrium swelling
        let nu = self.crosslink_density();
        // Simplified: Q^(5/3) ~ 1 / (2*chi * nu * Vc)
        let vc = self.mc / self.density; // molar volume
        (1.0 / (2.0 * chi * nu * vc)).powf(3.0 / 5.0)
    }

    /// Effective functionality of junctions.
    pub fn functionality(&self) -> f64 {
        if self.n_junctions == 0 {
            return 0.0;
        }
        2.0 * self.n_strands as f64 / self.n_junctions as f64
    }
}

/// Flory-Huggins free energy for polymer blends / block copolymers.
///
/// F/(nkT) = φ ln(φ)/N_A + (1-φ) ln(1-φ)/N_B + χ φ(1-φ)
#[derive(Debug, Clone)]
pub struct FloryHuggins {
    /// Degree of polymerization of component A.
    pub n_a: usize,
    /// Degree of polymerization of component B.
    pub n_b: usize,
    /// Flory-Huggins interaction parameter χ.
    pub chi: f64,
}

impl FloryHuggins {
    /// Create a new Flory-Huggins model.
    pub fn new(n_a: usize, n_b: usize, chi: f64) -> Self {
        FloryHuggins { n_a, n_b, chi }
    }

    /// Free energy of mixing per kT per monomer.
    pub fn free_energy(&self, phi: f64) -> f64 {
        let phi = phi.clamp(1e-10, 1.0 - 1e-10);
        let na = self.n_a as f64;
        let nb = self.n_b as f64;
        phi * phi.ln() / na + (1.0 - phi) * (1.0 - phi).ln() / nb + self.chi * phi * (1.0 - phi)
    }

    /// Second derivative d²F/dφ² (for spinodal condition).
    pub fn d2f_dphi2(&self, phi: f64) -> f64 {
        let phi = phi.clamp(1e-10, 1.0 - 1e-10);
        let na = self.n_a as f64;
        let nb = self.n_b as f64;
        1.0 / (na * phi) + 1.0 / (nb * (1.0 - phi)) - 2.0 * self.chi
    }

    /// Spinodal condition: d²F/dφ² = 0 → χ_s = (1/√N_A + 1/√N_B)² / 2.
    pub fn chi_spinodal(&self) -> f64 {
        let na = self.n_a as f64;
        let nb = self.n_b as f64;
        let val = 1.0 / na.sqrt() + 1.0 / nb.sqrt();
        val * val / 2.0
    }

    /// Critical composition φ_c = √N_B / (√N_A + √N_B).
    pub fn phi_critical(&self) -> f64 {
        let na = self.n_a as f64;
        let nb = self.n_b as f64;
        nb.sqrt() / (na.sqrt() + nb.sqrt())
    }

    /// Critical χ parameter.
    pub fn chi_critical(&self) -> f64 {
        self.chi_spinodal()
    }

    /// Is the system phase-separated (above spinodal)?
    pub fn is_phase_separated(&self, phi: f64) -> bool {
        self.d2f_dphi2(phi) < 0.0
    }

    /// Order parameter m = φ - φ_c.
    pub fn order_parameter(&self, phi: f64) -> f64 {
        phi - self.phi_critical()
    }
}

/// Polymer brush simulation (Alexander-de Gennes theory).
///
/// A brush is a layer of polymer chains end-grafted to a surface.
#[derive(Debug, Clone)]
pub struct PolymerBrush {
    /// Grafting density σ (chains/area).
    pub grafting_density: f64,
    /// Degree of polymerization N.
    pub n: usize,
    /// Statistical segment length b.
    pub b: f64,
    /// Thermal energy kT.
    pub kt: f64,
    /// Excluded volume parameter v.
    pub v_excluded: f64,
}

impl PolymerBrush {
    /// Create a new polymer brush.
    pub fn new(n: usize, b: f64, kt: f64, grafting_density: f64, v_excluded: f64) -> Self {
        PolymerBrush {
            grafting_density,
            n,
            b,
            kt,
            v_excluded,
        }
    }

    /// Brush height H (Alexander-de Gennes scaling).
    ///
    /// H ~ N * b * (σ b²)^(1/3).
    pub fn brush_height(&self) -> f64 {
        let n = self.n as f64;
        n * self.b * (self.grafting_density * self.b * self.b).powf(1.0 / 3.0)
    }

    /// Osmotic pressure in the brush (de Gennes).
    ///
    /// Π ~ kT * σ^(9/4) / b^3 at height z.
    pub fn osmotic_pressure(&self, _z: f64) -> f64 {
        let s = self.grafting_density;
        self.kt * s.powf(9.0 / 4.0) / (self.b * self.b * self.b)
    }

    /// Free energy per chain in brush (stretching + excluded volume).
    pub fn free_energy_per_chain(&self) -> f64 {
        let n = self.n as f64;
        let h = self.brush_height();
        let stretch = 3.0 * h * h / (2.0 * n * self.b * self.b);
        let ev = self.v_excluded * n * n / h;
        (stretch + ev) * self.kt
    }

    /// Compression force per chain when brush is compressed to height d.
    pub fn compression_force(&self, d: f64) -> f64 {
        let h = self.brush_height();
        if d >= h {
            return 0.0;
        }
        let x = d / h;
        // Alexander-de Gennes: f ~ kT/b * (H/d)^(9/4) - (d/H)^(3/4)
        self.kt / self.b * ((h / d).powf(5.0 / 4.0) - x.powf(7.0 / 4.0))
    }

    /// Overlap density threshold σ* = 1 / (π * Rf²) where Rf ~ b * N^(3/5).
    pub fn overlap_density(&self) -> f64 {
        let n = self.n as f64;
        let rf = self.b * n.powf(0.6);
        1.0 / (PI * rf * rf)
    }

    /// Is the brush in the stretched regime (σ > σ*)?
    pub fn is_stretched(&self) -> bool {
        self.grafting_density > self.overlap_density()
    }
}

/// Manning condensation model for polyelectrolyte chains.
///
/// Counterions condense onto the chain backbone when the Manning parameter
/// ξ = l_B / b > 1/|z|, reducing the effective charge.
#[derive(Debug, Clone)]
pub struct Polyelectrolyte {
    /// Bjerrum length l_B (nm), where electrostatic and thermal energies equal.
    pub bjerrum_length: f64,
    /// Monomer spacing b (nm).
    pub b: f64,
    /// Monomer charge valence z.
    pub z: i32,
    /// Counterion valence z_c.
    pub z_c: i32,
    /// Number of monomers N.
    pub n: usize,
    /// Salt concentration c_s (M).
    pub c_salt: f64,
    /// Temperature in kT units.
    pub kt: f64,
}

impl Polyelectrolyte {
    /// Create a new polyelectrolyte model.
    pub fn new(
        n: usize,
        b: f64,
        bjerrum_length: f64,
        z: i32,
        z_c: i32,
        c_salt: f64,
        kt: f64,
    ) -> Self {
        Polyelectrolyte {
            bjerrum_length,
            b,
            z,
            z_c,
            n,
            c_salt,
            kt,
        }
    }

    /// Manning parameter ξ = l_B / b.
    pub fn manning_parameter(&self) -> f64 {
        self.bjerrum_length / self.b
    }

    /// Critical Manning parameter for condensation: ξ_c = 1 / |z * z_c|.
    pub fn critical_manning_parameter(&self) -> f64 {
        let zz = (self.z.abs() * self.z_c.abs()) as f64;
        1.0 / zz
    }

    /// Fraction of condensed counterions θ.
    ///
    /// θ = 1 - 1/(|z| * ξ) when ξ > ξ_c, else 0.
    pub fn condensed_fraction(&self) -> f64 {
        let xi = self.manning_parameter();
        let xi_c = self.critical_manning_parameter();
        if xi <= xi_c {
            0.0
        } else {
            1.0 - 1.0 / (self.z.abs() as f64 * xi)
        }
    }

    /// Effective linear charge density after condensation.
    pub fn effective_charge_density(&self) -> f64 {
        let theta = self.condensed_fraction();
        let z = self.z.abs() as f64;
        z * (1.0 - theta) / self.b
    }

    /// Debye screening length κ⁻¹ = 1 / sqrt(8π l_B N_A c_s).
    ///
    /// Here c_s is in mol/L, N_A = 6.022e23, result in nm.
    pub fn debye_length(&self) -> f64 {
        let c_m3 = self.c_salt * 6.022e23 * 1000.0; // convert mol/L to 1/m³
        // In nm units: κ⁻¹ = 1/sqrt(8π l_B_m c_m3)
        let l_b_m = self.bjerrum_length * 1e-9;
        (8.0 * PI * l_b_m * c_m3).sqrt().recip() * 1e9 // back to nm
    }

    /// Electrostatic persistence length l_e (OSF theory).
    ///
    /// l_e = l_B / (4 κ² b²) in nm.
    pub fn electrostatic_persistence_length(&self) -> f64 {
        let kappa = 1.0 / self.debye_length();
        self.bjerrum_length / (4.0 * kappa * kappa * self.b * self.b)
    }

    /// Total persistence length l_p = l_p0 + l_e.
    pub fn total_persistence_length(&self, l_p0: f64) -> f64 {
        l_p0 + self.electrostatic_persistence_length()
    }
}

/// Debye-Hückel screened interaction between charges in solution.
///
/// U(r) = (z₁ z₂ e²) / (4πε₀ε_r) * exp(-κr) / r = kT * l_B * z₁ z₂ * exp(-κr) / r
#[derive(Debug, Clone)]
pub struct DebyeHuckel {
    /// Bjerrum length l_B (nm).
    pub bjerrum_length: f64,
    /// Inverse Debye length κ (1/nm).
    pub kappa: f64,
}

impl DebyeHuckel {
    /// Create a new Debye-Hückel model.
    pub fn new(bjerrum_length: f64, kappa: f64) -> Self {
        DebyeHuckel {
            bjerrum_length,
            kappa,
        }
    }

    /// Screened Coulomb potential U(r) in kT units.
    pub fn potential(&self, r: f64, z1: f64, z2: f64) -> f64 {
        if r < 1e-15 {
            return f64::INFINITY;
        }
        self.bjerrum_length * z1 * z2 * (-self.kappa * r).exp() / r
    }

    /// Force magnitude F(r) = -dU/dr.
    pub fn force_magnitude(&self, r: f64, z1: f64, z2: f64) -> f64 {
        if r < 1e-15 {
            return f64::INFINITY;
        }
        self.bjerrum_length * z1 * z2 * (-self.kappa * r).exp() / (r * r) * (1.0 + self.kappa * r)
    }

    /// Second virial coefficient B₂ from DH potential (numerical integration).
    pub fn second_virial(&self, n_points: usize) -> f64 {
        let r_max = 10.0 / self.kappa.max(0.1);
        let dr = r_max / n_points as f64;
        let mut b2 = 0.0f64;
        for i in 1..=n_points {
            let r = i as f64 * dr;
            let mayer_f = (-self.potential(r, 1.0, 1.0)).exp() - 1.0;
            b2 -= 2.0 * PI * r * r * mayer_f * dr;
        }
        b2
    }
}
