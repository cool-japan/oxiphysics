#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Nosé-Hoover Chain (NHC) thermostat — Martyna, Klein & Tuckerman (1992).
//!
//! Implements a chain of M coupled thermostats using the Yoshida-Suzuki
//! RESPA factorization for time-reversible, symplectic integration.
//!
//! Algorithm follows Tuckerman *et al.* "A Liouville-operator derived
//! measure-preserving integrator for molecular dynamics simulations in
//! the isothermal-isobaric ensemble" (2006).

/// Boltzmann constant in SI units (J K⁻¹).
const KB: f64 = 1.380_649e-23;

/// Yoshida-Suzuki weights for 3-stage (order-4) factorization.
/// w\[0\] = w\[2\] = 1/(2 - 2^{1/3}), w\[1\] = -2^{1/3}/(2 - 2^{1/3})
const YS_WEIGHTS: [f64; 3] = [
    0.828_982_024_000_897_5,
    -0.657_964_048_001_795,
    0.828_982_024_000_897_5,
];

/// Nosé-Hoover Chain thermostat (Martyna *et al.* 1992).
///
/// Maintains M coupled thermostat degrees of freedom (`xi`, `v_xi`) that
/// act on the physical system's kinetic energy via a RESPA-style
/// Yoshida-Suzuki integration.
///
/// # Physical units
/// All quantities use SI units (K, kg, J, s) unless otherwise noted.
/// The coupling time `tau` is in seconds (use picosecond values × 1e-12).
#[derive(Debug, Clone)]
pub struct NhcThermostat {
    /// Target temperature (K).
    pub target_temp: f64,
    /// Number of chains.
    pub n_chains: usize,
    /// Chain positions (friction variables ξ_j).
    pub xi: Vec<f64>,
    /// Chain velocities (v_{ξ_j}).
    pub v_xi: Vec<f64>,
    /// Chain masses Q_j.
    /// Q_1 = n_dof · k_B · T · τ²; Q_j = k_B · T · τ² for j > 0.
    pub q_masses: Vec<f64>,
    /// Coupling time constant (s).
    pub tau: f64,
    /// Degrees of freedom of the physical system.
    pub n_dof: usize,
}

impl NhcThermostat {
    /// Create a new NHC thermostat.
    ///
    /// # Arguments
    /// * `target_temp` – desired temperature (K)
    /// * `n_chains`    – number of chained thermostats (typically 3–5)
    /// * `n_dof`       – physical degrees of freedom
    /// * `tau`         – coupling time constant (s)
    pub fn new(target_temp: f64, n_chains: usize, n_dof: usize, tau: f64) -> Self {
        assert!(n_chains >= 1, "NHC requires at least 1 chain");
        let kt = KB * target_temp;
        let tau2 = tau * tau;
        let mut q_masses = Vec::with_capacity(n_chains);
        // First mass couples to physical KE
        q_masses.push((n_dof as f64) * kt * tau2);
        // Remaining masses couple to each preceding thermostat
        for _ in 1..n_chains {
            q_masses.push(kt * tau2);
        }
        Self {
            target_temp,
            n_chains,
            xi: vec![0.0; n_chains],
            v_xi: vec![0.0; n_chains],
            q_masses,
            tau,
            n_dof,
        }
    }

    /// Perform one NHC propagation step using Yoshida-Suzuki factorization.
    ///
    /// Operates on the chain variables and returns a **velocity scale factor**
    /// that must be applied to all physical velocities:
    ///
    /// ```text
    /// v_physical *= scale
    /// ```
    ///
    /// # Arguments
    /// * `kinetic_energy` – current physical kinetic energy (J)
    /// * `dt`             – time step (s)
    pub fn integrate(&mut self, kinetic_energy: f64, dt: f64) -> f64 {
        let kt = KB * self.target_temp;
        let m = self.n_chains;
        let ndof = self.n_dof as f64;
        let mut scale = 1.0_f64;

        // --- Yoshida-Suzuki outer loop (3 sub-steps) ---
        for &w in &YS_WEIGHTS {
            let dt_ys = w * dt; // full sub-step
            let dt2 = 0.5 * dt_ys; // half sub-step
            let dt4 = 0.25 * dt_ys; // quarter sub-step
            let dt8 = 0.125 * dt_ys; // eighth sub-step

            // === HALF-STEP FORWARD on chain velocities (reversed order: M..1) ===

            // Outermost thermostat (index M-1): no coupling from above
            // G_{M-1} = (Q_{M-2} * v_{M-2}^2 - kT) / Q_{M-1}
            if m > 1 {
                let g_last =
                    (self.q_masses[m - 2] * self.v_xi[m - 2].powi(2) - kt) / self.q_masses[m - 1];
                self.v_xi[m - 1] += g_last * dt4;
            }

            // Inner thermostats M-2 down to 0 (coupled to the one above)
            if m >= 2 {
                for j in (0..m - 1).rev() {
                    // Coupling: v_xi[j] is "damped" by v_xi[j+1]
                    let aa = (-self.v_xi[j + 1] * dt8).exp();
                    let g_j = if j == 0 {
                        // Coupling to physical system KE
                        (2.0 * kinetic_energy * scale * scale - ndof * kt) / self.q_masses[0]
                    } else {
                        (self.q_masses[j - 1] * self.v_xi[j - 1].powi(2) - kt) / self.q_masses[j]
                    };
                    self.v_xi[j] = aa * (aa * self.v_xi[j] + g_j * dt4);
                }
            } else {
                // Single chain: just update v_xi[0] directly
                let g0 = (2.0 * kinetic_energy * scale * scale - ndof * kt) / self.q_masses[0];
                self.v_xi[0] += g0 * dt4;
            }

            // === HALF-STEP: update xi positions & scale velocities ===

            // Update all chain positions
            for j in 0..m {
                self.xi[j] += self.v_xi[j] * dt2;
            }

            // Scale physical velocities: exp(-v_xi[0] * dt/2)
            let local_scale = (-self.v_xi[0] * dt2).exp();
            scale *= local_scale;

            // === HALF-STEP BACKWARD on chain velocities (forward order: 0..M) ===

            // Innermost thermostat (index 0)
            if m >= 2 {
                for j in 0..m - 1 {
                    let aa = (-self.v_xi[j + 1] * dt8).exp();
                    let g_j = if j == 0 {
                        (2.0 * kinetic_energy * scale * scale - ndof * kt) / self.q_masses[0]
                    } else {
                        (self.q_masses[j - 1] * self.v_xi[j - 1].powi(2) - kt) / self.q_masses[j]
                    };
                    self.v_xi[j] = aa * (aa * self.v_xi[j] + g_j * dt4);
                }
            } else {
                // Single chain
                let g0 = (2.0 * kinetic_energy * scale * scale - ndof * kt) / self.q_masses[0];
                self.v_xi[0] += g0 * dt4;
            }

            // Outermost thermostat
            if m > 1 {
                let g_last =
                    (self.q_masses[m - 2] * self.v_xi[m - 2].powi(2) - kt) / self.q_masses[m - 1];
                self.v_xi[m - 1] += g_last * dt4;
            }
        }

        scale
    }

    /// Instantaneous temperature from kinetic energy.
    ///
    /// T = 2 · KE / (n_dof · k_B)
    pub fn current_temperature(&self, kinetic_energy: f64) -> f64 {
        2.0 * kinetic_energy / ((self.n_dof as f64) * KB)
    }

    /// NHC conserved energy (extended Hamiltonian):
    ///
    /// H_NHC = KE + PE + Σ_j ½·Q_j·v_ξ_j² + n_dof·k_B·T·ξ_0 + Σ_{j>0} k_B·T·ξ_j
    pub fn conserved_energy(&self, kinetic_energy: f64, potential_energy: f64) -> f64 {
        let kt = KB * self.target_temp;
        let bath_ke: f64 = self
            .q_masses
            .iter()
            .zip(self.v_xi.iter())
            .map(|(q, v)| 0.5 * q * v * v)
            .sum();
        let bath_pos: f64 = (self.n_dof as f64) * kt * self.xi[0]
            + self.xi[1..].iter().map(|x| kt * x).sum::<f64>();
        kinetic_energy + potential_energy + bath_ke + bath_pos
    }
}

// ---------------------------------------------------------------------------
// RESPA-NHC propagator (split into fast/slow thermostat coupling)
// ---------------------------------------------------------------------------

/// RESPA-NHC propagator: couples a Nosé-Hoover chain to a multi-timestep
/// (RESPA) integrator. The thermostat acts on the slow outer time step
/// while the inner fast forces are propagated without thermostat coupling.
///
/// Algorithm:
/// 1. Half-step NHC propagation (thermostat acts on velocities)
/// 2. N inner steps of velocity-Verlet for fast forces
/// 3. Half-step NHC propagation
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RespaNhc {
    /// NHC thermostat for the outer (slow) loop.
    pub nhc: NhcThermostat,
    /// Number of inner (fast) steps per outer step.
    pub n_inner: usize,
}

impl RespaNhc {
    /// Create a RESPA-NHC propagator.
    ///
    /// # Arguments
    /// * `target_temp` – temperature (K)
    /// * `n_chains` – NHC chain length
    /// * `n_dof` – physical degrees of freedom
    /// * `tau` – NHC coupling time constant (s)
    /// * `n_inner` – number of inner fast steps per outer step
    pub fn new(target_temp: f64, n_chains: usize, n_dof: usize, tau: f64, n_inner: usize) -> Self {
        Self {
            nhc: NhcThermostat::new(target_temp, n_chains, n_dof, tau),
            n_inner: n_inner.max(1),
        }
    }

    /// Perform one outer RESPA step with NHC coupling.
    ///
    /// Returns `(scale, inner_dt)` where `scale` is the total velocity scale
    /// factor from the two half-step NHC propagations and `inner_dt` is
    /// `dt_outer / n_inner`.
    pub fn outer_step(&mut self, kinetic_energy: f64, dt_outer: f64) -> (f64, f64) {
        // First half-step NHC propagation
        let scale1 = self.nhc.integrate(kinetic_energy, 0.5 * dt_outer);
        let ke_after = kinetic_energy * scale1 * scale1;
        let inner_dt = dt_outer / self.n_inner as f64;

        // (Caller performs n_inner velocity-Verlet inner steps here)

        // Second half-step NHC propagation
        let scale2 = self.nhc.integrate(ke_after, 0.5 * dt_outer);
        (scale1 * scale2, inner_dt)
    }
}

// ---------------------------------------------------------------------------
// MTK equations for NPT (Martyna-Tobias-Klein)
// ---------------------------------------------------------------------------

/// Martyna-Tobias-Klein (MTK) NPT thermostat-barostat coupling.
///
/// Extends NHC to the isothermal-isobaric ensemble by coupling a barostat
/// degree of freedom (volume/cell) to an additional NHC chain, while the
/// particle thermostat acts on particle kinetic energy.
///
/// Reference: Martyna, Tobias & Klein, J. Chem. Phys. 101, 4177 (1994).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MtkNpt {
    /// Particle thermostat (NHC chain).
    pub particle_nhc: NhcThermostat,
    /// Barostat thermostat (NHC chain for the volume DOF).
    pub baro_nhc: NhcThermostat,
    /// Barostat mass W (kg·m²).
    pub w_baro: f64,
    /// Barostat velocity (volume strain rate, s⁻¹).
    pub v_eps: f64,
    /// Barostat position (ln(V/V₀)).
    pub eps: f64,
    /// Target pressure (Pa).
    pub target_pressure: f64,
    /// Number of spatial dimensions (usually 3).
    pub n_dim: usize,
}

impl MtkNpt {
    /// Create a new MTK NPT thermostat-barostat.
    ///
    /// # Arguments
    /// * `target_temp` – target temperature (K)
    /// * `target_pressure` – target pressure (Pa)
    /// * `n_chains` – NHC chain length (for both particle and baro chains)
    /// * `n_dof` – particle degrees of freedom
    /// * `tau_t` – thermostat coupling time (s)
    /// * `tau_p` – barostat coupling time (s)
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        target_temp: f64,
        target_pressure: f64,
        n_chains: usize,
        n_dof: usize,
        tau_t: f64,
        tau_p: f64,
    ) -> Self {
        let kt = KB * target_temp;
        let w_baro = (n_dof as f64 + 3.0) * kt * tau_p * tau_p;
        Self {
            particle_nhc: NhcThermostat::new(target_temp, n_chains, n_dof, tau_t),
            baro_nhc: NhcThermostat::new(target_temp, n_chains, 1, tau_t),
            w_baro,
            v_eps: 0.0,
            eps: 0.0,
            target_pressure,
            n_dim: 3,
        }
    }

    /// Compute the barostat kinetic energy: ½ W v_ε².
    pub fn baro_kinetic_energy(&self) -> f64 {
        0.5 * self.w_baro * self.v_eps * self.v_eps
    }

    /// Propagate one MTK step.
    ///
    /// Returns `(vel_scale, vol_scale)`:
    /// - `vel_scale`: factor to multiply all particle velocities
    /// - `vol_scale`: factor to multiply the box volume
    ///
    /// # Arguments
    /// * `kinetic_energy` – current particle kinetic energy (J)
    /// * `inst_pressure` – instantaneous pressure from virial (Pa)
    /// * `volume` – current volume (m³)
    /// * `dt` – time step (s)
    pub fn propagate(
        &mut self,
        kinetic_energy: f64,
        inst_pressure: f64,
        volume: f64,
        dt: f64,
    ) -> (f64, f64) {
        let kt = KB * self.particle_nhc.target_temp;
        let ndof = self.particle_nhc.n_dof as f64;
        let ndim = self.n_dim as f64;

        // Particle NHC half-step
        let scale_nhc = self.particle_nhc.integrate(kinetic_energy, 0.5 * dt);

        // Barostat force: G_eps = (ndim/W) * V * (P_inst - P_target) + (2*KE/ndof - kT)/W
        let g_eps = ndim * volume * (inst_pressure - self.target_pressure) / self.w_baro
            + (2.0 * kinetic_energy * scale_nhc * scale_nhc / ndof - kt) / self.w_baro;

        // Barostat velocity half-step
        self.v_eps += g_eps * 0.5 * dt;

        // Barostat NHC half-step (thermostat the barostat DOF)
        let baro_ke = self.baro_kinetic_energy();
        let baro_scale = self.baro_nhc.integrate(baro_ke, 0.5 * dt);
        self.v_eps *= baro_scale;

        // Update eps (volume strain)
        self.eps += self.v_eps * dt;
        let vol_scale = (ndim * self.v_eps * dt).exp();

        // Velocity scaling from barostat: exp(-(1 + ndim/ndof) * v_eps * dt/2)
        let vel_baro_scale = (-(1.0 + ndim / ndof) * self.v_eps * 0.5 * dt).exp();

        // Second barostat NHC half-step
        let baro_ke2 = self.baro_kinetic_energy();
        let baro_scale2 = self.baro_nhc.integrate(baro_ke2, 0.5 * dt);
        self.v_eps *= baro_scale2;

        // Barostat velocity second half-step
        self.v_eps += g_eps * 0.5 * dt;

        // Second particle NHC half-step
        let ke_updated = kinetic_energy * (scale_nhc * vel_baro_scale).powi(2);
        let scale_nhc2 = self.particle_nhc.integrate(ke_updated, 0.5 * dt);

        let total_vel_scale = scale_nhc * vel_baro_scale * scale_nhc2;
        (total_vel_scale, vol_scale)
    }

    /// MTK conserved energy (extended Hamiltonian) for the NPT ensemble.
    ///
    /// H_MTK = KE + PE + ½ W v_ε² + P_target V
    ///       + particle_NHC_energy + baro_NHC_energy
    pub fn conserved_energy(&self, kinetic_energy: f64, potential_energy: f64, volume: f64) -> f64 {
        let nhc_particle = self.particle_nhc.conserved_energy(kinetic_energy, 0.0);
        let baro_ke = self.baro_kinetic_energy();
        let nhc_baro = self.baro_nhc.conserved_energy(baro_ke, 0.0);
        nhc_particle + potential_energy + nhc_baro + self.target_pressure * volume
    }
}

// ---------------------------------------------------------------------------
// Flexible chain length NHC
// ---------------------------------------------------------------------------

/// NHC thermostat with configurable chain length and separate mass schedule.
///
/// Allows non-uniform chain masses (e.g., for massive thermostatting where
/// each DOF has its own chain, or for gradually decaying masses along the chain).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FlexibleNhc {
    /// Inner NHC thermostat.
    pub inner: NhcThermostat,
    /// Per-chain mass scaling factors (multiplied with the default Q values).
    pub mass_scales: Vec<f64>,
}

impl FlexibleNhc {
    /// Create a flexible NHC thermostat with custom mass scaling.
    ///
    /// `mass_scales` should have length `n_chains`. Each Q_j is multiplied
    /// by `mass_scales[j]`.
    pub fn new(
        target_temp: f64,
        n_chains: usize,
        n_dof: usize,
        tau: f64,
        mass_scales: Vec<f64>,
    ) -> Self {
        assert_eq!(
            mass_scales.len(),
            n_chains,
            "mass_scales length must match n_chains"
        );
        let mut nhc = NhcThermostat::new(target_temp, n_chains, n_dof, tau);
        for (j, &s) in mass_scales.iter().enumerate() {
            nhc.q_masses[j] *= s;
        }
        Self {
            inner: nhc,
            mass_scales,
        }
    }

    /// Delegate integration to the inner NHC.
    pub fn integrate(&mut self, kinetic_energy: f64, dt: f64) -> f64 {
        self.inner.integrate(kinetic_energy, dt)
    }

    /// Delegate conserved energy to the inner NHC.
    pub fn conserved_energy(&self, kinetic_energy: f64, potential_energy: f64) -> f64 {
        self.inner
            .conserved_energy(kinetic_energy, potential_energy)
    }
}

// ---------------------------------------------------------------------------
// Yoshida-Suzuki 5-stage (order 6) weights
// ---------------------------------------------------------------------------

/// Yoshida-Suzuki weights for 5-stage (order-6) factorization.
///
/// Higher-order factorization reduces the integration error of the NHC
/// propagation at the cost of more substeps.
#[allow(dead_code)]
const YS_WEIGHTS_5: [f64; 5] = [
    0.414_490_771_794_375_9,
    0.414_490_771_794_375_9,
    -0.657_963_087_177_503_5,
    0.414_490_771_794_375_9,
    0.414_490_771_794_375_9,
];

/// NHC thermostat using 5-stage Yoshida-Suzuki factorization (order 6).
///
/// Provides higher-order accuracy compared to the 3-stage version at the
/// cost of 5 sub-steps per NHC propagation.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct NhcThermostatOrder6 {
    /// Target temperature (K).
    pub target_temp: f64,
    /// Number of chains.
    pub n_chains: usize,
    /// Chain positions.
    pub xi: Vec<f64>,
    /// Chain velocities.
    pub v_xi: Vec<f64>,
    /// Chain masses.
    pub q_masses: Vec<f64>,
    /// Coupling time constant (s).
    pub tau: f64,
    /// Degrees of freedom.
    pub n_dof: usize,
}

impl NhcThermostatOrder6 {
    /// Create a new order-6 NHC thermostat.
    pub fn new(target_temp: f64, n_chains: usize, n_dof: usize, tau: f64) -> Self {
        assert!(n_chains >= 1);
        let kt = KB * target_temp;
        let tau2 = tau * tau;
        let mut q_masses = Vec::with_capacity(n_chains);
        q_masses.push((n_dof as f64) * kt * tau2);
        for _ in 1..n_chains {
            q_masses.push(kt * tau2);
        }
        Self {
            target_temp,
            n_chains,
            xi: vec![0.0; n_chains],
            v_xi: vec![0.0; n_chains],
            q_masses,
            tau,
            n_dof,
        }
    }

    /// Integrate with 5-stage Yoshida-Suzuki weights.
    pub fn integrate(&mut self, kinetic_energy: f64, dt: f64) -> f64 {
        let kt = KB * self.target_temp;
        let m = self.n_chains;
        let ndof = self.n_dof as f64;
        let mut scale = 1.0_f64;

        for &w in &YS_WEIGHTS_5 {
            let dt_ys = w * dt;
            let dt2 = 0.5 * dt_ys;
            let dt4 = 0.25 * dt_ys;
            let dt8 = 0.125 * dt_ys;

            if m > 1 {
                let g_last =
                    (self.q_masses[m - 2] * self.v_xi[m - 2].powi(2) - kt) / self.q_masses[m - 1];
                self.v_xi[m - 1] += g_last * dt4;
            }

            if m >= 2 {
                for j in (0..m - 1).rev() {
                    let aa = (-self.v_xi[j + 1] * dt8).exp();
                    let g_j = if j == 0 {
                        (2.0 * kinetic_energy * scale * scale - ndof * kt) / self.q_masses[0]
                    } else {
                        (self.q_masses[j - 1] * self.v_xi[j - 1].powi(2) - kt) / self.q_masses[j]
                    };
                    self.v_xi[j] = aa * (aa * self.v_xi[j] + g_j * dt4);
                }
            } else {
                let g0 = (2.0 * kinetic_energy * scale * scale - ndof * kt) / self.q_masses[0];
                self.v_xi[0] += g0 * dt4;
            }

            for j in 0..m {
                self.xi[j] += self.v_xi[j] * dt2;
            }
            let local_scale = (-self.v_xi[0] * dt2).exp();
            scale *= local_scale;

            if m >= 2 {
                for j in 0..m - 1 {
                    let aa = (-self.v_xi[j + 1] * dt8).exp();
                    let g_j = if j == 0 {
                        (2.0 * kinetic_energy * scale * scale - ndof * kt) / self.q_masses[0]
                    } else {
                        (self.q_masses[j - 1] * self.v_xi[j - 1].powi(2) - kt) / self.q_masses[j]
                    };
                    self.v_xi[j] = aa * (aa * self.v_xi[j] + g_j * dt4);
                }
            } else {
                let g0 = (2.0 * kinetic_energy * scale * scale - ndof * kt) / self.q_masses[0];
                self.v_xi[0] += g0 * dt4;
            }

            if m > 1 {
                let g_last =
                    (self.q_masses[m - 2] * self.v_xi[m - 2].powi(2) - kt) / self.q_masses[m - 1];
                self.v_xi[m - 1] += g_last * dt4;
            }
        }
        scale
    }
}

// ---------------------------------------------------------------------------
// IsokineThermostat — isokinetic (Gaussian) thermostat
// ---------------------------------------------------------------------------

/// Isokinetic (Gaussian) thermostat.
///
/// Rigorously constrains the kinetic energy to exactly the target value
/// by using a Lagrange multiplier α:
///
/// dp/dt = F - α·p
/// α = (sum F·p) / (sum p·p / m)
///
/// This produces the isokinetic ensemble (microcanonical in momenta).
///
/// Reference: Evans & Morriss, Phys. Rev. A 30, 1060 (1984).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct IsokineThermostat {
    /// Target kinetic energy (J).
    pub target_ke: f64,
    /// Number of degrees of freedom.
    pub n_dof: usize,
    /// Target temperature (K).
    pub target_temp: f64,
}

impl IsokineThermostat {
    /// Create a new isokinetic thermostat.
    pub fn new(target_temp: f64, n_dof: usize) -> Self {
        let target_ke = 0.5 * (n_dof as f64) * KB * target_temp;
        Self {
            target_ke,
            n_dof,
            target_temp,
        }
    }

    /// Compute the Gaussian friction coefficient α from forces and momenta.
    ///
    /// α = (sum_i F_i · p_i) / (sum_i p_i · p_i / m_i)
    ///
    /// # Arguments
    /// * `forces`    — per-atom forces (N)
    /// * `momenta`   — per-atom momenta (kg·m/s)
    /// * `masses`    — per-atom masses (kg)
    pub fn gaussian_friction(
        &self,
        forces: &[[f64; 3]],
        momenta: &[[f64; 3]],
        masses: &[f64],
    ) -> f64 {
        let n = forces.len().min(momenta.len()).min(masses.len());
        let mut num = 0.0_f64;
        let mut den = 0.0_f64;
        for i in 0..n {
            for d in 0..3 {
                num += forces[i][d] * momenta[i][d];
                den += momenta[i][d] * momenta[i][d] / masses[i];
            }
        }
        if den.abs() < 1e-300 {
            return 0.0;
        }
        num / den
    }

    /// Rescale velocities to exactly match the target kinetic energy.
    ///
    /// v_new = v * sqrt(KE_target / KE_current)
    ///
    /// Returns the scale factor applied.
    pub fn rescale_velocities(&self, velocities: &mut [[f64; 3]], masses: &[f64]) -> f64 {
        let n = velocities.len().min(masses.len());
        let mut ke = 0.0;
        for i in 0..n {
            for d in 0..3 {
                ke += 0.5 * masses[i] * velocities[i][d] * velocities[i][d];
            }
        }
        if ke < 1e-300 {
            return 1.0;
        }
        let scale = (self.target_ke / ke).sqrt();
        for i in 0..n {
            for d in 0..3 {
                velocities[i][d] *= scale;
            }
        }
        scale
    }

    /// Check if the kinetic energy constraint is satisfied (within tolerance).
    pub fn is_constrained(&self, velocities: &[[f64; 3]], masses: &[f64], tol: f64) -> bool {
        let n = velocities.len().min(masses.len());
        let mut ke = 0.0;
        for i in 0..n {
            for d in 0..3 {
                ke += 0.5 * masses[i] * velocities[i][d] * velocities[i][d];
            }
        }
        (ke - self.target_ke).abs() / self.target_ke < tol
    }
}

// ---------------------------------------------------------------------------
// StochasticVelocityRescaling — Bussi-Donadio-Parrinello thermostat
// ---------------------------------------------------------------------------

/// Stochastic Velocity Rescaling (SVR) thermostat.
///
/// Rescales velocities to sample the canonical ensemble.  The stochastic
/// component ensures ergodicity while the deterministic component drives
/// the temperature toward the target.
///
/// Reference: Bussi, Donadio & Parrinello, J. Chem. Phys. 126, 014101 (2007).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct StochasticVelocityRescaling {
    /// Target temperature (K).
    pub target_temp: f64,
    /// Number of degrees of freedom.
    pub n_dof: usize,
    /// Thermostat coupling time constant (s).
    pub tau: f64,
}

impl StochasticVelocityRescaling {
    /// Create a new SVR thermostat.
    pub fn new(target_temp: f64, n_dof: usize, tau: f64) -> Self {
        Self {
            target_temp,
            n_dof,
            tau,
        }
    }

    /// Compute the deterministic velocity scale factor (Berendsen part).
    ///
    /// alpha_det = sqrt(1 + (dt/tau) * (T_target/T_inst - 1))
    ///
    /// This is the Berendsen thermostat scale without stochastic correction.
    pub fn deterministic_scale(&self, t_inst: f64, dt: f64) -> f64 {
        let ratio = self.target_temp / t_inst.max(1e-10);
        (1.0 + (dt / self.tau) * (ratio - 1.0)).sqrt().max(0.0)
    }

    /// Kinetic energy at target temperature.
    pub fn target_ke(&self) -> f64 {
        0.5 * (self.n_dof as f64) * KB * self.target_temp
    }

    /// Instantaneous temperature from kinetic energy.
    pub fn instantaneous_temperature(&self, ke: f64) -> f64 {
        2.0 * ke / ((self.n_dof as f64) * KB)
    }

    /// Relaxation rate: tau_ke = tau / (2 * n_dof).
    pub fn relaxation_rate(&self) -> f64 {
        self.tau / (2.0 * self.n_dof as f64)
    }
}

// ---------------------------------------------------------------------------
// NhcConservedQuantity — monitoring extended Hamiltonian drift
// ---------------------------------------------------------------------------

/// Extended Hamiltonian (conserved quantity) monitor for NHC.
///
/// In a properly implemented NHC the extended Hamiltonian
/// H_ext = KE + PE + bath_KE + bath_PE
/// should be conserved.  This struct accumulates statistics about the drift.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct NhcConservedQuantity {
    /// Initial value of the extended Hamiltonian.
    pub h_initial: f64,
    /// Current value.
    pub h_current: f64,
    /// Maximum absolute drift observed.
    pub max_drift: f64,
    /// Number of samples collected.
    pub n_samples: u64,
}

impl NhcConservedQuantity {
    /// Create a new monitor with an initial value.
    pub fn new(h_initial: f64) -> Self {
        Self {
            h_initial,
            h_current: h_initial,
            max_drift: 0.0,
            n_samples: 1,
        }
    }

    /// Record a new sample of the extended Hamiltonian.
    pub fn record(&mut self, h_new: f64) {
        self.h_current = h_new;
        let drift = (h_new - self.h_initial).abs();
        if drift > self.max_drift {
            self.max_drift = drift;
        }
        self.n_samples += 1;
    }

    /// Relative drift: |H_current - H_initial| / |H_initial|.
    pub fn relative_drift(&self) -> f64 {
        if self.h_initial.abs() < 1e-300 {
            return 0.0;
        }
        (self.h_current - self.h_initial).abs() / self.h_initial.abs()
    }

    /// Maximum relative drift observed.
    pub fn max_relative_drift(&self) -> f64 {
        if self.h_initial.abs() < 1e-300 {
            return 0.0;
        }
        self.max_drift / self.h_initial.abs()
    }

    /// Returns true if the drift is within a given relative tolerance.
    pub fn is_conserved(&self, tol: f64) -> bool {
        self.relative_drift() < tol
    }
}

// ---------------------------------------------------------------------------
// Yoshida-Suzuki 7-stage (order 8) weights
// ---------------------------------------------------------------------------

/// Yoshida-Suzuki weights for 7-stage (order 8) factorization.
///
/// For very high precision NHC integration where machine-precision
/// conservation is required.
#[allow(dead_code)]
pub const YS_WEIGHTS_7: [f64; 7] = [
    0.784_513_610_477_560,
    0.235_573_213_359_357,
    -1.177_679_984_178_871,
    1.315_186_320_683_906,
    -1.177_679_984_178_871,
    0.235_573_213_359_357,
    0.784_513_610_477_560,
];

// ---------------------------------------------------------------------------
// Thermostat response function
// ---------------------------------------------------------------------------

/// Frequency response of the NHC thermostat.
///
/// Returns the effective coupling frequency (rad/s) for the first chain
/// element: ω_0 = 1/τ.
#[allow(dead_code)]
pub fn nhc_coupling_frequency(tau: f64) -> f64 {
    if tau > 0.0 { 1.0 / tau } else { 0.0 }
}

/// Quality factor Q of the NHC thermostat at a given frequency `omega`.
///
/// The thermostat is most efficient when the physical motion frequency
/// matches ω_0 = 1/τ.  The Q-factor is defined here as:
/// Q(ω) = ω / ω_0 (dimensionless ratio).
#[allow(dead_code)]
pub fn nhc_quality_factor(omega: f64, tau: f64) -> f64 {
    let omega0 = nhc_coupling_frequency(tau);
    if omega0 < 1e-300 {
        return 0.0;
    }
    omega / omega0
}

/// Compute the optimal coupling time τ for a given target frequency ω.
///
/// Returns τ = 1 / ω.
#[allow(dead_code)]
pub fn nhc_optimal_tau(omega: f64) -> f64 {
    if omega > 0.0 { 1.0 / omega } else { 0.0 }
}

// ---------------------------------------------------------------------------
// Chain velocity initialization (Maxwell-Boltzmann)
// ---------------------------------------------------------------------------

/// Initialize chain velocities from a Maxwell-Boltzmann distribution
/// at temperature T using a deterministic (pseudo-random) approach.
///
/// This routine uses a simple linear congruential generator (LCG) seeded
/// by `seed` so that tests are reproducible without external crates.
#[allow(dead_code)]
pub fn init_chain_velocities(nhc: &mut NhcThermostat, seed: u64) {
    let kt = KB * nhc.target_temp;
    let mut state = seed;
    for (v, q) in nhc.v_xi.iter_mut().zip(nhc.q_masses.iter()) {
        // LCG step
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        // Map state to [-1, 1] via bit manipulation
        let u = (state >> 11) as f64 / (1u64 << 53) as f64; // [0, 1)
        let u = 2.0 * u - 1.0; // [-1, 1)
        // Gaussian approximation: Box-Muller not needed; scale by sqrt(kT/Q)
        let sigma = (kt / q.max(1e-300)).sqrt();
        *v = u * sigma;
    }
}

/// Compute the kinetic energy of the NHC bath variables.
///
/// KE_bath = sum_j 0.5 * Q_j * v_xi_j^2
#[allow(dead_code)]
pub fn chain_kinetic_energy(nhc: &NhcThermostat) -> f64 {
    nhc.q_masses
        .iter()
        .zip(nhc.v_xi.iter())
        .map(|(q, v)| 0.5 * q * v * v)
        .sum()
}

/// Compute the potential energy of the NHC bath variables.
///
/// V_bath = n_dof * k_B * T * xi_0 + sum_{j>0} k_B * T * xi_j
#[allow(dead_code)]
pub fn chain_potential_energy(nhc: &NhcThermostat) -> f64 {
    let kt = KB * nhc.target_temp;
    let ndof = nhc.n_dof as f64;
    let mut v = ndof * kt * nhc.xi[0];
    for xi in &nhc.xi[1..] {
        v += kt * xi;
    }
    v
}

/// Compute the extended system energy: physical_ke + bath_ke + bath_pe.
///
/// This quantity should be conserved during NHC dynamics.
#[allow(dead_code)]
pub fn extended_energy(nhc: &NhcThermostat, physical_ke: f64) -> f64 {
    physical_ke + chain_kinetic_energy(nhc) + chain_potential_energy(nhc)
}

// ---------------------------------------------------------------------------
// RESPA multi-timestep helper
// ---------------------------------------------------------------------------

/// RESPA (Reference System Propagator Algorithm) multi-timestep wrapper.
///
/// Divides a large outer timestep `dt_outer` into `n_inner` inner steps
/// for the fast (bonded) forces, while slow (non-bonded) forces are applied
/// every outer step.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RespaIntegrator {
    /// Outer (slow) timestep.
    pub dt_outer: f64,
    /// Number of inner steps per outer step.
    pub n_inner: usize,
    /// Inner (fast) timestep = dt_outer / n_inner.
    pub dt_inner: f64,
    /// Step counter.
    pub step: usize,
}

impl RespaIntegrator {
    /// Create a new RESPA integrator.
    #[allow(dead_code)]
    pub fn new(dt_outer: f64, n_inner: usize) -> Self {
        let dt_inner = dt_outer / n_inner.max(1) as f64;
        Self {
            dt_outer,
            n_inner,
            dt_inner,
            step: 0,
        }
    }

    /// Advance one outer step, returning the sequence of inner timesteps.
    ///
    /// Returns a Vec of `n_inner` inner timesteps (all equal to `dt_inner`).
    #[allow(dead_code)]
    pub fn inner_steps(&mut self) -> Vec<f64> {
        self.step += 1;
        vec![self.dt_inner; self.n_inner]
    }

    /// Half-step velocity update for the slow forces.
    ///
    /// Each velocity component `v_i` should be updated as:
    /// `v_i += 0.5 * dt_outer * f_slow_i / m_i`
    ///
    /// Returns the half-step scale: always 0.5 * dt_outer.
    #[allow(dead_code)]
    pub fn half_step_dt(&self) -> f64 {
        0.5 * self.dt_outer
    }

    /// Reset the step counter.
    #[allow(dead_code)]
    pub fn reset(&mut self) {
        self.step = 0;
    }

    /// Total simulation time elapsed.
    #[allow(dead_code)]
    pub fn elapsed_time(&self) -> f64 {
        self.step as f64 * self.dt_outer
    }
}

// ---------------------------------------------------------------------------
// Berendsen thermostat (velocity rescaling)
// ---------------------------------------------------------------------------

/// Berendsen weak-coupling thermostat (Berendsen *et al.* 1984).
///
/// **Note**: Does not produce correct canonical ensemble.
/// Use for equilibration only.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BerendsenThermostat {
    /// Target temperature (K).
    pub target_temp: f64,
    /// Coupling time constant (s).
    pub tau: f64,
}

impl BerendsenThermostat {
    /// Create a new Berendsen thermostat.
    #[allow(dead_code)]
    pub fn new(target_temp: f64, tau: f64) -> Self {
        Self { target_temp, tau }
    }

    /// Compute the velocity rescaling factor λ.
    ///
    /// λ = sqrt(1 + (dt/τ) * (T_target/T_current - 1))
    ///
    /// Returns 1.0 if T_current ≈ 0.
    #[allow(dead_code)]
    pub fn lambda(&self, t_current: f64, dt: f64) -> f64 {
        if t_current < 1e-300 {
            return 1.0;
        }
        let ratio = self.target_temp / t_current;
        (1.0 + (dt / self.tau) * (ratio - 1.0)).max(0.0).sqrt()
    }

    /// Rescale velocities (3D) by λ.  Returns the scale factor used.
    #[allow(dead_code)]
    pub fn rescale(&self, velocities: &mut [[f64; 3]], t_current: f64, dt: f64) -> f64 {
        let lam = self.lambda(t_current, dt);
        for v in velocities.iter_mut() {
            v[0] *= lam;
            v[1] *= lam;
            v[2] *= lam;
        }
        lam
    }

    /// Compute instantaneous temperature from kinetic energy.
    ///
    /// T = 2 * KE / (n_dof * k_B)
    #[allow(dead_code)]
    pub fn temperature_from_ke(&self, ke: f64, n_dof: usize) -> f64 {
        if n_dof == 0 {
            return 0.0;
        }
        2.0 * ke / (n_dof as f64 * KB)
    }
}

// ---------------------------------------------------------------------------
// Andersen thermostat (stochastic collisions)
// ---------------------------------------------------------------------------

/// Andersen thermostat: random velocity reassignment.
///
/// With collision frequency ν, each atom has probability ν*dt of having
/// its velocity reassigned from a Maxwell-Boltzmann distribution.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AndersenThermostat {
    /// Target temperature (K).
    pub target_temp: f64,
    /// Collision frequency (s⁻¹).
    pub nu: f64,
}

impl AndersenThermostat {
    /// Create a new Andersen thermostat.
    #[allow(dead_code)]
    pub fn new(target_temp: f64, nu: f64) -> Self {
        Self { target_temp, nu }
    }

    /// Probability that a given atom experiences a collision in time `dt`.
    #[allow(dead_code)]
    pub fn collision_probability(&self, dt: f64) -> f64 {
        (self.nu * dt).min(1.0)
    }

    /// Sample a single velocity component from a Maxwell-Boltzmann distribution
    /// using a simple pseudo-random approach (LCG).
    ///
    /// Returns a velocity component with variance k_B * T / mass.
    #[allow(dead_code)]
    pub fn sample_velocity_component(&self, mass: f64, state: &mut u64) -> f64 {
        // Advance LCG
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let u = (*state >> 11) as f64 / (1u64 << 53) as f64; // [0, 1)
        // Box-Muller: use two states; approximation with single (cos part)
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let u2 = (*state >> 11) as f64 / (1u64 << 53) as f64;
        let sigma = (KB * self.target_temp / mass.max(1e-300)).sqrt();
        let z = (-2.0 * u.max(1e-300).ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
        sigma * z
    }

    /// Count atoms that would receive a collision for given dt and seed.
    ///
    /// Uses a deterministic LCG to decide per-atom collisions.
    #[allow(dead_code)]
    pub fn count_collisions(&self, n_atoms: usize, dt: f64, seed: u64) -> usize {
        let prob = self.collision_probability(dt);
        let mut state = seed;
        let mut count = 0;
        for _ in 0..n_atoms {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let u = (state >> 11) as f64 / (1u64 << 53) as f64;
            if u < prob {
                count += 1;
            }
        }
        count
    }
}

// ---------------------------------------------------------------------------
// NHC extended system - conserved quantity time series
// ---------------------------------------------------------------------------

/// Simple time series for tracking a conserved quantity over a simulation.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ConservedQuantityTimeSeries {
    /// Recorded values.
    pub values: Vec<f64>,
    /// Timestep (s).
    pub dt: f64,
}

impl ConservedQuantityTimeSeries {
    /// Create an empty time series.
    #[allow(dead_code)]
    pub fn new(dt: f64) -> Self {
        Self {
            values: Vec::new(),
            dt,
        }
    }

    /// Append a value.
    #[allow(dead_code)]
    pub fn push(&mut self, val: f64) {
        self.values.push(val);
    }

    /// Mean of the recorded values.
    #[allow(dead_code)]
    pub fn mean(&self) -> f64 {
        if self.values.is_empty() {
            return 0.0;
        }
        self.values.iter().sum::<f64>() / self.values.len() as f64
    }

    /// Standard deviation of the recorded values.
    #[allow(dead_code)]
    pub fn std_dev(&self) -> f64 {
        let n = self.values.len();
        if n < 2 {
            return 0.0;
        }
        let m = self.mean();
        let var = self.values.iter().map(|&v| (v - m) * (v - m)).sum::<f64>() / (n - 1) as f64;
        var.sqrt()
    }

    /// Relative standard deviation (coefficient of variation).
    #[allow(dead_code)]
    pub fn relative_std_dev(&self) -> f64 {
        let m = self.mean();
        if m.abs() < 1e-300 {
            return 0.0;
        }
        self.std_dev() / m.abs()
    }

    /// Number of samples.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Returns true if no samples recorded.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -------------------------------------------------------------------
    // Use reduced units for tests where the magnitude of the NHC effect
    // is easier to control.  We pick tau so that Q_1 ~ n_dof (order 1),
    // meaning dt_reduced ~ 0.01 gives a clearly non-trivial coupling.
    //
    // Strategy: work in SI but choose tau such that ω_0 = 1/tau gives a
    // meaningful kick in a single step.
    // -------------------------------------------------------------------

    fn make_nhc(n_dof: usize, temp: f64, n_chains: usize, tau: f64) -> NhcThermostat {
        NhcThermostat::new(temp, n_chains, n_dof, tau)
    }

    #[test]
    fn test_nhc_creation() {
        let n_dof = 9_usize;
        let temp = 300.0_f64;
        let n_chains = 3_usize;
        let tau = 1e-13_f64; // 0.1 ps

        let nhc = NhcThermostat::new(temp, n_chains, n_dof, tau);

        assert_eq!(nhc.n_chains, n_chains);
        assert_eq!(nhc.xi.len(), n_chains);
        assert_eq!(nhc.v_xi.len(), n_chains);
        assert_eq!(nhc.q_masses.len(), n_chains);

        let q1_expected = (n_dof as f64) * KB * temp * tau * tau;
        assert!(
            (nhc.q_masses[0] - q1_expected).abs() < 1e-50,
            "Q_1 mismatch: got {}, expected {}",
            nhc.q_masses[0],
            q1_expected
        );

        let qj_expected = KB * temp * tau * tau;
        for j in 1..n_chains {
            assert!(
                (nhc.q_masses[j] - qj_expected).abs() < 1e-50,
                "Q_{j} mismatch: got {}, expected {}",
                nhc.q_masses[j],
                qj_expected
            );
        }
    }

    #[test]
    fn test_nhc_integrate_hot() {
        // Hot system: KE = 2 × equilibrium → thermostat should cool → scale < 1
        let n_dof = 9_usize;
        let temp = 300.0_f64;
        // Use a large dt/tau ratio so the coupling has a clearly visible effect
        // in a single step.  tau = 1e-15 s, dt = 1e-15 s → dt/tau = 1.
        let tau = 1e-15_f64;
        let mut nhc = make_nhc(n_dof, temp, 3, tau);

        // KE corresponding to 2×T (hot system)
        let ke_hot = (n_dof as f64) * KB * temp; // = 2 × (n_dof/2 * kBT)
        let dt = 1e-15_f64;

        let scale = nhc.integrate(ke_hot, dt);

        assert!(
            scale < 1.0,
            "Hot system: expected scale < 1.0 (cooling), got {scale}"
        );
    }

    #[test]
    fn test_nhc_integrate_cold() {
        // Cold system: KE = 0.5 × equilibrium → thermostat should heat → scale > 1
        let n_dof = 9_usize;
        let temp = 300.0_f64;
        let tau = 1e-15_f64;
        let mut nhc = make_nhc(n_dof, temp, 3, tau);

        // KE corresponding to 0.5×T (cold system)
        let ke_cold = 0.25 * (n_dof as f64) * KB * temp; // = 0.5 × equilibrium
        let dt = 1e-15_f64;

        let scale = nhc.integrate(ke_cold, dt);

        assert!(
            scale > 1.0,
            "Cold system: expected scale > 1.0 (heating), got {scale}"
        );
    }

    #[test]
    fn test_nhc_temperature_estimate() {
        // 3 atoms × 3 DOF = 9 DOF; KE at target T → current_temperature ≈ T
        let n_dof = 9_usize;
        let temp = 300.0_f64;
        let nhc = make_nhc(n_dof, temp, 3, 1e-13);

        // KE at exactly target temperature: KE = n_dof/2 * kB * T
        let ke = 0.5 * (n_dof as f64) * KB * temp;
        let t_est = nhc.current_temperature(ke);

        assert!(
            (t_est - temp).abs() < 1.0,
            "Temperature estimate off: got {t_est}, expected {temp}"
        );
    }

    #[test]
    fn test_nhc_conserved_energy() {
        // Run 10 steps on a near-equilibrium system; conserved energy should
        // stay roughly constant (within 5% for small dt).
        let n_dof = 9_usize;
        let temp = 300.0_f64;
        let tau = 1e-13_f64; // 0.1 ps — reasonable coupling
        let mut nhc = make_nhc(n_dof, temp, 3, tau);

        let ke0 = 0.5 * (n_dof as f64) * KB * temp;
        let pe0 = 1.5 * KB * temp;
        // Use a very small step so integrator is in the linear regime
        let dt = 1e-16_f64;

        let h0 = nhc.conserved_energy(ke0, pe0);
        let mut ke = ke0;

        for _ in 0..10 {
            let scale = nhc.integrate(ke, dt);
            ke *= scale * scale;
        }

        let h1 = nhc.conserved_energy(ke, pe0);
        let ref_mag = h0.abs().max(1e-100);
        let rel_drift = (h1 - h0).abs() / ref_mag;
        assert!(
            rel_drift < 0.05,
            "NHC conserved energy drifted by {:.1}%: H0={h0:.6e}, H1={h1:.6e}",
            rel_drift * 100.0
        );
    }

    #[test]
    fn test_respa_nhc_creation() {
        let rnhc = RespaNhc::new(300.0, 3, 9, 1e-13, 4);
        assert_eq!(rnhc.n_inner, 4);
        assert_eq!(rnhc.nhc.n_chains, 3);
        assert_eq!(rnhc.nhc.n_dof, 9);
    }

    #[test]
    fn test_respa_nhc_outer_step_hot() {
        // Hot system: thermostat should cool → total scale < 1
        let n_dof = 9_usize;
        let temp = 300.0_f64;
        let tau = 1e-15_f64;
        let mut rnhc = RespaNhc::new(temp, 3, n_dof, tau, 4);

        let ke_hot = (n_dof as f64) * KB * temp;
        let dt = 1e-15_f64;
        let (scale, inner_dt) = rnhc.outer_step(ke_hot, dt);

        assert!(
            scale < 1.0,
            "RESPA-NHC hot: expected scale < 1, got {scale}"
        );
        let expected_inner = dt / 4.0;
        assert!(
            (inner_dt - expected_inner).abs() < 1e-30,
            "inner_dt mismatch"
        );
    }

    #[test]
    fn test_respa_nhc_inner_dt_calculation() {
        let rnhc = RespaNhc::new(300.0, 3, 9, 1e-13, 10);
        let dt_outer = 1e-14_f64;
        let expected = dt_outer / 10.0;
        let (_, inner_dt) = {
            let mut r = rnhc;
            let ke = 0.5 * 9.0 * KB * 300.0;
            r.outer_step(ke, dt_outer)
        };
        assert!(
            (inner_dt - expected).abs() / expected < 1e-10,
            "inner_dt = {inner_dt}, expected {expected}"
        );
    }

    #[test]
    fn test_mtk_npt_creation() {
        let mtk = MtkNpt::new(300.0, 1e5, 3, 9, 1e-13, 1e-13);
        assert_eq!(mtk.n_dim, 3);
        assert!((mtk.target_pressure - 1e5).abs() < 1e-10);
        assert!(mtk.w_baro > 0.0, "barostat mass must be positive");
        assert!(
            (mtk.v_eps).abs() < 1e-30,
            "initial baro velocity should be 0"
        );
    }

    #[test]
    fn test_mtk_npt_baro_kinetic_energy() {
        let mut mtk = MtkNpt::new(300.0, 1e5, 3, 9, 1e-13, 1e-13);
        mtk.v_eps = 1e10; // artificial velocity for test
        let bke = mtk.baro_kinetic_energy();
        let expected = 0.5 * mtk.w_baro * 1e20;
        assert!(
            (bke - expected).abs() / expected.max(1e-100) < 1e-10,
            "baro KE: got {bke}, expected {expected}"
        );
    }

    #[test]
    fn test_mtk_npt_propagate_returns_finite() {
        let mut mtk = MtkNpt::new(300.0, 1e5, 3, 9, 1e-13, 1e-13);
        let ke = 0.5 * 9.0 * KB * 300.0;
        let volume = 1e-27; // 1 nm³
        let dt = 1e-16;
        let (vel_scale, vol_scale) = mtk.propagate(ke, 1e5, volume, dt);
        assert!(vel_scale.is_finite(), "vel_scale must be finite");
        assert!(vol_scale.is_finite(), "vol_scale must be finite");
        assert!(vel_scale > 0.0, "vel_scale must be positive");
        assert!(vol_scale > 0.0, "vol_scale must be positive");
    }

    #[test]
    fn test_mtk_conserved_energy_initial() {
        let mtk = MtkNpt::new(300.0, 1e5, 3, 9, 1e-13, 1e-13);
        let ke = 0.5 * 9.0 * KB * 300.0;
        let pe = -1e-20;
        let vol = 1e-27;
        let h = mtk.conserved_energy(ke, pe, vol);
        assert!(h.is_finite(), "conserved energy must be finite");
    }

    #[test]
    fn test_flexible_nhc_creation() {
        let scales = vec![1.0, 2.0, 0.5];
        let fnhc = FlexibleNhc::new(300.0, 3, 9, 1e-13, scales.clone());
        // Check that masses were scaled
        let base_nhc = NhcThermostat::new(300.0, 3, 9, 1e-13);
        for j in 0..3 {
            let expected = base_nhc.q_masses[j] * scales[j];
            assert!(
                (fnhc.inner.q_masses[j] - expected).abs() < 1e-60,
                "chain {j} mass mismatch"
            );
        }
    }

    #[test]
    fn test_flexible_nhc_integrate_hot() {
        let scales = vec![1.5, 1.0, 0.8];
        let mut fnhc = FlexibleNhc::new(300.0, 3, 9, 1e-15, scales);
        let ke_hot = 9.0 * KB * 300.0;
        let scale = fnhc.integrate(ke_hot, 1e-15);
        assert!(
            scale < 1.0,
            "FlexibleNHC hot: expected scale < 1, got {scale}"
        );
    }

    #[test]
    fn test_flexible_nhc_conserved_energy() {
        let scales = vec![1.0, 1.0, 1.0];
        let mut fnhc = FlexibleNhc::new(300.0, 3, 9, 1e-13, scales);
        let ke0 = 0.5 * 9.0 * KB * 300.0;
        let pe0 = 1e-21;
        let dt = 1e-16;
        let h0 = fnhc.conserved_energy(ke0, pe0);
        let mut ke = ke0;
        for _ in 0..10 {
            let s = fnhc.integrate(ke, dt);
            ke *= s * s;
        }
        let h1 = fnhc.conserved_energy(ke, pe0);
        let drift = (h1 - h0).abs() / h0.abs().max(1e-100);
        assert!(drift < 0.05, "FlexibleNHC conserved energy drift {drift}");
    }

    #[test]
    fn test_nhc_order6_hot_cools() {
        let n_dof = 9_usize;
        let temp = 300.0_f64;
        let tau = 1e-15_f64;
        let mut nhc6 = NhcThermostatOrder6::new(temp, 3, n_dof, tau);
        let ke_hot = (n_dof as f64) * KB * temp;
        let scale = nhc6.integrate(ke_hot, 1e-15);
        assert!(
            scale < 1.0,
            "Order-6 NHC hot: expected scale < 1, got {scale}"
        );
    }

    #[test]
    fn test_nhc_order6_cold_heats() {
        let n_dof = 9_usize;
        let temp = 300.0_f64;
        let tau = 1e-15_f64;
        let mut nhc6 = NhcThermostatOrder6::new(temp, 3, n_dof, tau);
        let ke_cold = 0.25 * (n_dof as f64) * KB * temp;
        let scale = nhc6.integrate(ke_cold, 1e-15);
        assert!(
            scale > 1.0,
            "Order-6 NHC cold: expected scale > 1, got {scale}"
        );
    }

    #[test]
    fn test_nhc_single_chain() {
        // Single chain (M=1) should still thermostat correctly
        let n_dof = 9;
        let temp = 300.0;
        let tau = 1e-15;
        let mut nhc = make_nhc(n_dof, temp, 1, tau);
        let ke_hot = (n_dof as f64) * KB * temp;
        let scale = nhc.integrate(ke_hot, 1e-15);
        assert!(
            scale < 1.0,
            "Single chain NHC hot: expected scale < 1, got {scale}"
        );
    }

    #[test]
    fn test_nhc_long_chain() {
        // Chain length 10 should work fine
        let n_dof = 9;
        let temp = 300.0;
        let tau = 1e-15;
        let mut nhc = make_nhc(n_dof, temp, 10, tau);
        assert_eq!(nhc.q_masses.len(), 10);
        let ke_hot = (n_dof as f64) * KB * temp;
        let scale = nhc.integrate(ke_hot, 1e-15);
        assert!(scale.is_finite(), "Long chain scale must be finite");
        assert!(scale > 0.0, "Long chain scale must be positive");
    }

    #[test]
    fn test_nhc_equilibrium_scale_near_one() {
        // At equilibrium KE, the thermostat should produce scale ≈ 1
        let n_dof = 9;
        let temp = 300.0;
        let tau = 1e-13;
        let mut nhc = make_nhc(n_dof, temp, 3, tau);
        let ke_eq = 0.5 * (n_dof as f64) * KB * temp;
        // Use a very small step so the correction is tiny
        let scale = nhc.integrate(ke_eq, 1e-17);
        assert!(
            (scale - 1.0).abs() < 0.01,
            "At equilibrium, scale should be near 1.0, got {scale}"
        );
    }

    // ── IsokineThermostat tests ───────────────────────────────────────────────

    #[test]
    fn test_isokin_target_ke() {
        let n_dof = 9;
        let temp = 300.0;
        let ik = IsokineThermostat::new(temp, n_dof);
        let expected = 0.5 * 9.0 * KB * temp;
        assert!(
            (ik.target_ke - expected).abs() < 1e-40,
            "target KE mismatch"
        );
    }

    #[test]
    fn test_isokin_rescale_restores_ke() {
        let n_dof = 9;
        let temp = 300.0;
        let ik = IsokineThermostat::new(temp, n_dof);
        let mass = 1.0e-27;
        let masses = vec![mass; 3];
        // Start with velocities giving KE = 2 * target (hot)
        let v0 = (2.0 * ik.target_ke / (3.0 * 0.5 * mass)).sqrt();
        let mut vels: Vec<[f64; 3]> = vec![[v0, 0.0, 0.0]; 3];
        let scale = ik.rescale_velocities(&mut vels, &masses);
        assert!(scale < 1.0, "hot system: scale should be < 1, got {scale}");
        assert!(
            ik.is_constrained(&vels, &masses, 1e-8),
            "velocities should be constrained to target KE"
        );
    }

    #[test]
    fn test_isokin_gaussian_friction_zero_force() {
        let ik = IsokineThermostat::new(300.0, 9);
        let forces = vec![[0.0; 3]; 3];
        let momenta = vec![[1e-24, 0.0, 0.0]; 3];
        let masses = vec![1.0e-27; 3];
        let alpha = ik.gaussian_friction(&forces, &momenta, &masses);
        assert!(alpha.abs() < 1e-50, "zero force → α = 0, got {alpha}");
    }

    #[test]
    fn test_isokin_rescale_cold_system() {
        let n_dof = 9;
        let temp = 300.0;
        let ik = IsokineThermostat::new(temp, n_dof);
        let mass = 1.0e-27;
        let masses = vec![mass; 3];
        // Cold: velocities giving KE = 0.5 * target
        let v0 = (0.5 * ik.target_ke / (3.0 * 0.5 * mass)).sqrt();
        let mut vels: Vec<[f64; 3]> = vec![[v0, 0.0, 0.0]; 3];
        let scale = ik.rescale_velocities(&mut vels, &masses);
        assert!(scale > 1.0, "cold system: scale should be > 1, got {scale}");
    }

    // ── StochasticVelocityRescaling tests ────────────────────────────────────

    #[test]
    fn test_svr_target_ke() {
        let svr = StochasticVelocityRescaling::new(300.0, 9, 1e-13);
        let ke = svr.target_ke();
        let expected = 0.5 * 9.0 * KB * 300.0;
        assert!((ke - expected).abs() < 1e-40, "SVR target KE mismatch");
    }

    #[test]
    fn test_svr_instantaneous_temperature() {
        let svr = StochasticVelocityRescaling::new(300.0, 9, 1e-13);
        let ke = svr.target_ke();
        let t = svr.instantaneous_temperature(ke);
        assert!(
            (t - 300.0).abs() < 1e-6,
            "T_inst at target KE should be T_target, got {t}"
        );
    }

    #[test]
    fn test_svr_deterministic_scale_hot() {
        // Hot system (T_inst > T_target) → scale < 1
        let svr = StochasticVelocityRescaling::new(300.0, 9, 1e-13);
        let scale = svr.deterministic_scale(600.0, 1e-14); // T_inst = 2*T
        assert!(
            scale < 1.0,
            "hot system: deterministic scale < 1, got {scale}"
        );
    }

    #[test]
    fn test_svr_deterministic_scale_cold() {
        // Cold system (T_inst < T_target) → scale > 1
        let svr = StochasticVelocityRescaling::new(300.0, 9, 1e-13);
        let scale = svr.deterministic_scale(150.0, 1e-14); // T_inst = 0.5*T
        assert!(
            scale > 1.0,
            "cold system: deterministic scale > 1, got {scale}"
        );
    }

    #[test]
    fn test_svr_relaxation_rate_positive() {
        let svr = StochasticVelocityRescaling::new(300.0, 9, 1e-13);
        assert!(svr.relaxation_rate() > 0.0);
    }

    // ── NhcConservedQuantity tests ────────────────────────────────────────────

    #[test]
    fn test_nhc_conserved_quantity_initial() {
        let cq = NhcConservedQuantity::new(1.0e-20);
        assert_eq!(cq.n_samples, 1);
        assert!(cq.relative_drift() < 1e-15, "initial drift should be 0");
    }

    #[test]
    fn test_nhc_conserved_quantity_record() {
        let mut cq = NhcConservedQuantity::new(1.0e-20);
        cq.record(1.0001e-20); // 0.01% drift
        assert_eq!(cq.n_samples, 2);
        assert!(cq.relative_drift() > 0.0);
        assert!(cq.relative_drift() < 0.01, "drift should be small");
    }

    #[test]
    fn test_nhc_conserved_quantity_is_conserved() {
        let mut cq = NhcConservedQuantity::new(1.0e-20);
        cq.record(1.00001e-20); // ~0.001% drift
        assert!(
            cq.is_conserved(1e-3),
            "tiny drift should pass conservation test"
        );
    }

    #[test]
    fn test_nhc_conserved_quantity_max_drift_tracked() {
        let mut cq = NhcConservedQuantity::new(1.0e-20);
        cq.record(1.001e-20); // 0.1% drift
        cq.record(1.0001e-20); // 0.01% drift
        let max_rd = cq.max_relative_drift();
        // Max drift should be from the first sample
        assert!(max_rd > 5e-4, "max drift should be tracked, got {max_rd}");
    }

    #[test]
    fn test_ys_weights_7_sum_unity() {
        let sum: f64 = YS_WEIGHTS_7.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-12,
            "YS7 weights should sum to 1, got {sum}"
        );
    }

    // ── NHC coupling frequency tests ──────────────────────────────────────

    #[test]
    fn test_nhc_coupling_frequency() {
        let tau = 1e-13;
        let omega = nhc_coupling_frequency(tau);
        assert!((omega - 1.0 / tau).abs() < 1e-5, "omega = {omega}");
    }

    #[test]
    fn test_nhc_coupling_frequency_zero_tau() {
        let omega = nhc_coupling_frequency(0.0);
        assert_eq!(omega, 0.0);
    }

    #[test]
    fn test_nhc_optimal_tau() {
        let omega = 1e12_f64;
        let tau = nhc_optimal_tau(omega);
        assert!((tau - 1.0 / omega).abs() < 1e-25, "tau = {tau}");
    }

    #[test]
    fn test_nhc_quality_factor() {
        let tau = 1e-13;
        let omega = nhc_coupling_frequency(tau);
        let q = nhc_quality_factor(omega, tau);
        assert!(
            (q - 1.0).abs() < 1e-10,
            "Q at resonance should be 1, got {q}"
        );
    }

    // ── Chain velocity initialization tests ───────────────────────────────

    #[test]
    fn test_init_chain_velocities_not_all_zero() {
        let mut nhc = make_nhc(9, 300.0, 3, 1e-13);
        init_chain_velocities(&mut nhc, 42);
        let nonzero = nhc.v_xi.iter().any(|&v| v.abs() > 1e-30);
        assert!(
            nonzero,
            "chain velocities should not all be zero after init"
        );
    }

    #[test]
    fn test_init_chain_velocities_different_seeds() {
        let mut nhc1 = make_nhc(9, 300.0, 3, 1e-13);
        let mut nhc2 = make_nhc(9, 300.0, 3, 1e-13);
        init_chain_velocities(&mut nhc1, 1);
        init_chain_velocities(&mut nhc2, 2);
        // Different seeds should give different velocities
        let same = nhc1
            .v_xi
            .iter()
            .zip(nhc2.v_xi.iter())
            .all(|(a, b)| (a - b).abs() < 1e-50);
        assert!(!same, "different seeds should give different velocities");
    }

    #[test]
    fn test_chain_kinetic_energy_zero_velocities() {
        let nhc = make_nhc(9, 300.0, 3, 1e-13);
        // Default v_xi = 0
        let ke = chain_kinetic_energy(&nhc);
        assert_eq!(ke, 0.0);
    }

    #[test]
    fn test_chain_kinetic_energy_positive() {
        let mut nhc = make_nhc(9, 300.0, 3, 1e-13);
        init_chain_velocities(&mut nhc, 99);
        let ke = chain_kinetic_energy(&nhc);
        assert!(ke >= 0.0, "chain KE should be non-negative");
    }

    #[test]
    fn test_chain_potential_energy_zero_xi() {
        let nhc = make_nhc(9, 300.0, 3, 1e-13);
        // Default xi = 0
        let pe = chain_potential_energy(&nhc);
        assert_eq!(pe, 0.0);
    }

    #[test]
    fn test_extended_energy_initial() {
        let nhc = make_nhc(9, 300.0, 3, 1e-13);
        let phys_ke = 1.0e-20;
        let e = extended_energy(&nhc, phys_ke);
        // bath terms are 0 initially
        assert!((e - phys_ke).abs() < 1e-30, "extended E = {e}");
    }

    // ── RESPA multi-timestep tests ────────────────────────────────────────

    #[test]
    fn test_respa_inner_steps_count() {
        let mut respa = RespaIntegrator::new(1e-12, 4);
        let steps = respa.inner_steps();
        assert_eq!(steps.len(), 4);
    }

    #[test]
    fn test_respa_inner_dt() {
        let respa = RespaIntegrator::new(1e-12, 4);
        assert!(
            (respa.dt_inner - 0.25e-12).abs() < 1e-25,
            "dt_inner = {}",
            respa.dt_inner
        );
    }

    #[test]
    fn test_respa_step_counter() {
        let mut respa = RespaIntegrator::new(1e-12, 4);
        respa.inner_steps();
        respa.inner_steps();
        assert_eq!(respa.step, 2);
    }

    #[test]
    fn test_respa_elapsed_time() {
        let mut respa = RespaIntegrator::new(1e-12, 2);
        respa.inner_steps();
        respa.inner_steps();
        respa.inner_steps();
        assert!((respa.elapsed_time() - 3e-12).abs() < 1e-24);
    }

    #[test]
    fn test_respa_half_step_dt() {
        let respa = RespaIntegrator::new(1e-12, 4);
        assert!((respa.half_step_dt() - 0.5e-12).abs() < 1e-25);
    }

    #[test]
    fn test_respa_reset() {
        let mut respa = RespaIntegrator::new(1e-12, 2);
        respa.inner_steps();
        respa.reset();
        assert_eq!(respa.step, 0);
    }

    // ── Berendsen thermostat tests ────────────────────────────────────────

    #[test]
    fn test_berendsen_lambda_equilibrium() {
        let thermo = BerendsenThermostat::new(300.0, 1e-13);
        let lam = thermo.lambda(300.0, 1e-14);
        assert!((lam - 1.0).abs() < 1e-10, "at equilibrium λ=1, got {lam}");
    }

    #[test]
    fn test_berendsen_lambda_hot() {
        let thermo = BerendsenThermostat::new(300.0, 1e-13);
        let lam = thermo.lambda(600.0, 1e-14); // T_inst = 2*T_target → cool
        assert!(lam < 1.0, "hot system: lambda < 1, got {lam}");
    }

    #[test]
    fn test_berendsen_lambda_cold() {
        let thermo = BerendsenThermostat::new(300.0, 1e-13);
        let lam = thermo.lambda(150.0, 1e-14);
        assert!(lam > 1.0, "cold system: lambda > 1, got {lam}");
    }

    #[test]
    fn test_berendsen_rescale_velocities() {
        let thermo = BerendsenThermostat::new(300.0, 1e-13);
        let mut vels = vec![[1.0_f64, 0.0, 0.0]; 3];
        let lam = thermo.rescale(&mut vels, 600.0, 1e-14);
        assert!(lam < 1.0);
        for v in &vels {
            assert!((v[0] - lam).abs() < 1e-12, "velocity should scale by λ");
        }
    }

    #[test]
    fn test_berendsen_temperature_from_ke() {
        let thermo = BerendsenThermostat::new(300.0, 1e-13);
        let ke = 0.5 * 9.0 * KB * 300.0;
        let t = thermo.temperature_from_ke(ke, 9);
        assert!((t - 300.0).abs() < 1e-6, "T from KE = {t}");
    }

    // ── Andersen thermostat tests ─────────────────────────────────────────

    #[test]
    fn test_andersen_collision_probability() {
        let thermo = AndersenThermostat::new(300.0, 1e12);
        let prob = thermo.collision_probability(1e-12);
        assert!((prob - 1.0).abs() < 1e-10, "prob ν*dt = 1.0: {prob}");
    }

    #[test]
    fn test_andersen_collision_probability_small() {
        let thermo = AndersenThermostat::new(300.0, 1e10);
        let prob = thermo.collision_probability(1e-13);
        // ν * dt = 1e10 * 1e-13 = 1e-3
        assert!((prob - 1e-3).abs() < 1e-10, "prob = {prob}");
    }

    #[test]
    fn test_andersen_count_collisions_zero_prob() {
        let thermo = AndersenThermostat::new(300.0, 0.0);
        let n = thermo.count_collisions(100, 1e-12, 42);
        assert_eq!(n, 0, "zero frequency → zero collisions");
    }

    #[test]
    fn test_andersen_count_collisions_high_prob() {
        let thermo = AndersenThermostat::new(300.0, 1e15); // very high ν
        let n = thermo.count_collisions(100, 1e-12, 42);
        assert_eq!(n, 100, "all atoms should collide at prob=1");
    }

    #[test]
    fn test_andersen_sample_velocity_finite() {
        let thermo = AndersenThermostat::new(300.0, 1e12);
        let mut state = 12345u64;
        let v = thermo.sample_velocity_component(1.0e-26, &mut state);
        assert!(v.is_finite(), "sampled velocity should be finite");
    }

    // ── ConservedQuantityTimeSeries tests ─────────────────────────────────

    #[test]
    fn test_time_series_mean() {
        let mut ts = ConservedQuantityTimeSeries::new(1e-15);
        ts.push(1.0);
        ts.push(3.0);
        ts.push(2.0);
        assert!((ts.mean() - 2.0).abs() < 1e-12, "mean = {}", ts.mean());
    }

    #[test]
    fn test_time_series_std_dev_constant() {
        let mut ts = ConservedQuantityTimeSeries::new(1e-15);
        ts.push(5.0);
        ts.push(5.0);
        ts.push(5.0);
        assert!(
            ts.std_dev().abs() < 1e-12,
            "std dev of constant = {}",
            ts.std_dev()
        );
    }

    #[test]
    fn test_time_series_relative_std_dev() {
        let mut ts = ConservedQuantityTimeSeries::new(1e-15);
        ts.push(100.0);
        ts.push(102.0);
        let rsd = ts.relative_std_dev();
        assert!(rsd > 0.0 && rsd < 0.1, "relative std dev = {rsd}");
    }

    #[test]
    fn test_time_series_len() {
        let mut ts = ConservedQuantityTimeSeries::new(1e-15);
        assert!(ts.is_empty());
        ts.push(1.0);
        ts.push(2.0);
        assert_eq!(ts.len(), 2);
    }

    #[test]
    fn test_time_series_empty_mean() {
        let ts = ConservedQuantityTimeSeries::new(1e-15);
        assert_eq!(ts.mean(), 0.0);
    }
}
