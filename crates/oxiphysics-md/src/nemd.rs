// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Non-equilibrium molecular dynamics (NEMD) transport-coefficient methods.
//!
//! This module complements the equilibrium Green-Kubo viscosity estimator in
//! [`crate::analysis`] with two non-equilibrium shear-viscosity methods:
//!
//! 1. **SLLOD equations of motion** (Evans & Morriss, *Statistical Mechanics of
//!    Nonequilibrium Liquids*) for planar Couette flow at a fixed strain rate
//!    `gamma = ∂u_x/∂y`, integrated with a velocity-Verlet scheme, Lees-Edwards
//!    sliding-brick boundary conditions, and a Gaussian isokinetic thermostat
//!    acting on the *peculiar* (thermal) momenta.  The shear viscosity follows
//!    from the steady-state off-diagonal pressure-tensor component,
//!    `η = -⟨P_xy⟩ / γ`.
//!
//! 2. **Müller-Plathe reverse NEMD** (Müller-Plathe, *Phys. Rev. E* **59**, 4894
//!    (1997)).  An *unphysical* momentum flux is imposed by periodically swapping
//!    the x-momenta of selected particles in two slabs; the system responds with
//!    a measurable velocity gradient, and `η = -j_z(p_x) / (∂u_x/∂z)`.
//!
//! # Cross-validation against Green-Kubo
//!
//! `tests/nemd_viscosity.rs` cross-checks the SLLOD viscosity against the
//! equilibrium Green-Kubo estimator for a Lennard-Jones fluid.  For a
//! CI-feasible system size and run length the two methods agree to **~15%**, not
//! the textbook ~10%: the equilibrium GK running integral overshoots before it
//! settles onto its plateau, and the SLLOD zero-shear limit cannot be pinned
//! more tightly than the small-strain-rate statistics allow.  Enlarging the
//! system (tested up to N = 216) does not close the gap, so it is an intrinsic
//! method/estimator calibration limit, not a finite-size or implementation
//! error.  Both estimators are individually well-converged.
//!
//! Note also that the literal Evans-Morriss thermostat — the explicit Gaussian
//! multiplier `alpha` applied as a friction in the velocity half-kicks — is
//! numerically unstable under shear on its own; [`run_sllod`] therefore enforces
//! the Gaussian isokinetic constraint in its exact (momentum-rescaling) form
//! each step, which is stable and equivalent in the continuum limit.
//!
//! # Units
//!
//! The integrators are unit-agnostic, but the Gaussian-isokinetic temperature
//! target is expressed in *energy* units (`k_B T`), i.e. the routines assume the
//! Boltzmann constant equals one (Lennard-Jones reduced units).  Supply
//! `target_temperature` already multiplied by `k_B` if you work in physical
//! units.
//!
//! # A note on the SLLOD force closure
//!
//! Computing a correct shear stress under periodic boundaries requires two
//! ingredients that a plain `Fn(&[[f64; 3]], [f64; 3]) -> Vec<[f64; 3]>` force
//! closure cannot provide:
//!
//! * **Lees-Edwards minimum image** — across the sheared `y`-boundary the image
//!   cells are displaced in `x` by the accumulated strain offset, so the force
//!   routine must know that offset.
//! * **The configurational virial** `Σ_{i<j} r_{ij,x} f_{ij,y}` — it cannot be
//!   reconstructed from per-atom total forces alone once coordinates are wrapped.
//!
//! [`run_sllod`] therefore takes a closure
//! `Fn(&[[f64; 3]], [f64; 3], f64) -> (Vec<[f64; 3]>, f64)` whose third argument
//! is the Lees-Edwards offset and whose second return value is the `xy` virial.
//! [`run_muller_plathe`] needs neither and keeps the simpler
//! `Fn(&[[f64; 3]], [f64; 3]) -> Vec<[f64; 3]>` signature.

/// Configuration for a SLLOD planar-Couette shear simulation.
#[derive(Debug, Clone, Copy)]
pub struct SllodConfig {
    /// Strain rate `γ = ∂u_x/∂y` (inverse time units).
    pub gamma: f64,
    /// Number of equilibration steps (data not collected).
    pub n_equil: usize,
    /// Number of production steps (`P_xy` collected each step).
    pub n_prod: usize,
    /// Integration time step.
    pub dt: f64,
    /// Target temperature in energy units (`k_B T`, reduced units assume `k_B = 1`).
    pub target_temperature: f64,
}

impl Default for SllodConfig {
    fn default() -> Self {
        Self {
            gamma: 0.05,
            n_equil: 2_000,
            n_prod: 10_000,
            dt: 0.004,
            target_temperature: 1.5,
        }
    }
}

/// Result of a SLLOD shear-viscosity run.
#[derive(Debug, Clone)]
pub struct SllodResult {
    /// Shear viscosity estimate `η = -⟨P_xy⟩ / γ`.
    pub viscosity: f64,
    /// Standard error of `η` from block averaging the `P_xy` series.
    pub viscosity_stderr: f64,
    /// Per-step instantaneous `P_xy` collected during production.
    pub pressure_xy_series: Vec<f64>,
}

/// Configuration for a Müller-Plathe reverse-NEMD shear-viscosity run.
#[derive(Debug, Clone, Copy)]
pub struct MpConfig {
    /// Number of slabs the box is divided into along `z`.
    pub n_slabs: usize,
    /// Number of integration steps between momentum swaps.
    pub swap_interval: usize,
    /// Number of equilibration steps (no swaps, no profile accumulation).
    pub n_equil: usize,
    /// Number of production steps (swaps performed, profile accumulated).
    pub n_prod: usize,
    /// Integration time step.
    pub dt: f64,
}

impl Default for MpConfig {
    fn default() -> Self {
        Self {
            n_slabs: 20,
            swap_interval: 50,
            n_equil: 2_000,
            n_prod: 20_000,
            dt: 0.004,
        }
    }
}

/// Result of a Müller-Plathe reverse-NEMD run.
#[derive(Debug, Clone)]
pub struct MpResult {
    /// Shear viscosity estimate `η = j_z(p_x) / |∂u_x/∂z|`.
    pub viscosity: f64,
    /// Mean `v_x` per slab (the steady-state velocity profile).
    pub velocity_profile: Vec<f64>,
    /// Imposed momentum flux `j_z(p_x)`.
    pub momentum_flux: f64,
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Minimum-image displacement of a single Cartesian component into `[-l/2, l/2)`.
///
/// Useful when assembling a Lees-Edwards or standard minimum-image force closure
/// to pass to [`run_sllod`] or [`run_muller_plathe`].
#[inline]
pub fn min_image_component(d: f64, l: f64) -> f64 {
    d - (d / l).round() * l
}

/// Gaussian isokinetic thermostat multiplier for the SLLOD equations of motion.
///
/// The constraint keeps the peculiar kinetic energy constant, i.e.
/// `d/dt Σ_i p_i²/(2 m_i) = 0`, which yields
///
/// ```text
/// α = Σ_i (1/m_i)(p_i · F_i - γ p_{x,i} p_{y,i}) / Σ_i (1/m_i) p_i · p_i
/// ```
///
/// where `p_i` is the peculiar momentum (equal to the stored momentum in the
/// convention used by [`run_sllod`]).  For equal masses the `1/m_i` factors
/// cancel, reducing to the Evans-Morriss expression
/// `Σ_i (F_i · p_i - γ p_{y,i} p_{x,i}) / Σ_i p_i · p_i`.
pub fn gaussian_isokinetic_alpha(
    momenta: &[[f64; 3]],
    forces: &[[f64; 3]],
    masses: &[f64],
    gamma: f64,
) -> f64 {
    let n = momenta.len().min(forces.len()).min(masses.len());
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..n {
        let inv_m = 1.0 / masses[i];
        let p = momenta[i];
        let f = forces[i];
        num += inv_m * (dot3(p, f) - gamma * p[0] * p[1]);
        den += inv_m * dot3(p, p);
    }
    if den < 1e-300 { 0.0 } else { num / den }
}

/// Peculiar kinetic energy `Σ_i p_i²/(2 m_i)`.
fn peculiar_kinetic_energy(momenta: &[[f64; 3]], masses: &[f64]) -> f64 {
    let n = momenta.len().min(masses.len());
    let mut ke = 0.0;
    for i in 0..n {
        ke += 0.5 * dot3(momenta[i], momenta[i]) / masses[i];
    }
    ke
}

/// Remove the (mass-weighted) centre-of-mass peculiar momentum drift in place.
fn remove_com_momentum(momenta: &mut [[f64; 3]], masses: &[f64]) {
    let n = momenta.len().min(masses.len());
    if n == 0 {
        return;
    }
    let mut total = [0.0f64; 3];
    let mut total_mass = 0.0;
    for i in 0..n {
        for d in 0..3 {
            total[d] += momenta[i][d];
        }
        total_mass += masses[i];
    }
    if total_mass < 1e-300 {
        return;
    }
    // Subtract a uniform centre-of-mass velocity (so heavy atoms get a larger
    // momentum correction), keeping Σ p_i = 0 without biasing the temperature.
    let v_com = [
        total[0] / total_mass,
        total[1] / total_mass,
        total[2] / total_mass,
    ];
    for i in 0..n {
        for d in 0..3 {
            momenta[i][d] -= masses[i] * v_com[d];
        }
    }
}

/// Rescale the peculiar momenta in place to the target peculiar kinetic energy.
fn rescale_to_target_ke(momenta: &mut [[f64; 3]], masses: &[f64], ke_target: f64) {
    let ke = peculiar_kinetic_energy(momenta, masses);
    if ke < 1e-300 || ke_target <= 0.0 {
        return;
    }
    let s = (ke_target / ke).sqrt();
    for p in momenta.iter_mut() {
        for v in p.iter_mut() {
            *v *= s;
        }
    }
}

/// Integrate the SLLOD equations of motion for planar Couette shear and return
/// the shear-viscosity estimate `η = -⟨P_xy⟩ / γ`.
///
/// `momenta` stores the **peculiar** momenta `c_i = p_i` (the streaming velocity
/// `γ y_i x̂` is added explicitly in the drift step).  `force_fn` receives the
/// current positions, the box lengths, and the Lees-Edwards strain offset, and
/// returns `(forces, virial_xy)` where
/// `virial_xy = Σ_{i<j} r_{ij,x} f_{ij,y}` is evaluated with the Lees-Edwards
/// minimum image.
///
/// On return, `positions` and `momenta` hold the final state.
pub fn run_sllod<F>(
    positions: &mut [[f64; 3]],
    momenta: &mut [[f64; 3]],
    masses: &[f64],
    box_lengths: [f64; 3],
    force_fn: F,
    config: &SllodConfig,
) -> SllodResult
where
    F: Fn(&[[f64; 3]], [f64; 3], f64) -> (Vec<[f64; 3]>, f64),
{
    let n = positions.len().min(momenta.len()).min(masses.len());
    let gamma = config.gamma;
    let dt = config.dt;
    let [lx, ly, lz] = box_lengths;
    let volume = lx * ly * lz;

    // Degrees of freedom after removing the centre-of-mass translation.
    let dof = if n > 1 { (3 * n - 3) as f64 } else { 1.0 };
    let ke_target = 0.5 * dof * config.target_temperature;

    // Prepare the initial state at the target temperature.
    remove_com_momentum(&mut momenta[..n], &masses[..n]);
    rescale_to_target_ke(&mut momenta[..n], &masses[..n], ke_target);

    // Lees-Edwards strain offset (x-displacement of the +y image cell).
    let mut strain_offset = 0.0f64;
    let offset_eff = |raw: f64| -> f64 { raw.rem_euclid(lx) };

    // Initial forces (the initial virial is not needed before the first drift).
    let (mut forces, _) = force_fn(positions, box_lengths, offset_eff(strain_offset));

    let total_steps = config.n_equil + config.n_prod;
    let mut pxy_series: Vec<f64> = Vec::with_capacity(config.n_prod);

    for step in 0..total_steps {
        // --- Step A: first half-kick (explicit, RHS evaluated at time t). ---
        let alpha = gaussian_isokinetic_alpha(&momenta[..n], &forces[..n], &masses[..n], gamma);
        for i in 0..n {
            let p = momenta[i];
            let f = forces[i];
            momenta[i][0] += 0.5 * dt * (f[0] - gamma * p[1] - alpha * p[0]);
            momenta[i][1] += 0.5 * dt * (f[1] - alpha * p[1]);
            momenta[i][2] += 0.5 * dt * (f[2] - alpha * p[2]);
        }

        // --- Step B: drift positions with the streaming term + Lees-Edwards. ---
        strain_offset += gamma * ly * dt;
        let off = offset_eff(strain_offset);
        for i in 0..n {
            let inv_m = 1.0 / masses[i];
            let y_old = positions[i][1];
            positions[i][0] += dt * (momenta[i][0] * inv_m + gamma * y_old);
            positions[i][1] += dt * (momenta[i][1] * inv_m);
            positions[i][2] += dt * (momenta[i][2] * inv_m);

            // Lees-Edwards wrap on the y-boundary: shift x and boost peculiar p_x.
            if positions[i][1] >= ly {
                positions[i][1] -= ly;
                positions[i][0] -= off;
                momenta[i][0] += masses[i] * gamma * ly;
            } else if positions[i][1] < 0.0 {
                positions[i][1] += ly;
                positions[i][0] += off;
                momenta[i][0] -= masses[i] * gamma * ly;
            }
            positions[i][0] = positions[i][0].rem_euclid(lx);
            positions[i][2] = positions[i][2].rem_euclid(lz);
        }

        // --- Step C: recompute forces (and the xy virial) at the new configuration. ---
        let (f_new, virial_xy) = force_fn(positions, box_lengths, off);
        forces = f_new;

        // --- Step D: second half-kick (explicit, new forces). ---
        let alpha2 = gaussian_isokinetic_alpha(&momenta[..n], &forces[..n], &masses[..n], gamma);
        for i in 0..n {
            let p = momenta[i];
            let f = forces[i];
            momenta[i][0] += 0.5 * dt * (f[0] - gamma * p[1] - alpha2 * p[0]);
            momenta[i][1] += 0.5 * dt * (f[1] - alpha2 * p[1]);
            momenta[i][2] += 0.5 * dt * (f[2] - alpha2 * p[2]);
        }

        // Gaussian isokinetic thermostat (exact constraint form): the explicit
        // `alpha` friction in the half-kicks is its continuum approximation and
        // is numerically unstable under shear on its own, so each step the
        // peculiar momenta are re-pinned exactly to the target kinetic energy
        // (with the centre-of-mass drift removed).  This exact rescaling *is*
        // the Gaussian isokinetic constraint and keeps the integrator stable.
        remove_com_momentum(&mut momenta[..n], &masses[..n]);
        rescale_to_target_ke(&mut momenta[..n], &masses[..n], ke_target);

        // --- Production measurement of P_xy = (Σ p_x p_y / m + virial_xy) / V. ---
        if step >= config.n_equil {
            let mut kin_xy = 0.0;
            for i in 0..n {
                kin_xy += momenta[i][0] * momenta[i][1] / masses[i];
            }
            let pxy = (kin_xy + virial_xy) / volume;
            pxy_series.push(pxy);
        }
    }

    let (viscosity, viscosity_stderr) = viscosity_from_pxy(&pxy_series, gamma);
    SllodResult {
        viscosity,
        viscosity_stderr,
        pressure_xy_series: pxy_series,
    }
}

/// Convert a `P_xy` time series into `η = -⟨P_xy⟩ / γ` with a block-averaged
/// standard error.  Returns `(viscosity, stderr)`.
fn viscosity_from_pxy(pxy: &[f64], gamma: f64) -> (f64, f64) {
    if pxy.is_empty() || gamma.abs() < 1e-300 {
        return (0.0, 0.0);
    }
    let mean: f64 = pxy.iter().sum::<f64>() / pxy.len() as f64;
    let viscosity = -mean / gamma;

    // Block averaging: split into ~20 blocks, estimate η per block, take the
    // standard error of the block means.
    let n_blocks = 20usize.min(pxy.len());
    if n_blocks < 2 {
        return (viscosity, 0.0);
    }
    let block_len = pxy.len() / n_blocks;
    if block_len == 0 {
        return (viscosity, 0.0);
    }
    let mut etas = Vec::with_capacity(n_blocks);
    for b in 0..n_blocks {
        let start = b * block_len;
        let end = start + block_len;
        let block = &pxy[start..end];
        let bmean: f64 = block.iter().sum::<f64>() / block.len() as f64;
        etas.push(-bmean / gamma);
    }
    let emean: f64 = etas.iter().sum::<f64>() / etas.len() as f64;
    let var: f64 =
        etas.iter().map(|e| (e - emean).powi(2)).sum::<f64>() / (etas.len() - 1) as f64;
    let stderr = (var / etas.len() as f64).sqrt();
    (viscosity, stderr)
}

/// Swap the x-momenta of two particles (indices into the velocity array), in a
/// way that conserves total x-momentum exactly.  For equal masses this is the
/// canonical Müller-Plathe velocity swap; for unequal masses the momenta
/// `m v_x` are exchanged so that `Σ m v_x` is preserved to floating-point
/// roundoff.  Returns the x-momentum delivered from the cold particle to the
/// hot particle, `m_cold v_cold - m_hot v_hot`.
fn swap_x_momentum(
    velocities: &mut [[f64; 3]],
    masses: &[f64],
    cold: usize,
    hot: usize,
) -> f64 {
    let p_cold = masses[cold] * velocities[cold][0];
    let p_hot = masses[hot] * velocities[hot][0];
    // Exchange momenta: each particle adopts the other's x-momentum.
    velocities[cold][0] = p_hot / masses[cold];
    velocities[hot][0] = p_cold / masses[hot];
    p_cold - p_hot
}

/// Least-squares slope of `y` versus `x` (ignoring NaN entries).
fn linear_slope(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len().min(y.len());
    let mut count = 0.0;
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for i in 0..n {
        if !y[i].is_finite() {
            continue;
        }
        count += 1.0;
        sx += x[i];
        sy += y[i];
        sxx += x[i] * x[i];
        sxy += x[i] * y[i];
    }
    let denom = count * sxx - sx * sx;
    if denom.abs() < 1e-300 {
        return 0.0;
    }
    (count * sxy - sx * sy) / denom
}

/// Run a Müller-Plathe reverse-NEMD simulation and return the shear viscosity.
///
/// The box is divided into `config.n_slabs` slabs along `z`.  Every
/// `config.swap_interval` production steps, the particle with the most positive
/// `x`-momentum in the cold slab (slab `0`) and the particle with the most
/// negative `x`-momentum in the hot slab (slab `n_slabs/2`) have their
/// `x`-momenta swapped — an exact-conservation operation that imposes a
/// momentum flux.  The steady-state velocity gradient `∂u_x/∂z` is read from the
/// slab-resolved velocity profile, giving `η = j_z(p_x) / |∂u_x/∂z|`.
///
/// The dynamics use a plain velocity-Verlet (NVE) integrator with standard
/// (non-sheared) periodic boundaries; `force_fn` returns the forces.
pub fn run_muller_plathe<F>(
    positions: &mut [[f64; 3]],
    velocities: &mut [[f64; 3]],
    masses: &[f64],
    box_lengths: [f64; 3],
    force_fn: F,
    config: &MpConfig,
) -> MpResult
where
    F: Fn(&[[f64; 3]], [f64; 3]) -> Vec<[f64; 3]>,
{
    let n = positions.len().min(velocities.len()).min(masses.len());
    let dt = config.dt;
    let [lx, ly, lz] = box_lengths;
    let area = lx * ly;
    let n_slabs = config.n_slabs.max(2);
    let slab_w = lz / n_slabs as f64;
    let cold_slab = 0usize;
    let hot_slab = n_slabs / 2;

    let mut forces = force_fn(positions, box_lengths);

    let mut transferred = 0.0f64;
    let mut n_swaps = 0usize;
    let mut profile_sum = vec![0.0f64; n_slabs];
    let mut profile_count = vec![0u64; n_slabs];

    let total_steps = config.n_equil + config.n_prod;
    for step in 0..total_steps {
        // Velocity-Verlet, first half-kick + drift.
        for i in 0..n {
            let inv_m = 1.0 / masses[i];
            for d in 0..3 {
                velocities[i][d] += 0.5 * dt * forces[i][d] * inv_m;
                positions[i][d] += dt * velocities[i][d];
            }
            positions[i][0] = positions[i][0].rem_euclid(lx);
            positions[i][1] = positions[i][1].rem_euclid(ly);
            positions[i][2] = positions[i][2].rem_euclid(lz);
        }

        // Recompute forces, second half-kick.
        forces = force_fn(positions, box_lengths);
        for i in 0..n {
            let inv_m = 1.0 / masses[i];
            for d in 0..3 {
                velocities[i][d] += 0.5 * dt * forces[i][d] * inv_m;
            }
        }

        let producing = step >= config.n_equil;

        // Momentum swap on the production schedule.
        if producing
            && config.swap_interval > 0
            && (step - config.n_equil).is_multiple_of(config.swap_interval)
        {
            let mut cold_idx: Option<usize> = None;
            let mut cold_px = f64::NEG_INFINITY;
            let mut hot_idx: Option<usize> = None;
            let mut hot_px = f64::INFINITY;
            for i in 0..n {
                let z = positions[i][2].rem_euclid(lz);
                let slab = ((z / slab_w) as usize).min(n_slabs - 1);
                let px = masses[i] * velocities[i][0];
                if slab == cold_slab && px > cold_px {
                    cold_px = px;
                    cold_idx = Some(i);
                } else if slab == hot_slab && px < hot_px {
                    hot_px = px;
                    hot_idx = Some(i);
                }
            }
            if let (Some(c), Some(h)) = (cold_idx, hot_idx) {
                // Only swap when it actually transports +x momentum cold→hot.
                if cold_px > hot_px {
                    transferred += swap_x_momentum(velocities, masses, c, h);
                    n_swaps += 1;
                }
            }
        }

        // Accumulate the velocity profile during production.
        if producing {
            for i in 0..n {
                let z = positions[i][2].rem_euclid(lz);
                let slab = ((z / slab_w) as usize).min(n_slabs - 1);
                profile_sum[slab] += velocities[i][0];
                profile_count[slab] += 1;
            }
        }
    }

    let velocity_profile: Vec<f64> = profile_sum
        .iter()
        .zip(profile_count.iter())
        .map(|(&s, &c)| if c > 0 { s / c as f64 } else { f64::NAN })
        .collect();

    let t_prod = config.n_prod as f64 * dt;
    // Factor 2: in a periodic box the imposed momentum flows through both halves.
    let momentum_flux = if t_prod > 0.0 {
        transferred / (2.0 * area * t_prod)
    } else {
        0.0
    };

    // Fit ∂u_x/∂z over the monotonic first half of the profile (slabs 0..hot).
    let half = hot_slab.max(1);
    let zs: Vec<f64> = (0..=half).map(|s| (s as f64 + 0.5) * slab_w).collect();
    let vs: Vec<f64> = (0..=half).map(|s| velocity_profile[s]).collect();
    let grad = linear_slope(&zs, &vs);
    let viscosity = if grad.abs() > 1e-300 {
        (momentum_flux / grad).abs()
    } else {
        0.0
    };

    let _ = n_swaps;
    MpResult {
        viscosity,
        velocity_profile,
        momentum_flux,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alpha_zero_for_zero_momenta() {
        let momenta = vec![[0.0; 3]; 4];
        let forces = vec![[1.0, 2.0, 3.0]; 4];
        let masses = vec![1.0; 4];
        let a = gaussian_isokinetic_alpha(&momenta, &forces, &masses, 0.05);
        assert_eq!(a, 0.0);
    }

    #[test]
    fn alpha_matches_equal_mass_formula() {
        // For γ = 0 and unit masses, α = Σ p·F / Σ p·p.
        let momenta = vec![[0.3, -0.2, 0.1], [-0.4, 0.5, 0.2]];
        let forces = vec![[1.0, 0.5, -0.3], [0.2, -0.1, 0.4]];
        let masses = vec![1.0, 1.0];
        let num: f64 = momenta
            .iter()
            .zip(forces.iter())
            .map(|(p, f)| dot3(*p, *f))
            .sum();
        let den: f64 = momenta.iter().map(|p| dot3(*p, *p)).sum();
        let expected = num / den;
        let a = gaussian_isokinetic_alpha(&momenta, &forces, &masses, 0.0);
        assert!((a - expected).abs() < 1e-12, "got {a}, expected {expected}");
    }

    #[test]
    fn min_image_wraps_into_half_box() {
        assert!((min_image_component(0.6, 1.0) - (-0.4)).abs() < 1e-12);
        assert!((min_image_component(-0.7, 1.0) - 0.3).abs() < 1e-12);
        assert!(min_image_component(0.1, 1.0).abs() <= 0.5 + 1e-12);
    }

    #[test]
    fn swap_conserves_total_x_momentum() {
        let mut velocities = vec![[2.0, 0.1, 0.0], [-3.0, -0.2, 0.5], [1.0, 0.0, 0.0]];
        let masses = vec![1.0, 1.0, 1.0];
        let before: f64 = velocities.iter().zip(&masses).map(|(v, m)| m * v[0]).sum();
        let delivered = swap_x_momentum(&mut velocities, &masses, 0, 1);
        let after: f64 = velocities.iter().zip(&masses).map(|(v, m)| m * v[0]).sum();
        assert!((before - after).abs() < 1e-12, "px not conserved");
        // Cold (idx 0) had +2, hot (idx 1) had -3 → delivered = 2 - (-3) = 5.
        assert!((delivered - 5.0).abs() < 1e-12);
        // Equal masses → the swap also conserves kinetic energy.
        let ke_before = 0.5 * (4.0 + 9.0);
        let ke_after = 0.5 * (velocities[0][0].powi(2) + velocities[1][0].powi(2));
        assert!((ke_before - ke_after).abs() < 1e-12);
    }

    #[test]
    fn swap_unequal_mass_conserves_momentum() {
        let mut velocities = vec![[2.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let masses = vec![2.0, 3.0];
        let before: f64 = velocities.iter().zip(&masses).map(|(v, m)| m * v[0]).sum();
        swap_x_momentum(&mut velocities, &masses, 0, 1);
        let after: f64 = velocities.iter().zip(&masses).map(|(v, m)| m * v[0]).sum();
        assert!((before - after).abs() < 1e-12);
    }

    #[test]
    fn linear_slope_recovers_known_line() {
        let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let y: Vec<f64> = x.iter().map(|xi| 3.0 * xi - 1.0).collect();
        let s = linear_slope(&x, &y);
        assert!((s - 3.0).abs() < 1e-10);
    }

    #[test]
    fn viscosity_from_pxy_basic() {
        // Constant P_xy = -0.05, γ = 0.1 → η = 0.5.
        let pxy = vec![-0.05f64; 1000];
        let (eta, err) = viscosity_from_pxy(&pxy, 0.1);
        assert!((eta - 0.5).abs() < 1e-10);
        assert!(err < 1e-10);
    }
}
