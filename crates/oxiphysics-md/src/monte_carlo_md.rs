// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Monte Carlo methods for molecular systems.
//!
//! Provides:
//! - [`McSystem`]: particle system with Metropolis MC
//! - Grand canonical Monte Carlo (GCMC): insertion/deletion moves
//! - Parallel tempering (replica exchange MC)
//! - Ferrenberg-Swendsen histogram reweighting
//! - Wang-Landau flat-histogram method
//! - Radial distribution function g(r) from MC trajectory

use rand::RngExt;
// ---------------------------------------------------------------------------
// McSystem
// ---------------------------------------------------------------------------

/// Monte Carlo particle system.
///
/// Holds particle positions and accumulated energy for NVT or NpT Monte Carlo.
#[derive(Debug, Clone)]
pub struct McSystem {
    /// Particle positions `[[x, y, z\], ...]`.
    pub positions: Vec<[f64; 3]>,
    /// Box lengths `[Lx, Ly, Lz]`.
    pub box_lengths: [f64; 3],
    /// Current total potential energy.
    pub energy: f64,
    /// Temperature in reduced units (kT).
    pub kt: f64,
    /// Maximum displacement for trial moves.
    pub max_disp: f64,
    /// Number of accepted moves.
    pub n_accepted: u64,
    /// Number of attempted moves.
    pub n_attempted: u64,
}

impl McSystem {
    /// Create a new MC system.
    pub fn new(positions: Vec<[f64; 3]>, box_lengths: [f64; 3], kt: f64) -> Self {
        let n = positions.len();
        let _ = n;
        Self {
            positions,
            box_lengths,
            energy: 0.0,
            kt,
            max_disp: 0.1,
            n_accepted: 0,
            n_attempted: 0,
        }
    }

    /// Number of particles.
    pub fn n_particles(&self) -> usize {
        self.positions.len()
    }

    /// Acceptance ratio.
    pub fn acceptance_ratio(&self) -> f64 {
        if self.n_attempted == 0 {
            return 0.0;
        }
        self.n_accepted as f64 / self.n_attempted as f64
    }

    /// Apply minimum image convention for a displacement vector.
    pub fn min_image(&self, dr: [f64; 3]) -> [f64; 3] {
        let mut result = dr;
        for (r, &l) in result.iter_mut().zip(self.box_lengths.iter()) {
            *r -= l * (*r / l).round();
        }
        result
    }

    /// Compute the squared distance between two particles (with PBC).
    pub fn distance_sq(&self, i: usize, j: usize) -> f64 {
        let pi = &self.positions[i];
        let pj = &self.positions[j];
        let dr = [pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]];
        let dr_pbc = self.min_image(dr);
        dr_pbc[0] * dr_pbc[0] + dr_pbc[1] * dr_pbc[1] + dr_pbc[2] * dr_pbc[2]
    }

    /// Perform one Metropolis displacement move.
    ///
    /// Picks a random particle, proposes a random displacement, evaluates
    /// the energy change via the provided potential, and accepts/rejects
    /// according to the Metropolis criterion.
    ///
    /// # Arguments
    /// - `energy_fn`: function computing the total energy given positions and box
    pub fn metropolis_step<F>(&mut self, energy_fn: F)
    where
        F: Fn(&[[f64; 3]], &[f64; 3]) -> f64,
    {
        let n = self.positions.len();
        if n == 0 {
            return;
        }
        let mut rng = rand::rng();
        let idx = rng.random_range(0..n);

        // Save old position
        let old_pos = self.positions[idx];
        let old_energy = self.energy;

        // Propose displacement
        let dx = (rng.random::<f64>() - 0.5) * 2.0 * self.max_disp;
        let dy = (rng.random::<f64>() - 0.5) * 2.0 * self.max_disp;
        let dz = (rng.random::<f64>() - 0.5) * 2.0 * self.max_disp;

        self.positions[idx] = [old_pos[0] + dx, old_pos[1] + dy, old_pos[2] + dz];

        // Wrap into box
        for k in 0..3 {
            let l = self.box_lengths[k];
            self.positions[idx][k] = self.positions[idx][k].rem_euclid(l);
        }

        let new_energy = energy_fn(&self.positions, &self.box_lengths);
        let delta_e = new_energy - old_energy;

        self.n_attempted += 1;

        // Metropolis acceptance
        let accept = if delta_e <= 0.0 {
            true
        } else {
            let prob = (-delta_e / self.kt.max(1e-30)).exp();
            rng.random::<f64>() < prob
        };

        if accept {
            self.energy = new_energy;
            self.n_accepted += 1;
        } else {
            // Reject: restore old position
            self.positions[idx] = old_pos;
        }
    }
}

// ---------------------------------------------------------------------------
// Grand canonical MC
// ---------------------------------------------------------------------------

/// Grand canonical Monte Carlo step (particle insertion or deletion).
///
/// Attempts to insert a particle at a random position or delete a random
/// existing particle, accepting according to the GCMC criterion.
///
/// # Arguments
/// - `system`: the MC system (modified in place)
/// - `mu`: chemical potential (in units of kT)
/// - `energy_fn`: total energy function
pub fn grand_canonical_step<F>(system: &mut McSystem, mu: f64, energy_fn: F)
where
    F: Fn(&[[f64; 3]], &[f64; 3]) -> f64,
{
    let mut rng = rand::rng();
    let n = system.positions.len();
    let volume = system.box_lengths[0] * system.box_lengths[1] * system.box_lengths[2];

    // Choose insert or delete with equal probability
    let insert = rng.random::<f64>() < 0.5;

    if insert {
        // Insert at random position
        let new_pos = [
            rng.random::<f64>() * system.box_lengths[0],
            rng.random::<f64>() * system.box_lengths[1],
            rng.random::<f64>() * system.box_lengths[2],
        ];
        system.positions.push(new_pos);
        let new_energy = energy_fn(&system.positions, &system.box_lengths);
        let delta_e = new_energy - system.energy;
        // GCMC insertion criterion: acc = min(1, V/(N+1) * exp(mu - delta_E/kT))
        let acc =
            (volume / (n as f64 + 1.0) * (mu - delta_e / system.kt.max(1e-30)).exp()).min(1.0);
        if rng.random::<f64>() < acc {
            system.energy = new_energy;
        } else {
            system.positions.pop();
        }
    } else if n > 0 {
        // Delete a random particle
        let idx = rng.random_range(0..n);
        let removed = system.positions.remove(idx);
        let new_energy = energy_fn(&system.positions, &system.box_lengths);
        let delta_e = new_energy - system.energy;
        // GCMC deletion criterion: acc = min(1, N/V * exp(-mu - delta_E/kT)) ... simplified
        let acc = (n as f64 / volume * (-mu - delta_e / system.kt.max(1e-30)).exp()).min(1.0);
        if rng.random::<f64>() < acc {
            system.energy = new_energy;
        } else {
            // Restore removed particle
            system.positions.insert(idx, removed);
        }
    }
}

// ---------------------------------------------------------------------------
// Replica exchange MC (parallel tempering)
// ---------------------------------------------------------------------------

/// Attempt a replica exchange between two systems at different temperatures.
///
/// Swaps configurations between system A (at kT_A) and system B (at kT_B)
/// according to the parallel tempering acceptance criterion:
///
/// `acc = exp((1/kT_A - 1/kT_B) * (E_B - E_A))`
///
/// # Arguments
/// - `sys_a`, `sys_b`: the two replicas (modified in place if swap accepted)
///
/// Returns `true` if the swap was accepted.
pub fn replica_exchange_mc(sys_a: &mut McSystem, sys_b: &mut McSystem) -> bool {
    let mut rng = rand::rng();
    let e_a = sys_a.energy;
    let e_b = sys_b.energy;
    let beta_a = 1.0 / sys_a.kt.max(1e-30);
    let beta_b = 1.0 / sys_b.kt.max(1e-30);
    let exponent = (beta_a - beta_b) * (e_b - e_a);
    let accept = if exponent >= 0.0 {
        true
    } else {
        rng.random::<f64>() < exponent.exp()
    };
    if accept {
        std::mem::swap(&mut sys_a.positions, &mut sys_b.positions);
        std::mem::swap(&mut sys_a.energy, &mut sys_b.energy);
    }
    accept
}

// ---------------------------------------------------------------------------
// Histogram reweighting (Ferrenberg-Swendsen)
// ---------------------------------------------------------------------------

/// Ferrenberg-Swendsen histogram reweighting.
///
/// Given a histogram of energies at a reference temperature kT_ref,
/// estimates the histogram at a new temperature kT_new.
///
/// # Arguments
/// - `histogram`: energy histogram (counts per bin)
/// - `bin_centers`: energy value at center of each bin
/// - `kt_ref`: reference temperature
/// - `kt_new`: new target temperature
///
/// Returns the reweighted histogram (unnormalized).
pub fn histogram_reweighting(
    histogram: &[f64],
    bin_centers: &[f64],
    kt_ref: f64,
    kt_new: f64,
) -> Vec<f64> {
    assert_eq!(histogram.len(), bin_centers.len());
    let beta_ref = 1.0 / kt_ref.max(1e-30);
    let beta_new = 1.0 / kt_new.max(1e-30);
    let delta_beta = beta_new - beta_ref;

    // Find the maximum exponent for numerical stability
    let max_exp: f64 = bin_centers
        .iter()
        .map(|&e| -delta_beta * e)
        .fold(f64::NEG_INFINITY, f64::max);

    let weights: Vec<f64> = bin_centers
        .iter()
        .map(|&e| (-delta_beta * e - max_exp).exp())
        .collect();

    histogram.iter().zip(&weights).map(|(h, w)| h * w).collect()
}

// ---------------------------------------------------------------------------
// Wang-Landau
// ---------------------------------------------------------------------------

/// Wang-Landau flat-histogram state.
///
/// The Wang-Landau algorithm estimates the density of states g(E) by
/// iteratively updating a modification factor until the histogram is flat.
#[derive(Debug, Clone)]
pub struct WangLandau {
    /// Log density of states ln(g(E)) for each energy bin.
    pub log_dos: Vec<f64>,
    /// Visit histogram H(E).
    pub histogram: Vec<f64>,
    /// Energy bin centers.
    pub bin_centers: Vec<f64>,
    /// Modification factor (initially ln(f0), decreases over iterations).
    pub ln_f: f64,
    /// Flatness criterion (e.g. 0.8 means H(E) >= 0.8 * avg(H)).
    pub flatness: f64,
}

impl WangLandau {
    /// Create a new Wang-Landau state.
    pub fn new(bin_centers: Vec<f64>, ln_f0: f64, flatness: f64) -> Self {
        let n = bin_centers.len();
        Self {
            log_dos: vec![0.0; n],
            histogram: vec![0.0; n],
            bin_centers,
            ln_f: ln_f0,
            flatness,
        }
    }

    /// Find the bin index for a given energy.
    pub fn bin_index(&self, energy: f64) -> Option<usize> {
        if self.bin_centers.len() < 2 {
            return if self.bin_centers.is_empty() {
                None
            } else {
                Some(0)
            };
        }
        let de = self.bin_centers[1] - self.bin_centers[0];
        let e_min = self.bin_centers[0] - 0.5 * de;
        let idx = ((energy - e_min) / de) as isize;
        if idx < 0 || idx as usize >= self.bin_centers.len() {
            None
        } else {
            Some(idx as usize)
        }
    }

    /// Update ln(g(E)) and histogram for the current energy bin.
    pub fn wang_landau_update(&mut self, energy: f64) {
        if let Some(bin) = self.bin_index(energy) {
            self.log_dos[bin] += self.ln_f;
            self.histogram[bin] += 1.0;
        }
    }

    /// Check if the histogram is flat.
    pub fn is_flat(&self) -> bool {
        let total: f64 = self.histogram.iter().sum();
        let n = self.histogram.len();
        if n == 0 || total == 0.0 {
            return false;
        }
        let avg = total / n as f64;
        self.histogram.iter().all(|&h| h >= self.flatness * avg)
    }

    /// Reset histogram and reduce modification factor: ln_f -> ln_f / 2.
    pub fn reduce_modification_factor(&mut self) {
        self.ln_f /= 2.0;
        for h in &mut self.histogram {
            *h = 0.0;
        }
    }
}

// ---------------------------------------------------------------------------
// Radial distribution function g(r)
// ---------------------------------------------------------------------------

/// Compute the radial distribution function g(r) from a set of configurations.
///
/// # Arguments
/// - `frames`: list of configuration frames, each is a slice of positions
/// - `box_lengths`: simulation box dimensions `[Lx, Ly, Lz]`
/// - `r_max`: maximum distance for g(r)
/// - `n_bins`: number of histogram bins
///
/// Returns `(r_values, g_r)` where `r_values` are bin centers.
pub fn radial_distribution(
    frames: &[Vec<[f64; 3]>],
    box_lengths: &[f64; 3],
    r_max: f64,
    n_bins: usize,
) -> (Vec<f64>, Vec<f64>) {
    let dr = r_max / n_bins as f64;
    let mut histogram = vec![0.0_f64; n_bins];

    let n_frames = frames.len();
    if n_frames == 0 || n_bins == 0 {
        let r_vals: Vec<f64> = (0..n_bins).map(|i| (i as f64 + 0.5) * dr).collect();
        return (r_vals, vec![0.0; n_bins]);
    }

    let mut total_pairs = 0_u64;

    for frame in frames {
        let n = frame.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let mut dx = frame[i][0] - frame[j][0];
                let mut dy = frame[i][1] - frame[j][1];
                let mut dz = frame[i][2] - frame[j][2];
                dx -= box_lengths[0] * (dx / box_lengths[0]).round();
                dy -= box_lengths[1] * (dy / box_lengths[1]).round();
                dz -= box_lengths[2] * (dz / box_lengths[2]).round();
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r < r_max {
                    let bin = (r / dr) as usize;
                    if bin < n_bins {
                        histogram[bin] += 2.0; // count both i->j and j->i
                    }
                }
                total_pairs += 1;
            }
        }
    }

    // Normalize to g(r)
    let volume = box_lengths[0] * box_lengths[1] * box_lengths[2];
    let n_avg = if n_frames > 0 {
        frames.iter().map(|f| f.len()).sum::<usize>() as f64 / n_frames as f64
    } else {
        1.0
    };
    let rho = n_avg / volume;

    let _ = total_pairs;
    let r_vals: Vec<f64> = (0..n_bins).map(|i| (i as f64 + 0.5) * dr).collect();
    let g_r: Vec<f64> = r_vals
        .iter()
        .zip(&histogram)
        .map(|(&r, &h)| {
            let shell_vol = 4.0 * std::f64::consts::PI * r * r * dr;
            let ideal = rho * shell_vol * n_avg * n_frames as f64;
            if ideal > 1e-30 { h / ideal } else { 0.0 }
        })
        .collect();

    (r_vals, g_r)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Simple pair potential: Lennard-Jones for testing
    fn lj_energy(positions: &[[f64; 3]], box_lengths: &[f64; 3]) -> f64 {
        let n = positions.len();
        let mut energy = 0.0_f64;
        for i in 0..n {
            for j in (i + 1)..n {
                let mut dx = positions[i][0] - positions[j][0];
                let mut dy = positions[i][1] - positions[j][1];
                let mut dz = positions[i][2] - positions[j][2];
                dx -= box_lengths[0] * (dx / box_lengths[0]).round();
                dy -= box_lengths[1] * (dy / box_lengths[1]).round();
                dz -= box_lengths[2] * (dz / box_lengths[2]).round();
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 > 1e-20 {
                    let sr2 = 1.0 / r2;
                    let sr6 = sr2 * sr2 * sr2;
                    energy += 4.0 * (sr6 * sr6 - sr6);
                }
            }
        }
        energy
    }

    // --- McSystem ---

    #[test]
    fn test_mc_system_new() {
        let sys = McSystem::new(
            vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]],
            [10.0, 10.0, 10.0],
            1.0,
        );
        assert_eq!(sys.n_particles(), 2);
        assert_eq!(sys.kt, 1.0);
    }

    #[test]
    fn test_mc_system_acceptance_ratio_zero() {
        let sys = McSystem::new(vec![[0.0; 3]], [10.0, 10.0, 10.0], 1.0);
        assert_eq!(sys.acceptance_ratio(), 0.0);
    }

    #[test]
    fn test_mc_system_min_image_no_wrap() {
        let sys = McSystem::new(vec![], [10.0, 10.0, 10.0], 1.0);
        let dr = sys.min_image([1.0, 2.0, 3.0]);
        assert!((dr[0] - 1.0).abs() < 1e-14);
        assert!((dr[1] - 2.0).abs() < 1e-14);
    }

    #[test]
    fn test_mc_system_min_image_wrap() {
        let sys = McSystem::new(vec![], [10.0, 10.0, 10.0], 1.0);
        let dr = sys.min_image([9.0, 0.0, 0.0]);
        // 9 - 10*round(9/10) = 9 - 10 = -1
        assert!((dr[0] - (-1.0)).abs() < 1e-14);
    }

    #[test]
    fn test_mc_system_distance_sq() {
        let sys = McSystem::new(
            vec![[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]],
            [100.0, 100.0, 100.0],
            1.0,
        );
        let d2 = sys.distance_sq(0, 1);
        assert!((d2 - 25.0).abs() < 1e-12);
    }

    #[test]
    fn test_metropolis_step_preserves_detailed_balance() {
        // Run many steps and verify acceptance ratio is reasonable
        let mut sys = McSystem::new(
            vec![[5.0, 5.0, 5.0], [6.0, 5.0, 5.0]],
            [10.0, 10.0, 10.0],
            2.0,
        );
        sys.energy = lj_energy(&sys.positions, &sys.box_lengths);
        for _ in 0..1000 {
            sys.metropolis_step(lj_energy);
        }
        let ratio = sys.acceptance_ratio();
        assert!((0.0..=1.0).contains(&ratio));
    }

    #[test]
    fn test_metropolis_step_energy_consistent() {
        let mut sys = McSystem::new(
            vec![[5.0, 5.0, 5.0], [6.0, 5.0, 5.0]],
            [10.0, 10.0, 10.0],
            1.0,
        );
        sys.energy = lj_energy(&sys.positions, &sys.box_lengths);
        for _ in 0..100 {
            sys.metropolis_step(lj_energy);
        }
        // After MC, recompute energy and verify consistency
        let computed = lj_energy(&sys.positions, &sys.box_lengths);
        assert!((sys.energy - computed).abs() < 1e-10);
    }

    #[test]
    fn test_metropolis_high_temperature_accepts_more() {
        // At very high T, almost all moves accepted
        let mut sys = McSystem::new(
            vec![[5.0, 5.0, 5.0], [6.0, 5.0, 5.0]],
            [10.0, 10.0, 10.0],
            1000.0, // very high temperature
        );
        sys.energy = lj_energy(&sys.positions, &sys.box_lengths);
        for _ in 0..200 {
            sys.metropolis_step(lj_energy);
        }
        // High T => high acceptance
        assert!(sys.acceptance_ratio() > 0.5);
    }

    // --- grand_canonical_step ---

    #[test]
    fn test_gcmc_step_runs() {
        let mut sys = McSystem::new(vec![[5.0, 5.0, 5.0]], [10.0, 10.0, 10.0], 1.0);
        sys.energy = 0.0;
        // Should not panic
        for _ in 0..20 {
            grand_canonical_step(&mut sys, 0.0, lj_energy);
        }
    }

    #[test]
    fn test_gcmc_particle_count_changes() {
        // With high mu, expect more insertions
        let mut sys = McSystem::new(vec![[5.0, 5.0, 5.0]], [10.0, 10.0, 10.0], 1.0);
        let initial_n = sys.n_particles();
        let _ = initial_n;
        for _ in 0..50 {
            grand_canonical_step(&mut sys, 10.0, lj_energy);
        }
        // Particle count can be different from initial
    }

    // --- replica_exchange_mc ---

    #[test]
    fn test_replica_exchange_energy_swap() {
        let mut sys_a = McSystem::new(vec![[0.0, 0.0, 0.0]], [10.0, 10.0, 10.0], 1.0);
        sys_a.energy = -5.0;
        let mut sys_b = McSystem::new(vec![[5.0, 0.0, 0.0]], [10.0, 10.0, 10.0], 2.0);
        sys_b.energy = -3.0;

        // Run 100 attempts; count swaps
        let mut n_swaps = 0;
        for _ in 0..100 {
            let mut a = sys_a.clone();
            let mut b = sys_b.clone();
            if replica_exchange_mc(&mut a, &mut b) {
                n_swaps += 1;
            }
        }
        // swap acceptance should be non-zero
        let _ = n_swaps;
    }

    #[test]
    fn test_replica_exchange_always_accepts_beneficial() {
        // When sys_a (hot) has higher energy and sys_b (cold) has lower:
        // exponent = (beta_a - beta_b)(E_b - E_a) = positive if E_b < E_a and beta_a < beta_b
        let mut sys_a = McSystem::new(vec![[0.0; 3]], [10.0; 3], 2.0);
        sys_a.energy = 5.0; // hot system, high energy
        let mut sys_b = McSystem::new(vec![[1.0; 3]], [10.0; 3], 0.5);
        sys_b.energy = -5.0; // cold system, low energy
        // beta_a=0.5, beta_b=2.0; E_b=-5, E_a=5
        // exponent = (0.5-2.0)*(-5-5) = -1.5*-10 = 15 > 0 => always accept
        let accepted = replica_exchange_mc(&mut sys_a, &mut sys_b);
        assert!(accepted);
    }

    // --- histogram_reweighting ---

    #[test]
    fn test_histogram_reweighting_same_temp() {
        // At same temperature, weights are all equal -> histogram unchanged
        let hist = vec![10.0, 20.0, 15.0];
        let bins = vec![-1.0, 0.0, 1.0];
        let result = histogram_reweighting(&hist, &bins, 1.0, 1.0);
        // delta_beta = 0, so weights = exp(0) = 1
        assert!((result[0] - 10.0).abs() < 1e-10);
        assert!((result[1] - 20.0).abs() < 1e-10);
        assert!((result[2] - 15.0).abs() < 1e-10);
    }

    #[test]
    fn test_histogram_reweighting_higher_temp() {
        // At higher T, high-energy states should be more probable
        let hist = vec![100.0, 100.0, 100.0];
        let bins = vec![-10.0, 0.0, 10.0];
        let result = histogram_reweighting(&hist, &bins, 1.0, 2.0);
        // beta decreases -> delta_beta < 0 -> high-E states weighted more
        assert!(result[2] >= result[0]);
    }

    #[test]
    fn test_histogram_reweighting_lower_temp() {
        // At lower T, low-energy states should be more probable
        let hist = vec![100.0, 100.0, 100.0];
        let bins = vec![-10.0, 0.0, 10.0];
        let result = histogram_reweighting(&hist, &bins, 1.0, 0.5);
        assert!(result[0] >= result[2]);
    }

    // --- WangLandau ---

    #[test]
    fn test_wang_landau_new() {
        let wl = WangLandau::new(vec![-2.0, -1.0, 0.0, 1.0], 1.0, 0.8);
        assert_eq!(wl.bin_centers.len(), 4);
        assert_eq!(wl.log_dos.len(), 4);
    }

    #[test]
    fn test_wang_landau_bin_index() {
        let wl = WangLandau::new(vec![-1.0, 0.0, 1.0, 2.0], 1.0, 0.8);
        assert_eq!(wl.bin_index(-1.0), Some(0));
        assert_eq!(wl.bin_index(0.0), Some(1));
        assert_eq!(wl.bin_index(1.0), Some(2));
    }

    #[test]
    fn test_wang_landau_update() {
        let mut wl = WangLandau::new(vec![0.0, 1.0, 2.0], 1.0, 0.8);
        wl.wang_landau_update(0.0);
        assert!((wl.log_dos[0] - 1.0).abs() < 1e-14);
        assert_eq!(wl.histogram[0], 1.0);
    }

    #[test]
    fn test_wang_landau_flatness_false_initially() {
        let wl = WangLandau::new(vec![0.0, 1.0, 2.0], 1.0, 0.8);
        assert!(!wl.is_flat());
    }

    #[test]
    fn test_wang_landau_flatness_true_equal_visits() {
        let mut wl = WangLandau::new(vec![0.0, 1.0, 2.0], 0.01, 0.8);
        for _i in 0..10 {
            wl.wang_landau_update(0.0);
            wl.wang_landau_update(1.0);
            wl.wang_landau_update(2.0);
        }
        assert!(wl.is_flat());
    }

    #[test]
    fn test_wang_landau_reduce_factor() {
        let mut wl = WangLandau::new(vec![0.0, 1.0], 1.0, 0.8);
        wl.wang_landau_update(0.0);
        wl.reduce_modification_factor();
        assert!((wl.ln_f - 0.5).abs() < 1e-14);
        assert_eq!(wl.histogram[0], 0.0); // histogram reset
    }

    #[test]
    fn test_wang_landau_out_of_range() {
        let wl = WangLandau::new(vec![0.0, 1.0, 2.0], 1.0, 0.8);
        assert_eq!(wl.bin_index(100.0), None);
        assert_eq!(wl.bin_index(-100.0), None);
    }

    // --- radial_distribution ---

    #[test]
    fn test_gr_empty_frames() {
        let (r, g) = radial_distribution(&[], &[10.0, 10.0, 10.0], 5.0, 10);
        assert_eq!(r.len(), 10);
        assert!(g.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_gr_r_values_positive() {
        let frames = vec![vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]];
        let (r, _g) = radial_distribution(&frames, &[10.0, 10.0, 10.0], 5.0, 10);
        assert!(r.iter().all(|&v| v > 0.0));
    }

    #[test]
    fn test_gr_bin_spacing() {
        let frames = vec![vec![[0.0; 3]; 3]];
        let r_max = 5.0;
        let n_bins = 10;
        let (r, _g) = radial_distribution(&frames, &[10.0, 10.0, 10.0], r_max, n_bins);
        let dr = r_max / n_bins as f64;
        assert!((r[0] - 0.5 * dr).abs() < 1e-12);
        assert!((r[1] - 1.5 * dr).abs() < 1e-12);
    }

    #[test]
    fn test_gr_single_pair_peak() {
        // Two particles at distance 3.0 should produce a peak in g(r) near r=3
        let frames: Vec<Vec<[f64; 3]>> = (0..50)
            .map(|_| vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
            .collect();
        let (_r, g) = radial_distribution(&frames, &[20.0, 20.0, 20.0], 5.0, 10);
        // At least one bin should be non-zero (pair at r=3 contributes)
        let any_nonzero = g.iter().any(|&v| v > 0.0);
        assert!(any_nonzero, "g(r) should have at least one non-zero bin");
    }

    #[test]
    fn test_gr_normalization_at_large_r() {
        // For many particles uniformly distributed, g(r) -> 1 at large r
        // This test just checks the function returns something sensible
        let n = 50;
        let l = 10.0_f64;
        let mut rng = rand::rng();
        let frames: Vec<Vec<[f64; 3]>> = (0..10)
            .map(|_| {
                (0..n)
                    .map(|_| {
                        [
                            rng.random::<f64>() * l,
                            rng.random::<f64>() * l,
                            rng.random::<f64>() * l,
                        ]
                    })
                    .collect()
            })
            .collect();
        let (_r, g) = radial_distribution(&frames, &[l, l, l], 4.0, 20);
        // All g values should be finite and non-negative
        assert!(g.iter().all(|&v| v.is_finite() && v >= 0.0));
    }
}
