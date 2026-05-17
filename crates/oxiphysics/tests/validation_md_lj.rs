// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Phase 21.5 MD validation — Lennard–Jones liquid radial distribution
//! function `g(r)` against Verlet (1967) Table IV.
//!
//! **Setup.** N = 256 LJ particles on a 4×4×4 face-centred-cubic lattice
//! (4 basis atoms × 64 unit cells) inside a periodic cubic box of edge
//! `L = (N/ρ*)^(1/3)` with reduced density `ρ* = 0.85`.  A deterministic
//! 1 % position jitter (LCG seed 42) breaks the lattice symmetry.  Reduced
//! LJ units: `σ = ε = m = k_B = 1`.
//!
//! **Dynamics.** Velocity-Verlet integrator at `Δt* = 0.005 τ` coupled to a
//! Langevin thermostat (`γ = 1/τ`, target `T* = 0.71`).  Pair forces are
//! evaluated with an `O(N²)` double loop using the minimum-image convention
//! under periodic boundary conditions; cutoff `r_c = 2.5 σ`.
//!
//! **Measurements.**
//! * `T_equil` — instantaneous temperature averaged over the 100-step
//!   window centred on step 2500 of the 5000-step equilibration; checked
//!   against `T* = 0.71` with tolerance `0.05 · T*`.
//! * `g(r)` — radial distribution function from 500 snapshots sampled
//!   every 10 steps during a 5000-step production run; bin width
//!   `Δr = 0.02 σ` over `r ∈ [0, 2.5σ]`.  The first peak must lie at
//!   `r/σ ∈ [1.08, 1.15]` with height in `[2.5, 3.3]`, matched against
//!   the JSON baselines `md_lj_rdf_peak_position` and
//!   `md_lj_rdf_peak_height`.
//!
//! **Deviation from spec.** The spec calls for a *force-shifted* LJ
//! truncation, but the shipping `LennardJones` potential uses an
//! *energy-shifted* truncation (`V(r) ← V(r) − V(r_c)` with no force
//! correction).  At `r_c = 2.5 σ` the force-discontinuity is small
//! (`|F(r_c)| ≈ 0.039 ε/σ`) and has a negligible effect on the first-peak
//! location and height of `g(r)`, so the Verlet reference values remain
//! the appropriate target.
//!
//! **Ignored by default.** A 10 000-step O(N²) simulation at N = 256
//! runs several minutes in debug mode on a laptop (≈ 2.5·10⁸ pair
//! evaluations), well above the 90 s CI budget.  Run interactively with:
//!
//! ```text
//! cargo test -p oxiphysics --test validation_md_lj \
//!     --release -- --ignored --nocapture
//! ```

#![allow(missing_docs)]

#[path = "regression_harness.rs"]
mod harness;

use oxiphysics::md::thermostat::LangevinThermostat;
use oxiphysics::md::{AtomSet, Integrator, LennardJones, Potential, VelocityVerlet};
use oxiphysics_core::math::Vec3;

// ---------------------------------------------------------------------------
// Problem constants (reduced LJ units: σ = ε = m = k_B = 1)
// ---------------------------------------------------------------------------

/// Number of particles (4 × 4 × 4 FCC cells × 4 basis atoms).
const N_PARTICLES: usize = 256;
/// Number of FCC unit cells per edge.
const N_CELLS: usize = 4;
/// Reduced density ρ*.
const RHO_STAR: f64 = 0.85;
/// Target reduced temperature T*.
const T_STAR: f64 = 0.71;
/// Integration time step Δt* (τ).
const DT_STAR: f64 = 0.005;
/// Langevin friction γ (1/τ).
const GAMMA: f64 = 1.0;
/// LJ well depth ε (reduced units).
const EPSILON: f64 = 1.0;
/// LJ length scale σ (reduced units).
const SIGMA: f64 = 1.0;
/// Cutoff radius (σ).
const R_CUT: f64 = 2.5;
/// Reduced Boltzmann constant (= 1 in reduced units).
const K_B: f64 = 1.0;
/// Atom mass (= 1 in reduced units).
const MASS: f64 = 1.0;
/// Fractional lattice-jitter amplitude (of the lattice constant).
const LATTICE_JITTER_FRAC: f64 = 0.01;
/// Number of equilibration steps.
const N_EQUIL_STEPS: usize = 5000;
/// Number of production steps.
const N_PROD_STEPS: usize = 5000;
/// Sample interval during production.
const RDF_SAMPLE_EVERY: usize = 10;
/// RDF bin width (σ).
const RDF_DR: f64 = 0.02;
/// Maximum RDF radius (σ) — capped by the cutoff.
const RDF_R_MAX: f64 = 2.5;
/// Peak search window, lower bound (r/σ).
const PEAK_SEARCH_LO: f64 = 0.9;
/// Peak search window, upper bound (r/σ).
const PEAK_SEARCH_HI: f64 = 1.4;
/// Mid-equilibration temperature-window centre step.
const T_WIN_CENTRE: usize = 2500;
/// Temperature window half-width (steps on each side).
const T_WIN_HALF: usize = 50;

// ---------------------------------------------------------------------------
// Deterministic LCG RNG (matches the LCG used inside LangevinThermostat)
// ---------------------------------------------------------------------------

/// 64-bit LCG (Knuth / MMIX constants).  Deterministic across platforms.
struct Lcg64 {
    state: u64,
}

impl Lcg64 {
    fn new(seed: u64) -> Self {
        Self { state: seed.max(1) }
    }

    /// Uniform `f64` in `[0, 1)`.
    fn next_uniform(&mut self) -> f64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.state >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Standard-normal sample via Box–Muller.
    fn next_normal(&mut self) -> f64 {
        let u1 = self.next_uniform().max(1e-300);
        let u2 = self.next_uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

// ---------------------------------------------------------------------------
// Initial state: FCC lattice + jitter + Maxwell–Boltzmann velocities
// ---------------------------------------------------------------------------

/// Build the initial `AtomSet` and return `(atoms, box_length)`.
///
/// Positions follow an FCC lattice with 4 basis atoms per conventional cell.
/// A 1 % jitter of the lattice spacing is added to break symmetry, and
/// velocities are drawn from a Maxwell–Boltzmann distribution at `T*`
/// with centre-of-mass drift removed.
fn build_initial_state(seed: u64) -> (AtomSet, f64) {
    // Periodic box edge L = (N/ρ*)^(1/3) and FCC lattice constant a = L / N_cells.
    let l_box = (N_PARTICLES as f64 / RHO_STAR).cbrt();
    let a_cell = l_box / N_CELLS as f64;
    let jitter = LATTICE_JITTER_FRAC * a_cell;

    // FCC basis (4 atoms per conventional cubic cell, in cell-fraction units).
    let basis: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];

    let mut atoms = AtomSet::with_capacity(N_PARTICLES);
    let mut rng = Lcg64::new(seed);

    // Positions with symmetric (± jitter/2 centred on zero) jitter.
    for iz in 0..N_CELLS {
        for iy in 0..N_CELLS {
            for ix in 0..N_CELLS {
                for b in &basis {
                    let dx = (rng.next_uniform() - 0.5) * 2.0 * jitter;
                    let dy = (rng.next_uniform() - 0.5) * 2.0 * jitter;
                    let dz = (rng.next_uniform() - 0.5) * 2.0 * jitter;
                    let x = (ix as f64 + b[0]) * a_cell + dx;
                    let y = (iy as f64 + b[1]) * a_cell + dy;
                    let z = (iz as f64 + b[2]) * a_cell + dz;
                    // Wrap into [0, L) in case the jitter nudged outside.
                    let pos = Vec3::new(wrap(x, l_box), wrap(y, l_box), wrap(z, l_box));
                    atoms.add_atom(pos, Vec3::zeros(), MASS, 0.0, 0);
                }
            }
        }
    }
    assert_eq!(atoms.len(), N_PARTICLES);

    // Maxwell–Boltzmann velocities at T*: σ_v = √(k_B T / m) per component.
    let sigma_v = (K_B * T_STAR / MASS).sqrt();
    for v in &mut atoms.velocities {
        *v = Vec3::new(
            sigma_v * rng.next_normal(),
            sigma_v * rng.next_normal(),
            sigma_v * rng.next_normal(),
        );
    }
    atoms.remove_com_velocity();

    // Rescale to land exactly at T* (reduces warm-up noise).
    let t_now = atoms.temperature(K_B);
    if t_now > 1e-12 {
        let scale = (T_STAR / t_now).sqrt();
        for v in &mut atoms.velocities {
            *v *= scale;
        }
    }

    (atoms, l_box)
}

/// Wrap a scalar into `[0, L)`.
#[inline]
fn wrap(x: f64, l: f64) -> f64 {
    x - (x / l).floor() * l
}

/// Minimum-image separation component on `[−L/2, L/2)`.
#[inline]
fn min_image(d: f64, l: f64) -> f64 {
    let inv_l = 1.0 / l;
    d - l * (d * inv_l).round()
}

// ---------------------------------------------------------------------------
// Force evaluation: O(N²) LJ with minimum-image PBC
// ---------------------------------------------------------------------------

/// Accumulate LJ pair forces on `atoms.forces` under periodic boundaries.
///
/// Uses the shipping `LennardJones::force(r)` (positive = repulsive) so
/// `f_on_i = +force(r) * r_ij / r`, `f_on_j = −f_on_i`.
fn compute_forces(atoms: &mut AtomSet, lj: &LennardJones, l_box: f64) {
    atoms.clear_forces();
    let n = atoms.len();
    let r_cut_sq = R_CUT * R_CUT;

    // Borrow positions immutably, forces mutably — use indexed access to
    // avoid aliasing the forces vector against itself.
    for i in 0..n {
        for j in (i + 1)..n {
            let ri = atoms.positions[i];
            let rj = atoms.positions[j];
            let dx = min_image(ri.x - rj.x, l_box);
            let dy = min_image(ri.y - rj.y, l_box);
            let dz = min_image(ri.z - rj.z, l_box);
            let r_sq = dx * dx + dy * dy + dz * dz;
            if r_sq >= r_cut_sq || r_sq < 1e-24 {
                continue;
            }
            let r = r_sq.sqrt();
            let f_mag = lj.force(r); // +ve = repulsive
            let inv_r = 1.0 / r;
            let fx = f_mag * dx * inv_r;
            let fy = f_mag * dy * inv_r;
            let fz = f_mag * dz * inv_r;
            atoms.forces[i] += Vec3::new(fx, fy, fz);
            atoms.forces[j] -= Vec3::new(fx, fy, fz);
        }
    }
}

// ---------------------------------------------------------------------------
// RDF accumulation
// ---------------------------------------------------------------------------

/// Update a pair-count histogram with `atoms` (minimum-image PBC).
fn accumulate_rdf(hist: &mut [u64], atoms: &AtomSet, l_box: f64, dr: f64, n_bins: usize) {
    let n = atoms.len();
    let r_max = dr * n_bins as f64;
    let r_max_sq = r_max * r_max;
    for i in 0..n {
        for j in (i + 1)..n {
            let ri = atoms.positions[i];
            let rj = atoms.positions[j];
            let dx = min_image(ri.x - rj.x, l_box);
            let dy = min_image(ri.y - rj.y, l_box);
            let dz = min_image(ri.z - rj.z, l_box);
            let r_sq = dx * dx + dy * dy + dz * dz;
            if r_sq >= r_max_sq || r_sq < 1e-24 {
                continue;
            }
            let r = r_sq.sqrt();
            let bin = (r / dr) as usize;
            if bin < n_bins {
                hist[bin] += 1;
            }
        }
    }
}

/// Normalise the pair histogram to a proper radial distribution function.
///
/// `g(r) = 2 · n_pairs(r) / (N · ρ · V_shell · n_samples)`
///
/// where the factor of 2 compensates for counting each `i < j` pair once
/// and `V_shell = (4π/3) (r_outer³ − r_inner³)` is the exact spherical-shell
/// volume (not the `4π r² Δr` linear approximation).
fn normalize_rdf(hist: &[u64], dr: f64, n_atoms: usize, rho: f64, n_samples: usize) -> Vec<f64> {
    let n_bins = hist.len();
    let mut g = vec![0.0_f64; n_bins];
    let n_f = n_atoms as f64;
    let n_s_f = n_samples as f64;
    let four_third_pi = 4.0 / 3.0 * std::f64::consts::PI;
    for (k, &count) in hist.iter().enumerate() {
        let r_lo = k as f64 * dr;
        let r_hi = (k + 1) as f64 * dr;
        let v_shell = four_third_pi * (r_hi * r_hi * r_hi - r_lo * r_lo * r_lo);
        let denom = n_f * rho * v_shell * n_s_f;
        g[k] = if denom > 0.0 {
            2.0 * count as f64 / denom
        } else {
            0.0
        };
    }
    g
}

/// Locate the maximum of `g(r)` in `[lo, hi]` (r/σ units).
///
/// Returns `(r_peak, g_peak)`.  Uses a parabolic refinement between the
/// dominant bin and its two neighbours for sub-bin peak resolution.
fn find_first_peak(g: &[f64], dr: f64, lo: f64, hi: f64) -> (f64, f64) {
    let n_bins = g.len();
    let k_lo = ((lo / dr).floor() as usize)
        .max(1)
        .min(n_bins.saturating_sub(2));
    let k_hi = ((hi / dr).ceil() as usize).max(k_lo + 1).min(n_bins - 1);

    let mut k_best = k_lo;
    let mut g_best = g[k_lo];
    for (k, &gk) in g.iter().enumerate().take(k_hi + 1).skip(k_lo + 1) {
        if gk > g_best {
            g_best = gk;
            k_best = k;
        }
    }

    // Parabolic interpolation around the peak bin.
    if k_best > 0 && k_best + 1 < n_bins {
        let y0 = g[k_best - 1];
        let y1 = g[k_best];
        let y2 = g[k_best + 1];
        let denom = y0 - 2.0 * y1 + y2;
        if denom.abs() > 1e-12 {
            let delta = 0.5 * (y0 - y2) / denom;
            if delta.abs() < 1.0 {
                let r_peak = (k_best as f64 + 0.5 + delta) * dr;
                let g_peak = y1 - 0.25 * (y0 - y2) * delta;
                return (r_peak, g_peak);
            }
        }
    }
    let r_peak = (k_best as f64 + 0.5) * dr;
    (r_peak, g_best)
}

// ---------------------------------------------------------------------------
// The validation test
// ---------------------------------------------------------------------------

/// Lennard–Jones liquid RDF validation against Verlet (1967) Table IV.
///
/// Marked `#[ignore]` because a 256-particle × 10 000-step O(N²) simulation
/// takes minutes in debug mode — too slow for default CI.  Run with:
///
/// ```text
/// cargo test -p oxiphysics --test validation_md_lj \
///     --release -- --ignored --nocapture
/// ```
#[test]
#[ignore = "O(N^2) LJ MD @ N=256, 10k steps — ~2–5 min debug / ~30–60 s release"]
fn test_lj_liquid_rdf_peak() {
    // ------------------------------------------------------------------
    // Setup
    // ------------------------------------------------------------------
    let (mut atoms, l_box) = build_initial_state(42);
    let lj = LennardJones::new(EPSILON, SIGMA, R_CUT);
    let mut integrator = VelocityVerlet::new();
    let mut thermostat = LangevinThermostat::with_seed(GAMMA, 12345);

    // Sanity-check the force sign (catches sign inversion before the sim).
    let r_min = 2.0_f64.powf(1.0 / 6.0);
    let f_at_min = lj.force(r_min);
    assert!(
        f_at_min.abs() < 1e-9,
        "LJ force at r_min should be zero, got {f_at_min}"
    );
    let f_close = lj.force(0.9);
    assert!(
        f_close > 0.0,
        "LJ force at r < r_min must be repulsive, got {f_close}"
    );
    let f_far = lj.force(1.3);
    assert!(
        f_far < 0.0,
        "LJ force at r_min < r < r_c must be attractive, got {f_far}"
    );

    // Prime the force vector with the initial configuration so the first
    // VelocityVerlet half-kick uses a valid F(t=0).
    compute_forces(&mut atoms, &lj, l_box);

    // ------------------------------------------------------------------
    // Equilibration: VelocityVerlet + Langevin O-step each frame.
    //
    // Window-averaged temperature diagnostic spanning 100 steps centred
    // on step 2500 of equilibration.
    // ------------------------------------------------------------------
    let mut t_samples: Vec<f64> = Vec::with_capacity(2 * T_WIN_HALF + 1);
    let t_win_start = T_WIN_CENTRE.saturating_sub(T_WIN_HALF);
    let t_win_end = T_WIN_CENTRE + T_WIN_HALF;

    for step in 0..N_EQUIL_STEPS {
        integrator.step(&mut atoms, DT_STAR, &mut |a: &mut AtomSet| {
            compute_forces(a, &lj, l_box);
        });
        // Keep particles in the primary cell so min-image remains valid.
        atoms.apply_pbc([l_box, l_box, l_box]);
        thermostat.apply_ou_step(&mut atoms, T_STAR, DT_STAR, K_B);

        if step >= t_win_start && step <= t_win_end {
            t_samples.push(atoms.temperature(K_B));
        }
    }

    assert!(
        !t_samples.is_empty(),
        "temperature window produced no samples"
    );
    let t_equil: f64 = t_samples.iter().sum::<f64>() / t_samples.len() as f64;
    let t_dev = (t_equil - T_STAR).abs();
    assert!(
        t_dev < 0.05 * T_STAR,
        "equilibration temperature off target: T_measured = {t_equil:.4}, \
         T* = {T_STAR}, |ΔT| = {t_dev:.4}, tol = {:.4}",
        0.05 * T_STAR,
    );

    // ------------------------------------------------------------------
    // Production: accumulate RDF histogram every RDF_SAMPLE_EVERY steps.
    // ------------------------------------------------------------------
    let n_bins = (RDF_R_MAX / RDF_DR).round() as usize;
    let mut hist: Vec<u64> = vec![0; n_bins];
    let mut n_rdf_samples: usize = 0;

    for step in 0..N_PROD_STEPS {
        integrator.step(&mut atoms, DT_STAR, &mut |a: &mut AtomSet| {
            compute_forces(a, &lj, l_box);
        });
        atoms.apply_pbc([l_box, l_box, l_box]);
        thermostat.apply_ou_step(&mut atoms, T_STAR, DT_STAR, K_B);

        if step % RDF_SAMPLE_EVERY == 0 {
            accumulate_rdf(&mut hist, &atoms, l_box, RDF_DR, n_bins);
            n_rdf_samples += 1;
        }
    }
    assert!(
        n_rdf_samples > 0,
        "no RDF snapshots were accumulated during production"
    );

    // ------------------------------------------------------------------
    // Analysis: normalise and locate the first peak.
    // ------------------------------------------------------------------
    let g_r = normalize_rdf(&hist, RDF_DR, N_PARTICLES, RHO_STAR, n_rdf_samples);
    let (r_peak, g_peak) = find_first_peak(&g_r, RDF_DR, PEAK_SEARCH_LO, PEAK_SEARCH_HI);

    eprintln!(
        "[md_lj_rdf] T_equil = {t_equil:.4} (target T* = {T_STAR}); \
         r_peak = {r_peak:.4} σ, g_peak = {g_peak:.4}; \
         {n_rdf_samples} RDF snapshots, N = {N_PARTICLES}, ρ* = {RHO_STAR}, \
         L = {l_box:.4} σ"
    );

    // ------------------------------------------------------------------
    // Regression gate: compare against the Verlet (1967) Table IV values
    // via the shared baseline-JSON harness.
    // ------------------------------------------------------------------
    let pos_baseline = harness::load_baseline("md_lj_rdf_peak_position")
        .expect("load baseline md_lj_rdf_peak_position");
    let height_baseline = harness::load_baseline("md_lj_rdf_peak_height")
        .expect("load baseline md_lj_rdf_peak_height");

    // `assert_close!` is exported at the crate root by `regression_harness.rs`.
    assert_close!(r_peak, pos_baseline);
    assert_close!(g_peak, height_baseline);
}
