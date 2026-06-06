//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{IlSimParams, IonModel, IonState, IonType, PairType, RdfResult, WaldenPoint};

/// Boltzmann constant (J K⁻¹).
pub(super) const K_B: f64 = 1.380_649e-23;
/// Universal gas constant (J mol⁻¹ K⁻¹).
pub(super) const R_GAS: f64 = 8.314_462_618;
/// Dot product of two 3-vectors.
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Euclidean norm of a 3-vector.
#[inline]
pub(super) fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
/// Subtract two 3-vectors (a − b).
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Add two 3-vectors.
#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Scale a 3-vector by a scalar.
#[inline]
pub(super) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Minimum image convention for periodic boundaries.
#[inline]
pub(super) fn min_image(dr: [f64; 3], box_len: f64) -> [f64; 3] {
    [
        dr[0] - box_len * (dr[0] / box_len).round(),
        dr[1] - box_len * (dr[1] / box_len).round(),
        dr[2] - box_len * (dr[2] / box_len).round(),
    ]
}
/// Approximate complementary error function using Abramowitz & Stegun 7.1.26.
pub(super) fn erfc_approx(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.327_591_1 * x.abs());
    let poly = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    let result = poly * (-x * x).exp();
    if x >= 0.0 { result } else { 2.0 - result }
}
/// Compute LJ sigma from Lorentz combining rule.
pub fn mix_sigma(s1: f64, s2: f64) -> f64 {
    0.5 * (s1 + s2)
}
/// Compute LJ epsilon from Berthelot combining rule.
pub fn mix_epsilon(e1: f64, e2: f64) -> f64 {
    (e1 * e2).sqrt()
}
/// Compute LJ energy at distance r with given sigma and epsilon.
pub fn lj_energy(r: f64, sigma: f64, epsilon: f64) -> f64 {
    if r < 0.1 {
        return 1e10;
    }
    let sr6 = (sigma / r).powi(6);
    4.0 * epsilon * (sr6 * sr6 - sr6)
}
/// Compute LJ force magnitude at distance r (positive = repulsive at short r).
pub fn lj_force(r: f64, sigma: f64, epsilon: f64) -> f64 {
    if r < 0.1 {
        return 1e10;
    }
    let sr6 = (sigma / r).powi(6);
    24.0 * epsilon / r * (2.0 * sr6 * sr6 - sr6)
}
/// Compute all pairwise forces for an ionic liquid system.
pub fn compute_forces(state: &mut IonState, params: &IlSimParams) {
    state.zero_forces();
    let n = state.n;
    let bl = params.box_len;
    let rc2_lj = params.lj_cutoff * params.lj_cutoff;
    for i in 0..n {
        for j in (i + 1)..n {
            let dr = min_image(sub3(state.positions[i], state.positions[j]), bl);
            let r2 = dot3(dr, dr);
            let r = r2.sqrt();
            if r2 < rc2_lj && r > 0.1 {
                let sig = mix_sigma(state.lj_sigmas[i], state.lj_sigmas[j]);
                let eps = mix_epsilon(state.lj_epsilons[i], state.lj_epsilons[j]);
                let f_lj = lj_force(r, sig, eps);
                let fvec = scale3(dr, f_lj / r);
                state.forces[i] = add3(state.forces[i], fvec);
                state.forces[j] = sub3(state.forces[j], fvec);
            }
            let f_coul = params.wolf.force(state.charges[i], state.charges[j], r);
            if f_coul.abs() > 1e-30 {
                let fvec = scale3(dr, f_coul / r);
                state.forces[i] = add3(state.forces[i], fvec);
                state.forces[j] = sub3(state.forces[j], fvec);
            }
        }
    }
}
/// Total potential energy of the system.
pub fn compute_potential_energy(state: &IonState, params: &IlSimParams) -> f64 {
    let n = state.n;
    let bl = params.box_len;
    let rc2_lj = params.lj_cutoff * params.lj_cutoff;
    let mut energy = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let dr = min_image(sub3(state.positions[i], state.positions[j]), bl);
            let r2 = dot3(dr, dr);
            let r = r2.sqrt();
            if r2 < rc2_lj && r > 0.1 {
                let sig = mix_sigma(state.lj_sigmas[i], state.lj_sigmas[j]);
                let eps = mix_epsilon(state.lj_epsilons[i], state.lj_epsilons[j]);
                energy += lj_energy(r, sig, eps);
            }
            energy += params.wolf.energy(state.charges[i], state.charges[j], r);
        }
    }
    energy
}
/// One velocity Verlet step.
pub fn velocity_verlet_step(state: &mut IonState, params: &IlSimParams) {
    let dt = params.dt;
    let half_dt = 0.5 * dt;
    for i in 0..state.n {
        let inv_m = 1.0 / state.masses[i];
        for d in 0..3 {
            state.velocities[i][d] += half_dt * state.forces[i][d] * inv_m;
            state.positions[i][d] += dt * state.velocities[i][d];
        }
        for d in 0..3 {
            while state.positions[i][d] < 0.0 {
                state.positions[i][d] += params.box_len;
            }
            while state.positions[i][d] >= params.box_len {
                state.positions[i][d] -= params.box_len;
            }
        }
    }
    compute_forces(state, params);
    for i in 0..state.n {
        let inv_m = 1.0 / state.masses[i];
        for d in 0..3 {
            state.velocities[i][d] += half_dt * state.forces[i][d] * inv_m;
        }
    }
}
/// Apply Berendsen thermostat velocity scaling.
pub fn berendsen_thermostat(state: &mut IonState, params: &IlSimParams) {
    let current_t = state.temperature();
    if current_t < 1e-10 {
        return;
    }
    let lambda = (1.0 + params.dt / params.tau_t * (params.temperature / current_t - 1.0)).sqrt();
    for v in &mut state.velocities {
        *v = scale3(*v, lambda);
    }
}
/// Compute radial distribution function for specified pair type.
pub fn compute_rdf(
    state: &IonState,
    params: &IlSimParams,
    pair_type: PairType,
    n_bins: usize,
    r_max: f64,
) -> RdfResult {
    let dr = r_max / n_bins as f64;
    let mut histogram = vec![0.0_f64; n_bins];
    let bl = params.box_len;
    let mut n_pairs = 0_usize;
    let n_type_i;
    let n_type_j;
    for i in 0..state.n {
        for j in (i + 1)..state.n {
            let include = match pair_type {
                PairType::CationAnion => {
                    (state.ion_types[i] == IonType::Cation && state.ion_types[j] == IonType::Anion)
                        || (state.ion_types[i] == IonType::Anion
                            && state.ion_types[j] == IonType::Cation)
                }
                PairType::CationCation => {
                    state.ion_types[i] == IonType::Cation && state.ion_types[j] == IonType::Cation
                }
                PairType::AnionAnion => {
                    state.ion_types[i] == IonType::Anion && state.ion_types[j] == IonType::Anion
                }
                PairType::All => true,
            };
            if !include {
                continue;
            }
            let rij = min_image(sub3(state.positions[i], state.positions[j]), bl);
            let r = norm3(rij);
            let bin = (r / dr) as usize;
            if bin < n_bins {
                histogram[bin] += 2.0;
                n_pairs += 1;
            }
        }
    }
    match pair_type {
        PairType::CationAnion => {
            n_type_i = state.n_cations();
            n_type_j = state.n_anions();
        }
        PairType::CationCation => {
            n_type_i = state.n_cations();
            n_type_j = state.n_cations();
        }
        PairType::AnionAnion => {
            n_type_i = state.n_anions();
            n_type_j = state.n_anions();
        }
        PairType::All => {
            n_type_i = state.n;
            n_type_j = state.n;
        }
    }
    let volume = bl * bl * bl;
    let density_j = if pair_type == PairType::CationCation
        || pair_type == PairType::AnionAnion
        || pair_type == PairType::All
    {
        if n_type_j > 1 {
            (n_type_j - 1) as f64 / volume
        } else {
            1.0 / volume
        }
    } else {
        n_type_j as f64 / volume
    };
    let mut r_vals = Vec::with_capacity(n_bins);
    let mut g_vals = Vec::with_capacity(n_bins);
    for (k, &h_k) in histogram.iter().enumerate() {
        let r_lo = k as f64 * dr;
        let r_hi = r_lo + dr;
        let r_mid = 0.5 * (r_lo + r_hi);
        let shell_vol = (4.0 / 3.0) * PI * (r_hi.powi(3) - r_lo.powi(3));
        let ideal = n_type_i as f64 * density_j * shell_vol;
        let g = if ideal > 0.0 { h_k / ideal } else { 0.0 };
        r_vals.push(r_mid);
        g_vals.push(g);
    }
    let _ = n_pairs;
    RdfResult {
        r: r_vals,
        g_r: g_vals,
        n_bins,
        dr,
    }
}
/// Compute running coordination number from an RDF result.
///
/// n(r) = 4π ρ ∫₀ʳ g(r′) r′² dr′
pub fn coordination_number(rdf: &RdfResult, density: f64) -> Vec<f64> {
    let mut cn = Vec::with_capacity(rdf.n_bins);
    let mut running = 0.0;
    for (k, &r) in rdf.r.iter().enumerate() {
        running += 4.0 * PI * density * rdf.g_r[k] * r * r * rdf.dr;
        cn.push(running);
    }
    cn
}
/// Find the first minimum in g(r) to determine first shell cutoff.
pub fn first_shell_cutoff(rdf: &RdfResult) -> Option<f64> {
    let mut found_peak = false;
    let mut _peak_val = 0.0;
    for k in 1..rdf.n_bins.saturating_sub(1) {
        if rdf.g_r[k] > rdf.g_r[k - 1] && rdf.g_r[k] > rdf.g_r[k + 1] && !found_peak {
            found_peak = true;
            _peak_val = rdf.g_r[k];
        }
        if found_peak && rdf.g_r[k] < rdf.g_r[k - 1] && rdf.g_r[k] < rdf.g_r[k + 1] {
            return Some(rdf.r[k]);
        }
    }
    None
}
/// Approximate Voronoi coordination number using nearest-neighbor distances.
///
/// For each ion of type `center_type`, counts how many ions of type
/// `neighbor_type` are within a distance cutoff.
pub fn voronoi_coordination(
    state: &IonState,
    box_len: f64,
    center_type: IonType,
    neighbor_type: IonType,
    cutoff: f64,
) -> Vec<usize> {
    let mut counts = Vec::new();
    for (i, &it) in state.ion_types.iter().enumerate() {
        if it != center_type {
            continue;
        }
        let mut count = 0_usize;
        for (j, &jt) in state.ion_types.iter().enumerate() {
            if i == j || jt != neighbor_type {
                continue;
            }
            let dr = min_image(sub3(state.positions[i], state.positions[j]), box_len);
            let r = norm3(dr);
            if r < cutoff {
                count += 1;
            }
        }
        counts.push(count);
    }
    counts
}
/// Average Voronoi coordination number.
pub fn avg_voronoi_coordination(counts: &[usize]) -> f64 {
    if counts.is_empty() {
        return 0.0;
    }
    counts.iter().sum::<usize>() as f64 / counts.len() as f64
}
/// Compute mean squared displacement between reference and current positions.
pub fn compute_msd(
    ref_positions: &[[f64; 3]],
    cur_positions: &[[f64; 3]],
    ion_types: &[IonType],
    target_type: Option<IonType>,
) -> f64 {
    let n = ref_positions.len().min(cur_positions.len());
    let mut sum = 0.0;
    let mut count = 0_usize;
    for i in 0..n {
        if let Some(tt) = target_type
            && ion_types[i] != tt
        {
            continue;
        }
        let dr = sub3(cur_positions[i], ref_positions[i]);
        sum += dot3(dr, dr);
        count += 1;
    }
    if count > 0 { sum / count as f64 } else { 0.0 }
}
/// Compute diffusion coefficient from MSD using Einstein relation.
///
/// D = MSD / (6 · t)
pub fn diffusion_from_msd(msd: f64, time: f64) -> f64 {
    if time <= 0.0 {
        return 0.0;
    }
    msd / (6.0 * time)
}
/// Compute ion pair correlation function C(t).
///
/// C(t) = <h(0)·h(t)> / <h(0)²>
/// where h(t)=1 if a specific cation-anion pair is within `cutoff` at time t.
pub fn ion_pair_lifetime(
    trajectory: &[Vec<[f64; 3]>],
    ion_types: &[IonType],
    box_len: f64,
    cutoff: f64,
    max_lag: usize,
) -> Vec<f64> {
    let n_frames = trajectory.len();
    let n_ions = ion_types.len();
    if n_frames < 2 || n_ions < 2 {
        return vec![1.0];
    }
    let lag_max = max_lag.min(n_frames);
    let mut c_t = vec![0.0_f64; lag_max];
    let mut norm = 0.0_f64;
    let mut pairs = Vec::new();
    for i in 0..n_ions {
        for j in (i + 1)..n_ions {
            if (ion_types[i] == IonType::Cation && ion_types[j] == IonType::Anion)
                || (ion_types[i] == IonType::Anion && ion_types[j] == IonType::Cation)
            {
                pairs.push((i, j));
            }
        }
    }
    for &(pi, pj) in &pairs {
        for t0 in 0..n_frames {
            let dr0 = min_image(sub3(trajectory[t0][pi], trajectory[t0][pj]), box_len);
            let r0 = norm3(dr0);
            let h0 = if r0 < cutoff { 1.0 } else { 0.0 };
            if h0 < 0.5 {
                continue;
            }
            norm += 1.0;
            for (lag, c_lag) in c_t[..lag_max].iter_mut().enumerate() {
                let t1 = t0 + lag;
                if t1 >= n_frames {
                    break;
                }
                let dr1 = min_image(sub3(trajectory[t1][pi], trajectory[t1][pj]), box_len);
                let r1 = norm3(dr1);
                let h1 = if r1 < cutoff { 1.0 } else { 0.0 };
                *c_lag += h0 * h1;
            }
        }
    }
    if norm > 0.0 {
        for val in &mut c_t {
            *val /= norm;
        }
    }
    c_t
}
/// Estimate ion pair lifetime from the correlation function (1/e decay time).
pub fn pair_lifetime_from_correlation(c_t: &[f64], dt: f64) -> f64 {
    let threshold = 1.0 / std::f64::consts::E;
    for (i, &val) in c_t.iter().enumerate() {
        if val < threshold {
            return i as f64 * dt;
        }
    }
    c_t.len() as f64 * dt
}
/// Compute electric current vector J(t) = Σ q_i · v_i for a single frame.
pub fn electric_current(charges: &[f64], velocities: &[[f64; 3]]) -> [f64; 3] {
    let mut j = [0.0; 3];
    for (i, &q) in charges.iter().enumerate() {
        j[0] += q * velocities[i][0];
        j[1] += q * velocities[i][1];
        j[2] += q * velocities[i][2];
    }
    j
}
/// Compute current autocorrelation function <J(0)·J(t)>.
pub fn current_autocorrelation(current_trajectory: &[[f64; 3]], max_lag: usize) -> Vec<f64> {
    let n = current_trajectory.len();
    let lag_max = max_lag.min(n);
    let mut cac = vec![0.0; lag_max];
    for t0 in 0..n {
        for (lag, c_lag) in cac[..lag_max].iter_mut().enumerate() {
            let t1 = t0 + lag;
            if t1 >= n {
                break;
            }
            *c_lag += dot3(current_trajectory[t0], current_trajectory[t1]);
        }
    }
    for (lag, c_lag) in cac[..lag_max].iter_mut().enumerate() {
        let count = (n - lag) as f64;
        if count > 0.0 {
            *c_lag /= count;
        }
    }
    cac
}
/// Compute conductivity via Green-Kubo from current autocorrelation.
///
/// σ = (1 / 3VkT) ∫₀^∞ <J(0)·J(t)> dt
pub fn green_kubo_conductivity(cac: &[f64], dt: f64, volume: f64, temperature: f64) -> f64 {
    let integral: f64 = cac.iter().sum::<f64>() * dt;
    integral / (3.0 * volume * K_B * temperature)
}
/// Compute Nernst-Einstein conductivity from diffusion coefficients.
///
/// σ_NE = (n_+ q_+² D_+ + n_- q_-² D_-) / (V k_B T)
pub fn nernst_einstein_conductivity(
    n_cations: usize,
    n_anions: usize,
    q_cat: f64,
    q_an: f64,
    d_cat: f64,
    d_an: f64,
    volume: f64,
    temperature: f64,
) -> f64 {
    let num = n_cations as f64 * q_cat * q_cat * d_cat + n_anions as f64 * q_an * q_an * d_an;
    num / (volume * K_B * temperature)
}
/// Compute Haven ratio.
///
/// H_R = σ_NE / σ_GK
/// Values < 1 indicate correlated ion motion.
pub fn haven_ratio(sigma_ne: f64, sigma_gk: f64) -> f64 {
    if sigma_gk.abs() < 1e-30 {
        return 1.0;
    }
    sigma_ne / sigma_gk
}
/// Compute off-diagonal stress tensor element from positions and forces.
///
/// σ_αβ = (1/V) Σ_i \[ m_i v_iα v_iβ + r_iα f_iβ \]
pub fn stress_tensor_element(state: &IonState, volume: f64, alpha: usize, beta: usize) -> f64 {
    let mut sigma = 0.0;
    for (i, &m) in state.masses.iter().enumerate() {
        sigma += m * state.velocities[i][alpha] * state.velocities[i][beta];
        sigma += state.positions[i][alpha] * state.forces[i][beta];
    }
    sigma / volume
}
/// Compute stress autocorrelation function.
pub fn stress_autocorrelation(stress_trajectory: &[f64], max_lag: usize) -> Vec<f64> {
    let n = stress_trajectory.len();
    let lag_max = max_lag.min(n);
    let mut sac = vec![0.0; lag_max];
    for t0 in 0..n {
        for (lag, s_lag) in sac[..lag_max].iter_mut().enumerate() {
            let t1 = t0 + lag;
            if t1 >= n {
                break;
            }
            *s_lag += stress_trajectory[t0] * stress_trajectory[t1];
        }
    }
    for (lag, s_lag) in sac[..lag_max].iter_mut().enumerate() {
        let count = (n - lag) as f64;
        if count > 0.0 {
            *s_lag /= count;
        }
    }
    sac
}
/// Compute viscosity from stress autocorrelation via Green-Kubo.
///
/// η = (V / k_B T) ∫₀^∞ <σ_xy(0)·σ_xy(t)> dt
pub fn green_kubo_viscosity(sac: &[f64], dt: f64, volume: f64, temperature: f64) -> f64 {
    let integral: f64 = sac.iter().sum::<f64>() * dt;
    volume * integral / (K_B * temperature)
}
/// Generate Walden plot data from conductivity and viscosity at various temperatures.
pub fn walden_plot(
    conductivities: &[f64],
    viscosities: &[f64],
    temperatures: &[f64],
    _molar_mass: f64,
    density: f64,
) -> Vec<WaldenPoint> {
    let n = conductivities
        .len()
        .min(viscosities.len())
        .min(temperatures.len());
    let mut points = Vec::with_capacity(n);
    for i in 0..n {
        let sigma = conductivities[i];
        let eta = viscosities[i];
        if sigma <= 0.0 || eta <= 0.0 || density <= 0.0 {
            continue;
        }
        let log_lambda = sigma.log10();
        let log_inv_eta = (1.0 / eta).log10();
        points.push(WaldenPoint {
            log_lambda,
            log_inv_eta,
            temperature: temperatures[i],
        });
    }
    points
}
/// Compute Walden product (Λ · η).
pub fn walden_product(conductivity: f64, viscosity: f64) -> f64 {
    conductivity * viscosity
}
/// Classify ionicity from Walden plot deviation.
///
/// Returns a classification string based on how far below the ideal KCl line.
pub fn walden_classification(log_lambda: f64, log_inv_eta: f64) -> &'static str {
    let deviation = log_lambda - log_inv_eta;
    if deviation > -0.3 {
        "good ionic liquid"
    } else if deviation > -0.7 {
        "poor ionic liquid"
    } else {
        "non-ionic / associated"
    }
}
/// Compute cage radius for a given ion.
///
/// The cage is defined by the nearest `n_neighbors` counter-ions.
pub fn cage_radius(state: &IonState, ion_idx: usize, box_len: f64, n_neighbors: usize) -> f64 {
    let target_type = match state.ion_types[ion_idx] {
        IonType::Cation => IonType::Anion,
        IonType::Anion => IonType::Cation,
    };
    let mut distances: Vec<f64> = Vec::new();
    for (j, &jt) in state.ion_types.iter().enumerate() {
        if j == ion_idx || jt != target_type {
            continue;
        }
        let dr = min_image(sub3(state.positions[ion_idx], state.positions[j]), box_len);
        distances.push(norm3(dr));
    }
    distances.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let k = n_neighbors.min(distances.len());
    if k == 0 {
        return 0.0;
    }
    distances[..k].iter().sum::<f64>() / k as f64
}
/// Compute average cage dynamics for all ions of a given type.
pub fn average_cage_dynamics(
    state: &IonState,
    box_len: f64,
    target_type: IonType,
    n_neighbors: usize,
) -> f64 {
    let mut sum = 0.0;
    let mut count = 0;
    for (i, &it) in state.ion_types.iter().enumerate() {
        if it != target_type {
            continue;
        }
        sum += cage_radius(state, i, box_len, n_neighbors);
        count += 1;
    }
    if count > 0 { sum / count as f64 } else { 0.0 }
}
/// Compute cage rattling amplitude from short-time MSD.
///
/// Rattling amplitude ≈ sqrt(MSD at plateau) before diffusive regime.
pub fn rattling_amplitude(msd_values: &[f64], _dt: f64) -> f64 {
    if msd_values.len() < 3 {
        return 0.0;
    }
    let mut min_slope = f64::MAX;
    let mut plateau_msd = msd_values[0];
    for i in 1..msd_values.len() {
        let slope = msd_values[i] - msd_values[i - 1];
        if slope < min_slope {
            min_slope = slope;
            plateau_msd = msd_values[i];
        }
    }
    plateau_msd.sqrt()
}
/// Detect glass transition temperature from density vs temperature data.
///
/// Uses a two-line intersection method: fits linear segments above and
/// below T_g and finds their crossing.
pub fn detect_glass_transition(temperatures: &[f64], densities: &[f64]) -> Option<f64> {
    let n = temperatures.len().min(densities.len());
    if n < 6 {
        return None;
    }
    let mut best_tg = 0.0;
    let mut best_error = f64::MAX;
    for split in 3..(n - 3) {
        let (a1, b1) = linear_fit(&temperatures[..split], &densities[..split]);
        let (a2, b2) = linear_fit(&temperatures[split..], &densities[split..]);
        let mut error = 0.0;
        for (i, &temp) in temperatures[..split].iter().enumerate() {
            let predicted = a1 * temp + b1;
            error += (predicted - densities[i]).powi(2);
        }
        for i in split..n {
            let predicted = a2 * temperatures[i] + b2;
            error += (predicted - densities[i]).powi(2);
        }
        if error < best_error {
            best_error = error;
            let denom = a1 - a2;
            if denom.abs() > 1e-15 {
                best_tg = (b2 - b1) / denom;
            }
        }
    }
    if best_tg > 0.0 { Some(best_tg) } else { None }
}
/// Simple linear least squares fit: y = a·x + b.
pub(super) fn linear_fit(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len() as f64;
    if n < 2.0 {
        return (0.0, y.first().copied().unwrap_or(0.0));
    }
    let sx: f64 = x.iter().sum();
    let sy: f64 = y.iter().sum();
    let sxx: f64 = x.iter().map(|&v| v * v).sum();
    let sxy: f64 = x.iter().zip(y.iter()).map(|(&xi, &yi)| xi * yi).sum();
    let denom = n * sxx - sx * sx;
    if denom.abs() < 1e-30 {
        return (0.0, sy / n);
    }
    let a = (n * sxy - sx * sy) / denom;
    let b = (sy - a * sx) / n;
    (a, b)
}
/// Initialize a simple cubic lattice of alternating cations and anions.
pub fn init_cubic_lattice(
    cation: &IonModel,
    anion: &IonModel,
    n_pairs: usize,
    box_len: f64,
) -> IonState {
    let mut state = IonState::new();
    let n_side = ((2 * n_pairs) as f64).cbrt().ceil() as usize;
    let spacing = box_len / n_side as f64;
    let mut count = 0_usize;
    'outer: for ix in 0..n_side {
        for iy in 0..n_side {
            for iz in 0..n_side {
                if count >= 2 * n_pairs {
                    break 'outer;
                }
                let pos = [
                    (ix as f64 + 0.5) * spacing,
                    (iy as f64 + 0.5) * spacing,
                    (iz as f64 + 0.5) * spacing,
                ];
                let model = if count.is_multiple_of(2) {
                    cation
                } else {
                    anion
                };
                state.add_ion(model, pos, [0.0; 3]);
                count += 1;
            }
        }
    }
    state
}
/// Assign Maxwell-Boltzmann velocities to all ions.
pub fn assign_velocities(state: &mut IonState, temperature: f64) {
    use rand::RngExt;
    let mut rng = rand::rng();
    for (i, vel) in state.velocities.iter_mut().enumerate() {
        let sigma = (R_GAS * temperature / (state.masses[i] * 1000.0)).sqrt();
        for v in vel.iter_mut().take(3) {
            let u1: f64 = rng.random_range(1e-10_f64..1.0_f64);
            let u2: f64 = rng.random_range(0.0_f64..2.0 * PI);
            *v = sigma * (-2.0 * u1.ln()).sqrt() * u2.cos();
        }
    }
    let mut v_com = [0.0; 3];
    let mut total_mass = 0.0;
    for (i, &vel) in state.velocities.iter().enumerate() {
        v_com = add3(v_com, scale3(vel, state.masses[i]));
        total_mass += state.masses[i];
    }
    if total_mass > 0.0 {
        v_com = scale3(v_com, 1.0 / total_mass);
        for v in &mut state.velocities {
            *v = sub3(*v, v_com);
        }
    }
}
/// Vogel-Fulcher-Tammann viscosity at temperature T.
///
/// η(T) = η₀ · exp(B / (T - T₀))
pub fn vft_viscosity(eta_0: f64, b: f64, t_0: f64, t: f64) -> f64 {
    if (t - t_0).abs() < 1e-10 {
        return f64::MAX;
    }
    eta_0 * (b / (t - t_0)).exp()
}
/// Compute fragility index from VFT parameters.
///
/// m = B · T_g / (T_g - T₀)²  (at T = T_g)
pub fn fragility_index(b: f64, t_g: f64, t_0: f64) -> f64 {
    let denom = (t_g - t_0).powi(2);
    if denom.abs() < 1e-10 {
        return 0.0;
    }
    b * t_g / denom
}
/// Classify fragility as "strong" or "fragile".
pub fn fragility_class(m: f64) -> &'static str {
    if m < 40.0 { "strong" } else { "fragile" }
}
/// Compute velocity autocorrelation function for a single ion type.
pub fn velocity_autocorrelation(
    velocity_trajectory: &[Vec<[f64; 3]>],
    ion_types: &[IonType],
    target_type: IonType,
    max_lag: usize,
) -> Vec<f64> {
    let n_frames = velocity_trajectory.len();
    let lag_max = max_lag.min(n_frames);
    let mut vac = vec![0.0; lag_max];
    let mut count_per_lag = vec![0_usize; lag_max];
    for (i, &ion_type) in ion_types.iter().enumerate() {
        if ion_type != target_type {
            continue;
        }
        for t0 in 0..n_frames {
            for lag in 0..lag_max {
                let t1 = t0 + lag;
                if t1 >= n_frames {
                    break;
                }
                vac[lag] += dot3(velocity_trajectory[t0][i], velocity_trajectory[t1][i]);
                count_per_lag[lag] += 1;
            }
        }
    }
    for lag in 0..lag_max {
        if count_per_lag[lag] > 0 {
            vac[lag] /= count_per_lag[lag] as f64;
        }
    }
    vac
}
/// Compute self-diffusion coefficient from velocity autocorrelation (Green-Kubo).
///
/// D = (1/3) ∫₀^∞ <v(0)·v(t)> dt
pub fn diffusion_from_vac(vac: &[f64], dt: f64) -> f64 {
    let integral: f64 = vac.iter().sum::<f64>() * dt;
    integral / 3.0
}
/// Structure factor S(q) from positions (for a single q-vector magnitude).
pub fn structure_factor(
    positions: &[[f64; 3]],
    ion_types: &[IonType],
    target_type: IonType,
    q_mag: f64,
) -> f64 {
    let mut cos_sum = 0.0;
    let mut sin_sum = 0.0;
    let mut count = 0_usize;
    for (i, pos) in positions.iter().enumerate() {
        if ion_types[i] != target_type {
            continue;
        }
        let qr = q_mag * pos[0];
        cos_sum += qr.cos();
        sin_sum += qr.sin();
        count += 1;
    }
    if count == 0 {
        return 0.0;
    }
    (cos_sum * cos_sum + sin_sum * sin_sum) / count as f64
}
/// Non-Gaussian parameter α₂(t) = (3/5)·<r⁴>/`r²`² − 1.
pub fn non_gaussian_parameter(displacements: &[[f64; 3]]) -> f64 {
    if displacements.is_empty() {
        return 0.0;
    }
    let mut r2_sum = 0.0;
    let mut r4_sum = 0.0;
    for dr in displacements {
        let r2 = dot3(*dr, *dr);
        r2_sum += r2;
        r4_sum += r2 * r2;
    }
    let n = displacements.len() as f64;
    let r2_avg = r2_sum / n;
    let r4_avg = r4_sum / n;
    if r2_avg.abs() < 1e-30 {
        return 0.0;
    }
    (3.0 / 5.0) * r4_avg / (r2_avg * r2_avg) - 1.0
}
/// Compute density of the system (g/cm³).
pub fn compute_density(state: &IonState, box_len: f64) -> f64 {
    let total_mass: f64 = state.masses.iter().sum();
    let volume_angstrom3 = box_len * box_len * box_len;
    total_mass * 1.660_54e-24 / (volume_angstrom3 * 1e-24)
}
/// Compute charge density profile along one axis.
pub fn charge_density_profile(
    state: &IonState,
    box_len: f64,
    axis: usize,
    n_bins: usize,
) -> (Vec<f64>, Vec<f64>) {
    let dz = box_len / n_bins as f64;
    let mut profile = vec![0.0; n_bins];
    let volume_slice = box_len * box_len * dz;
    for (i, pos) in state.positions.iter().enumerate() {
        let z = pos[axis];
        let bin = ((z / dz) as usize).min(n_bins - 1);
        profile[bin] += state.charges[i];
    }
    for val in &mut profile {
        *val /= volume_slice;
    }
    let positions: Vec<f64> = (0..n_bins).map(|k| (k as f64 + 0.5) * dz).collect();
    (positions, profile)
}
/// Compute number density profile along one axis for a specific ion type.
pub fn number_density_profile(
    state: &IonState,
    box_len: f64,
    axis: usize,
    n_bins: usize,
    target_type: IonType,
) -> (Vec<f64>, Vec<f64>) {
    let dz = box_len / n_bins as f64;
    let mut profile = vec![0.0; n_bins];
    let volume_slice = box_len * box_len * dz;
    for (i, &it) in state.ion_types.iter().enumerate() {
        if it != target_type {
            continue;
        }
        let z = state.positions[i][axis];
        let bin = ((z / dz) as usize).min(n_bins - 1);
        profile[bin] += 1.0;
    }
    for val in &mut profile {
        *val /= volume_slice;
    }
    let positions: Vec<f64> = (0..n_bins).map(|k| (k as f64 + 0.5) * dz).collect();
    (positions, profile)
}
