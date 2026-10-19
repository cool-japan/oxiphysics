//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{
    HeatFluxAcf, MeanSquaredDisplacement, RadialDistributionFunction, VelocityAutocorrelation,
};

/// Boltzmann constant in kJ/mol/K
pub(super) const KB: f64 = 8.314e-3;
/// Compute the radial distribution function g(r) for a set of atomic positions
/// inside a cubic periodic box.
///
/// # Arguments
/// * `positions` – Atomic positions (Cartesian, same units as `box_size`).
/// * `box_size`  – Side length of the cubic simulation box.
/// * `n_bins`    – Number of histogram bins.
/// * `cutoff`    – Maximum distance to include (must be ≤ box_size / 2 for
///   correct PBC treatment).
///
/// # Returns
/// A [`RadialDistributionFunction`] with `n_bins` bins spanning `[0, cutoff)`.
pub fn compute_rdf(
    positions: &[[f64; 3]],
    box_size: f64,
    n_bins: usize,
    cutoff: f64,
) -> RadialDistributionFunction {
    let n = positions.len();
    if n < 2 || n_bins == 0 || box_size <= 0.0 || cutoff <= 0.0 {
        return RadialDistributionFunction {
            bins: vec![0.0; n_bins],
            g_r: vec![0.0; n_bins],
        };
    }
    let bin_width = cutoff / n_bins as f64;
    let mut hist = vec![0u64; n_bins];
    for i in 0..n {
        for j in (i + 1)..n {
            let r = dist_pbc(positions[i], positions[j], box_size);
            if r < cutoff {
                let bin = (r / bin_width) as usize;
                let bin = bin.min(n_bins - 1);
                hist[bin] += 1;
            }
        }
    }
    let volume = box_size * box_size * box_size;
    let n_pairs = (n * (n - 1)) as f64 / 2.0;
    let bins: Vec<f64> = (0..n_bins).map(|i| i as f64 * bin_width).collect();
    let g_r: Vec<f64> = hist
        .iter()
        .enumerate()
        .map(|(i, &count)| {
            let r_low = i as f64 * bin_width;
            let r_high = r_low + bin_width;
            let v_shell = 4.0 / 3.0 * PI * (r_high.powi(3) - r_low.powi(3));
            let expected = n_pairs * v_shell / volume;
            if expected > 0.0 {
                count as f64 / expected
            } else {
                0.0
            }
        })
        .collect();
    RadialDistributionFunction { bins, g_r }
}
/// Compute the mean squared displacement from a set of particle trajectories.
///
/// # Arguments
/// * `trajectories` – One trajectory per particle; each trajectory is a
///   `Vec<[f64; 3]>` of positions at successive time steps.
///   All trajectories must have the same length.
/// * `dt`           – Time step between consecutive frames.
///
/// # Returns
/// A [`MeanSquaredDisplacement`] with one entry per frame (lag time).
pub fn compute_msd(trajectories: &[Vec<[f64; 3]>], dt: f64) -> MeanSquaredDisplacement {
    if trajectories.is_empty() {
        return MeanSquaredDisplacement {
            times: vec![],
            msd: vec![],
        };
    }
    let n_frames = trajectories[0].len();
    let _n_particles = trajectories.len();
    if n_frames == 0 {
        return MeanSquaredDisplacement {
            times: vec![],
            msd: vec![],
        };
    }
    let mut times = Vec::with_capacity(n_frames);
    let mut msd_vals = Vec::with_capacity(n_frames);
    for lag in 0..n_frames {
        times.push(lag as f64 * dt);
        let n_origins = n_frames - lag;
        let mut sum_sq = 0.0;
        let mut count = 0u64;
        for traj in trajectories {
            if traj.len() < n_frames {
                continue;
            }
            for t0 in 0..n_origins {
                let r0 = traj[t0];
                let r1 = traj[t0 + lag];
                let dx = r1[0] - r0[0];
                let dy = r1[1] - r0[1];
                let dz = r1[2] - r0[2];
                sum_sq += dx * dx + dy * dy + dz * dz;
                count += 1;
            }
        }
        msd_vals.push(if count > 0 {
            sum_sq / count as f64
        } else {
            0.0
        });
    }
    MeanSquaredDisplacement {
        times,
        msd: msd_vals,
    }
}
/// Compute the normalised velocity autocorrelation function.
///
/// # Arguments
/// * `velocities` – One velocity time series per particle.
///   `velocities[particle][frame]` = `[vx, vy, vz]`.
///   All time series must have the same length.
/// * `dt`         – Time step between frames.
///
/// # Returns
/// A [`VelocityAutocorrelation`] normalised so that `C(0) = 1`.
pub fn compute_vacf(velocities: &[Vec<[f64; 3]>], dt: f64) -> VelocityAutocorrelation {
    if velocities.is_empty() {
        return VelocityAutocorrelation {
            times: vec![],
            vacf: vec![],
        };
    }
    let n_frames = velocities[0].len();
    let _n_particles = velocities.len();
    if n_frames == 0 {
        return VelocityAutocorrelation {
            times: vec![],
            vacf: vec![],
        };
    }
    let mut c_raw = vec![0.0f64; n_frames];
    let mut counts = vec![0u64; n_frames];
    for v_traj in velocities {
        if v_traj.len() < n_frames {
            continue;
        }
        for lag in 0..n_frames {
            let n_origins = n_frames - lag;
            for t0 in 0..n_origins {
                let v0 = v_traj[t0];
                let v1 = v_traj[t0 + lag];
                c_raw[lag] += v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
                counts[lag] += 1;
            }
        }
    }
    for (c, cnt) in c_raw.iter_mut().zip(counts.iter()) {
        if *cnt > 0 {
            *c /= *cnt as f64;
        }
    }
    let c0 = c_raw[0];
    let vacf: Vec<f64> = if c0.abs() > 1e-30 {
        c_raw.iter().map(|&c| c / c0).collect()
    } else {
        vec![0.0; n_frames]
    };
    let times: Vec<f64> = (0..n_frames).map(|i| i as f64 * dt).collect();
    VelocityAutocorrelation { times, vacf }
}
/// Compute the instantaneous temperature from velocities and masses.
///
/// T = (2 · E_kin) / (N_dof · k_B)
///
/// where N_dof = 3 · N (no constraints assumed) and
/// E_kin = ½ Σ_i m_i |v_i|².
///
/// # Arguments
/// * `velocities` – Velocities in Å/ps (or any consistent unit system).
/// * `masses`     – Masses in atomic mass units (amu = g/mol).
///
/// # Returns
/// Temperature in Kelvin (using kB = 8.314×10⁻³ kJ/mol/K,
/// assuming velocities in Å/ps and masses in amu).
pub fn compute_temperature(velocities: &[[f64; 3]], masses: &[f64]) -> f64 {
    if velocities.is_empty() || masses.is_empty() {
        return 0.0;
    }
    let n = velocities.len().min(masses.len());
    if n == 0 {
        return 0.0;
    }
    let e_kin = compute_kinetic_energy(velocities, masses);
    let n_dof = 3 * n;
    2.0 * e_kin / (n_dof as f64 * KB)
}
/// Compute the total kinetic energy.
///
/// E_kin = ½ Σ_i m_i (v_ix² + v_iy² + v_iz²)
///
/// Units: if masses in amu and velocities in Å/ps, the result is in
/// kJ/mol (using the conversion 1 amu·Å²/ps² = 0.01 kJ/mol, so
/// this function returns kJ/mol when the caller uses MD-standard units).
///
/// For simplicity this implementation returns the raw value in
/// consistent units; the caller is responsible for unit interpretation.
pub fn compute_kinetic_energy(velocities: &[[f64; 3]], masses: &[f64]) -> f64 {
    let n = velocities.len().min(masses.len());
    let mut e_kin = 0.0;
    for i in 0..n {
        let v = velocities[i];
        let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        e_kin += 0.5 * masses[i] * v2;
    }
    e_kin
}
/// Compute the pressure tensor from virial and kinetic energy contributions.
///
/// P_αβ = (N·kB·T·δ_αβ + W_αβ) / V
///
/// where W is the virial tensor.  Here we only compute the diagonal
/// (isotropic) pressure contribution from kinetic energy:
///
/// P = (2·E_kin) / (3·V) + P_virial_correction
///
/// # Arguments
/// * `e_kin`  – Total kinetic energy (kJ/mol).
/// * `virial` – Scalar virial sum Σ r_ij · f_ij (kJ/mol, isotropic).
/// * `volume` – Box volume (nm³).
///
/// # Returns
/// Pressure in bar (1 kJ/mol/nm³ = 16.6054 bar).
pub fn compute_pressure(e_kin: f64, virial: f64, volume: f64) -> f64 {
    pub(super) const KJMOL_NM3_TO_BAR: f64 = 16.6054;
    if volume <= 0.0 {
        return 0.0;
    }
    (2.0 * e_kin + virial) / (3.0 * volume) * KJMOL_NM3_TO_BAR
}
/// Compute the centre of mass of a set of atoms.
///
/// # Arguments
/// * `positions` – Atomic positions.
/// * `masses`    – Atomic masses.
///
/// # Returns
/// Centre of mass position `[x, y, z]`, or `[0,0,0]` if empty.
pub fn centre_of_mass(positions: &[[f64; 3]], masses: &[f64]) -> [f64; 3] {
    let n = positions.len().min(masses.len());
    if n == 0 {
        return [0.0; 3];
    }
    let total_mass: f64 = masses[..n].iter().sum();
    if total_mass < 1e-30 {
        return [0.0; 3];
    }
    let mut com = [0.0f64; 3];
    for i in 0..n {
        let m = masses[i];
        com[0] += m * positions[i][0];
        com[1] += m * positions[i][1];
        com[2] += m * positions[i][2];
    }
    [
        com[0] / total_mass,
        com[1] / total_mass,
        com[2] / total_mass,
    ]
}
/// Compute the root-mean-square deviation (RMSD) between two sets of positions.
///
/// RMSD = sqrt\[ (1/N) Σ_i |r_i - r_i_ref|² \]
///
/// Does NOT perform any alignment (rotation/translation).
pub fn compute_rmsd(positions: &[[f64; 3]], reference: &[[f64; 3]]) -> f64 {
    let n = positions.len().min(reference.len());
    if n == 0 {
        return 0.0;
    }
    let sum_sq: f64 = (0..n)
        .map(|i| {
            let dx = positions[i][0] - reference[i][0];
            let dy = positions[i][1] - reference[i][1];
            let dz = positions[i][2] - reference[i][2];
            dx * dx + dy * dy + dz * dz
        })
        .sum();
    (sum_sq / n as f64).sqrt()
}
/// Compute the radius of gyration for a set of atoms.
///
/// Rg = sqrt\[ Σ_i m_i |r_i - r_com|² / Σ_i m_i \]
pub fn compute_radius_of_gyration(positions: &[[f64; 3]], masses: &[f64]) -> f64 {
    let com = centre_of_mass(positions, masses);
    let n = positions.len().min(masses.len());
    let total_mass: f64 = masses[..n].iter().sum();
    if total_mass < 1e-30 {
        return 0.0;
    }
    let sum_sq: f64 = (0..n)
        .map(|i| {
            let dx = positions[i][0] - com[0];
            let dy = positions[i][1] - com[1];
            let dz = positions[i][2] - com[2];
            masses[i] * (dx * dx + dy * dy + dz * dz)
        })
        .sum();
    (sum_sq / total_mass).sqrt()
}
/// Compute the diffusion coefficient directly from a velocity trajectory using
/// the Green-Kubo relation:
///
/// D = (1/3) ∫₀^∞ ⟨v(t)·v(0)⟩ dt
///
/// Integrates the (non-normalised) VACF using the trapezoidal rule.
///
/// # Arguments
/// * `velocities` – `velocities[particle][frame] = [vx, vy, vz]`.
/// * `dt`         – Time step between frames.
///
/// # Returns
/// Diffusion coefficient in length² / time.
pub fn compute_diffusion_coefficient_green_kubo(velocities: &[Vec<[f64; 3]>], dt: f64) -> f64 {
    if velocities.is_empty() {
        return 0.0;
    }
    let n_frames = velocities[0].len();
    if n_frames < 2 {
        return 0.0;
    }
    let mut c_raw = vec![0.0f64; n_frames];
    let mut counts = vec![0u64; n_frames];
    for v_traj in velocities {
        if v_traj.len() < n_frames {
            continue;
        }
        for lag in 0..n_frames {
            let n_origins = n_frames - lag;
            for t0 in 0..n_origins {
                let v0 = v_traj[t0];
                let v1 = v_traj[t0 + lag];
                c_raw[lag] += v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
                counts[lag] += 1;
            }
        }
    }
    for (c, cnt) in c_raw.iter_mut().zip(counts.iter()) {
        if *cnt > 0 {
            *c /= *cnt as f64;
        }
    }
    let integral: f64 = c_raw.windows(2).map(|w| (w[0] + w[1]) * 0.5 * dt).sum();
    integral / 3.0
}
/// Compute the shear viscosity via the Green-Kubo relation from the
/// off-diagonal pressure tensor autocorrelation:
///
/// η = (V / (kB T)) ∫₀^∞ ⟨σ_xy(t) σ_xy(0)⟩ dt
///
/// # Arguments
/// * `stress_acf`  – Time series of one off-diagonal stress component σ_αβ(t).
/// * `dt`          – Time step.
/// * `volume`      – System volume (nm³ or Å³; same units as stress).
/// * `temperature` – Temperature (K).
///
/// # Returns
/// Viscosity in units consistent with stress × time × volume / (kB T).
pub fn compute_viscosity_green_kubo(
    stress_acf: &[f64],
    dt: f64,
    volume: f64,
    temperature: f64,
) -> f64 {
    if stress_acf.is_empty() || temperature < 1e-10 || volume < 1e-30 {
        return 0.0;
    }
    let n = stress_acf.len();
    let mut acf = vec![0.0f64; n];
    for lag in 0..n {
        let n_origins = n - lag;
        for t0 in 0..n_origins {
            acf[lag] += stress_acf[t0] * stress_acf[t0 + lag];
        }
        acf[lag] /= n_origins as f64;
    }
    let integral: f64 = acf.windows(2).map(|w| (w[0] + w[1]) * 0.5 * dt).sum();
    volume / (KB * temperature) * integral
}
/// Compute the full 3×3 pressure tensor from kinetic and virial contributions.
///
/// P_αβ = (Σ_i m_i v_iα v_iβ + W_αβ) / V
///
/// # Arguments
/// * `velocities`    – Atom velocities \[vx, vy, vz\].
/// * `masses`        – Atom masses.
/// * `virial_tensor` – Full virial tensor W_αβ (3×3 matrix, row-major).
/// * `volume`        – Box volume (Å³).
///
/// # Returns
/// Pressure tensor P_αβ in bar (using 1 kJ/mol/Å³ = 16605.4 bar).
pub fn compute_pressure_tensor(
    velocities: &[[f64; 3]],
    masses: &[f64],
    virial_tensor: &[[f64; 3]; 3],
    volume: f64,
) -> [[f64; 3]; 3] {
    pub(super) const KJ_ANG3_TO_BAR: f64 = 16_605.4;
    if volume < 1e-30 {
        return [[0.0; 3]; 3];
    }
    let n = velocities.len().min(masses.len());
    let mut kinetic = [[0.0f64; 3]; 3];
    for i in 0..n {
        let m = masses[i];
        let v = velocities[i];
        for a in 0..3 {
            for b in 0..3 {
                kinetic[a][b] += m * v[a] * v[b];
            }
        }
    }
    let mut pt = [[0.0f64; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            pt[a][b] = (kinetic[a][b] + virial_tensor[a][b]) / volume * KJ_ANG3_TO_BAR;
        }
    }
    pt
}
/// Compute the scalar pressure from the pressure tensor (trace / 3).
pub fn pressure_from_tensor(pt: &[[f64; 3]; 3]) -> f64 {
    (pt[0][0] + pt[1][1] + pt[2][2]) / 3.0
}
/// Compute the normalised heat flux autocorrelation function.
///
/// # Arguments
/// * `heat_flux` – Time series of heat flux vectors J(t) = \[Jx, Jy, Jz\].
/// * `dt`        – Time step between frames.
///
/// # Returns
/// A [`HeatFluxAcf`] normalised so that ACF(0) = 1.
pub fn compute_heat_flux_acf(heat_flux: &[[f64; 3]], dt: f64) -> HeatFluxAcf {
    let n = heat_flux.len();
    if n == 0 {
        return HeatFluxAcf {
            times: vec![],
            acf: vec![],
            j0_sq: 0.0,
        };
    }
    let mut acf_raw = vec![0.0f64; n];
    for lag in 0..n {
        let n_origins = n - lag;
        for t0 in 0..n_origins {
            let j0 = heat_flux[t0];
            let jt = heat_flux[t0 + lag];
            acf_raw[lag] += j0[0] * jt[0] + j0[1] * jt[1] + j0[2] * jt[2];
        }
        acf_raw[lag] /= n_origins as f64;
    }
    let j0_sq = acf_raw[0];
    let times: Vec<f64> = (0..n).map(|i| i as f64 * dt).collect();
    let acf: Vec<f64> = if j0_sq.abs() > 1e-30 {
        acf_raw.iter().map(|&c| c / j0_sq).collect()
    } else {
        vec![0.0; n]
    };
    HeatFluxAcf { times, acf, j0_sq }
}
/// Fully normalised associated Legendre function of the polar angle,
///
/// `P̄_l^m(cos θ) = √((2l+1)/(4π) · (l−m)!/(l+m)!) · P_l^m(cos θ)`,
///
/// including the Condon-Shortley phase `(−1)^m`, so that the complex
/// spherical harmonic is `Y_l^m(θ, φ) = P̄_l^m(cos θ) · e^{imφ}` for
/// `0 ≤ m ≤ l`.  Returns `0.0` for `m > l`.
///
/// `cos_theta` and `sin_theta` (`≥ 0`) are passed separately so that
/// directions close to the poles keep their full precision.  The stable
/// three-term recurrence in `l` at fixed `m` is seeded with
/// `P̄_m^m = −√((2m+1)/(2m)) · sin θ · P̄_{m−1}^{m−1}` (`P̄_0^0 = 1/√(4π)`) and
/// `P̄_{m+1}^m = √(2m+3) · cos θ · P̄_m^m`.
fn normalised_legendre(l: usize, m: usize, cos_theta: f64, sin_theta: f64) -> f64 {
    if m > l {
        return 0.0;
    }
    let mut p_mm = (1.0 / (4.0 * PI)).sqrt();
    for k in 1..=m {
        let kf = k as f64;
        p_mm *= -((2.0 * kf + 1.0) / (2.0 * kf)).sqrt() * sin_theta;
    }
    if l == m {
        return p_mm;
    }
    let mf = m as f64;
    let mut p_lower = p_mm;
    let mut p_upper = (2.0 * mf + 3.0).sqrt() * cos_theta * p_mm;
    for k in (m + 2)..=l {
        let kf = k as f64;
        let scale = ((4.0 * kf * kf - 1.0) / (kf * kf - mf * mf)).sqrt();
        let shrink =
            (((kf - 1.0) * (kf - 1.0) - mf * mf) / (4.0 * (kf - 1.0) * (kf - 1.0) - 1.0)).sqrt();
        let p_next = scale * (cos_theta * p_upper - shrink * p_lower);
        p_lower = p_upper;
        p_upper = p_next;
    }
    p_upper
}

/// Add the real, orthonormal spherical harmonics of degree `l` evaluated in
/// the direction `(θ, φ)` to the accumulator `acc`.
///
/// The direction is given by the cosine and sine (`≥ 0`) of the polar angle and
/// by the azimuth `phi`.  `acc` must hold at least `2l + 1` entries (a shorter
/// slice is left untouched); the components are
///
/// ```text
/// [Y_l0, √2·Re Y_l1, √2·Im Y_l1, …, √2·Re Y_ll, √2·Im Y_ll]
/// ```
///
/// The `√2` makes the basis orthonormal on the sphere: since
/// `|Q_{l,−m}| = |Q_{l,m}|` for every real bond distribution, the pair
/// `(Re, Im)` of order `m > 0` carries the weight of both `±m` and
/// `Σ_m |Q_lm|²` equals the sum of squares of the bond-averaged components.
/// For a single direction the addition theorem gives
/// `Σ_k y_k² = (2l + 1)/(4π)`.
pub(super) fn accumulate_real_harmonics(
    acc: &mut [f64],
    l: usize,
    cos_theta: f64,
    sin_theta: f64,
    phi: f64,
) {
    let Some(acc) = acc.get_mut(..2 * l + 1) else {
        return;
    };
    acc[0] += normalised_legendre(l, 0, cos_theta, sin_theta);
    for m in 1..=l {
        let amplitude = std::f64::consts::SQRT_2 * normalised_legendre(l, m, cos_theta, sin_theta);
        let angle = m as f64 * phi;
        acc[2 * m - 1] += amplitude * angle.cos();
        acc[2 * m] += amplitude * angle.sin();
    }
}

/// Steinhardt order parameter `Q_l = √(4π/(2l+1) · Σ_m |Q_lm|²)` from the
/// bond sums of the real harmonic components of degree `l`.
fn steinhardt_order_parameter(component_sums: &[f64], n_bonds: f64, l: usize) -> f64 {
    let sum_sq: f64 = component_sums
        .iter()
        .map(|&q| {
            let mean = q / n_bonds;
            mean * mean
        })
        .sum();
    (4.0 * PI / (2 * l + 1) as f64 * sum_sq).sqrt()
}

/// Compute global Steinhardt bond-orientational order parameters Q4 and Q6.
///
/// `Q_l = sqrt( (4π/(2l+1)) · Σ_{m=−l}^{l} |Q_lm|² )` with
/// `Q_lm = ⟨Y_lm(r̂_ij)⟩`, where the average is over all bonds `(i, j)` with
/// `|r_ij| < cutoff` (each pair counted once).  The sum over `m` is evaluated
/// with the real orthonormal spherical harmonics of degree `l` (the `√2`-weighted
/// `Re`/`Im` pairs of the complex `Y_lm`, orthonormal on the sphere).
///
/// `Q_l` is invariant under rotations and translations of the configuration.
/// Reference values for the bulk lattices (nearest-neighbour bonds) are
/// FCC `Q4 = 0.190941`, `Q6 = 0.574524`; a single bond has `Q4 = Q6 = 1`.
///
/// # Arguments
/// * `positions` – Atomic positions.
/// * `cutoff`    – Neighbour cutoff distance.
///
/// # Returns
/// `(Q4, Q6)` — both lie in \[0, 1\] (up to floating-point rounding); `(0, 0)`
/// if there are fewer than two atoms or no bond.
pub fn compute_bond_order_params(positions: &[[f64; 3]], cutoff: f64) -> (f64, f64) {
    let n = positions.len();
    if n < 2 {
        return (0.0, 0.0);
    }
    let mut qlm4 = [0.0f64; 9];
    let mut qlm6 = [0.0f64; 13];
    let mut n_bonds = 0u64;
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[j][0] - positions[i][0];
            let dy = positions[j][1] - positions[i][1];
            let dz = positions[j][2] - positions[i][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r >= cutoff || r < 1e-14 {
                continue;
            }
            let cos_theta = (dz / r).clamp(-1.0, 1.0);
            let sin_theta = dx.hypot(dy) / r;
            let phi = dy.atan2(dx);
            accumulate_real_harmonics(&mut qlm4, 4, cos_theta, sin_theta, phi);
            accumulate_real_harmonics(&mut qlm6, 6, cos_theta, sin_theta, phi);
            n_bonds += 1;
        }
    }
    if n_bonds == 0 {
        return (0.0, 0.0);
    }
    let nb = n_bonds as f64;
    (
        steinhardt_order_parameter(&qlm4, nb, 4),
        steinhardt_order_parameter(&qlm6, nb, 6),
    )
}
/// Compute the running coordination number N(r) from a radial distribution function.
///
/// N(r) = ρ_0 * ∫₀ʳ g(r') * 4πr'² dr'
///
/// where ρ_0 = N / V is the average number density.
///
/// # Arguments
/// * `rdf`      – Pre-computed RDF.
/// * `n_atoms`  – Total number of atoms.
/// * `box_size` – Cubic box length.
///
/// # Returns
/// A vector of length `rdf.n_bins()` with the running coordination number.
pub fn running_coordination_number(
    rdf: &RadialDistributionFunction,
    n_atoms: usize,
    box_size: f64,
) -> Vec<f64> {
    let n_bins = rdf.n_bins();
    if n_bins == 0 || box_size <= 0.0 || n_atoms == 0 {
        return vec![0.0; n_bins];
    }
    let bw = rdf.bin_width();
    let volume = box_size * box_size * box_size;
    let rho = n_atoms as f64 / volume;
    let mut coord = vec![0.0; n_bins];
    let mut running = 0.0;
    for (i, (&_bin, &g)) in rdf.bins.iter().zip(rdf.g_r.iter()).enumerate() {
        let r_lo = i as f64 * bw;
        let r_hi = r_lo + bw;
        let v_shell = 4.0 / 3.0 * PI * (r_hi.powi(3) - r_lo.powi(3));
        running += rho * g * v_shell;
        coord[i] = running;
    }
    coord
}
/// Compute the diffusion coefficient directly from MSD data using a linear fit.
///
/// D = slope / 6  (3D Einstein relation: MSD = 6 D t)
///
/// Uses the second half of the MSD curve to avoid ballistic regime.
pub fn diffusion_from_msd_slope(msd: &MeanSquaredDisplacement) -> f64 {
    msd.diffusion_coefficient()
}
/// Compute the power spectrum of the velocity autocorrelation function via DFT.
///
/// The phonon density of states g(ω) ∝ ∫ VACF(t) * exp(-i ω t) dt.
///
/// We compute the single-sided (ω ≥ 0) spectrum using the discrete cosine
/// transform (real DFT of even signal):
///
/// ```text
/// S(k) = Σ_{n=0}^{N-1} VACF(n) * cos(2π k n / N)
/// ```
///
/// # Arguments
/// * `vacf`  – A [`VelocityAutocorrelation`] (normalised or raw, both work).
/// * `dt`    – Time step (ps).
///
/// # Returns
/// `(frequencies, spectrum)` where `frequencies[k]` is in THz (= 1/ps) and
/// `spectrum[k]` is the (real, positive) spectral amplitude.
pub fn power_spectrum_from_vacf(vacf: &VelocityAutocorrelation, dt: f64) -> (Vec<f64>, Vec<f64>) {
    let n = vacf.vacf.len();
    if n == 0 || dt <= 0.0 {
        return (vec![], vec![]);
    }
    let n_half = n / 2 + 1;
    let mut freqs = Vec::with_capacity(n_half);
    let mut spectrum = Vec::with_capacity(n_half);
    for k in 0..n_half {
        let freq = k as f64 / (n as f64 * dt);
        let mut re = 0.0_f64;
        for m in 0..n {
            re += vacf.vacf[m] * (2.0 * PI * k as f64 * m as f64 / n as f64).cos();
        }
        freqs.push(freq);
        spectrum.push(re.abs() * dt);
    }
    (freqs, spectrum)
}
/// Compute a simple pair-distance histogram (unnormalized count per bin).
///
/// Unlike `compute_rdf`, this returns raw counts without normalisation.
///
/// # Arguments
/// * `positions` – Atom positions.
/// * `cutoff`    – Maximum distance to histogram.
/// * `n_bins`    – Number of histogram bins.
///
/// # Returns
/// `(bin_edges, counts)` — bin_edges has length `n_bins`, counts same length.
pub fn pair_distance_histogram(
    positions: &[[f64; 3]],
    cutoff: f64,
    n_bins: usize,
) -> (Vec<f64>, Vec<u64>) {
    if n_bins == 0 || cutoff <= 0.0 {
        return (vec![], vec![]);
    }
    let bw = cutoff / n_bins as f64;
    let bins: Vec<f64> = (0..n_bins).map(|i| i as f64 * bw).collect();
    let mut counts = vec![0u64; n_bins];
    let n = positions.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[j][0] - positions[i][0];
            let dy = positions[j][1] - positions[i][1];
            let dz = positions[j][2] - positions[i][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < cutoff {
                let b = ((r / bw) as usize).min(n_bins - 1);
                counts[b] += 1;
            }
        }
    }
    (bins, counts)
}
/// Lindemann melting criterion.
///
/// The Lindemann parameter is the ratio of the root-mean-square displacement
/// (from equilibrium) to the mean nearest-neighbour distance:
///
/// L = sqrt(⟨u²⟩) / d_nn
///
/// A value L > 0.1–0.15 is typically associated with melting.
///
/// # Arguments
/// * `positions_t0`   – Reference (equilibrium) positions.
/// * `positions_t`    – Current positions.
/// * `nearest_cutoff` – Cutoff for nearest-neighbour distance estimation.
///
/// # Returns
/// Lindemann parameter (dimensionless).
pub fn lindemann_parameter(
    positions_t0: &[[f64; 3]],
    positions_t: &[[f64; 3]],
    nearest_cutoff: f64,
) -> f64 {
    let n = positions_t0.len().min(positions_t.len());
    if n == 0 {
        return 0.0;
    }
    let msd: f64 = (0..n)
        .map(|i| {
            let dx = positions_t[i][0] - positions_t0[i][0];
            let dy = positions_t[i][1] - positions_t0[i][1];
            let dz = positions_t[i][2] - positions_t0[i][2];
            dx * dx + dy * dy + dz * dz
        })
        .sum::<f64>()
        / n as f64;
    let mut d_nn_sum = 0.0_f64;
    let mut d_nn_count = 0_u64;
    for i in 0..n {
        let mut min_r = f64::INFINITY;
        for j in 0..n {
            if i == j {
                continue;
            }
            let dx = positions_t0[j][0] - positions_t0[i][0];
            let dy = positions_t0[j][1] - positions_t0[i][1];
            let dz = positions_t0[j][2] - positions_t0[i][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < nearest_cutoff && r < min_r {
                min_r = r;
            }
        }
        if min_r.is_finite() {
            d_nn_sum += min_r;
            d_nn_count += 1;
        }
    }
    if d_nn_count == 0 {
        return 0.0;
    }
    let d_nn = d_nn_sum / d_nn_count as f64;
    if d_nn < 1e-30 {
        return 0.0;
    }
    msd.sqrt() / d_nn
}
/// Assign each atom to a cluster based on distance cutoff.
///
/// Two atoms belong to the same cluster if their distance is ≤ `cutoff`.
/// Uses union-find to identify connected components.
///
/// # Returns
/// A vector of length `n_atoms` with cluster labels (0-based, consecutive integers).
pub fn cluster_analysis(positions: &[[f64; 3]], cutoff: f64) -> Vec<usize> {
    let n = positions.len();
    if n == 0 {
        return vec![];
    }
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[j][0] - positions[i][0];
            let dy = positions[j][1] - positions[i][1];
            let dz = positions[j][2] - positions[i][2];
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 <= cutoff * cutoff {
                let ri = find(&mut parent, i);
                let rj = find(&mut parent, j);
                if ri != rj {
                    parent[rj] = ri;
                }
            }
        }
    }
    for k in 0..n {
        let r = find(&mut parent, k);
        parent[k] = r;
    }
    let mut map = std::collections::HashMap::new();
    let mut next_id = 0_usize;
    let mut labels = vec![0_usize; n];
    for (k, &p) in parent.iter().enumerate() {
        let id = *map.entry(p).or_insert_with(|| {
            let v = next_id;
            next_id += 1;
            v
        });
        labels[k] = id;
    }
    labels
}
/// Count the number of distinct clusters.
pub fn count_clusters(positions: &[[f64; 3]], cutoff: f64) -> usize {
    let labels = cluster_analysis(positions, cutoff);
    labels
        .iter()
        .cloned()
        .collect::<std::collections::HashSet<_>>()
        .len()
}
/// Return the size of the largest cluster.
pub fn largest_cluster_size(positions: &[[f64; 3]], cutoff: f64) -> usize {
    let labels = cluster_analysis(positions, cutoff);
    let mut counts = std::collections::HashMap::<usize, usize>::new();
    for &l in &labels {
        *counts.entry(l).or_insert(0) += 1;
    }
    counts.values().cloned().max().unwrap_or(0)
}
/// Minimum-image distance between two points in a cubic box.
pub(super) fn dist_pbc(a: [f64; 3], b: [f64; 3], box_size: f64) -> f64 {
    let mut dx = a[0] - b[0];
    let mut dy = a[1] - b[1];
    let mut dz = a[2] - b[2];
    dx -= box_size * (dx / box_size).round();
    dy -= box_size * (dy / box_size).round();
    dz -= box_size * (dz / box_size).round();
    (dx * dx + dy * dy + dz * dz).sqrt()
}
