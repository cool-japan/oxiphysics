//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::{
    BootstrapResult, DihedralFit, HarmonicAngleFit, HarmonicBondFit, LjFitResult, RespFitResult,
};

/// Lennard-Jones 12-6 energy: `4·ε·[(σ/r)¹² − (σ/r)⁶]`
#[allow(dead_code)]
pub fn lj_energy(r: f64, epsilon: f64, sigma: f64) -> f64 {
    let sr = sigma / r;
    let sr6 = sr.powi(6);
    let sr12 = sr6 * sr6;
    4.0 * epsilon * (sr12 - sr6)
}
/// Harmonic energy: `0.5·k·(x − x₀)²`
#[allow(dead_code)]
pub fn harmonic_energy(x: f64, x0: f64, k: f64) -> f64 {
    let dx = x - x0;
    0.5 * k * dx * dx
}
/// OPLS-style single-term torsion energy: `k·(1 + cos(n·φ − δ))`
#[allow(dead_code)]
pub fn torsion_energy(phi: f64, k: f64, n: u32, delta: f64) -> f64 {
    k * (1.0 + (n as f64 * phi - delta).cos())
}
pub(super) fn rmse(predicted: &[f64], target: &[f64]) -> f64 {
    assert_eq!(predicted.len(), target.len());
    let n = predicted.len() as f64;
    let sse: f64 = predicted
        .iter()
        .zip(target.iter())
        .map(|(p, t)| (p - t).powi(2))
        .sum();
    (sse / n).sqrt()
}
/// Fit Lennard-Jones 12-6 parameters (ε, σ) to (distance, energy) data using
/// gradient-descent with backtracking line search (Levenberg-Marquardt style).
///
/// Optimises in log-parameter space (log ε, log σ) for numerical stability.
/// Initial guesses: σ = median(distances), ε from the minimum energy.
/// Returns [`LjFitResult`] with the best-fit parameters and the RMSE.
#[allow(dead_code)]
pub fn fit_lj_parameters(distances: &[f64], energies: &[f64]) -> LjFitResult {
    assert_eq!(distances.len(), energies.len());
    assert!(!distances.is_empty());
    let min_idx = energies
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);
    let r_min_est = distances[min_idx];
    let sig0 = (r_min_est / 2.0_f64.powf(1.0 / 6.0)).max(0.1);
    let eps0 = energies[min_idx].abs().max(0.01);
    let mut u = eps0.ln();
    let mut v = sig0.ln();
    let loss = |u_: f64, v_: f64| -> f64 {
        let e = u_.exp();
        let s = v_.exp();
        let n = distances.len() as f64;
        distances
            .iter()
            .zip(energies.iter())
            .map(|(&r, &y)| (lj_energy(r, e, s) - y).powi(2))
            .sum::<f64>()
            / n
    };
    let max_iter = 20_000_usize;
    let tol = 1e-10_f64;
    let h = 1e-5_f64;
    for _ in 0..max_iter {
        let l0 = loss(u, v);
        let grad_u = (loss(u + h, v) - l0) / h;
        let grad_v = (loss(u, v + h) - l0) / h;
        let mut step = 0.01_f64;
        let mut attempts = 0;
        loop {
            let nu = u - step * grad_u;
            let nv = v - step * grad_v;
            if loss(nu, nv) < l0 || attempts >= 30 {
                break;
            }
            step *= 0.5;
            attempts += 1;
        }
        let new_u = u - step * grad_u;
        let new_v = v - step * grad_v;
        let du = (new_u - u).abs();
        let dv = (new_v - v).abs();
        u = new_u;
        v = new_v;
        if du < tol && dv < tol {
            break;
        }
    }
    let eps = u.exp();
    let sig = v.exp();
    let pred: Vec<f64> = distances.iter().map(|&r| lj_energy(r, eps, sig)).collect();
    let rmse_val = rmse(&pred, energies);
    LjFitResult {
        epsilon: eps,
        sigma: sig,
        rmse: rmse_val,
    }
}
/// Fit a harmonic bond potential `V = 0.5·k·(r − r₀)²` using an analytic
/// parabolic (least-squares) fit.
///
/// The fit solves the normal equations for the parabola `a·r² + b·r + c`
/// and recovers `r0 = -b/(2a)` and `k = 2a`.
#[allow(dead_code)]
pub fn fit_harmonic_bond(distances: &[f64], energies: &[f64]) -> HarmonicBondFit {
    fit_harmonic_generic(distances, energies)
}
/// Fit a harmonic angle potential `V = 0.5·k·(θ − θ₀)²` using an analytic
/// parabolic (least-squares) fit.
#[allow(dead_code)]
pub fn fit_harmonic_angle(angles: &[f64], energies: &[f64]) -> HarmonicAngleFit {
    let r = fit_harmonic_generic(angles, energies);
    HarmonicAngleFit {
        theta0: r.r0,
        k: r.k,
    }
}
/// Internal: analytic least-squares parabolic fit `a·x² + b·x + c`.
/// Returns (r0, k) where r0 = -b/(2a), k = 2a.
pub(super) fn fit_harmonic_generic(xs: &[f64], ys: &[f64]) -> HarmonicBondFit {
    assert_eq!(xs.len(), ys.len());
    assert!(xs.len() >= 3, "Need at least 3 points for parabolic fit");
    let n = xs.len() as f64;
    let s1: f64 = xs.iter().sum();
    let s2: f64 = xs.iter().map(|&x| x * x).sum();
    let s3: f64 = xs.iter().map(|&x| x.powi(3)).sum();
    let s4: f64 = xs.iter().map(|&x| x.powi(4)).sum();
    let t0: f64 = ys.iter().sum();
    let t1: f64 = xs.iter().zip(ys.iter()).map(|(&x, &y)| x * y).sum();
    let t2: f64 = xs.iter().zip(ys.iter()).map(|(&x, &y)| x * x * y).sum();
    let det = s4 * (s2 * n - s1 * s1) - s3 * (s3 * n - s1 * s2) + s2 * (s3 * s1 - s2 * s2);
    let (a, b) = if det.abs() < 1e-30 {
        let r0 = xs.iter().cloned().fold(f64::INFINITY, f64::min);
        (1.0, -2.0 * r0)
    } else {
        let det_a = t2 * (s2 * n - s1 * s1) - s3 * (t1 * n - s1 * t0) + s2 * (t1 * s1 - s2 * t0);
        let det_b = s4 * (t1 * n - s1 * t0) - t2 * (s3 * n - s1 * s2) + s2 * (s3 * t0 - s2 * t1);
        (det_a / det, det_b / det)
    };
    let k = (2.0 * a).max(0.0);
    let r0 = if a.abs() > 1e-30 {
        -b / (2.0 * a)
    } else {
        xs[0]
    };
    HarmonicBondFit { r0, k }
}
/// Fit a single-term OPLS dihedral `V = k·(1 + cos(n·φ − δ))` to energy data.
///
/// Strategy: for each candidate periodicity n ∈ {1,2,3,4,6}, perform a linear
/// least-squares fit in the Fourier basis `{cos(n·φ), sin(n·φ), 1}` and keep
/// the best (lowest RMSE) solution.  Then recover k and δ from the cosine and
/// sine coefficients.
#[allow(dead_code)]
pub fn fit_dihedral(angles: &[f64], energies: &[f64]) -> DihedralFit {
    assert_eq!(angles.len(), energies.len());
    assert!(!angles.is_empty());
    let candidates: &[u32] = &[1, 2, 3, 4, 6];
    let mut best_rmse = f64::INFINITY;
    let mut best = DihedralFit {
        k: 0.0,
        n: 1,
        delta: 0.0,
    };
    for &n in candidates {
        let m = angles.len() as f64;
        let c: Vec<f64> = angles.iter().map(|&p| (n as f64 * p).cos()).collect();
        let s: Vec<f64> = angles.iter().map(|&p| (n as f64 * p).sin()).collect();
        let sc2: f64 = c.iter().map(|&x| x * x).sum();
        let ss2: f64 = s.iter().map(|&x| x * x).sum();
        let scs: f64 = c.iter().zip(s.iter()).map(|(&x, &y)| x * y).sum();
        let sc: f64 = c.iter().sum();
        let ss: f64 = s.iter().sum();
        let scy: f64 = c.iter().zip(energies.iter()).map(|(&x, &y)| x * y).sum();
        let ssy: f64 = s.iter().zip(energies.iter()).map(|(&x, &y)| x * y).sum();
        let sy: f64 = energies.iter().sum();
        let mat = [[sc2, scs, sc], [scs, ss2, ss], [sc, ss, m]];
        let rhs = [scy, ssy, sy];
        let det3 = det3x3(&mat);
        let (coeff_a, coeff_b) = if det3.abs() < 1e-30 {
            (0.0_f64, 0.0_f64)
        } else {
            let ma = [
                [rhs[0], mat[0][1], mat[0][2]],
                [rhs[1], mat[1][1], mat[1][2]],
                [rhs[2], mat[2][1], mat[2][2]],
            ];
            let mb = [
                [mat[0][0], rhs[0], mat[0][2]],
                [mat[1][0], rhs[1], mat[1][2]],
                [mat[2][0], rhs[2], mat[2][2]],
            ];
            (det3x3(&ma) / det3, det3x3(&mb) / det3)
        };
        let k_fit = (coeff_a.powi(2) + coeff_b.powi(2)).sqrt();
        let delta_fit = coeff_b.atan2(coeff_a);
        let pred: Vec<f64> = angles
            .iter()
            .map(|&p| torsion_energy(p, k_fit, n, delta_fit))
            .collect();
        let r = rmse(&pred, energies);
        if r < best_rmse {
            best_rmse = r;
            best = DihedralFit {
                k: k_fit,
                n,
                delta: delta_fit,
            };
        }
    }
    best
}
pub(super) fn det3x3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
/// Fit restrained partial charges to an electrostatic potential map.
///
/// Minimises: Σ_k (V_esp(r_k) - Σ_i q_i/|r_k - R_i|)² + λ Σ_i q_i²
/// subject to Σ_i q_i = Q_total.
///
/// Uses a simple iterative steepest-descent approach with the total-charge
/// constraint enforced via projection.
///
/// # Arguments
/// * `esp_points`    – Grid points where ESP is known: `(r, V)` pairs.
/// * `atom_positions`– Positions of atoms.
/// * `q_total`       – Required total charge.
/// * `restraint`     – RESP hyperbolic restraint strength a (kJ/mol/e²).
/// * `max_iter`      – Maximum number of iterations.
#[allow(dead_code)]
pub fn fit_resp_charges(
    esp_points: &[([f64; 3], f64)],
    atom_positions: &[[f64; 3]],
    q_total: f64,
    restraint: f64,
    max_iter: usize,
) -> RespFitResult {
    let n_atoms = atom_positions.len();
    let n_pts = esp_points.len();
    if n_atoms == 0 || n_pts == 0 {
        return RespFitResult {
            charges: vec![],
            rmse: 0.0,
            net_charge_error: 0.0,
        };
    }
    let inv_r: Vec<Vec<f64>> = esp_points
        .iter()
        .map(|&(r, _)| {
            atom_positions
                .iter()
                .map(|&ri| {
                    let dr = [r[0] - ri[0], r[1] - ri[1], r[2] - ri[2]];
                    let d = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
                    if d < 1e-10 { 0.0 } else { 1.0 / d }
                })
                .collect()
        })
        .collect();
    let v_target: Vec<f64> = esp_points.iter().map(|&(_, v)| v).collect();
    let mut q = vec![q_total / n_atoms as f64; n_atoms];
    let step = 0.001;
    for _ in 0..max_iter {
        let mut grad = vec![0.0f64; n_atoms];
        for k in 0..n_pts {
            let v_pred: f64 = (0..n_atoms).map(|i| q[i] * inv_r[k][i]).sum();
            let res = v_target[k] - v_pred;
            for i in 0..n_atoms {
                grad[i] -= 2.0 * res * inv_r[k][i];
            }
        }
        for i in 0..n_atoms {
            grad[i] += 2.0 * restraint * q[i];
        }
        for i in 0..n_atoms {
            q[i] -= step * grad[i];
        }
        let q_sum: f64 = q.iter().sum();
        let correction = (q_total - q_sum) / n_atoms as f64;
        for qi in q.iter_mut() {
            *qi += correction;
        }
    }
    let mut sse = 0.0f64;
    for k in 0..n_pts {
        let v_pred: f64 = (0..n_atoms).map(|i| q[i] * inv_r[k][i]).sum();
        sse += (v_target[k] - v_pred).powi(2);
    }
    let rmse_val = (sse / n_pts as f64).sqrt();
    let net_charge_error = (q.iter().sum::<f64>() - q_total).abs();
    RespFitResult {
        charges: q,
        rmse: rmse_val,
        net_charge_error,
    }
}
/// Solve a k×k linear system A x = b via Gaussian elimination.
pub(super) fn solve_linear_system(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let k = b.len();
    let mut aug: Vec<Vec<f64>> = (0..k)
        .map(|i| {
            let mut row = a[i].clone();
            row.push(b[i]);
            row
        })
        .collect();
    for col in 0..k {
        let mut max_row = col;
        for row in (col + 1)..k {
            if aug[row][col].abs() > aug[max_row][col].abs() {
                max_row = row;
            }
        }
        aug.swap(col, max_row);
        let pivot = aug[col][col];
        if pivot.abs() < 1e-15 {
            continue;
        }
        for row in (col + 1)..k {
            let factor = aug[row][col] / pivot;
            for j in col..=k {
                let v = aug[col][j];
                aug[row][j] -= factor * v;
            }
        }
    }
    let mut x = vec![0.0f64; k];
    for i in (0..k).rev() {
        let mut sum = aug[i][k];
        for j in (i + 1)..k {
            sum -= aug[i][j] * x[j];
        }
        if aug[i][i].abs() > 1e-15 {
            x[i] = sum / aug[i][i];
        }
    }
    x
}
/// Mean absolute error.
#[allow(dead_code)]
pub fn mae(predicted: &[f64], target: &[f64]) -> f64 {
    assert_eq!(predicted.len(), target.len());
    let n = predicted.len() as f64;
    predicted
        .iter()
        .zip(target.iter())
        .map(|(p, t)| (p - t).abs())
        .sum::<f64>()
        / n
}
/// R² (coefficient of determination).
#[allow(dead_code)]
pub fn r_squared(predicted: &[f64], target: &[f64]) -> f64 {
    assert_eq!(predicted.len(), target.len());
    let n = target.len() as f64;
    let mean_t = target.iter().sum::<f64>() / n;
    let ss_tot: f64 = target.iter().map(|&t| (t - mean_t).powi(2)).sum();
    let ss_res: f64 = predicted
        .iter()
        .zip(target.iter())
        .map(|(p, t)| (p - t).powi(2))
        .sum();
    if ss_tot.abs() < 1e-30 {
        return 1.0;
    }
    1.0 - ss_res / ss_tot
}
/// Maximum absolute error.
#[allow(dead_code)]
pub fn max_error(predicted: &[f64], target: &[f64]) -> f64 {
    assert_eq!(predicted.len(), target.len());
    predicted
        .iter()
        .zip(target.iter())
        .map(|(p, t)| (p - t).abs())
        .fold(0.0_f64, f64::max)
}
/// Perform bootstrap estimation of harmonic bond fit uncertainty.
#[allow(dead_code)]
pub fn bootstrap_harmonic_bond(
    data: &[(f64, f64)],
    n_bootstrap: usize,
    seed: u64,
) -> BootstrapResult {
    let n = data.len();
    if n < 3 {
        return BootstrapResult {
            std_k: 0.0,
            std_r0: 0.0,
            n_bootstrap: 0,
        };
    }
    let mut state = seed ^ 6364136223846793005;
    let mut lcg_next = move || -> usize {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as usize
    };
    let mut ks = Vec::with_capacity(n_bootstrap);
    let mut r0s = Vec::with_capacity(n_bootstrap);
    for _ in 0..n_bootstrap {
        let xs: Vec<f64> = (0..n).map(|_| data[lcg_next() % n].0).collect();
        let ys: Vec<f64> = (0..n).map(|_| data[lcg_next() % n].1).collect();
        if xs.len() < 3 {
            continue;
        }
        let fit = fit_harmonic_bond(&xs, &ys);
        ks.push(fit.k);
        r0s.push(fit.r0);
    }
    let std_k = std_dev(&ks);
    let std_r0 = std_dev(&r0s);
    BootstrapResult {
        std_k,
        std_r0,
        n_bootstrap,
    }
}
pub(super) fn std_dev(vals: &[f64]) -> f64 {
    if vals.len() < 2 {
        return 0.0;
    }
    let n = vals.len() as f64;
    let mean = vals.iter().sum::<f64>() / n;
    let var = vals.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / n;
    var.sqrt()
}
/// Compute k-fold cross-validation RMSE for a harmonic bond model.
///
/// The data are split into `k_folds` contiguous folds; each fold is used once
/// as a validation set while the remaining data train a harmonic fit.  Returns
/// the average RMSE across all folds.
#[allow(dead_code)]
pub fn cross_validate(data: &[(f64, f64)], k_folds: usize) -> f64 {
    assert!(k_folds >= 2, "Need at least 2 folds");
    assert!(data.len() >= k_folds, "Need at least k_folds data points");
    let n = data.len();
    let fold_size = n / k_folds;
    let mut total_rmse = 0.0_f64;
    for fold in 0..k_folds {
        let val_start = fold * fold_size;
        let val_end = if fold == k_folds - 1 {
            n
        } else {
            val_start + fold_size
        };
        let train_xs: Vec<f64> = data[..val_start]
            .iter()
            .chain(&data[val_end..])
            .map(|&(x, _)| x)
            .collect();
        let train_ys: Vec<f64> = data[..val_start]
            .iter()
            .chain(&data[val_end..])
            .map(|&(_, y)| y)
            .collect();
        let val_xs: Vec<f64> = data[val_start..val_end].iter().map(|&(x, _)| x).collect();
        let val_ys: Vec<f64> = data[val_start..val_end].iter().map(|&(_, y)| y).collect();
        if train_xs.len() < 3 {
            continue;
        }
        let fit = fit_harmonic_generic(&train_xs, &train_ys);
        let pred: Vec<f64> = val_xs
            .iter()
            .map(|&x| harmonic_energy(x, fit.r0, fit.k))
            .collect();
        total_rmse += rmse(&pred, &val_ys);
    }
    total_rmse / k_folds as f64
}
/// Compute per-atom RMSD between two coordinate sets.
///
/// Both slices must have the same length.  Returns 0.0 for empty inputs.
#[allow(dead_code)]
pub fn atom_rmsd(coords_a: &[[f64; 3]], coords_b: &[[f64; 3]]) -> f64 {
    let n = coords_a.len().min(coords_b.len());
    if n == 0 {
        return 0.0;
    }
    let sum_sq: f64 = (0..n)
        .map(|i| {
            let dx = coords_a[i][0] - coords_b[i][0];
            let dy = coords_a[i][1] - coords_b[i][1];
            let dz = coords_a[i][2] - coords_b[i][2];
            dx * dx + dy * dy + dz * dz
        })
        .sum();
    (sum_sq / n as f64).sqrt()
}
/// Root-Mean-Square Fluctuation (RMSF) per atom over a trajectory.
///
/// Given a list of frames (each frame is a slice of \[f64; 3\] atom coords),
/// returns the per-atom RMSF relative to the mean position.
#[allow(dead_code)]
pub fn rmsf_per_atom(trajectory: &[Vec<[f64; 3]>]) -> Vec<f64> {
    let n_frames = trajectory.len();
    if n_frames == 0 {
        return vec![];
    }
    let n_atoms = trajectory[0].len();
    if n_atoms == 0 {
        return vec![];
    }
    let mut mean = vec![[0.0f64; 3]; n_atoms];
    for frame in trajectory {
        for (i, pos) in frame.iter().enumerate() {
            mean[i][0] += pos[0];
            mean[i][1] += pos[1];
            mean[i][2] += pos[2];
        }
    }
    let nf = n_frames as f64;
    for m in mean.iter_mut() {
        m[0] /= nf;
        m[1] /= nf;
        m[2] /= nf;
    }
    (0..n_atoms)
        .map(|i| {
            let sum_sq: f64 = trajectory
                .iter()
                .map(|frame| {
                    let dx = frame[i][0] - mean[i][0];
                    let dy = frame[i][1] - mean[i][1];
                    let dz = frame[i][2] - mean[i][2];
                    dx * dx + dy * dy + dz * dz
                })
                .sum();
            (sum_sq / nf).sqrt()
        })
        .collect()
}
/// Mean RMSF over all atoms.
#[allow(dead_code)]
pub fn mean_rmsf(trajectory: &[Vec<[f64; 3]>]) -> f64 {
    let rmsf_vals = rmsf_per_atom(trajectory);
    if rmsf_vals.is_empty() {
        return 0.0;
    }
    rmsf_vals.iter().sum::<f64>() / rmsf_vals.len() as f64
}
/// Compute the displacement auto-correlation function (VACF-like).
///
/// For each time lag τ = 0, 1, ..., max_lag-1, computes:
///   C(τ) = <Δr(t) · Δr(t+τ)> / <Δr(t) · Δr(t)>
///
/// where Δr(t) = r(t) - r_mean.
/// Normalized so C(0) = 1.
#[allow(dead_code)]
pub fn displacement_autocorrelation(
    trajectory: &[Vec<[f64; 3]>],
    atom_idx: usize,
    max_lag: usize,
) -> Vec<f64> {
    let n_frames = trajectory.len();
    if n_frames == 0 || max_lag == 0 {
        return vec![];
    }
    let traj_atom: Vec<[f64; 3]> = trajectory
        .iter()
        .map(|f| {
            if atom_idx < f.len() {
                f[atom_idx]
            } else {
                [0.0; 3]
            }
        })
        .collect();
    let nf = n_frames as f64;
    let mean: [f64; 3] = {
        let s = traj_atom.iter().fold([0.0; 3], |acc, p| {
            [acc[0] + p[0], acc[1] + p[1], acc[2] + p[2]]
        });
        [s[0] / nf, s[1] / nf, s[2] / nf]
    };
    let disp: Vec<[f64; 3]> = traj_atom
        .iter()
        .map(|p| [p[0] - mean[0], p[1] - mean[1], p[2] - mean[2]])
        .collect();
    let c0: f64 = disp
        .iter()
        .map(|d| d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
        .sum::<f64>()
        / nf;
    if c0 < 1e-30 {
        return vec![0.0; max_lag.min(n_frames)];
    }
    let actual_max = max_lag.min(n_frames);
    (0..actual_max)
        .map(|lag| {
            let n_pairs = n_frames - lag;
            if n_pairs == 0 {
                return 0.0;
            }
            let sum: f64 = (0..n_pairs)
                .map(|t| {
                    let a = disp[t];
                    let b = disp[t + lag];
                    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
                })
                .sum();
            sum / (n_pairs as f64 * c0)
        })
        .collect()
}
/// Lipid tail order parameter S₂ = ½⟨3cos²θ - 1⟩.
///
/// `cos_theta_vals` are cos(θ) values where θ is the angle between each
/// bond vector and the director (e.g., membrane normal).
#[allow(dead_code)]
pub fn s2_order_parameter(cos_theta_vals: &[f64]) -> f64 {
    if cos_theta_vals.is_empty() {
        return 0.0;
    }
    let n = cos_theta_vals.len() as f64;
    let mean: f64 = cos_theta_vals
        .iter()
        .map(|&c| 3.0 * c * c - 1.0)
        .sum::<f64>()
        / n;
    0.5 * mean
}
/// Nematic order parameter Q for a set of unit vectors.
///
/// Q = ½⟨3cos²θ - 1⟩ where θ is angle with the average director.
#[allow(dead_code)]
pub fn nematic_order_parameter(vectors: &[[f64; 3]]) -> f64 {
    if vectors.is_empty() {
        return 0.0;
    }
    let n = vectors.len() as f64;
    let mean_v: [f64; 3] = {
        let s = vectors.iter().fold([0.0; 3], |acc, v| {
            [acc[0] + v[0], acc[1] + v[1], acc[2] + v[2]]
        });
        [s[0] / n, s[1] / n, s[2] / n]
    };
    let norm = (mean_v[0] * mean_v[0] + mean_v[1] * mean_v[1] + mean_v[2] * mean_v[2]).sqrt();
    if norm < 1e-30 {
        return 0.0;
    }
    let director = [mean_v[0] / norm, mean_v[1] / norm, mean_v[2] / norm];
    let cos_vals: Vec<f64> = vectors
        .iter()
        .map(|v| v[0] * director[0] + v[1] * director[1] + v[2] * director[2])
        .collect();
    s2_order_parameter(&cos_vals)
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    use std::f64::consts::PI;
    #[test]
    fn test_lj_energy_zero_crossing() {
        assert!((lj_energy(1.0, 1.0, 1.0)).abs() < 1e-12);
    }
    #[test]
    fn test_lj_energy_minimum() {
        let sigma = 1.0;
        let eps = 2.0;
        let r_min = 2.0_f64.powf(1.0 / 6.0) * sigma;
        assert!((lj_energy(r_min, eps, sigma) - (-eps)).abs() < 1e-10);
    }
    #[test]
    fn test_lj_energy_repulsive_at_short_range() {
        assert!(lj_energy(0.5, 1.0, 1.0) > 0.0);
    }
    #[test]
    fn test_lj_energy_sign_far() {
        assert!(lj_energy(2.0, 1.0, 1.0) < 0.0);
    }
    #[test]
    fn test_harmonic_energy_at_equilibrium() {
        assert_eq!(harmonic_energy(1.5, 1.5, 100.0), 0.0);
    }
    #[test]
    fn test_harmonic_energy_displaced() {
        assert!((harmonic_energy(2.0, 1.5, 100.0) - 12.5).abs() < 1e-12);
    }
    #[test]
    fn test_harmonic_energy_symmetric() {
        let k = 50.0;
        let x0 = 1.0;
        let dx = 0.3;
        assert!((harmonic_energy(x0 + dx, x0, k) - harmonic_energy(x0 - dx, x0, k)).abs() < 1e-12);
    }
    #[test]
    fn test_torsion_energy_at_minimum() {
        let k = 2.0;
        let n = 2_u32;
        let delta = 0.0_f64;
        let phi_min = PI / n as f64;
        let e = torsion_energy(phi_min, k, n, delta);
        assert!((e).abs() < 1e-12);
    }
    #[test]
    fn test_torsion_energy_at_maximum() {
        let k = 3.0;
        let n = 1_u32;
        let delta = 0.0_f64;
        let phi_max = 0.0;
        let e = torsion_energy(phi_max, k, n, delta);
        assert!((e - 2.0 * k).abs() < 1e-12);
    }
    #[test]
    fn test_torsion_energy_periodicity() {
        let k = 1.5;
        let n = 3_u32;
        let delta = 0.5;
        let phi = 0.7;
        let e1 = torsion_energy(phi, k, n, delta);
        let e2 = torsion_energy(phi + 2.0 * PI / n as f64, k, n, delta);
        assert!((e1 - e2).abs() < 1e-10);
    }
    #[test]
    fn test_fit_lj_recovers_parameters() {
        let eps_true = 1.0;
        let sig_true = 1.0;
        let rs: Vec<f64> = (0..15).map(|i| 1.0 + i as f64 * 0.15).collect();
        let es: Vec<f64> = rs
            .iter()
            .map(|&r| lj_energy(r, eps_true, sig_true))
            .collect();
        let result = fit_lj_parameters(&rs, &es);
        assert!(result.rmse < 0.15, "RMSE too large: {}", result.rmse);
    }
    #[test]
    fn test_fit_lj_result_has_finite_rmse() {
        let rs = vec![1.0, 1.5, 2.0, 2.5, 3.0];
        let es: Vec<f64> = rs.iter().map(|&r| lj_energy(r, 1.0, 1.0)).collect();
        let result = fit_lj_parameters(&rs, &es);
        assert!(result.rmse.is_finite());
    }
    #[test]
    fn test_fit_lj_single_point() {
        let result = fit_lj_parameters(&[2.0], &[lj_energy(2.0, 1.0, 1.0)]);
        assert!(result.rmse.is_finite());
    }
    #[test]
    fn test_fit_harmonic_bond_exact() {
        let k_true = 100.0;
        let r0_true = 1.5;
        let rs: Vec<f64> = (0..10).map(|i| 1.0 + i as f64 * 0.1).collect();
        let es: Vec<f64> = rs
            .iter()
            .map(|&r| harmonic_energy(r, r0_true, k_true))
            .collect();
        let fit = fit_harmonic_bond(&rs, &es);
        assert!((fit.r0 - r0_true).abs() < 1e-6, "r0 off: {}", fit.r0);
        assert!((fit.k - k_true).abs() < 1e-4, "k off: {}", fit.k);
    }
    #[test]
    fn test_fit_harmonic_bond_three_points() {
        let rs = vec![1.3, 1.5, 1.7];
        let es: Vec<f64> = rs.iter().map(|&r| harmonic_energy(r, 1.5, 200.0)).collect();
        let fit = fit_harmonic_bond(&rs, &es);
        assert!((fit.r0 - 1.5).abs() < 1e-6);
        assert!((fit.k - 200.0).abs() < 1e-3);
    }
    #[test]
    fn test_fit_harmonic_bond_returns_positive_k() {
        let rs: Vec<f64> = (0..5).map(|i| 1.0 + i as f64 * 0.2).collect();
        let es: Vec<f64> = rs.iter().map(|&r| harmonic_energy(r, 1.4, 50.0)).collect();
        let fit = fit_harmonic_bond(&rs, &es);
        assert!(fit.k >= 0.0);
    }
    #[test]
    fn test_fit_harmonic_angle_exact() {
        let k_true = 80.0;
        let theta0_true = PI / 2.0;
        let thetas: Vec<f64> = (0..8).map(|i| 1.0 + i as f64 * 0.1).collect();
        let es: Vec<f64> = thetas
            .iter()
            .map(|&t| harmonic_energy(t, theta0_true, k_true))
            .collect();
        let fit = fit_harmonic_angle(&thetas, &es);
        assert!(
            (fit.theta0 - theta0_true).abs() < 1e-5,
            "theta0 off: {}",
            fit.theta0
        );
        assert!((fit.k - k_true).abs() < 1e-3, "k off: {}", fit.k);
    }
    #[test]
    fn test_fit_harmonic_angle_returns_positive_k() {
        let thetas: Vec<f64> = (0..6).map(|i| 1.2 + i as f64 * 0.05).collect();
        let es: Vec<f64> = thetas
            .iter()
            .map(|&t| harmonic_energy(t, 1.4, 60.0))
            .collect();
        let fit = fit_harmonic_angle(&thetas, &es);
        assert!(fit.k >= 0.0);
    }
    #[test]
    fn test_fit_dihedral_recovers_periodicity() {
        let n_true = 2_u32;
        let k_true = 1.5;
        let delta_true = 0.0_f64;
        let phis: Vec<f64> = (0..20).map(|i| i as f64 * PI / 10.0).collect();
        let es: Vec<f64> = phis
            .iter()
            .map(|&p| torsion_energy(p, k_true, n_true, delta_true))
            .collect();
        let fit = fit_dihedral(&phis, &es);
        assert_eq!(fit.n, n_true, "Wrong periodicity: got {}", fit.n);
    }
    #[test]
    fn test_fit_dihedral_rmse_small() {
        let n_true = 1_u32;
        let k_true = 2.0;
        let delta_true = PI / 4.0;
        let phis: Vec<f64> = (0..24).map(|i| i as f64 * PI / 12.0).collect();
        let es: Vec<f64> = phis
            .iter()
            .map(|&p| torsion_energy(p, k_true, n_true, delta_true))
            .collect();
        let fit = fit_dihedral(&phis, &es);
        let pred: Vec<f64> = phis
            .iter()
            .map(|&p| torsion_energy(p, fit.k, fit.n, fit.delta))
            .collect();
        let r = rmse(&pred, &es);
        assert!(r < 0.05, "RMSE too large: {r}");
    }
    #[test]
    fn test_fit_dihedral_k_non_negative() {
        let phis: Vec<f64> = (0..12).map(|i| i as f64 * PI / 6.0).collect();
        let es: Vec<f64> = phis
            .iter()
            .map(|&p| torsion_energy(p, 1.0, 3, 0.0))
            .collect();
        let fit = fit_dihedral(&phis, &es);
        assert!(fit.k >= 0.0);
    }
    #[test]
    fn test_cross_validate_perfect_fit() {
        let k = 100.0;
        let r0 = 1.5;
        let data: Vec<(f64, f64)> = (0..20)
            .map(|i| {
                let x = 1.0 + i as f64 * 0.05;
                (x, harmonic_energy(x, r0, k))
            })
            .collect();
        let cv_rmse = cross_validate(&data, 5);
        assert!(cv_rmse < 1.0, "CV RMSE unexpectedly large: {cv_rmse}");
    }
    #[test]
    fn test_cross_validate_two_folds() {
        let data: Vec<(f64, f64)> = (0..10)
            .map(|i| {
                let x = 1.0 + i as f64 * 0.1;
                (x, harmonic_energy(x, 1.5, 50.0))
            })
            .collect();
        let cv_rmse = cross_validate(&data, 2);
        assert!(cv_rmse.is_finite());
        assert!(cv_rmse >= 0.0);
    }
    #[test]
    fn test_cross_validate_returns_finite() {
        let data: Vec<(f64, f64)> = (0..12)
            .map(|i| {
                let x = 0.8 + i as f64 * 0.1;
                (x, harmonic_energy(x, 1.2, 80.0))
            })
            .collect();
        let result = cross_validate(&data, 4);
        assert!(result.is_finite());
    }
    #[test]
    fn test_det3x3_identity() {
        let m = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!((det3x3(&m) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_det3x3_known() {
        let m = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 10.0]];
        assert!((det3x3(&m) - (-3.0)).abs() < 1e-12);
    }
    #[test]
    fn test_rmse_identical() {
        let v = vec![1.0, 2.0, 3.0];
        assert_eq!(rmse(&v, &v.clone()), 0.0);
    }
    #[test]
    fn test_rmse_known_value() {
        let pred = vec![1.0, 3.0];
        let target = vec![2.0, 2.0];
        assert!((rmse(&pred, &target) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_mae_identical() {
        let v = vec![1.0, 2.0, 3.0];
        assert!((mae(&v, &v.clone())).abs() < 1e-12);
    }
    #[test]
    fn test_mae_known_value() {
        let pred = vec![1.0, 3.0];
        let target = vec![2.0, 2.0];
        assert!((mae(&pred, &target) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_r_squared_perfect() {
        let v = vec![1.0, 2.0, 3.0];
        assert!((r_squared(&v, &v.clone()) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_r_squared_mean_prediction() {
        let target = vec![1.0, 2.0, 3.0];
        let pred = vec![2.0, 2.0, 2.0];
        let r2 = r_squared(&pred, &target);
        assert!(r2.abs() < 1e-10, "R² for mean prediction = {r2}");
    }
    #[test]
    fn test_max_error_known() {
        let pred = vec![0.0, 5.0, 2.0];
        let target = vec![1.0, 2.0, 2.0];
        let me = max_error(&pred, &target);
        assert!((me - 3.0).abs() < 1e-12, "max_error = {me}");
    }
    #[test]
    fn test_max_error_zero_for_perfect() {
        let v = vec![1.0, 2.0];
        assert!(max_error(&v, &v.clone()).abs() < 1e-12);
    }
    #[test]
    fn test_resp_charges_net_charge_satisfied() {
        let atom_positions = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let esp_points: Vec<([f64; 3], f64)> = vec![
            ([0.5, 0.5, 0.0], 0.5),
            ([0.5, -0.5, 0.0], 0.5),
            ([0.0, 0.5, 0.5], 0.3),
        ];
        let result = fit_resp_charges(&esp_points, &atom_positions, 0.0, 0.01, 50);
        assert!(
            result.net_charge_error < 0.01,
            "net charge error = {}",
            result.net_charge_error
        );
        assert_eq!(result.charges.len(), 2);
    }
    #[test]
    fn test_resp_charges_sum_to_total() {
        let atom_positions = vec![[0.0f64; 3], [2.0, 0.0, 0.0]];
        let esp_points = vec![([1.0f64, 0.0, 0.0], 0.0)];
        let result = fit_resp_charges(&esp_points, &atom_positions, -1.0, 0.1, 100);
        let total: f64 = result.charges.iter().sum();
        assert!((total - (-1.0)).abs() < 0.01, "total charge = {total}");
    }
    #[test]
    fn test_resp_empty_input() {
        let result = fit_resp_charges(&[], &[], 0.0, 0.1, 10);
        assert!(result.charges.is_empty());
    }
    #[test]
    fn test_parametric_potential_constant_basis() {
        let pot = ParametricPotential::new(vec![
            Box::new(|_x| 1.0),
            Box::new(|x| x),
            Box::new(|x| x * x),
        ]);
        let xs: Vec<f64> = (0..10).map(|i| i as f64 * 0.1).collect();
        let ys: Vec<f64> = xs.iter().map(|&x| 1.0 + 2.0 * x + 3.0 * x * x).collect();
        let coeffs = pot.fit(&xs, &ys);
        assert_eq!(coeffs.len(), 3);
        assert!((coeffs[0] - 1.0).abs() < 1e-4, "c0 = {}", coeffs[0]);
        assert!((coeffs[1] - 2.0).abs() < 1e-4, "c1 = {}", coeffs[1]);
        assert!((coeffs[2] - 3.0).abs() < 1e-4, "c2 = {}", coeffs[2]);
    }
    #[test]
    fn test_parametric_potential_predict() {
        let pot = ParametricPotential::new(vec![Box::new(|x: f64| x * x)]);
        let xs: Vec<f64> = vec![1.0, 2.0, 3.0];
        let ys: Vec<f64> = vec![2.0, 8.0, 18.0];
        let coeffs = pot.fit(&xs, &ys);
        let pred = pot.predict(4.0, &coeffs);
        assert!((pred - 32.0).abs() < 0.5, "pred = {pred}");
    }
    #[test]
    fn test_bootstrap_harmonic_bond_returns_result() {
        let data: Vec<(f64, f64)> = (0..20)
            .map(|i| {
                let x = 1.0 + i as f64 * 0.05;
                (x, harmonic_energy(x, 1.5, 100.0))
            })
            .collect();
        let result = bootstrap_harmonic_bond(&data, 20, 42);
        assert_eq!(result.n_bootstrap, 20);
        assert!(result.std_k >= 0.0);
        assert!(result.std_r0 >= 0.0);
    }
    #[test]
    fn test_bootstrap_harmonic_bond_small_data() {
        let data = vec![(1.0, 0.0), (1.5, 0.0)];
        let result = bootstrap_harmonic_bond(&data, 10, 1);
        assert_eq!(result.n_bootstrap, 0);
    }
    #[test]
    fn test_solve_linear_system_identity() {
        let a = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let b = vec![3.0, 7.0];
        let x = solve_linear_system(&a, &b);
        assert!((x[0] - 3.0).abs() < 1e-10);
        assert!((x[1] - 7.0).abs() < 1e-10);
    }
    #[test]
    fn test_solve_linear_system_2x2() {
        let a = vec![vec![3.0, 1.0], vec![1.0, 2.0]];
        let b = vec![9.0, 8.0];
        let x = solve_linear_system(&a, &b);
        assert!((x[0] - 2.0).abs() < 1e-8, "x[0] = {}", x[0]);
        assert!((x[1] - 3.0).abs() < 1e-8, "x[1] = {}", x[1]);
    }
}
