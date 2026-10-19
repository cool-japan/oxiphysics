// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Harmonic (steady-state frequency-domain) analysis by direct solution.
//!
//! A harmonic load `F e^{iωt}` on the damped system `M ü + C u̇ + K u = F e^{iωt}`
//! produces the steady state `u(t) = û e^{iωt}` whose complex amplitude solves
//! the *dynamic stiffness* system
//!
//! ```text
//! (K − ω² M + iω C) û = F,        C = α M + β K   (Rayleigh damping).
//! ```
//!
//! With `û = u_r + i u_i` and `A = K − ω² M` the real and imaginary parts give the
//! real `2n x 2n` block system
//!
//! ```text
//! [  A   −ωC ] [u_r]   [F]
//! [ ωC     A ] [u_i] = [0].
//! ```
//!
//! `A` is indefinite above the first natural frequency, so a conjugate-gradient
//! solve (or a fixed number of block Gauss–Seidel sweeps, which diverge once
//! `ω‖C‖` is not small) is not valid.  The system is therefore solved directly:
//! the unknowns are interleaved (`2i` = real part, `2i + 1` = imaginary part of
//! DOF `i`), the DOFs are renumbered by reverse Cuthill–McKee to keep the
//! band narrow, and a banded LU factorisation with partial pivoting is used.  If
//! the band storage would exceed the memory budget, ILU(0)-preconditioned GMRES
//! is used instead and its *true* residual is verified.

use super::harmonic_solve::{BAND_STORAGE_LIMIT, BandLu, gmres_ilu0, reverse_cuthill_mckee};
use crate::error::{Error, Result};
use crate::sparse::CsrMatrix;

/// Relative residual at which the iterative fallback stops.
const GMRES_TOL: f64 = 1.0e-12;
/// Largest true relative residual accepted from the iterative fallback.
const GMRES_ACCEPT: f64 = 1.0e-8;
/// Krylov dimension of the iterative fallback.
const GMRES_RESTART: usize = 80;
/// Maximum number of GMRES restart cycles of the iterative fallback.
const GMRES_MAX_CYCLES: usize = 100;

/// Result of harmonic analysis at a single frequency.
///
/// The steady-state displacement is `u(t) = Re[(u_r + i u_i) e^{iωt}]
/// = u_r cos ωt − u_i sin ωt`, so a damped system responds with a negative
/// imaginary part (the response lags the load).
#[derive(Debug, Clone)]
pub struct HarmonicResult {
    /// Excitation frequency ω (rad/s).
    pub omega: f64,
    /// Real part of displacement response.
    pub displacement_real: Vec<f64>,
    /// Imaginary part of displacement response.
    pub displacement_imag: Vec<f64>,
    /// Peak displacement magnitude at each DOF.
    pub displacement_amplitude: Vec<f64>,
}

/// Perform harmonic analysis using direct frequency response.
///
/// Solves `(K − ω² M + iω C) û = F` for the complex displacement amplitude
/// `û = u_r + i u_i` (time convention `e^{+iωt}`), where the damping matrix is
/// Rayleigh damping `C = α M + β K`.
///
/// # Arguments
/// * `stiffness` - Global stiffness matrix K
/// * `mass` - Global mass matrix M (its sparsity pattern may differ from K's)
/// * `omega` - Excitation frequency (rad/s); `0` gives the static solution
/// * `alpha` - Mass-proportional damping coefficient
/// * `beta` - Stiffness-proportional damping coefficient
/// * `force_vector` - Applied harmonic force amplitude (real)
///
/// # Errors
///
/// * [`Error::InvalidInput`] if the matrices are not square and of equal size, if
///   `force_vector` has the wrong length, or if any input is non-finite.
/// * [`Error::SingularSystem`] if the dynamic stiffness matrix is singular: an
///   undamped system excited exactly at a natural frequency, or an unconstrained
///   structure analysed statically.
/// * [`Error::NotConverged`] only when the direct solver cannot be used for memory
///   reasons and the iterative fallback does not reach a relative residual of
///   `1e-8`.
pub fn harmonic_analysis(
    stiffness: &CsrMatrix,
    mass: &CsrMatrix,
    omega: f64,
    alpha: f64,
    beta: f64,
    force_vector: &[f64],
) -> Result<HarmonicResult> {
    harmonic_analysis_with(stiffness, mass, omega, alpha, beta, force_vector, false)
}

/// Check the inputs of [`harmonic_analysis`] and return the number of DOFs.
fn validate_inputs(
    stiffness: &CsrMatrix,
    mass: &CsrMatrix,
    omega: f64,
    alpha: f64,
    beta: f64,
    force_vector: &[f64],
) -> Result<usize> {
    let n = stiffness.nrows;
    if stiffness.ncols != n {
        return Err(Error::InvalidInput(format!(
            "the stiffness matrix must be square, got {} x {}",
            stiffness.nrows, stiffness.ncols
        )));
    }
    if mass.nrows != n || mass.ncols != n {
        return Err(Error::InvalidInput(format!(
            "the mass matrix must be {n} x {n} like the stiffness matrix, got {} x {}",
            mass.nrows, mass.ncols
        )));
    }
    if force_vector.len() != n {
        return Err(Error::InvalidInput(format!(
            "the force vector has {} entries but the system has {n} DOFs",
            force_vector.len()
        )));
    }
    if !omega.is_finite() || !alpha.is_finite() || !beta.is_finite() {
        return Err(Error::InvalidInput(
            "omega, alpha and beta must be finite".to_string(),
        ));
    }
    if force_vector.iter().any(|f| !f.is_finite()) {
        return Err(Error::InvalidInput(
            "the force vector holds a non-finite entry".to_string(),
        ));
    }
    for (name, matrix) in [("stiffness", stiffness), ("mass", mass)] {
        let well_formed = matrix.row_ptr.len() == n + 1
            && matrix.row_ptr.first().copied() == Some(0)
            && matrix.row_ptr.windows(2).all(|pair| pair[0] <= pair[1])
            && matrix.col_indices.len() == matrix.values.len()
            && matrix.row_ptr.last().copied() == Some(matrix.values.len())
            && matrix.col_indices.iter().all(|&c| c < n);
        if !well_formed {
            return Err(Error::InvalidInput(format!(
                "the {name} matrix is not a well-formed CSR matrix"
            )));
        }
        if matrix.values.iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput(format!(
                "the {name} matrix holds a non-finite entry"
            )));
        }
    }
    Ok(n)
}

/// Assemble the interleaved real block system of the dynamic stiffness.
///
/// DOF `i` is placed at renumbered position `position[i]`; its real part is
/// unknown `2 position[i]` and its imaginary part `2 position[i] + 1`.  Every
/// stored entry of `K` and of `M` contributes (duplicates are summed by
/// [`CsrMatrix::from_triplets`]), so the two sparsity patterns need not agree.
fn assemble_block_system(
    stiffness: &CsrMatrix,
    mass: &CsrMatrix,
    omega: f64,
    alpha: f64,
    beta: f64,
    position: &[usize],
) -> CsrMatrix {
    let n = stiffness.nrows;
    let omega_sq = omega * omega;
    let mut triplets: Vec<(usize, usize, f64)> =
        Vec::with_capacity(4 * (stiffness.nnz() + mass.nnz()));
    // Block of DOF pair (row, col): [[a, -ωc], [ωc, a]] with a = K − ω² M entry, c = C entry.
    let mut push_block = |row: usize, col: usize, a: f64, c: f64| {
        let r = 2 * position[row];
        let s = 2 * position[col];
        let coupling = omega * c;
        if a != 0.0 {
            triplets.push((r, s, a));
            triplets.push((r + 1, s + 1, a));
        }
        if coupling != 0.0 {
            triplets.push((r, s + 1, -coupling));
            triplets.push((r + 1, s, coupling));
        }
    };
    for row in 0..n {
        for idx in stiffness.row_ptr[row]..stiffness.row_ptr[row + 1] {
            let value = stiffness.values[idx];
            push_block(row, stiffness.col_indices[idx], value, beta * value);
        }
        for idx in mass.row_ptr[row]..mass.row_ptr[row + 1] {
            let value = mass.values[idx];
            push_block(row, mass.col_indices[idx], -omega_sq * value, alpha * value);
        }
    }
    CsrMatrix::from_triplets(2 * n, 2 * n, &triplets)
}

/// Solve the real block system: banded LU when it fits in memory (and unless
/// `prefer_iterative`), ILU(0)-GMRES with a verified residual otherwise.
fn solve_block_system(system: &CsrMatrix, rhs: &[f64], prefer_iterative: bool) -> Result<Vec<f64>> {
    let (kl, ku) = BandLu::bandwidths(system);
    let fits = BandLu::storage_len(rhs.len(), kl, ku) <= BAND_STORAGE_LIMIT;
    if fits && !prefer_iterative {
        let x = BandLu::factor(system)?.solve(rhs);
        return if x.iter().all(|v| v.is_finite()) {
            Ok(x)
        } else {
            Err(Error::SingularSystem)
        };
    }
    let outcome = gmres_ilu0(system, rhs, GMRES_RESTART, GMRES_MAX_CYCLES, GMRES_TOL);
    if outcome.relative_residual <= GMRES_ACCEPT && outcome.x.iter().all(|v| v.is_finite()) {
        Ok(outcome.x)
    } else {
        Err(Error::NotConverged {
            solver: "ILU(0)-GMRES".to_string(),
            iterations: outcome.iterations,
            residual: outcome.relative_residual,
        })
    }
}

/// [`harmonic_analysis`] with the choice of the linear solver exposed for tests.
fn harmonic_analysis_with(
    stiffness: &CsrMatrix,
    mass: &CsrMatrix,
    omega: f64,
    alpha: f64,
    beta: f64,
    force_vector: &[f64],
    prefer_iterative: bool,
) -> Result<HarmonicResult> {
    let n = validate_inputs(stiffness, mass, omega, alpha, beta, force_vector)?;
    if n == 0 {
        return Ok(HarmonicResult {
            omega,
            displacement_real: Vec::new(),
            displacement_imag: Vec::new(),
            displacement_amplitude: Vec::new(),
        });
    }

    // Renumber the DOFs to keep the bandwidth of the block system small.
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); n];
    for matrix in [stiffness, mass] {
        for (row, neighbours) in adjacency.iter_mut().enumerate() {
            let range = matrix.row_ptr[row]..matrix.row_ptr[row + 1];
            neighbours.extend(matrix.col_indices[range].iter().copied());
        }
    }
    let order = reverse_cuthill_mckee(&adjacency);
    let mut position = vec![0usize; n];
    for (new_index, &old_index) in order.iter().enumerate() {
        position[old_index] = new_index;
    }

    let system = assemble_block_system(stiffness, mass, omega, alpha, beta, &position);
    let mut rhs = vec![0.0_f64; 2 * n];
    for (dof, &load) in force_vector.iter().enumerate() {
        rhs[2 * position[dof]] = load;
    }
    let solution = solve_block_system(&system, &rhs, prefer_iterative)?;

    let displacement_real: Vec<f64> = position.iter().map(|&p| solution[2 * p]).collect();
    let displacement_imag: Vec<f64> = position.iter().map(|&p| solution[2 * p + 1]).collect();
    let displacement_amplitude: Vec<f64> = displacement_real
        .iter()
        .zip(displacement_imag.iter())
        .map(|(re, im)| re.hypot(*im))
        .collect();
    Ok(HarmonicResult {
        omega,
        displacement_real,
        displacement_imag,
        displacement_amplitude,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modal::{ModalResult, frequency_response_function};
    use std::f64::consts::PI;

    /// Complex number as `(re, im)`.
    type Cx = (f64, f64);

    fn cx_mul(a: Cx, b: Cx) -> Cx {
        (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
    }

    fn cx_div(a: Cx, b: Cx) -> Cx {
        let d = b.0 * b.0 + b.1 * b.1;
        ((a.0 * b.0 + a.1 * b.1) / d, (a.1 * b.0 - a.0 * b.1) / d)
    }

    /// Independent reference: complex Gaussian elimination with partial pivoting.
    fn solve_complex(mut a: Vec<Vec<Cx>>, mut b: Vec<Cx>) -> Vec<Cx> {
        let n = b.len();
        for col in 0..n {
            let mut p = col;
            for row in (col + 1)..n {
                let mag = |z: Cx| z.0.hypot(z.1);
                if mag(a[row][col]) > mag(a[p][col]) {
                    p = row;
                }
            }
            a.swap(col, p);
            b.swap(col, p);
            for row in (col + 1)..n {
                let f = cx_div(a[row][col], a[col][col]);
                let pivot_values = a[col].clone();
                for (entry, upper) in a[row].iter_mut().zip(pivot_values.iter()).skip(col) {
                    let t = cx_mul(f, *upper);
                    entry.0 -= t.0;
                    entry.1 -= t.1;
                }
                let t = cx_mul(f, b[col]);
                b[row].0 -= t.0;
                b[row].1 -= t.1;
            }
        }
        let mut x = vec![(0.0, 0.0); n];
        for row in (0..n).rev() {
            let mut s = b[row];
            for k in (row + 1)..n {
                let t = cx_mul(a[row][k], x[k]);
                s.0 -= t.0;
                s.1 -= t.1;
            }
            x[row] = cx_div(s, a[row][row]);
        }
        x
    }

    fn from_dense(rows: &[Vec<f64>]) -> CsrMatrix {
        let n = rows.len();
        let mut triplets = Vec::new();
        for (i, row) in rows.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                if v != 0.0 {
                    triplets.push((i, j, v));
                }
            }
        }
        CsrMatrix::from_triplets(n, n, &triplets)
    }

    /// Reference response `(K − ω²M + iω(αM + βK))⁻¹ F` from dense data.
    fn reference_response(
        k: &[Vec<f64>],
        m: &[Vec<f64>],
        omega: f64,
        alpha: f64,
        beta: f64,
        f: &[f64],
    ) -> Vec<Cx> {
        let n = f.len();
        let z: Vec<Vec<Cx>> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        (
                            k[i][j] - omega * omega * m[i][j],
                            omega * (alpha * m[i][j] + beta * k[i][j]),
                        )
                    })
                    .collect()
            })
            .collect();
        let rhs: Vec<Cx> = f.iter().map(|&v| (v, 0.0)).collect();
        solve_complex(z, rhs)
    }

    fn assert_matches_reference(result: &HarmonicResult, reference: &[Cx], tol: f64) {
        assert_eq!(result.displacement_real.len(), reference.len());
        for (i, z) in reference.iter().enumerate() {
            let scale = 1.0 + z.0.hypot(z.1);
            assert!(
                (result.displacement_real[i] - z.0).abs() <= tol * scale,
                "dof {i}: u_r {} vs {}",
                result.displacement_real[i],
                z.0
            );
            assert!(
                (result.displacement_imag[i] - z.1).abs() <= tol * scale,
                "dof {i}: u_i {} vs {}",
                result.displacement_imag[i],
                z.1
            );
            assert!(
                (result.displacement_amplitude[i] - z.0.hypot(z.1)).abs() <= tol * scale,
                "dof {i}: |u| {} vs {}",
                result.displacement_amplitude[i],
                z.0.hypot(z.1)
            );
        }
    }

    /// `k_s·tridiag(−1, 2, −1)` fixed-fixed chain stiffness and `m_s·I` mass.
    fn uniform_chain(n: usize, k_s: f64, m_s: f64) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
        let mut k = vec![vec![0.0; n]; n];
        let mut m = vec![vec![0.0; n]; n];
        for i in 0..n {
            k[i][i] = 2.0 * k_s;
            m[i][i] = m_s;
            if i + 1 < n {
                k[i][i + 1] = -k_s;
                k[i + 1][i] = -k_s;
            }
        }
        (k, m)
    }

    /// Largest residual of the two real block equations `A u_r − ωC u_i = F` and
    /// `ωC u_r + A u_i = 0` (`A = K − ω²M`, `C = αM + βK`) for a computed response.
    fn block_residual(
        k: &CsrMatrix,
        m: &CsrMatrix,
        response: &HarmonicResult,
        (omega, alpha, beta): (f64, f64, f64),
        force: &[f64],
    ) -> f64 {
        let apply = |a: f64, b: f64, v: &[f64]| -> Vec<f64> {
            // a K v + b M v
            k.mul_vec(v)
                .iter()
                .zip(m.mul_vec(v).iter())
                .map(|(kv, mv)| a * kv + b * mv)
                .collect()
        };
        let (u_r, u_i) = (&response.displacement_real, &response.displacement_imag);
        let a_ur = apply(1.0, -omega * omega, u_r);
        let a_ui = apply(1.0, -omega * omega, u_i);
        let c_ur = apply(beta, alpha, u_r);
        let c_ui = apply(beta, alpha, u_i);
        let mut worst = 0.0_f64;
        for i in 0..force.len() {
            worst = worst.max((a_ur[i] - omega * c_ui[i] - force[i]).abs());
            worst = worst.max((omega * c_ur[i] + a_ui[i]).abs());
        }
        worst
    }

    #[test]
    fn undamped_one_dof_above_resonance_is_in_antiphase() {
        // (k − ω²m) u = F with k = m = 1, ω = 2, F = 1  =>  u = −1/3, purely real.
        let k = from_dense(&[vec![1.0]]);
        let m = from_dense(&[vec![1.0]]);
        let r = harmonic_analysis(&k, &m, 2.0, 0.0, 0.0, &[1.0]).expect("regular system");
        assert!((r.displacement_real[0] + 1.0 / 3.0).abs() < 1e-14, "{r:?}");
        assert!(r.displacement_imag[0].abs() < 1e-14);
        assert!((r.displacement_amplitude[0] - 1.0 / 3.0).abs() < 1e-14);
        assert_eq!(r.omega, 2.0);
    }

    #[test]
    fn damped_one_dof_matches_the_complex_closed_form() {
        // u = F / (k − ω²m + iωβk) with k = m = 1, β = 0.02, ω = 0.9.
        let k = from_dense(&[vec![1.0]]);
        let m = from_dense(&[vec![1.0]]);
        let (omega, beta, f) = (0.9, 0.02, 1.0);
        let r = harmonic_analysis(&k, &m, omega, 0.0, beta, &[f]).expect("regular system");
        let expected = cx_div((f, 0.0), (1.0 - omega * omega, omega * beta));
        assert!((r.displacement_real[0] - expected.0).abs() < 1e-8);
        assert!((r.displacement_imag[0] - expected.1).abs() < 1e-8);
        // Time convention e^{+iωt}: the damped response lags the load.
        assert!(r.displacement_imag[0] < 0.0);
    }

    #[test]
    fn mass_proportional_damping_enters_with_the_mass() {
        // (k − ω²m + iω α m) u = F with m = 2.
        let k = from_dense(&[vec![5.0]]);
        let m = from_dense(&[vec![2.0]]);
        let (omega, alpha) = (1.3, 0.4);
        let r = harmonic_analysis(&k, &m, omega, alpha, 0.0, &[3.0]).expect("regular system");
        let expected = cx_div((3.0, 0.0), (5.0 - omega * omega * 2.0, omega * alpha * 2.0));
        assert!((r.displacement_real[0] - expected.0).abs() < 1e-12);
        assert!((r.displacement_imag[0] - expected.1).abs() < 1e-12);
    }

    #[test]
    fn damped_resonance_peak_equals_f_over_omega_beta_k() {
        // At ω² = k/m the real part of the dynamic stiffness vanishes: u = F/(iωβk).
        let k = from_dense(&[vec![1.0]]);
        let m = from_dense(&[vec![1.0]]);
        let r = harmonic_analysis(&k, &m, 1.0, 0.0, 0.02, &[1.0]).expect("damped resonance");
        assert!(r.displacement_real[0].abs() < 1e-12);
        assert!((r.displacement_imag[0] + 50.0).abs() < 1e-9);
        assert!((r.displacement_amplitude[0] - 50.0).abs() < 1e-9);
    }

    #[test]
    fn static_limit_returns_k_inverse_f() {
        let k = from_dense(&[vec![100.0]]);
        let m = from_dense(&[vec![1.0]]);
        let r = harmonic_analysis(&k, &m, 0.0, 0.3, 0.01, &[10.0]).expect("static solve");
        assert!((r.displacement_real[0] - 0.1).abs() < 1e-14);
        assert_eq!(r.displacement_imag[0], 0.0);
    }

    #[test]
    fn damped_two_dof_above_first_resonance_matches_a_direct_complex_solve() {
        // m1 = 2, m2 = 1; springs k1 = 40 (ground), k2 = 25 (between the masses).
        let k = vec![vec![65.0, -25.0], vec![-25.0, 25.0]];
        let m = vec![vec![2.0, 0.0], vec![0.0, 1.0]];
        // Natural frequencies: ω₁ = 3.268 and ω₂ = 6.843 rad/s, so ω = 4.5 lies between them.
        let (omega, alpha, beta) = (4.5, 0.15, 0.004);
        let f = [1.0, -0.5];
        let result = harmonic_analysis(&from_dense(&k), &from_dense(&m), omega, alpha, beta, &f)
            .expect("regular system");
        let reference = reference_response(&k, &m, omega, alpha, beta, &f);
        assert_matches_reference(&result, &reference, 1e-12);
        // Golden values of the 2x2 complex system (Cramer's rule, computed independently).
        let golden = [
            (0.015_200_477_490_805_02, -0.001_253_962_762_674_119_1),
            (-0.024_963_584_567_907_78, 0.000_752_669_145_874_849_6),
        ];
        for (i, (re, im)) in golden.iter().enumerate() {
            assert!(
                (result.displacement_real[i] - re).abs() < 1e-12,
                "dof {i} u_r"
            );
            assert!(
                (result.displacement_imag[i] - im).abs() < 1e-12,
                "dof {i} u_i"
            );
        }
    }

    #[test]
    fn three_dof_chain_matches_the_modal_frf_with_equivalent_modal_damping() {
        // Uniform chain: closed-form mass-normalised modes, so the full modal sum is exact.
        let n = 3usize;
        let (k_s, m_s) = (2.0, 0.5);
        let (alpha, beta) = (0.05, 0.002);
        let (k, m) = uniform_chain(n, k_s, m_s);
        let frequencies: Vec<f64> = (1..=n)
            .map(|j| ((2.0 * k_s / m_s) * (1.0 - (j as f64 * PI / (n as f64 + 1.0)).cos())).sqrt())
            .collect();
        let mode_shapes: Vec<Vec<f64>> = (1..=n)
            .map(|j| {
                (1..=n)
                    .map(|i| {
                        (2.0 / ((n as f64 + 1.0) * m_s)).sqrt()
                            * (i as f64 * j as f64 * PI / (n as f64 + 1.0)).sin()
                    })
                    .collect()
            })
            .collect();
        // C = αM + βK is diagonal in the modal basis: 2 ζ_r ω_r = α + β ω_r².
        let zeta: Vec<f64> = frequencies
            .iter()
            .map(|&w| alpha / (2.0 * w) + beta * w / 2.0)
            .collect();
        let modal = ModalResult {
            frequencies: frequencies.clone(),
            mode_shapes,
        };
        let f = [1.0, 0.5, -0.25];
        // Above the first resonance (ω_1 ≈ 1.53) and below the second (ω_2 ≈ 2.83).
        for &omega in &[2.2_f64, 1.0, 3.3] {
            let direct =
                harmonic_analysis(&from_dense(&k), &from_dense(&m), omega, alpha, beta, &f)
                    .expect("regular system");
            let frf = frequency_response_function(&modal, omega, &zeta, &f);
            for i in 0..n {
                let (re, im) = (
                    frf.magnitude[i] * frf.phase[i].cos(),
                    frf.magnitude[i] * frf.phase[i].sin(),
                );
                assert!(
                    (direct.displacement_real[i] - re).abs() < 1e-9,
                    "omega {omega} dof {i}: u_r {} vs FRF {re}",
                    direct.displacement_real[i]
                );
                assert!(
                    (direct.displacement_imag[i] - im).abs() < 1e-9,
                    "omega {omega} dof {i}: u_i {} vs FRF {im}",
                    direct.displacement_imag[i]
                );
            }
        }
    }

    #[test]
    fn mass_entries_outside_the_stiffness_pattern_are_not_dropped() {
        // K is diagonal; the consistent mass couples the DOFs.  The previous
        // implementation only looked up M on K's pattern and ignored this coupling.
        let k = vec![vec![4.0, 0.0], vec![0.0, 9.0]];
        let m = vec![vec![1.0, 0.3], vec![0.3, 1.0]];
        let (omega, alpha, beta) = (1.5, 0.02, 0.01);
        let f = [1.0, 2.0];
        let result = harmonic_analysis(&from_dense(&k), &from_dense(&m), omega, alpha, beta, &f)
            .expect("regular system");
        let reference = reference_response(&k, &m, omega, alpha, beta, &f);
        assert_matches_reference(&result, &reference, 1e-12);
    }

    #[test]
    fn dof_numbering_does_not_change_the_physical_response() {
        let n = 7usize;
        let (k, m) = uniform_chain(n, 3.0, 1.2);
        let f: Vec<f64> = (0..n).map(|i| 1.0 + 0.3 * i as f64).collect();
        let (omega, alpha, beta) = (2.4, 0.03, 0.002);
        let base = harmonic_analysis(&from_dense(&k), &from_dense(&m), omega, alpha, beta, &f)
            .expect("regular system");
        // Scramble the numbering: DOF i of the new system is DOF perm[i] of the old one.
        let perm: Vec<usize> = (0..n).map(|i| (i * 3 + 2) % n).collect();
        let permute_matrix = |a: &[Vec<f64>]| -> Vec<Vec<f64>> {
            (0..n)
                .map(|i| (0..n).map(|j| a[perm[i]][perm[j]]).collect())
                .collect()
        };
        let f_perm: Vec<f64> = (0..n).map(|i| f[perm[i]]).collect();
        let scrambled = harmonic_analysis(
            &from_dense(&permute_matrix(&k)),
            &from_dense(&permute_matrix(&m)),
            omega,
            alpha,
            beta,
            &f_perm,
        )
        .expect("regular system");
        for (i, &old) in perm.iter().enumerate() {
            assert!((scrambled.displacement_real[i] - base.displacement_real[old]).abs() < 1e-11);
            assert!((scrambled.displacement_imag[i] - base.displacement_imag[old]).abs() < 1e-11);
        }
    }

    #[test]
    fn long_damped_chain_satisfies_both_block_equations() {
        // 400 DOFs, far above many resonances: the old fixed-sweep CG could not do this.
        let n = 400usize;
        let (k, m) = uniform_chain(n, 50.0, 1.0);
        let k_csr = from_dense(&k);
        let m_csr = from_dense(&m);
        let (omega, alpha, beta) = (7.3, 0.02, 0.0005);
        let mut f = vec![0.0; n];
        f[n / 3] = 1.0;
        f[n - 5] = -2.0;
        let r = harmonic_analysis(&k_csr, &m_csr, omega, alpha, beta, &f).expect("regular");
        let worst = block_residual(&k_csr, &m_csr, &r, (omega, alpha, beta), &f);
        assert!(worst < 1e-9, "residual {worst}");
        assert!(
            r.displacement_imag.iter().any(|v| v.abs() > 1e-6),
            "response must be damped"
        );
    }

    #[test]
    fn scrambled_two_dimensional_grid_satisfies_both_block_equations() {
        // 24 x 18 DOFs with a grounded five-point stiffness stencil and a lumped mass,
        // numbered in a scrambled order so that the reordering has to do real work.  The
        // spectrum of K spans (0, 24) and ω² = 9.3 lies inside it: the dynamic stiffness
        // is indefinite and its band is far wider than a chain's.
        let (nx, ny) = (24usize, 18usize);
        let n = nx * ny;
        let label: Vec<usize> = (0..n).map(|i| (i * 131) % n).collect();
        let mut sorted = label.clone();
        sorted.sort_unstable();
        assert!(
            sorted.iter().copied().eq(0..n),
            "the scrambling must be a permutation"
        );
        let mut k_triplets: Vec<(usize, usize, f64)> = Vec::new();
        let mut m_triplets: Vec<(usize, usize, f64)> = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                let p = label[j * nx + i];
                k_triplets.push((p, p, 12.0));
                m_triplets.push((p, p, 1.0 + 0.01 * (p % 7) as f64));
                if i + 1 < nx {
                    let q = label[j * nx + i + 1];
                    k_triplets.push((p, q, -3.0));
                    k_triplets.push((q, p, -3.0));
                }
                if j + 1 < ny {
                    let q = label[(j + 1) * nx + i];
                    k_triplets.push((p, q, -3.0));
                    k_triplets.push((q, p, -3.0));
                }
            }
        }
        let k_csr = CsrMatrix::from_triplets(n, n, &k_triplets);
        let m_csr = CsrMatrix::from_triplets(n, n, &m_triplets);
        let (omega, alpha, beta) = (9.3_f64.sqrt(), 0.02, 0.001);
        let mut force = vec![0.0; n];
        force[label[5 * nx + 7]] = 1.0;
        force[label[12 * nx + 19]] = -0.5;
        let r = harmonic_analysis(&k_csr, &m_csr, omega, alpha, beta, &force).expect("regular");
        let worst = block_residual(&k_csr, &m_csr, &r, (omega, alpha, beta), &force);
        assert!(worst < 1e-10, "residual {worst}");
        // The response is not trivial: it is damped and reaches the loaded DOFs.
        assert!(r.displacement_imag.iter().any(|v| v.abs() > 1e-4));
        assert!(r.displacement_amplitude[label[5 * nx + 7]] > 1e-3);
    }

    #[test]
    fn iterative_fallback_agrees_with_the_direct_solution() {
        let n = 30usize;
        let (k, m) = uniform_chain(n, 20.0, 1.0);
        let k_csr = from_dense(&k);
        let m_csr = from_dense(&m);
        let f: Vec<f64> = (0..n).map(|i| ((i % 5) as f64) - 1.5).collect();
        let (omega, alpha, beta) = (2.0, 0.05, 0.002);
        let direct =
            harmonic_analysis_with(&k_csr, &m_csr, omega, alpha, beta, &f, false).expect("direct");
        let iterative = harmonic_analysis_with(&k_csr, &m_csr, omega, alpha, beta, &f, true)
            .expect("iterative");
        for i in 0..n {
            assert!(
                (direct.displacement_real[i] - iterative.displacement_real[i]).abs() < 1e-7,
                "dof {i}"
            );
            assert!(
                (direct.displacement_imag[i] - iterative.displacement_imag[i]).abs() < 1e-7,
                "dof {i}"
            );
        }
    }

    #[test]
    fn undamped_excitation_exactly_at_a_natural_frequency_is_singular() {
        let k = from_dense(&[vec![1.0]]);
        let m = from_dense(&[vec![1.0]]);
        assert!(matches!(
            harmonic_analysis(&k, &m, 1.0, 0.0, 0.0, &[1.0]),
            Err(Error::SingularSystem)
        ));
        // A free-free spring chain analysed statically has a rigid-body mode.
        let free = from_dense(&[vec![1.0, -1.0], vec![-1.0, 1.0]]);
        let masses = from_dense(&[vec![1.0, 0.0], vec![0.0, 1.0]]);
        assert!(matches!(
            harmonic_analysis(&free, &masses, 0.0, 0.0, 0.0, &[1.0, 0.0]),
            Err(Error::SingularSystem)
        ));
        // ... but is regular at any ω that is not a natural frequency.
        let r = harmonic_analysis(&free, &masses, 0.5, 0.0, 0.0, &[1.0, 0.0]).expect("regular");
        assert!(r.displacement_real.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn invalid_input_is_reported_not_panicked() {
        let k = from_dense(&[vec![1.0, 0.0], vec![0.0, 1.0]]);
        let m = from_dense(&[vec![1.0, 0.0], vec![0.0, 1.0]]);
        let m1 = from_dense(&[vec![1.0]]);
        let rect = CsrMatrix::from_triplets(2, 3, &[(0, 0, 1.0), (1, 1, 1.0)]);
        // Wrong force length.
        assert!(matches!(
            harmonic_analysis(&k, &m, 1.0, 0.0, 0.0, &[1.0]),
            Err(Error::InvalidInput(_))
        ));
        // Mismatched mass size.
        assert!(matches!(
            harmonic_analysis(&k, &m1, 1.0, 0.0, 0.0, &[1.0, 0.0]),
            Err(Error::InvalidInput(_))
        ));
        // Non-square stiffness.
        assert!(matches!(
            harmonic_analysis(&rect, &m, 1.0, 0.0, 0.0, &[1.0, 0.0]),
            Err(Error::InvalidInput(_))
        ));
        // Non-finite parameters and data.
        assert!(matches!(
            harmonic_analysis(&k, &m, f64::NAN, 0.0, 0.0, &[1.0, 0.0]),
            Err(Error::InvalidInput(_))
        ));
        assert!(matches!(
            harmonic_analysis(&k, &m, 1.0, f64::INFINITY, 0.0, &[1.0, 0.0]),
            Err(Error::InvalidInput(_))
        ));
        assert!(matches!(
            harmonic_analysis(&k, &m, 1.0, 0.0, 0.0, &[f64::NAN, 0.0]),
            Err(Error::InvalidInput(_))
        ));
        let bad = CsrMatrix::from_triplets(2, 2, &[(0, 0, f64::INFINITY), (1, 1, 1.0)]);
        assert!(matches!(
            harmonic_analysis(&bad, &m, 1.0, 0.0, 0.0, &[1.0, 0.0]),
            Err(Error::InvalidInput(_))
        ));
    }

    #[test]
    fn malformed_csr_structure_is_reported_not_panicked() {
        let m = from_dense(&[vec![1.0, 0.0], vec![0.0, 1.0]]);
        let force = [1.0, 0.0];
        let well_formed = from_dense(&[vec![2.0, -1.0], vec![-1.0, 2.0]]);
        let corrupt = |edit: &dyn Fn(&mut CsrMatrix)| {
            let mut k = well_formed.clone();
            edit(&mut k);
            k
        };
        // Row pointers that go backwards (4 values are stored, so only the monotonicity
        // test can catch this) would make a row slice panic.
        let backwards = corrupt(&|k| k.row_ptr = vec![0, 5, 4]);
        // The first row pointer must be zero.
        let offset_start = corrupt(&|k| k.row_ptr = vec![1, 3, 4]);
        // The last row pointer must equal the number of stored values.
        let short_pointer = corrupt(&|k| k.row_ptr = vec![0, 2, 3]);
        // Column indices must address an existing DOF.
        let column_out_of_range = corrupt(&|k| k.col_indices[0] = 2);
        // Index and value arrays must have the same length.
        let missing_value = corrupt(&|k| {
            k.values.pop();
        });
        for (name, k) in [
            ("backwards", backwards),
            ("offset_start", offset_start),
            ("short_pointer", short_pointer),
            ("column_out_of_range", column_out_of_range),
            ("missing_value", missing_value),
        ] {
            assert!(
                matches!(
                    harmonic_analysis(&k, &m, 1.0, 0.0, 0.0, &force),
                    Err(Error::InvalidInput(_))
                ),
                "{name}: a corrupt stiffness matrix must be rejected"
            );
            assert!(
                matches!(
                    harmonic_analysis(&m, &k, 1.0, 0.0, 0.0, &force),
                    Err(Error::InvalidInput(_))
                ),
                "{name}: a corrupt mass matrix must be rejected"
            );
        }
        // The untouched matrix is accepted.
        assert!(harmonic_analysis(&well_formed, &m, 0.5, 0.0, 0.0, &force).is_ok());
    }

    #[test]
    fn empty_system_gives_an_empty_result() {
        let k = CsrMatrix::new(0, 0);
        let m = CsrMatrix::new(0, 0);
        let r = harmonic_analysis(&k, &m, 3.0, 0.1, 0.1, &[]).expect("empty system");
        assert_eq!(r.omega, 3.0);
        assert!(r.displacement_real.is_empty() && r.displacement_amplitude.is_empty());
    }
}
