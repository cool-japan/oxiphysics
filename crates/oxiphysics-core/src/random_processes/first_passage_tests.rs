//! Cross-checks of [`ContinuousTimeMarkov::mean_first_passage_time`] on
//! pseudo-random generators.
//!
//! The reference solution deletes the target state from the generator and
//! solves the remaining dense system `-Q_rest m = 1` with the LU factorisation
//! of `nalgebra`; it shares no code with the elimination used by the
//! implementation.  Recurrence times are compared with Kac's formula
//! `1 / (pi_j q_j)`, which needs the stationary distribution instead.

use super::{ContinuousTimeMarkov, Lcg};
use nalgebra::{DMatrix, DVector};

/// Irreducible generator on `n` states: a ring `i -> i + 1` of positive random
/// rates plus random chords (every other ordered pair is connected with
/// probability one half).
fn random_irreducible_ctmc(n: usize, seed: u64) -> ContinuousTimeMarkov {
    let mut rng = Lcg::new(seed);
    let mut q = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let on_ring = j == (i + 1) % n;
            let chord = rng.uniform() < 0.5;
            if on_ring || chord {
                q[i * n + j] = 0.1 + 3.0 * rng.uniform();
            }
        }
        let exit_rate: f64 = (0..n).filter(|&j| j != i).map(|j| q[i * n + j]).sum();
        q[i * n + i] = -exit_rate;
    }
    ContinuousTimeMarkov::new(q, n)
}

/// Mean first-passage times into `target`, indexed by start state (the entry of
/// `target` itself is zero).
fn dense_first_passage_reference(ctmc: &ContinuousTimeMarkov, target: usize) -> Vec<f64> {
    let n = ctmc.n_states;
    let others: Vec<usize> = (0..n).filter(|&i| i != target).collect();
    let r = others.len();
    let minus_q_rest = DMatrix::from_fn(r, r, |row, col| -ctmc.q[others[row] * n + others[col]]);
    let ones = DVector::from_element(r, 1.0);
    let solution = minus_q_rest
        .lu()
        .solve(&ones)
        .expect("the reduced generator of an irreducible chain is non-singular");
    let mut times = vec![0.0f64; n];
    for (position, &state) in others.iter().enumerate() {
        times[state] = solution[position];
    }
    times
}

#[test]
fn first_passage_times_match_a_dense_solve_for_every_pair_of_states() {
    for (n, seed) in [(4usize, 1u64), (6, 2), (6, 3), (7, 4), (9, 5)] {
        let ctmc = random_irreducible_ctmc(n, seed);
        for target in 0..n {
            let reference = dense_first_passage_reference(&ctmc, target);
            for from in (0..n).filter(|&i| i != target) {
                let got = ctmc.mean_first_passage_time(from, target);
                let want = reference[from];
                assert!(
                    want > 0.0 && ((got - want) / want).abs() < 1e-10,
                    "n = {n}, seed = {seed}, {from} -> {target}: got {got}, dense solve {want}"
                );
            }
        }
    }
}

#[test]
fn first_passage_times_depend_on_the_start_state() {
    // The old implementation returned the same number for every `from`.
    let ctmc = random_irreducible_ctmc(6, 2);
    for target in 0..6 {
        let times: Vec<f64> = (0..6)
            .filter(|&from| from != target)
            .map(|from| ctmc.mean_first_passage_time(from, target))
            .collect();
        let spread = times.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            - times.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(spread > 1e-3, "target {target}: times = {times:?}");
    }
}

#[test]
fn recurrence_times_match_kacs_formula_for_random_generators() {
    for (n, seed) in [(4usize, 1u64), (6, 2), (7, 4)] {
        let ctmc = random_irreducible_ctmc(n, seed);
        let pi = ctmc.stationary_distribution(1e-14, 1_000_000);
        for (j, &pi_j) in pi.iter().enumerate() {
            let kac = 1.0 / (pi_j * -ctmc.q[j * n + j]);
            let got = ctmc.mean_first_passage_time(j, j);
            assert!(
                ((got - kac) / kac).abs() < 1e-6,
                "n = {n}, seed = {seed}, state {j}: recurrence {got}, Kac {kac}"
            );
        }
    }
}

#[test]
fn recurrence_time_is_one_over_exit_rate_plus_the_mean_return_trip() {
    // E[T_jj] = 1/q_j + sum_k (q_jk / q_j) m_kj  with m_kj the first-passage times into j.
    let ctmc = random_irreducible_ctmc(6, 3);
    let n = ctmc.n_states;
    for j in 0..n {
        let exit_rate = -ctmc.q[j * n + j];
        let reference = dense_first_passage_reference(&ctmc, j);
        let return_trip: f64 = (0..n)
            .filter(|&k| k != j)
            .map(|k| ctmc.q[j * n + k] / exit_rate * reference[k])
            .sum();
        let expected = 1.0 / exit_rate + return_trip;
        let got = ctmc.mean_first_passage_time(j, j);
        assert!(
            ((got - expected) / expected).abs() < 1e-10,
            "state {j}: {got} vs {expected}"
        );
    }
}

#[test]
fn first_passage_times_obey_the_triangle_inequality() {
    // Reaching `c` from `a` can be no faster than reaching `b` first and continuing:
    // m(a -> c) <= m(a -> b) + m(b -> c).
    let ctmc = random_irreducible_ctmc(7, 4);
    for a in 0..7 {
        for b in (0..7).filter(|&b| b != a) {
            for c in (0..7).filter(|&c| c != a && c != b) {
                let direct = ctmc.mean_first_passage_time(a, c);
                let via = ctmc.mean_first_passage_time(a, b) + ctmc.mean_first_passage_time(b, c);
                assert!(
                    direct <= via * (1.0 + 1e-12),
                    "{a} -> {c} = {direct} exceeds {a} -> {b} -> {c} = {via}"
                );
            }
        }
    }
}
