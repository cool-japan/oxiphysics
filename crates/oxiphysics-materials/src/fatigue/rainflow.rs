//! Rainflow cycle counting (ASTM E1049) and Palmgren-Miner damage from counted cycles.
//!
//! A load history is reduced to its reversals (peaks and valleys) and counted
//! with the four-point algorithm, which is equivalent to the three-point stack
//! algorithm of ASTM E1049-85 §5.4.4 including its start-point rule: closed
//! cycles count as one full cycle, and every range left over at the end of the
//! history (the residue) counts as a half cycle.  Half cycles are therefore never
//! lost and the counts add up to `(reversals − 1) / 2`.

use super::types::{BasquinCurve, RainflowCycle};

/// Reduce a sampled history to its reversals (peaks and valleys).
///
/// Non-finite samples are ignored.  Runs of equal consecutive samples are
/// collapsed to a single sample, so a flat peak or valley is kept as one
/// reversal and a monotone stretch containing a plateau contributes no interior
/// reversal.  The first and last retained samples are always kept.
fn extract_turning_points(signal: &[f64]) -> Vec<f64> {
    let mut samples: Vec<f64> = Vec::with_capacity(signal.len());
    for &sample in signal {
        if sample.is_finite() && samples.last() != Some(&sample) {
            samples.push(sample);
        }
    }
    let (Some(&first), Some(&last)) = (samples.first(), samples.last()) else {
        return samples;
    };
    if samples.len() < 3 {
        return samples;
    }
    let mut reversals: Vec<f64> = Vec::with_capacity(samples.len());
    reversals.push(first);
    for window in samples.windows(3) {
        // Consecutive samples differ, so a reversal is a change of direction.
        if (window[1] > window[0]) != (window[2] > window[1]) {
            reversals.push(window[1]);
        }
    }
    reversals.push(last);
    reversals
}
/// Perform rainflow cycle counting on a stress/strain time series.
///
/// Reduces the history to its reversals (plateaus collapse, non-finite samples
/// are ignored) and applies the four-point rainflow algorithm, which gives the
/// same counts as the ASTM E1049 three-point stack algorithm.  For four
/// consecutive reversals `S1..S4` with
/// `ΔS1 = |S2 − S1|`, `ΔS2 = |S3 − S2|` and `ΔS3 = |S4 − S3|`, the inner range
/// is a closed cycle when `ΔS2 ≤ ΔS1` and `ΔS2 ≤ ΔS3`: it is counted as one full
/// cycle (`count = 1.0`) and `S2`, `S3` are discarded.  The ranges remaining on
/// the stack at the end of the history are counted as half cycles
/// (`count = 0.5`).
///
/// Every reversal-to-reversal step is accounted for exactly once, so
/// `Σ 2·count·range` equals the total variation of the reversal sequence.  For
/// the ASTM E1049 history `[−2, 1, −3, 5, −1, 3, −4, 4, −2]` the ranges
/// 3, 4, 6, 8 and 9 receive the counts 0.5, 1.5, 0.5, 1.0 and 0.5.  A history
/// with fewer than two distinct finite samples has no cycles.
pub fn rainflow_count(signal: &[f64]) -> Vec<RainflowCycle> {
    let reversals = extract_turning_points(signal);
    if reversals.len() < 2 {
        return Vec::new();
    }
    let mut stack: Vec<f64> = Vec::with_capacity(reversals.len());
    let mut cycles: Vec<RainflowCycle> = Vec::new();
    for &point in &reversals {
        stack.push(point);
        while stack.len() >= 4 {
            let n = stack.len();
            let (s1, s2, s3, s4) = (stack[n - 4], stack[n - 3], stack[n - 2], stack[n - 1]);
            let first = (s2 - s1).abs();
            let inner = (s3 - s2).abs();
            let last = (s4 - s3).abs();
            if last >= inner && inner <= first {
                cycles.push(RainflowCycle {
                    range: inner,
                    mean: 0.5 * (s2 + s3),
                    count: 1.0,
                });
                stack.truncate(n - 3);
                stack.push(s4);
            } else {
                break;
            }
        }
    }
    for pair in stack.windows(2) {
        cycles.push(RainflowCycle {
            range: (pair[1] - pair[0]).abs(),
            mean: 0.5 * (pair[0] + pair[1]),
            count: 0.5,
        });
    }
    cycles
}
/// Rainflow counting returning `(range, mean, count)` triples.
///
/// A thin mapping of [`rainflow_count`]: `count` is 1.0 for a closed cycle and
/// 0.5 for a half cycle of the residue.  The counts must be carried into damage
/// sums, see [`miner_damage_from_rainflow`].
pub fn rain_flow_count_with_half_cycles(signal: &[f64]) -> Vec<(f64, f64, f64)> {
    rainflow_count(signal)
        .into_iter()
        .map(|cycle| (cycle.range, cycle.mean, cycle.count))
        .collect()
}
/// Rainflow counting returning `(range, mean, count)` triples.
///
/// Identical to [`rain_flow_count_with_half_cycles`]: the ASTM E1049 count
/// consists of full cycles (`count = 1.0`) and half cycles (`count = 0.5`), so a
/// bare `(range, mean)` pair cannot represent it without losing damage.  A
/// two-point history is one half cycle.
pub fn rain_flow_count(signal: &[f64]) -> Vec<(f64, f64, f64)> {
    rain_flow_count_with_half_cycles(signal)
}
/// Compute Palmgren-Miner damage from rainflow-counted cycles using an S-N curve.
///
/// Each entry is `(range, mean, count)` as returned by
/// [`rain_flow_count_with_half_cycles`].  The stress amplitude `range / 2` gives
/// the cycles to failure `N_f` from the S-N curve and the damage is
/// `D = Σ count / N_f`, so half cycles contribute half the damage of a full
/// cycle of the same range.  Cycles with infinite life (at or below the
/// endurance limit) contribute nothing.
pub fn miner_damage_from_rainflow(cycles: &[(f64, f64, f64)], sn: &BasquinCurve) -> f64 {
    let mut damage = 0.0;
    for &(range, _mean, count) in cycles {
        let stress_amp = range / 2.0;
        let nf = sn.cycles_to_failure(stress_amp);
        if nf.is_finite() && nf > 0.0 {
            damage += count / nf;
        }
    }
    damage
}
#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    /// The ASTM E1049 example history.
    const ASTM_HISTORY: [f64; 9] = [-2.0, 1.0, -3.0, 5.0, -1.0, 3.0, -4.0, 4.0, -2.0];

    /// Sum the counts of equal ranges (integer ranges only).
    fn counts_by_range(cycles: &[(f64, f64, f64)]) -> BTreeMap<i64, f64> {
        let mut map = BTreeMap::new();
        for &(range, _, count) in cycles {
            *map.entry(range.round() as i64).or_insert(0.0) += count;
        }
        map
    }

    /// Deterministic generator for the property tests.
    struct Lcg(u64);

    impl Lcg {
        fn next_u32(&mut self) -> u32 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 33) as u32
        }
        fn below(&mut self, bound: u32) -> u32 {
            self.next_u32() % bound
        }
    }

    /// Independent ASTM E1049-85 §5.4.4 three-point stack algorithm with the
    /// start-point rule, used as an oracle for the four-point implementation.
    fn astm_three_point(signal: &[f64]) -> Vec<(f64, f64, f64)> {
        let mut samples: Vec<f64> = Vec::new();
        for &s in signal {
            if samples.last() != Some(&s) {
                samples.push(s);
            }
        }
        let mut reversals: Vec<f64> = Vec::new();
        if samples.len() < 3 {
            reversals = samples.clone();
        } else {
            reversals.push(samples[0]);
            for window in samples.windows(3) {
                let rising_in = window[1] > window[0];
                let rising_out = window[2] > window[1];
                if rising_in != rising_out {
                    reversals.push(window[1]);
                }
            }
            reversals.push(samples[samples.len() - 1]);
        }
        if reversals.len() < 2 {
            return Vec::new();
        }
        let mut stack: Vec<f64> = Vec::new();
        let mut start = 0usize;
        let mut cycles = Vec::new();
        for &point in &reversals {
            stack.push(point);
            while stack.len() - start >= 3 {
                let n = stack.len();
                let (s0, s1, s2) = (stack[n - 3], stack[n - 2], stack[n - 1]);
                let y = (s1 - s0).abs();
                let x = (s2 - s1).abs();
                if x < y {
                    break;
                }
                if n - 3 == start {
                    cycles.push((y, 0.5 * (s0 + s1), 0.5));
                    start += 1;
                } else {
                    cycles.push((y, 0.5 * (s0 + s1), 1.0));
                    stack.truncate(n - 3);
                    stack.push(s2);
                }
            }
        }
        for pair in stack[start..].windows(2) {
            cycles.push(((pair[1] - pair[0]).abs(), 0.5 * (pair[0] + pair[1]), 0.5));
        }
        cycles
    }

    /// Aggregate (range, mean) -> doubled count, for histories with half-integer means.
    fn aggregate(cycles: &[(f64, f64, f64)]) -> BTreeMap<(i64, i64), i64> {
        let mut map = BTreeMap::new();
        for &(range, mean, count) in cycles {
            let key = ((2.0 * range).round() as i64, (2.0 * mean).round() as i64);
            *map.entry(key).or_insert(0) += (2.0 * count).round() as i64;
        }
        map
    }

    #[test]
    fn start_of_history_half_cycles_are_counted_not_dropped() {
        // [0, 5, -5]: ranges 5 and 10 are both half cycles.  The old three-point loop
        // reported one range-5 cycle and lost range 10.
        let cycles = rain_flow_count_with_half_cycles(&[0.0, 5.0, -5.0]);
        assert_eq!(cycles, vec![(5.0, 2.5, 0.5), (10.0, 0.0, 0.5)]);
    }

    #[test]
    fn test_rain_flow_count_simple() {
        // [0, 1, -1, 0]: ASTM gives 0.5 @ 1, 0.5 @ 2, 0.5 @ 1 (the old code returned two
        // full cycles of range 1 and no range 2).
        let cycles = rain_flow_count(&[0.0, 1.0, -1.0, 0.0]);
        assert_eq!(
            cycles,
            vec![(1.0, 0.5, 0.5), (2.0, 0.0, 0.5), (1.0, -0.5, 0.5)]
        );
    }

    #[test]
    fn astm_e1049_example_history() {
        let cycles = rain_flow_count_with_half_cycles(&ASTM_HISTORY);
        let by_range = counts_by_range(&cycles);
        let expected: BTreeMap<i64, f64> = [(3, 0.5), (4, 1.5), (6, 0.5), (8, 1.0), (9, 0.5)]
            .into_iter()
            .collect();
        assert_eq!(by_range, expected);
        let total: f64 = cycles.iter().map(|c| c.2).sum();
        assert_eq!(total, 4.0);
        // Exactly one closed cycle: range 4 between the reversals -1 and 3 (mean 1).
        let full: Vec<_> = cycles.iter().filter(|c| c.2 == 1.0).collect();
        assert_eq!(full, vec![&(4.0, 1.0, 1.0)]);
    }

    #[test]
    fn rainflow_count_structs_match_the_triples() {
        let structs = rainflow_count(&ASTM_HISTORY);
        let triples = rain_flow_count_with_half_cycles(&ASTM_HISTORY);
        assert_eq!(structs.len(), triples.len());
        for (cycle, triple) in structs.iter().zip(&triples) {
            assert_eq!((cycle.range, cycle.mean, cycle.count), *triple);
        }
        assert_eq!(rain_flow_count(&ASTM_HISTORY), triples);
    }

    #[test]
    fn test_rain_flow_count_empty() {
        assert!(rain_flow_count(&[]).is_empty());
        assert!(rain_flow_count(&[3.0]).is_empty());
        assert_eq!(rain_flow_count(&[0.0, 1.0]), vec![(1.0, 0.5, 0.5)]);
        assert_eq!(rain_flow_count(&[2.0, 2.0, 1.0]), vec![(1.0, 1.5, 0.5)]);
    }

    #[test]
    fn test_rain_flow_count_constant() {
        // A flat history used to produce zero-range cycles in the four-point extractor.
        assert!(rain_flow_count(&[1.0; 10]).is_empty());
        assert!(rainflow_count(&[1.0; 10]).is_empty());
        assert!(rain_flow_count(&[7.0, 7.0]).is_empty());
    }

    #[test]
    fn test_rainflow_with_half_cycles() {
        // Reversals 0, 2, -1, 3, -2, 1 never close a cycle: five half cycles.
        let cycles = rain_flow_count_with_half_cycles(&[0.0, 2.0, -1.0, 3.0, -2.0, 1.0]);
        let ranges: Vec<f64> = cycles.iter().map(|c| c.0).collect();
        assert_eq!(ranges, vec![2.0, 3.0, 4.0, 5.0, 3.0]);
        assert!(cycles.iter().all(|c| c.2 == 0.5));
        let total: f64 = cycles.iter().map(|c| c.2).sum();
        assert_eq!(total, 2.5);
    }

    #[test]
    fn plateaus_collapse_to_a_single_reversal() {
        // 0 -> 5 -> (flat) -> 0 is the history [0, 5, 0]: two half cycles of range 5.
        let flat_peak = rain_flow_count(&[0.0, 5.0, 5.0, 5.0, 0.0]);
        assert_eq!(flat_peak, rain_flow_count(&[0.0, 5.0, 0.0]));
        assert_eq!(flat_peak, vec![(5.0, 2.5, 0.5), (5.0, 2.5, 0.5)]);
        // A plateau inside a monotone ramp is not a reversal: one half cycle of range 2.
        assert_eq!(
            rain_flow_count(&[0.0, 1.0, 1.0, 2.0]),
            vec![(2.0, 1.0, 0.5)]
        );
        // A repeated start sample does not create a zero-range half cycle.
        assert_eq!(rain_flow_count(&[0.0, 0.0, 1.0]), vec![(1.0, 0.5, 0.5)]);
    }

    #[test]
    fn non_finite_samples_are_ignored() {
        let clean = rain_flow_count(&[0.0, 10.0, -4.0, 6.0]);
        let dirty = rain_flow_count(&[0.0, f64::NAN, 10.0, -4.0, f64::INFINITY, 6.0]);
        assert_eq!(clean, dirty);
        assert!(rain_flow_count(&[f64::NAN, 1.0, f64::INFINITY]).is_empty());
    }

    #[test]
    fn four_point_count_equals_the_astm_three_point_algorithm_on_random_histories() {
        let mut rng = Lcg(0x5eed_1234_abcd_ef01);
        for trial in 0..20_000 {
            let length = rng.below(17) as usize;
            let signal: Vec<f64> = (0..length).map(|_| f64::from(rng.below(9)) - 4.0).collect();
            let ours = aggregate(&rain_flow_count(&signal));
            let reference = aggregate(&astm_three_point(&signal));
            assert_eq!(ours, reference, "trial {trial}: signal {signal:?}");
        }
    }

    #[test]
    fn every_reversal_step_is_counted_exactly_once() {
        let mut rng = Lcg(42);
        for trial in 0..2_000 {
            let length = 2 + rng.below(40) as usize;
            let signal: Vec<f64> = (0..length)
                .map(|_| f64::from(rng.below(2001)) * 0.01 - 10.0)
                .collect();
            let cycles = rain_flow_count(&signal);
            let reversals = extract_turning_points(&signal);
            let variation: f64 = reversals.windows(2).map(|w| (w[1] - w[0]).abs()).sum();
            let counted: f64 = cycles.iter().map(|c| 2.0 * c.2 * c.0).sum();
            assert!(
                (counted - variation).abs() <= 1e-9 * variation.max(1.0),
                "trial {trial}: counted variation {counted} vs {variation}"
            );
            let total_count: f64 = cycles.iter().map(|c| c.2).sum();
            let expected_count = 0.5 * (reversals.len() as f64 - 1.0);
            assert!(
                (total_count - expected_count).abs() < 1e-12,
                "trial {trial}: count {total_count} vs {expected_count}"
            );
            assert!(
                cycles
                    .iter()
                    .all(|c| c.0 > 0.0 && (c.2 == 0.5 || c.2 == 1.0))
            );
        }
    }

    #[test]
    fn test_miner_damage_from_rainflow() {
        let sn = BasquinCurve::new(1500.0e6, -0.1, 100.0e6);
        // Basquin: N_f = (sigma_a / A)^(1/b) with b = -0.1.
        let life = |amplitude: f64| (amplitude / 1500.0e6).powf(-10.0);
        let cycles = [(400.0e6, 200.0e6, 1.0), (200.0e6, 100.0e6, 0.5)];
        let damage = miner_damage_from_rainflow(&cycles, &sn);
        // 200 MPa amplitude is a full cycle; 100 MPa is exactly the endurance limit.
        let expected = 1.0 / life(200.0e6);
        assert!(
            (damage - expected).abs() < 1e-12 * expected,
            "D = {damage}, expected {expected}"
        );
        let above_endurance = [(400.0e6, 0.0, 1.0), (300.0e6, 0.0, 0.5)];
        let expected = 1.0 / life(200.0e6) + 0.5 / life(150.0e6);
        let damage = miner_damage_from_rainflow(&above_endurance, &sn);
        assert!((damage - expected).abs() < 1e-12 * expected);
    }

    #[test]
    fn miner_damage_of_half_cycles_is_half_of_full_cycles() {
        let sn = BasquinCurve::new(1500.0e6, -0.1, 100.0e6);
        let full = miner_damage_from_rainflow(&[(800.0e6, 0.0, 1.0)], &sn);
        let half = miner_damage_from_rainflow(&[(800.0e6, 0.0, 0.5)], &sn);
        assert!(full > 0.0);
        assert!((half - 0.5 * full).abs() < 1e-15 * full);
    }

    #[test]
    fn miner_damage_from_a_history_includes_the_larger_residual_range() {
        // [0, +400, -400] MPa: half cycles of range 400 and 800 MPa.  The old code
        // counted one range-400 cycle and lost the damaging 800 MPa range.
        let sn = BasquinCurve::new(1500.0e6, -0.1, 100.0e6);
        let cycles = rain_flow_count_with_half_cycles(&[0.0, 400.0e6, -400.0e6]);
        assert_eq!(cycles.len(), 2);
        let life = |amplitude: f64| (amplitude / 1500.0e6).powf(-10.0);
        let expected = 0.5 / life(200.0e6) + 0.5 / life(400.0e6);
        let damage = miner_damage_from_rainflow(&cycles, &sn);
        assert!(
            (damage - expected).abs() < 1e-12 * expected,
            "D = {damage}, expected {expected}"
        );
        // The 800 MPa half cycle dominates the damage.
        assert!(0.5 / life(400.0e6) > 100.0 * (0.5 / life(200.0e6)));
    }
}
