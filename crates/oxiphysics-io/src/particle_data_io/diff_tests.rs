// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the particle delta encoder: sub-threshold motion must not be lost
//! and a reconstruction must carry the recorded mass, radius, step and time.

use super::{ParticleDataset, ParticleDelta, ParticleDiffWriter, make_test_dataset};

type Particle = (u64, [f64; 3], [f64; 3], f64, f64);

fn dataset(step: u64, time: f64, particles: &[Particle]) -> ParticleDataset {
    let mut ds = ParticleDataset::with_time(time, step);
    for &(id, position, velocity, mass, radius) in particles {
        ds.add_particle(id, position, velocity, mass, radius);
    }
    ds
}

fn assert_within(actual: [f64; 3], expected: [f64; 3], tolerance: f64, what: &str) {
    for (axis, (a, e)) in actual.iter().zip(&expected).enumerate() {
        assert!(
            (a - e).abs() <= tolerance,
            "{what}: component {axis} is {a}, expected {e} within {tolerance}"
        );
    }
}

fn make_ds(n: usize) -> ParticleDataset {
    make_test_dataset(n, 0.0)
}

#[test]
fn first_frame_marks_every_particle_new() {
    let mut w = ParticleDiffWriter::new();
    let n = w.write_frame(&make_ds(4));
    assert_eq!(n, 4);
    assert!(w.deltas[0].2.iter().all(|d| d.is_new && !d.is_removed));
}

#[test]
fn unchanged_frame_writes_nothing() {
    let mut w = ParticleDiffWriter::new();
    let ds = make_ds(3);
    w.write_frame(&ds);
    assert_eq!(w.write_frame(&ds), 0);
    assert!(w.deltas[1].2.is_empty());
}

#[test]
fn moved_particle_is_written_as_a_delta() {
    let mut w = ParticleDiffWriter::new();
    let ds1 = make_ds(2);
    let mut ds2 = ds1.clone();
    ds2.positions[0][0] += 1.0;
    w.write_frame(&ds1);
    assert_eq!(w.write_frame(&ds2), 1);
    let delta = &w.deltas[1].2[0];
    assert_eq!(delta.id, ds2.ids[0]);
    assert!(!delta.is_new);
    assert!((delta.delta_pos[0] - 1.0).abs() < 1e-12);
    assert_eq!(delta.delta_pos[1], 0.0);
}

#[test]
fn removed_particle_is_flagged_with_its_last_attributes() {
    let mut w = ParticleDiffWriter::new();
    let ds1 = dataset(
        0,
        0.0,
        &[
            (1, [0.0; 3], [0.0; 3], 1.0, 0.1),
            (2, [1.0; 3], [0.0; 3], 2.5, 0.3),
            (3, [2.0; 3], [0.0; 3], 1.0, 0.1),
        ],
    );
    let ds2 = dataset(
        1,
        0.1,
        &[
            (1, [0.0; 3], [0.0; 3], 1.0, 0.1),
            (3, [2.0; 3], [0.0; 3], 1.0, 0.1),
        ],
    );
    w.write_frame(&ds1);
    assert_eq!(w.write_frame(&ds2), 1);
    let removed = w.deltas[1].2.iter().find(|d| d.id == 2).unwrap();
    assert!(removed.is_removed && !removed.is_new);
    assert_eq!((removed.mass, removed.radius), (2.5, 0.3));
}

#[test]
fn removal_records_are_written_in_ascending_id_order() {
    let mut w = ParticleDiffWriter::new();
    let many = dataset(
        0,
        0.0,
        &[
            (9, [0.0; 3], [0.0; 3], 1.0, 0.1),
            (3, [0.0; 3], [0.0; 3], 1.0, 0.1),
            (7, [0.0; 3], [0.0; 3], 1.0, 0.1),
            (5, [0.0; 3], [0.0; 3], 1.0, 0.1),
        ],
    );
    let one = dataset(1, 0.1, &[(7, [0.0; 3], [0.0; 3], 1.0, 0.1)]);
    w.write_frame(&many);
    assert_eq!(w.write_frame(&one), 3);
    let order: Vec<u64> = w.deltas[1].2.iter().map(|d| d.id).collect();
    assert_eq!(order, vec![3, 5, 9]);
}

#[test]
fn added_particle_is_flagged_new() {
    let mut w = ParticleDiffWriter::new();
    let ds1 = make_ds(2);
    let mut ds2 = ds1.clone();
    ds2.add_particle(99, [10.0; 3], [0.0; 3], 1.0, 0.1);
    w.write_frame(&ds1);
    assert_eq!(w.write_frame(&ds2), 1);
    let added = &w.deltas[1].2[0];
    assert!(added.is_new);
    assert_eq!(added.id, 99);
    assert_eq!(added.delta_pos, [10.0; 3]);
}

#[test]
fn total_deltas_counts_all_records() {
    let mut w = ParticleDiffWriter::new();
    w.write_frame(&make_ds(5));
    assert_eq!(w.total_deltas(), 5);
    let mut moved = make_ds(5);
    moved.positions[2][1] += 3.0;
    w.write_frame(&moved);
    assert_eq!(w.total_deltas(), 6);
}

#[test]
fn single_sub_threshold_step_is_not_written() {
    let mut w = ParticleDiffWriter::with_thresholds(1.0, 1.0);
    let ds1 = make_ds(2);
    let mut ds2 = ds1.clone();
    ds2.positions[0][0] += 0.5;
    w.write_frame(&ds1);
    assert_eq!(w.write_frame(&ds2), 0);
}

#[test]
fn reconstruct_returns_every_recorded_particle() {
    let mut w = ParticleDiffWriter::new();
    let ds = make_ds(4);
    w.write_frame(&ds);
    let rec = w.reconstruct(0);
    assert_eq!(rec.len(), 4);
    assert_eq!(rec.ids, ds.ids);
    assert_eq!(rec.positions, ds.positions);
    assert_eq!(rec.velocities, ds.velocities);
}

#[test]
fn reconstruction_keeps_mass_radius_step_and_time() {
    let mut w = ParticleDiffWriter::new();
    w.write_frame(&dataset(
        40,
        1.25,
        &[
            (10, [0.0; 3], [0.0; 3], 2.5, 0.3),
            (11, [1.0, 2.0, 3.0], [0.5, 0.0, 0.0], 7.0, 0.05),
        ],
    ));
    w.write_frame(&dataset(
        50,
        1.75,
        &[
            (10, [0.5, 0.0, 0.0], [0.0; 3], 2.5, 0.3),
            (11, [1.0, 2.0, 3.0], [0.5, 0.0, 0.0], 7.0, 0.05),
        ],
    ));
    let rec = w.reconstruct(50);
    assert_eq!(rec.step, 50);
    assert_eq!(rec.time, 1.75);
    assert_eq!(rec.ids, vec![10, 11]);
    assert_eq!(rec.masses, vec![2.5, 7.0]);
    assert_eq!(rec.radii, vec![0.3, 0.05]);
    assert_eq!(rec.positions[0], [0.5, 0.0, 0.0]);
}

#[test]
fn reconstruct_stops_at_the_requested_step() {
    let mut w = ParticleDiffWriter::new();
    for (step, time, x) in [(10u64, 0.1, 0.0), (20, 0.2, 5.0), (30, 0.3, 9.0)] {
        w.write_frame(&dataset(
            step,
            time,
            &[(1, [x, 0.0, 0.0], [0.0; 3], 1.0, 0.1)],
        ));
    }
    let mid = w.reconstruct(25);
    assert_eq!((mid.step, mid.time), (20, 0.2));
    assert_eq!(mid.positions, vec![[5.0, 0.0, 0.0]]);
    let exact = w.reconstruct(10);
    assert_eq!((exact.step, exact.time), (10, 0.1));
    assert_eq!(exact.positions, vec![[0.0; 3]]);
    let late = w.reconstruct(1000);
    assert_eq!((late.step, late.time), (30, 0.3));
    assert_eq!(late.positions, vec![[9.0, 0.0, 0.0]]);
}

#[test]
fn reconstruct_before_the_first_frame_is_empty() {
    let mut w = ParticleDiffWriter::new();
    assert!(w.reconstruct(7).is_empty());
    w.write_frame(&dataset(10, 0.5, &[(1, [1.0; 3], [0.0; 3], 1.0, 0.1)]));
    let early = w.reconstruct(5);
    assert!(early.is_empty());
    assert_eq!(early.step, 5);
}

#[test]
fn slow_drift_is_not_lost_when_each_step_is_below_the_threshold() {
    let threshold = 1.0;
    let per_step = 0.4 * threshold;
    let mut w = ParticleDiffWriter::with_thresholds(threshold, threshold);
    let mut truth = Vec::new();
    for frame in 0..100u64 {
        let x = per_step * frame as f64;
        let position = [x, -0.5 * x, 3.0];
        truth.push(position);
        w.write_frame(&dataset(
            frame,
            0.01 * frame as f64,
            &[(1, position, [0.0; 3], 1.0, 0.1)],
        ));
    }
    // The reconstruction never strays further than the threshold from the
    // data that was recorded, at any step of the history.
    for (frame, expected) in truth.iter().enumerate() {
        let rec = w.reconstruct(frame as u64);
        assert_eq!(rec.len(), 1);
        assert_within(rec.positions[0], *expected, threshold, "drifting particle");
    }
    // 0.4 per frame crosses the threshold every third frame (1.2 > 1.0).
    for (frame, (_, _, records)) in w.deltas.iter().enumerate().skip(1) {
        let expected = usize::from(frame % 3 == 0);
        assert_eq!(records.len(), expected, "records in frame {frame}");
    }
    let last = w.reconstruct(99);
    assert_within(last.positions[0], truth[99], threshold, "final position");
    assert!(
        last.positions[0][0] > 30.0,
        "the particle must have travelled, got {:?}",
        last.positions[0]
    );
}

#[test]
fn half_threshold_steps_do_not_reconstruct_as_stationary() {
    // A particle that moves half the threshold every frame used to be
    // reconstructed at its starting point: no single step was ever written.
    let threshold = 0.02;
    let mut w = ParticleDiffWriter::with_thresholds(threshold, threshold);
    let mut final_position = [0.0; 3];
    for frame in 0..100u64 {
        final_position = [0.5 * threshold * frame as f64, 0.0, 1.0];
        w.write_frame(&dataset(
            frame,
            0.1 * frame as f64,
            &[(3, final_position, [0.0; 3], 1.0, 0.1)],
        ));
    }
    let rec = w.reconstruct(99);
    assert_within(
        rec.positions[0],
        final_position,
        threshold,
        "final position",
    );
    assert!(
        rec.positions[0][0] > 0.95,
        "the particle travelled 0.99 but was reconstructed at {:?}",
        rec.positions[0]
    );
    // One record roughly every three frames: far fewer than one per frame,
    // yet the reconstruction stays within the threshold at every step.
    assert!(w.total_deltas() < 50, "{} records", w.total_deltas());
    for frame in 0..100u64 {
        let truth = [0.5 * threshold * frame as f64, 0.0, 1.0];
        assert_within(
            w.reconstruct(frame).positions[0],
            truth,
            threshold,
            "history",
        );
    }
}

#[test]
fn slow_drift_of_the_velocity_is_not_lost_either() {
    let threshold = 0.5;
    let mut w = ParticleDiffWriter::with_thresholds(1e9, threshold);
    let mut truth = Vec::new();
    for frame in 0..60u64 {
        let v = [0.2 * frame as f64, 0.0, -0.1 * frame as f64];
        truth.push(v);
        w.write_frame(&dataset(frame, 0.0, &[(4, [0.0; 3], v, 1.0, 0.1)]));
    }
    for (frame, expected) in truth.iter().enumerate() {
        let rec = w.reconstruct(frame as u64);
        assert_within(rec.velocities[0], *expected, threshold, "velocity");
    }
}

#[test]
fn zero_thresholds_reproduce_the_history_to_rounding() {
    let mut w = ParticleDiffWriter::with_thresholds(0.0, 0.0);
    let mut truth = Vec::new();
    for frame in 0..200u64 {
        let t = frame as f64;
        let position = [100.0 + (0.1 * t).sin() * 37.1, 0.3 * t, 1e-3 * t * t];
        truth.push(position);
        w.write_frame(&dataset(frame, t, &[(1, position, [0.0; 3], 1.0, 0.1)]));
    }
    for (frame, expected) in truth.iter().enumerate() {
        let rec = w.reconstruct(frame as u64);
        assert_within(rec.positions[0], *expected, 1e-12, "closed loop");
    }
}

#[test]
fn mass_only_and_radius_only_changes_are_recorded() {
    let mut w = ParticleDiffWriter::new();
    let base = [(1u64, [1.0, 1.0, 1.0], [0.0; 3], 1.0, 0.1)];
    w.write_frame(&dataset(0, 0.0, &base));
    assert_eq!(
        w.write_frame(&dataset(1, 0.1, &[(1, [1.0; 3], [0.0; 3], 4.0, 0.1)])),
        1
    );
    assert_eq!(w.reconstruct(1).masses, vec![4.0]);
    assert_eq!(
        w.write_frame(&dataset(2, 0.2, &[(1, [1.0; 3], [0.0; 3], 4.0, 0.25)])),
        1
    );
    let rec = w.reconstruct(2);
    assert_eq!((rec.masses[0], rec.radii[0]), (4.0, 0.25));
    // An identical frame afterwards writes nothing.
    assert_eq!(
        w.write_frame(&dataset(3, 0.3, &[(1, [1.0; 3], [0.0; 3], 4.0, 0.25)])),
        0
    );
}

#[test]
fn velocity_change_above_the_threshold_is_recorded() {
    let mut w = ParticleDiffWriter::with_thresholds(1e-3, 1e-3);
    w.write_frame(&dataset(0, 0.0, &[(1, [0.0; 3], [0.0; 3], 1.0, 0.1)]));
    assert_eq!(
        w.write_frame(&dataset(
            1,
            0.1,
            &[(1, [0.0; 3], [5e-4, 0.0, 0.0], 1.0, 0.1)]
        )),
        0
    );
    assert_eq!(
        w.write_frame(&dataset(
            2,
            0.2,
            &[(1, [0.0; 3], [2e-3, 0.0, 0.0], 1.0, 0.1)]
        )),
        1
    );
    assert_eq!(w.reconstruct(2).velocities, vec![[2e-3, 0.0, 0.0]]);
}

#[test]
fn removed_and_readded_particle_restarts_from_absolute_values() {
    let mut w = ParticleDiffWriter::new();
    let present = |step: u64, x: f64| {
        dataset(
            step,
            0.1 * step as f64,
            &[(1, [x, 0.0, 0.0], [0.0; 3], 1.0, 0.1)],
        )
    };
    w.write_frame(&present(0, 1.0));
    w.write_frame(&dataset(1, 0.1, &[]));
    assert!(w.reconstruct(1).is_empty());
    w.write_frame(&present(2, 8.0));
    assert!(w.deltas[2].2[0].is_new);
    assert_eq!(w.reconstruct(2).positions, vec![[8.0, 0.0, 0.0]]);
    assert!(w.reconstruct(1).is_empty(), "history stays addressable");
}

#[test]
fn orphan_change_record_does_not_fabricate_a_particle() {
    let mut w = ParticleDiffWriter::new();
    w.deltas.push((
        0,
        0.0,
        vec![ParticleDelta {
            id: 77,
            delta_pos: [1.0, 2.0, 3.0],
            delta_vel: [0.0; 3],
            mass: 1.0,
            radius: 0.1,
            is_new: false,
            is_removed: false,
        }],
    ));
    assert!(w.reconstruct(0).is_empty());
}

#[test]
fn datasets_with_ragged_arrays_do_not_panic() {
    let mut ds = dataset(
        0,
        0.0,
        &[
            (1, [0.0; 3], [0.0; 3], 1.0, 0.1),
            (2, [1.0; 3], [0.0; 3], 1.0, 0.1),
            (3, [2.0; 3], [0.0; 3], 1.0, 0.1),
        ],
    );
    ds.masses.pop();
    let mut w = ParticleDiffWriter::new();
    assert_eq!(w.write_frame(&ds), 2);
    assert_eq!(w.reconstruct(0).ids, vec![1, 2]);
}

fn u64_at(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(bytes[offset..offset + 8].try_into().unwrap())
}

fn f64_at(bytes: &[u8], offset: usize) -> f64 {
    f64::from_le_bytes(bytes[offset..offset + 8].try_into().unwrap())
}

#[test]
fn wire_format_carries_time_mass_and_radius() {
    let mut w = ParticleDiffWriter::new();
    w.write_frame(&dataset(
        12,
        0.75,
        &[
            (5, [1.0, 2.0, 3.0], [4.0, 5.0, 6.0], 2.5, 0.3),
            (6, [7.0, 8.0, 9.0], [0.0; 3], 3.5, 0.4),
        ],
    ));
    w.write_frame(&dataset(
        13,
        0.80,
        &[
            (5, [1.5, 2.0, 3.0], [4.0, 5.0, 6.0], 2.5, 0.3),
            (6, [7.0, 8.0, 9.0], [0.0; 3], 3.5, 0.4),
        ],
    ));
    let bytes = w.to_bytes();
    assert_eq!(&bytes[..8], b"PDIFF\0\0\x02");
    assert_eq!(u64_at(&bytes, 8), 2, "frame count");

    const RECORD: usize = 8 + 1 + 6 * 8 + 2 * 8;
    // Frame 0 starts at 16: step, time, record count, then two records.
    assert_eq!(u64_at(&bytes, 16), 12);
    assert_eq!(f64_at(&bytes, 24), 0.75);
    assert_eq!(u64_at(&bytes, 32), 2);
    let rec0 = 40;
    assert_eq!(u64_at(&bytes, rec0), 5);
    assert_eq!(bytes[rec0 + 8], 0b01, "new flag");
    assert_eq!(f64_at(&bytes, rec0 + 9), 1.0);
    assert_eq!(f64_at(&bytes, rec0 + 9 + 3 * 8), 4.0);
    assert_eq!(f64_at(&bytes, rec0 + 9 + 6 * 8), 2.5, "mass");
    assert_eq!(f64_at(&bytes, rec0 + 9 + 7 * 8), 0.3, "radius");
    let rec1 = rec0 + RECORD;
    assert_eq!(u64_at(&bytes, rec1), 6);
    assert_eq!(f64_at(&bytes, rec1 + 9 + 6 * 8), 3.5);
    // Frame 1: one record, the 0.5 displacement of particle 5.
    let frame1 = rec1 + RECORD;
    assert_eq!(u64_at(&bytes, frame1), 13);
    assert_eq!(f64_at(&bytes, frame1 + 8), 0.80);
    assert_eq!(u64_at(&bytes, frame1 + 16), 1);
    let delta = frame1 + 24;
    assert_eq!(u64_at(&bytes, delta), 5);
    assert_eq!(bytes[delta + 8], 0);
    assert_eq!(f64_at(&bytes, delta + 9), 0.5);
    assert_eq!(bytes.len(), delta + RECORD);
}

#[test]
fn removed_flag_is_serialised_as_bit_one() {
    let mut w = ParticleDiffWriter::new();
    w.write_frame(&dataset(0, 0.0, &[(1, [0.0; 3], [0.0; 3], 1.0, 0.1)]));
    w.write_frame(&dataset(1, 0.1, &[]));
    let bytes = w.to_bytes();
    // header 16 + frame0 (24 + record) + frame1 header 24, then the record.
    const RECORD: usize = 8 + 1 + 6 * 8 + 2 * 8;
    let removed = 16 + 24 + RECORD + 24;
    assert_eq!(bytes[removed + 8], 0b10);
}
