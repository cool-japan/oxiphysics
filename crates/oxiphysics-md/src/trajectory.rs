// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! MD trajectory analysis: RMSD, radius of gyration, MSD, end-to-end distance,
//! contact maps, and principal axes via Jacobi iteration.

/// A single frame of an MD trajectory.
#[derive(Debug, Clone)]
pub struct TrajectoryFrame {
    /// Atom positions (Å or nm, consistent units).
    pub positions: Vec<[f64; 3]>,
    /// Atom velocities (optional).
    pub velocities: Option<Vec<[f64; 3]>>,
    /// Simulation box lengths along x, y, z.
    pub box_lengths: [f64; 3],
    /// Simulation time for this frame.
    pub time: f64,
}

// ── helpers ────────────────────────────────────────────────────────────────

fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

// ── RMSD ──────────────────────────────────────────────────────────────────

/// Root mean square displacement between two sets of positions.
///
/// Both slices must have the same length (number of atoms).
/// No mass-weighting; all atoms contribute equally.
pub fn rmsd(ref_pos: &[[f64; 3]], cur_pos: &[[f64; 3]]) -> f64 {
    assert_eq!(
        ref_pos.len(),
        cur_pos.len(),
        "rmsd: position slices must have equal length"
    );
    let n = ref_pos.len();
    if n == 0 {
        return 0.0;
    }
    let sum_sq: f64 = ref_pos
        .iter()
        .zip(cur_pos.iter())
        .map(|(r, c)| {
            let d = vec3_sub(*c, *r);
            vec3_dot(d, d)
        })
        .sum();
    (sum_sq / n as f64).sqrt()
}

// ── Radius of gyration ─────────────────────────────────────────────────────

/// Mass-weighted radius of gyration.
///
/// `positions` and `masses` must have the same length.
pub fn radius_of_gyration(positions: &[[f64; 3]], masses: &[f64]) -> f64 {
    assert_eq!(
        positions.len(),
        masses.len(),
        "radius_of_gyration: positions and masses must have equal length"
    );
    let n = positions.len();
    if n == 0 {
        return 0.0;
    }

    // Mass-weighted centroid
    let total_mass: f64 = masses.iter().sum();
    let mut com = [0.0f64; 3];
    for (pos, &m) in positions.iter().zip(masses.iter()) {
        com[0] += m * pos[0];
        com[1] += m * pos[1];
        com[2] += m * pos[2];
    }
    com[0] /= total_mass;
    com[1] /= total_mass;
    com[2] /= total_mass;

    // Weighted sum of squared distances from COM
    let sum: f64 = positions
        .iter()
        .zip(masses.iter())
        .map(|(pos, &m)| {
            let d = vec3_sub(*pos, com);
            m * vec3_dot(d, d)
        })
        .sum();
    (sum / total_mass).sqrt()
}

// ── Mean square displacement ───────────────────────────────────────────────

/// Mean square displacement of a single atom over trajectory frames.
///
/// Returns a `Vec`f64` of length `frames.len()` where entry `t` is the
/// squared displacement from frame 0 to frame `t`.
pub fn mean_square_displacement(frames: &[TrajectoryFrame], atom_idx: usize) -> Vec<f64> {
    if frames.is_empty() {
        return Vec::new();
    }
    let ref_pos = frames[0].positions[atom_idx];
    frames
        .iter()
        .map(|f| {
            let d = vec3_sub(f.positions[atom_idx], ref_pos);
            vec3_dot(d, d)
        })
        .collect()
}

// ── End-to-end distance ────────────────────────────────────────────────────

/// Euclidean distance between atom `first` and atom `last` in `positions`.
pub fn end_to_end_distance(positions: &[[f64; 3]], first: usize, last: usize) -> f64 {
    let d = vec3_sub(positions[last], positions[first]);
    vec3_norm(d)
}

// ── Contact map ────────────────────────────────────────────────────────────

/// Boolean contact matrix: `map\[i\]\[j\]` is `true` when the distance between
/// atoms `i` and `j` is ≤ `cutoff`.
pub fn contact_map(positions: &[[f64; 3]], cutoff: f64) -> Vec<Vec<bool>> {
    let n = positions.len();
    let cutoff2 = cutoff * cutoff;
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    let d = vec3_sub(positions[j], positions[i]);
                    vec3_dot(d, d) <= cutoff2
                })
                .collect()
        })
        .collect()
}

// ── Principal axes via Jacobi iteration ───────────────────────────────────

/// Compute principal axes and moments of inertia for a set of (optionally
/// mass-weighted) positions.
///
/// Returns `(axes, eigenvalues)` where `axes\[k\]` is the k-th principal axis
/// (unit vector) and `eigenvalues\[k\]` is the corresponding moment, sorted in
/// ascending order.
///
/// Uses the classical Jacobi iterative method for the symmetric 3×3 inertia
/// tensor.  Convergence is typically reached in < 100 sweeps.
pub fn principal_axes(positions: &[[f64; 3]], masses: &[f64]) -> ([[f64; 3]; 3], [f64; 3]) {
    assert_eq!(positions.len(), masses.len());

    // Mass-weighted centroid
    let total_mass: f64 = masses.iter().sum();
    let mut com = [0.0f64; 3];
    for (pos, &m) in positions.iter().zip(masses.iter()) {
        for k in 0..3 {
            com[k] += m * pos[k];
        }
    }
    for c in &mut com {
        *c /= total_mass;
    }

    // Build inertia tensor I (3×3, flattened row-major)
    let mut i_mat = [[0.0f64; 3]; 3];
    for (pos, &m) in positions.iter().zip(masses.iter()) {
        let r = [pos[0] - com[0], pos[1] - com[1], pos[2] - com[2]];
        let r2 = vec3_dot(r, r);
        for a in 0..3 {
            for b in 0..3 {
                let delta = if a == b { 1.0 } else { 0.0 };
                i_mat[a][b] += m * (r2 * delta - r[a] * r[b]);
            }
        }
    }

    // Jacobi diagonalisation of the symmetric 3×3 matrix
    jacobi3x3(i_mat)
}

/// Jacobi iterative diagonalisation for a real symmetric 3×3 matrix.
/// Returns `(eigenvectors_as_rows, eigenvalues)` sorted ascending by eigenvalue.
fn jacobi3x3(mut a: [[f64; 3]; 3]) -> ([[f64; 3]; 3], [f64; 3]) {
    // V accumulates the rotation matrix (columns = eigenvectors)
    let mut v = [[0.0f64; 3]; 3];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }

    for _ in 0..100 {
        // Find the largest off-diagonal element
        let mut p = 0usize;
        let mut q = 1usize;
        let mut max_val = a[0][1].abs();
        for (pp, qq, val) in [(0usize, 2usize, a[0][2].abs()), (1, 2, a[1][2].abs())] {
            if val > max_val {
                max_val = val;
                p = pp;
                q = qq;
            }
        }
        if max_val < 1e-12 {
            break;
        }

        // Compute Jacobi rotation angle
        let theta = if (a[q][q] - a[p][p]).abs() < 1e-14 {
            std::f64::consts::FRAC_PI_4
        } else {
            0.5 * ((2.0 * a[p][q]) / (a[q][q] - a[p][p])).atan()
        };
        let c = theta.cos();
        let s = theta.sin();

        // Apply rotation: A' = G^T A G  (in-place)
        let a_pp = c * c * a[p][p] - 2.0 * s * c * a[p][q] + s * s * a[q][q];
        let a_qq = s * s * a[p][p] + 2.0 * s * c * a[p][q] + c * c * a[q][q];
        let a_pq = 0.0f64; // by construction
        a[p][p] = a_pp;
        a[q][q] = a_qq;
        a[p][q] = a_pq;
        a[q][p] = a_pq;

        // Update remaining rows/columns
        let r = if p == 0 && q == 1 {
            2
        } else if p == 0 && q == 2 {
            1
        } else {
            0
        };
        let a_rp = c * a[r][p] - s * a[r][q];
        let a_rq = s * a[r][p] + c * a[r][q];
        a[r][p] = a_rp;
        a[p][r] = a_rp;
        a[r][q] = a_rq;
        a[q][r] = a_rq;

        // Accumulate eigenvectors
        for row in v.iter_mut() {
            let v_ip = c * row[p] - s * row[q];
            let v_iq = s * row[p] + c * row[q];
            row[p] = v_ip;
            row[q] = v_iq;
        }
    }

    // Eigenvalues are on the diagonal; eigenvectors are columns of v
    let mut pairs: [(f64, [f64; 3]); 3] = [
        (a[0][0], [v[0][0], v[1][0], v[2][0]]),
        (a[1][1], [v[0][1], v[1][1], v[2][1]]),
        (a[2][2], [v[0][2], v[1][2], v[2][2]]),
    ];
    // Sort ascending by eigenvalue
    pairs.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));

    let axes = [pairs[0].1, pairs[1].1, pairs[2].1];
    let evals = [pairs[0].0, pairs[1].0, pairs[2].0];
    (axes, evals)
}

// ── Frame interpolation ────────────────────────────────────────────────────

/// Linear interpolation between two trajectory frames.
///
/// Returns a new frame at fractional time `t` ∈ [0, 1] between `frame_a`
/// (at t=0) and `frame_b` (at t=1).
pub fn interpolate_frames(
    frame_a: &TrajectoryFrame,
    frame_b: &TrajectoryFrame,
    t: f64,
) -> TrajectoryFrame {
    assert_eq!(
        frame_a.positions.len(),
        frame_b.positions.len(),
        "interpolation requires same number of atoms"
    );
    let n = frame_a.positions.len();
    let t_clamped = t.clamp(0.0, 1.0);
    let one_minus_t = 1.0 - t_clamped;

    let positions: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            [
                one_minus_t * frame_a.positions[i][0] + t_clamped * frame_b.positions[i][0],
                one_minus_t * frame_a.positions[i][1] + t_clamped * frame_b.positions[i][1],
                one_minus_t * frame_a.positions[i][2] + t_clamped * frame_b.positions[i][2],
            ]
        })
        .collect();

    let velocities = match (&frame_a.velocities, &frame_b.velocities) {
        (Some(va), Some(vb)) => Some(
            (0..n)
                .map(|i| {
                    [
                        one_minus_t * va[i][0] + t_clamped * vb[i][0],
                        one_minus_t * va[i][1] + t_clamped * vb[i][1],
                        one_minus_t * va[i][2] + t_clamped * vb[i][2],
                    ]
                })
                .collect(),
        ),
        _ => None,
    };

    let box_lengths = [
        one_minus_t * frame_a.box_lengths[0] + t_clamped * frame_b.box_lengths[0],
        one_minus_t * frame_a.box_lengths[1] + t_clamped * frame_b.box_lengths[1],
        one_minus_t * frame_a.box_lengths[2] + t_clamped * frame_b.box_lengths[2],
    ];

    let time = one_minus_t * frame_a.time + t_clamped * frame_b.time;

    TrajectoryFrame {
        positions,
        velocities,
        box_lengths,
        time,
    }
}

// ── Trajectory alignment (centroid alignment) ─────────────────────────────

/// Translate positions so that their centroid is at the origin.
///
/// Returns the centroid that was subtracted.
pub fn center_positions(positions: &mut [[f64; 3]]) -> [f64; 3] {
    let n = positions.len();
    if n == 0 {
        return [0.0; 3];
    }
    let inv_n = 1.0 / n as f64;
    let mut centroid = [0.0f64; 3];
    for pos in positions.iter() {
        centroid[0] += pos[0];
        centroid[1] += pos[1];
        centroid[2] += pos[2];
    }
    centroid[0] *= inv_n;
    centroid[1] *= inv_n;
    centroid[2] *= inv_n;

    for pos in positions.iter_mut() {
        pos[0] -= centroid[0];
        pos[1] -= centroid[1];
        pos[2] -= centroid[2];
    }

    centroid
}

/// Align `mobile` positions to `reference` by translating centroids.
///
/// Modifies `mobile` in place. Returns the RMSD after alignment.
pub fn align_by_centroid(reference: &[[f64; 3]], mobile: &mut [[f64; 3]]) -> f64 {
    assert_eq!(reference.len(), mobile.len());
    let n = reference.len();
    if n == 0 {
        return 0.0;
    }

    // Compute centroids
    let inv_n = 1.0 / n as f64;
    let mut ref_com = [0.0f64; 3];
    let mut mob_com = [0.0f64; 3];
    for i in 0..n {
        for k in 0..3 {
            ref_com[k] += reference[i][k];
            mob_com[k] += mobile[i][k];
        }
    }
    for k in 0..3 {
        ref_com[k] *= inv_n;
        mob_com[k] *= inv_n;
    }

    // Translate mobile to match reference centroid
    let shift = [
        ref_com[0] - mob_com[0],
        ref_com[1] - mob_com[1],
        ref_com[2] - mob_com[2],
    ];
    for m in mobile.iter_mut().take(n) {
        for (mk, &s) in m.iter_mut().zip(shift.iter()) {
            *mk += s;
        }
    }

    rmsd(reference, mobile)
}

// ── PBC unwrapping ────────────────────────────────────────────────────────

/// Unwrap trajectory positions to remove PBC jumps.
///
/// When an atom crosses a periodic boundary between consecutive frames,
/// the position jump is larger than half the box length. This function
/// detects such jumps and applies the minimum-image correction to produce
/// a continuous trajectory.
///
/// Modifies `frames` in place.
pub fn unwrap_pbc(frames: &mut [TrajectoryFrame]) {
    if frames.len() < 2 {
        return;
    }
    let n_atoms = frames[0].positions.len();

    for t in 1..frames.len() {
        let box_l = frames[t].box_lengths;
        for i in 0..n_atoms {
            let prev_pos = frames[t - 1].positions[i];
            for (k, (cur, &bl)) in frames[t].positions[i]
                .iter_mut()
                .zip(box_l.iter())
                .enumerate()
            {
                let _ = k;
                let prev = prev_pos[k];
                let mut diff = *cur - prev;
                if bl > 0.0 {
                    // Apply minimum image convention
                    while diff > 0.5 * bl {
                        diff -= bl;
                    }
                    while diff < -0.5 * bl {
                        diff += bl;
                    }
                }
                *cur = prev + diff;
            }
        }
    }
}

// ── Bond angle computation ────────────────────────────────────────────────

/// Compute bond angle (in radians) between three atoms i-j-k.
///
/// The angle is at atom j: angle(r_ji, r_jk).
pub fn bond_angle(positions: &[[f64; 3]], i: usize, j: usize, k: usize) -> f64 {
    let rji = vec3_sub(positions[i], positions[j]);
    let rjk = vec3_sub(positions[k], positions[j]);
    let dot = vec3_dot(rji, rjk);
    let norm_ji = vec3_norm(rji);
    let norm_jk = vec3_norm(rjk);
    if norm_ji < 1e-15 || norm_jk < 1e-15 {
        return 0.0;
    }
    let cos_theta = (dot / (norm_ji * norm_jk)).clamp(-1.0, 1.0);
    cos_theta.acos()
}

// ── Dihedral angle computation ────────────────────────────────────────────

fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Compute dihedral angle (in radians) for atoms i-j-k-l.
///
/// Uses the atan2 convention to return values in (-π, π].
pub fn dihedral_angle(positions: &[[f64; 3]], i: usize, j: usize, k: usize, l: usize) -> f64 {
    let b1 = vec3_sub(positions[j], positions[i]);
    let b2 = vec3_sub(positions[k], positions[j]);
    let b3 = vec3_sub(positions[l], positions[k]);

    let n1 = vec3_cross(b1, b2);
    let n2 = vec3_cross(b2, b3);

    let m1 = vec3_cross(n1, b2);
    let b2_norm = vec3_norm(b2);
    if b2_norm < 1e-15 {
        return 0.0;
    }
    let m1_scaled = [m1[0] / b2_norm, m1[1] / b2_norm, m1[2] / b2_norm];

    let x = vec3_dot(n1, n2);
    let y = vec3_dot(m1_scaled, n2);
    (-y).atan2(x)
}

// ── Autocorrelation function ──────────────────────────────────────────────

/// Compute the normalized autocorrelation function of a 1D time series.
///
/// C(t) = ⟨(x(t') − ⟨x⟩)(x(t'+t) − ⟨x⟩)⟩ / Var(x)
///
/// Returns a vector of length `max_lag + 1` with C(0), C(1), ..., C(max_lag).
/// C(0) = 1 by definition.
pub fn autocorrelation(data: &[f64], max_lag: usize) -> Vec<f64> {
    let n = data.len();
    if n == 0 {
        return vec![];
    }
    let mean = data.iter().sum::<f64>() / n as f64;
    let var: f64 = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n as f64;
    if var.abs() < 1e-30 {
        return vec![1.0; (max_lag + 1).min(n)];
    }

    let max_t = max_lag.min(n - 1);
    let mut acf = Vec::with_capacity(max_t + 1);

    for lag in 0..=max_t {
        let count = n - lag;
        let c: f64 = (0..count)
            .map(|t| (data[t] - mean) * (data[t + lag] - mean))
            .sum::<f64>()
            / count as f64;
        acf.push(c / var);
    }

    acf
}

// ── RMSD-based clustering ─────────────────────────────────────────────────

/// Cluster trajectory frames by RMSD using a simple leader algorithm.
///
/// Assigns each frame to the first cluster whose representative is within
/// `cutoff` RMSD. If no cluster is within cutoff, a new cluster is created.
///
/// Returns a vector of cluster indices (length = number of frames) and
/// the list of representative frame indices.
pub fn rmsd_cluster(frames: &[TrajectoryFrame], cutoff: f64) -> (Vec<usize>, Vec<usize>) {
    if frames.is_empty() {
        return (vec![], vec![]);
    }
    let mut assignments = Vec::with_capacity(frames.len());
    let mut representatives: Vec<usize> = Vec::new();

    for (i, frame) in frames.iter().enumerate() {
        let mut assigned = false;
        for (c_idx, &rep_idx) in representatives.iter().enumerate() {
            let r = rmsd(&frames[rep_idx].positions, &frame.positions);
            if r <= cutoff {
                assignments.push(c_idx);
                assigned = true;
                break;
            }
        }
        if !assigned {
            let new_cluster_id = representatives.len();
            representatives.push(i);
            assignments.push(new_cluster_id);
        }
    }
    (assignments, representatives)
}

/// Cluster population counts.
///
/// Returns a vector of counts, one per cluster, in cluster index order.
pub fn cluster_populations(assignments: &[usize], n_clusters: usize) -> Vec<usize> {
    let mut counts = vec![0usize; n_clusters];
    for &a in assignments {
        if a < n_clusters {
            counts[a] += 1;
        }
    }
    counts
}

// ── PCA of trajectory ─────────────────────────────────────────────────────

/// Flatten a trajectory frame into a coordinate vector of length 3N.
pub fn frame_to_vector(frame: &TrajectoryFrame) -> Vec<f64> {
    frame
        .positions
        .iter()
        .flat_map(|&[x, y, z]| [x, y, z])
        .collect()
}

/// Compute the covariance matrix of a set of trajectory frames (after centering).
///
/// Each frame is converted to a 3N coordinate vector.  The covariance matrix
/// is of size 3N × 3N and is returned row-major as a `Vec<Vec`f64`>`.
///
/// This is an expensive O(n_frames × (3N)²) operation; use for small systems.
pub fn trajectory_covariance(frames: &[TrajectoryFrame]) -> Vec<Vec<f64>> {
    if frames.is_empty() {
        return vec![];
    }
    let n_atoms = frames[0].positions.len();
    let dim = 3 * n_atoms;

    // Convert to vectors
    let vecs: Vec<Vec<f64>> = frames.iter().map(frame_to_vector).collect();
    let n_frames = vecs.len();

    // Compute mean
    let mut mean = vec![0.0f64; dim];
    for v in &vecs {
        for (d, &x) in v.iter().enumerate() {
            mean[d] += x;
        }
    }
    for m in &mut mean {
        *m /= n_frames as f64;
    }

    // Build covariance matrix
    let mut cov = vec![vec![0.0f64; dim]; dim];
    for v in &vecs {
        let delta: Vec<f64> = v.iter().zip(mean.iter()).map(|(&x, &m)| x - m).collect();
        for i in 0..dim {
            for j in i..dim {
                cov[i][j] += delta[i] * delta[j];
            }
        }
    }
    for i in 0..dim {
        let (top, bot) = cov.split_at_mut(i + 1);
        top[i][i] /= n_frames as f64;
        for (jj, row_j) in bot.iter_mut().enumerate() {
            let j = i + 1 + jj;
            top[i][j] /= n_frames as f64;
            row_j[i] = top[i][j];
        }
    }
    cov
}

/// Project trajectory frames onto the first `n_pcs` principal components
/// (columns of `pc_matrix`, each of length `dim`).
///
/// Returns a Vec of length `n_frames`, each entry being a Vec of length `n_pcs`.
pub fn project_onto_pcs(frames: &[TrajectoryFrame], pcs: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if frames.is_empty() || pcs.is_empty() {
        return vec![];
    }
    let mean_pos: Vec<f64> = {
        let dim = 3 * frames[0].positions.len();
        let mut m = vec![0.0f64; dim];
        for f in frames {
            for (d, x) in frame_to_vector(f).iter().enumerate() {
                m[d] += x;
            }
        }
        let nf = frames.len() as f64;
        m.iter().map(|&s| s / nf).collect()
    };

    frames
        .iter()
        .map(|frame| {
            let v = frame_to_vector(frame);
            let delta: Vec<f64> = v
                .iter()
                .zip(mean_pos.iter())
                .map(|(&x, &m)| x - m)
                .collect();
            pcs.iter()
                .map(|pc| {
                    delta
                        .iter()
                        .zip(pc.iter())
                        .map(|(&d, &p)| d * p)
                        .sum::<f64>()
                })
                .collect()
        })
        .collect()
}

// ── Dynamic Cross-Correlation Matrix (DCCM) ───────────────────────────────

/// Compute the dynamic cross-correlation matrix (DCCM) of atomic positions.
///
/// C_ij = ⟨Δr_i · Δr_j⟩ / sqrt(⟨|Δr_i|²⟩ · ⟨|Δr_j|²⟩)
///
/// Returns an N×N correlation matrix where N is the number of atoms.
/// Values in [−1, 1] where +1 means fully correlated motion.
pub fn dynamic_cross_correlation(frames: &[TrajectoryFrame]) -> Vec<Vec<f64>> {
    let n_frames = frames.len();
    if n_frames < 2 {
        return vec![];
    }
    let n_atoms = frames[0].positions.len();
    if n_atoms == 0 {
        return vec![];
    }

    // Compute mean position for each atom
    let mut mean_pos = vec![[0.0f64; 3]; n_atoms];
    for frame in frames {
        for (i, &pos) in frame.positions.iter().enumerate() {
            mean_pos[i][0] += pos[0];
            mean_pos[i][1] += pos[1];
            mean_pos[i][2] += pos[2];
        }
    }
    let nf = n_frames as f64;
    for mp in &mut mean_pos {
        mp[0] /= nf;
        mp[1] /= nf;
        mp[2] /= nf;
    }

    // Cross-correlation: ⟨Δri · Δrj⟩ and ⟨|Δri|²⟩
    let mut cov = vec![vec![0.0f64; n_atoms]; n_atoms];
    let mut var = vec![0.0f64; n_atoms];
    for frame in frames {
        let deltas: Vec<[f64; 3]> = (0..n_atoms)
            .map(|i| {
                [
                    frame.positions[i][0] - mean_pos[i][0],
                    frame.positions[i][1] - mean_pos[i][1],
                    frame.positions[i][2] - mean_pos[i][2],
                ]
            })
            .collect();
        for (i, (&di, var_i)) in deltas.iter().zip(var.iter_mut()).enumerate() {
            for (j, &dj) in deltas.iter().enumerate().skip(i) {
                let dot = di[0] * dj[0] + di[1] * dj[1] + di[2] * dj[2];
                cov[i][j] += dot;
                if j > i {
                    cov[j][i] += dot;
                }
            }
            *var_i += di[0].powi(2) + di[1].powi(2) + di[2].powi(2);
        }
    }
    for (cov_row, var_i) in cov.iter_mut().zip(var.iter_mut()) {
        for c in cov_row.iter_mut() {
            *c /= nf;
        }
        *var_i /= nf;
    }

    // Normalise
    let mut dccm = vec![vec![0.0f64; n_atoms]; n_atoms];
    for i in 0..n_atoms {
        for j in 0..n_atoms {
            let denom = (var[i] * var[j]).sqrt();
            dccm[i][j] = if denom > 1e-30 {
                cov[i][j] / denom
            } else {
                0.0
            };
        }
    }
    dccm
}

// ── Contact frequency map ─────────────────────────────────────────────────

/// Compute the contact frequency map from a trajectory.
///
/// For each pair (i, j), counts the fraction of frames in which atoms i and j
/// are within `cutoff` of each other.
///
/// Returns an N×N matrix of contact frequencies ∈ \[0, 1\].
pub fn contact_frequency_map(frames: &[TrajectoryFrame], cutoff: f64) -> Vec<Vec<f64>> {
    if frames.is_empty() {
        return vec![];
    }
    let n = frames[0].positions.len();
    let cutoff2 = cutoff * cutoff;
    let mut freq = vec![vec![0.0f64; n]; n];

    for frame in frames {
        for (i, freq_row) in freq.iter_mut().enumerate() {
            for (j, f) in freq_row.iter_mut().enumerate() {
                let d = vec3_sub(frame.positions[j], frame.positions[i]);
                if vec3_dot(d, d) <= cutoff2 {
                    *f += 1.0;
                }
            }
        }
    }
    let nf = frames.len() as f64;
    for row in &mut freq {
        for v in row.iter_mut() {
            *v /= nf;
        }
    }
    freq
}

// ── Representative structure from trajectory ──────────────────────────────

/// Extract the structure (frame) closest to the trajectory mean.
///
/// Computes the coordinate-space mean position and returns the index of
/// the frame whose positions are closest (min RMSD) to the mean.
pub fn representative_structure(frames: &[TrajectoryFrame]) -> Option<usize> {
    if frames.is_empty() {
        return None;
    }
    let n_atoms = frames[0].positions.len();
    if n_atoms == 0 {
        return Some(0);
    }

    // Compute mean positions
    let mut mean_pos = vec![[0.0f64; 3]; n_atoms];
    let nf = frames.len() as f64;
    for frame in frames {
        for (i, &pos) in frame.positions.iter().enumerate() {
            mean_pos[i][0] += pos[0] / nf;
            mean_pos[i][1] += pos[1] / nf;
            mean_pos[i][2] += pos[2] / nf;
        }
    }

    // Find frame closest to mean
    let (best_idx, _) = frames
        .iter()
        .enumerate()
        .map(|(i, frame)| {
            let r = rmsd(&mean_pos, &frame.positions);
            (i, r)
        })
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))?;

    Some(best_idx)
}

/// Compute per-atom positional fluctuations (root mean square fluctuation, RMSF).
///
/// RMSF_i = sqrt(⟨|r_i − ⟨r_i⟩|²⟩)
///
/// Returns a vector of length N.
pub fn rmsf(frames: &[TrajectoryFrame]) -> Vec<f64> {
    if frames.is_empty() {
        return vec![];
    }
    let n_atoms = frames[0].positions.len();
    let nf = frames.len() as f64;

    let mut mean_pos = vec![[0.0f64; 3]; n_atoms];
    for frame in frames {
        for (i, &pos) in frame.positions.iter().enumerate() {
            mean_pos[i][0] += pos[0];
            mean_pos[i][1] += pos[1];
            mean_pos[i][2] += pos[2];
        }
    }
    for mp in &mut mean_pos {
        mp[0] /= nf;
        mp[1] /= nf;
        mp[2] /= nf;
    }

    let mut sum_sq = vec![0.0f64; n_atoms];
    for frame in frames {
        for (i, &pos) in frame.positions.iter().enumerate() {
            let d = vec3_sub(pos, mean_pos[i]);
            sum_sq[i] += vec3_dot(d, d);
        }
    }
    sum_sq.iter().map(|&s| (s / nf).sqrt()).collect()
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: build a regular n-atom ring of radius R in the xy-plane
    fn ring_positions(n: usize, radius: f64) -> Vec<[f64; 3]> {
        (0..n)
            .map(|i| {
                let angle = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
                [radius * angle.cos(), radius * angle.sin(), 0.0]
            })
            .collect()
    }

    #[test]
    fn test_rmsd_identical() {
        let pos: Vec<[f64; 3]> = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [0.0, 0.0, 0.0]];
        assert!(
            rmsd(&pos, &pos) < 1e-14,
            "RMSD of identical configs must be 0"
        );
    }

    #[test]
    fn test_rmsd_known() {
        let ref_pos: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0]];
        let cur_pos: Vec<[f64; 3]> = vec![[3.0, 4.0, 0.0]];
        let r = rmsd(&ref_pos, &cur_pos);
        assert!((r - 5.0).abs() < 1e-12, "expected 5.0, got {r}");
    }

    #[test]
    fn test_radius_of_gyration_uniform_sphere() {
        // 6 atoms at ±R along each axis with equal masses → Rg = R
        let r = 2.0f64;
        let positions: Vec<[f64; 3]> = vec![
            [r, 0.0, 0.0],
            [-r, 0.0, 0.0],
            [0.0, r, 0.0],
            [0.0, -r, 0.0],
            [0.0, 0.0, r],
            [0.0, 0.0, -r],
        ];
        let masses = vec![1.0f64; 6];
        let rg = radius_of_gyration(&positions, &masses);
        assert!(
            (rg - r).abs() < 1e-12,
            "Rg of uniform sphere: expected {r}, got {rg}"
        );
    }

    #[test]
    fn test_msd_linear_free_diffusion() {
        // Simulate free diffusion: position at time t is (D*t, 0, 0)
        // MSD should equal (D*t)^2 (ballistic here, but checks linearity)
        let d = 1.5f64;
        let frames: Vec<TrajectoryFrame> = (0..10)
            .map(|t| TrajectoryFrame {
                positions: vec![[d * t as f64, 0.0, 0.0]],
                velocities: None,
                box_lengths: [100.0, 100.0, 100.0],
                time: t as f64,
            })
            .collect();
        let msd = mean_square_displacement(&frames, 0);
        for (t, &val) in msd.iter().enumerate() {
            let expected = (d * t as f64).powi(2);
            assert!(
                (val - expected).abs() < 1e-12,
                "MSD[{t}] = {val}, expected {expected}"
            );
        }
    }

    #[test]
    fn test_end_to_end_distance() {
        let pos: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]];
        let d = end_to_end_distance(&pos, 0, 1);
        assert!((d - 5.0).abs() < 1e-12, "expected 5.0, got {d}");
    }

    #[test]
    fn test_contact_map_self_contact() {
        let pos = ring_positions(4, 1.0);
        let map = contact_map(&pos, 0.1);
        // Each atom is at distance 0 from itself → always in contact
        for (i, row) in map.iter().enumerate() {
            assert!(row[i], "atom {i} should be self-contact");
        }
    }

    #[test]
    fn test_contact_map_cutoff() {
        // Two atoms 2 Å apart
        let pos: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        assert!(
            contact_map(&pos, 2.5)[0][1],
            "should be in contact at cutoff 2.5"
        );
        assert!(
            !contact_map(&pos, 1.5)[0][1],
            "should not be in contact at cutoff 1.5"
        );
    }

    #[test]
    fn test_principal_axes_symmetric() {
        // 6 atoms on ±axis with equal mass → known degenerate inertia tensor
        let r = 1.0f64;
        let positions: Vec<[f64; 3]> = vec![
            [r, 0.0, 0.0],
            [-r, 0.0, 0.0],
            [0.0, r, 0.0],
            [0.0, -r, 0.0],
            [0.0, 0.0, r],
            [0.0, 0.0, -r],
        ];
        let masses = vec![1.0f64; 6];
        let (_, evals) = principal_axes(&positions, &masses);
        // All eigenvalues should equal 4.0 (sum of m*r^2 for 4 off-axis atoms × 2 axes)
        for &ev in &evals {
            assert!(
                (ev - 4.0).abs() < 1e-10,
                "expected eigenvalue 4.0, got {ev}"
            );
        }
    }

    // --- Frame interpolation tests ---

    #[test]
    fn test_interpolate_at_endpoints() {
        let fa = TrajectoryFrame {
            positions: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            velocities: None,
            box_lengths: [10.0, 10.0, 10.0],
            time: 0.0,
        };
        let fb = TrajectoryFrame {
            positions: vec![[2.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
            velocities: None,
            box_lengths: [12.0, 12.0, 12.0],
            time: 1.0,
        };

        let f0 = interpolate_frames(&fa, &fb, 0.0);
        assert!((f0.positions[0][0]).abs() < 1e-12);
        assert!((f0.positions[1][0] - 1.0).abs() < 1e-12);

        let f1 = interpolate_frames(&fa, &fb, 1.0);
        assert!((f1.positions[0][0] - 2.0).abs() < 1e-12);
        assert!((f1.positions[1][0] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_interpolate_midpoint() {
        let fa = TrajectoryFrame {
            positions: vec![[0.0, 0.0, 0.0]],
            velocities: None,
            box_lengths: [10.0, 10.0, 10.0],
            time: 0.0,
        };
        let fb = TrajectoryFrame {
            positions: vec![[4.0, 6.0, 8.0]],
            velocities: None,
            box_lengths: [20.0, 20.0, 20.0],
            time: 2.0,
        };

        let fmid = interpolate_frames(&fa, &fb, 0.5);
        assert!((fmid.positions[0][0] - 2.0).abs() < 1e-12);
        assert!((fmid.positions[0][1] - 3.0).abs() < 1e-12);
        assert!((fmid.positions[0][2] - 4.0).abs() < 1e-12);
        assert!((fmid.time - 1.0).abs() < 1e-12);
        assert!((fmid.box_lengths[0] - 15.0).abs() < 1e-12);
    }

    // --- Alignment tests ---

    #[test]
    fn test_center_positions() {
        let mut pos = vec![[1.0, 2.0, 3.0], [3.0, 4.0, 5.0]];
        let centroid = center_positions(&mut pos);
        assert!((centroid[0] - 2.0).abs() < 1e-12);
        assert!((centroid[1] - 3.0).abs() < 1e-12);
        assert!((centroid[2] - 4.0).abs() < 1e-12);

        // After centering, centroid of positions should be origin
        let sum_x: f64 = pos.iter().map(|p| p[0]).sum();
        assert!(
            sum_x.abs() < 1e-12,
            "centroid x should be 0 after centering"
        );
    }

    #[test]
    fn test_align_by_centroid_identical() {
        let reference = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut mobile = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let r = align_by_centroid(&reference, &mut mobile);
        assert!(r < 1e-12, "RMSD of identical after alignment should be 0");
    }

    #[test]
    fn test_align_by_centroid_translated() {
        let reference = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut mobile = vec![[10.0, 10.0, 10.0], [12.0, 10.0, 10.0]];
        let r = align_by_centroid(&reference, &mut mobile);
        assert!(
            r < 1e-10,
            "RMSD after centroid alignment of translation should be ~0, got {r}"
        );
    }

    // --- PBC unwrapping tests ---

    #[test]
    fn test_unwrap_pbc_no_jump() {
        let mut frames = vec![
            TrajectoryFrame {
                positions: vec![[1.0, 0.0, 0.0]],
                velocities: None,
                box_lengths: [10.0, 10.0, 10.0],
                time: 0.0,
            },
            TrajectoryFrame {
                positions: vec![[1.5, 0.0, 0.0]],
                velocities: None,
                box_lengths: [10.0, 10.0, 10.0],
                time: 1.0,
            },
        ];
        unwrap_pbc(&mut frames);
        assert!((frames[1].positions[0][0] - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_unwrap_pbc_with_jump() {
        // Atom at x=9.5 jumps to x=0.5 (box=10), should unwrap to 10.5
        let mut frames = vec![
            TrajectoryFrame {
                positions: vec![[9.5, 0.0, 0.0]],
                velocities: None,
                box_lengths: [10.0, 10.0, 10.0],
                time: 0.0,
            },
            TrajectoryFrame {
                positions: vec![[0.5, 0.0, 0.0]],
                velocities: None,
                box_lengths: [10.0, 10.0, 10.0],
                time: 1.0,
            },
        ];
        unwrap_pbc(&mut frames);
        assert!(
            (frames[1].positions[0][0] - 10.5).abs() < 1e-12,
            "Unwrapped x should be 10.5, got {}",
            frames[1].positions[0][0]
        );
    }

    #[test]
    fn test_unwrap_pbc_backward_jump() {
        // Atom at x=0.5 jumps to x=9.5 (backward crossing)
        let mut frames = vec![
            TrajectoryFrame {
                positions: vec![[0.5, 0.0, 0.0]],
                velocities: None,
                box_lengths: [10.0, 10.0, 10.0],
                time: 0.0,
            },
            TrajectoryFrame {
                positions: vec![[9.5, 0.0, 0.0]],
                velocities: None,
                box_lengths: [10.0, 10.0, 10.0],
                time: 1.0,
            },
        ];
        unwrap_pbc(&mut frames);
        assert!(
            (frames[1].positions[0][0] - (-0.5)).abs() < 1e-12,
            "Backward unwrap: expected -0.5, got {}",
            frames[1].positions[0][0]
        );
    }

    // --- Bond angle tests ---

    #[test]
    fn test_bond_angle_90_degrees() {
        let pos = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let angle = bond_angle(&pos, 0, 1, 2);
        assert!(
            (angle - std::f64::consts::FRAC_PI_2).abs() < 1e-10,
            "Bond angle should be pi/2, got {angle}"
        );
    }

    #[test]
    fn test_bond_angle_180_degrees() {
        let pos = vec![[-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let angle = bond_angle(&pos, 0, 1, 2);
        assert!(
            (angle - std::f64::consts::PI).abs() < 1e-10,
            "Bond angle should be pi, got {angle}"
        );
    }

    // --- Dihedral angle tests ---

    #[test]
    fn test_dihedral_angle_cis() {
        // cis configuration: dihedral ≈ 0
        let pos = vec![
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
        ];
        let phi = dihedral_angle(&pos, 0, 1, 2, 3);
        assert!(phi.abs() < 0.1, "Cis dihedral should be near 0, got {phi}");
    }

    #[test]
    fn test_dihedral_angle_trans() {
        // trans configuration: dihedral ≈ ±π
        let pos = vec![
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, -1.0, 0.0],
        ];
        let phi = dihedral_angle(&pos, 0, 1, 2, 3);
        assert!(
            (phi.abs() - std::f64::consts::PI).abs() < 0.1,
            "Trans dihedral should be near pi, got {phi}"
        );
    }

    // --- Autocorrelation tests ---

    #[test]
    fn test_autocorrelation_constant() {
        let data = vec![5.0; 100];
        let acf = autocorrelation(&data, 10);
        assert_eq!(acf.len(), 11);
        // For constant data, ACF should be 1 everywhere
        for &c in &acf {
            assert!(
                (c - 1.0).abs() < 1e-10,
                "ACF of constant should be 1, got {c}"
            );
        }
    }

    #[test]
    fn test_autocorrelation_zero_lag() {
        let data: Vec<f64> = (0..50).map(|i| (i as f64 * 0.1).sin()).collect();
        let acf = autocorrelation(&data, 20);
        assert!(
            (acf[0] - 1.0).abs() < 1e-10,
            "ACF(0) should be 1, got {}",
            acf[0]
        );
    }

    #[test]
    fn test_autocorrelation_decays() {
        // Noisy data: ACF should decay from 1
        let data: Vec<f64> = (0..200)
            .map(|i| {
                let t = i as f64 * 0.05;
                t.sin() + ((t * 7.3 + 1.0).sin() * 0.5)
            })
            .collect();
        let acf = autocorrelation(&data, 50);
        // ACF at large lag should be less than at lag 0
        assert!(
            acf[50].abs() < acf[0],
            "ACF should decay: ACF(0)={}, ACF(50)={}",
            acf[0],
            acf[50]
        );
    }

    // --- Mass-weighted RMSD test ---

    #[test]
    fn test_rmsd_two_atoms_known_displacement() {
        let ref_pos = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let cur_pos = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let r = rmsd(&ref_pos, &cur_pos);
        // Both displaced by 1.0 → RMSD = 1.0
        assert!((r - 1.0).abs() < 1e-12, "expected RMSD=1, got {r}");
    }

    #[test]
    fn test_radius_of_gyration_single_atom() {
        let pos = vec![[5.0, 3.0, 1.0]];
        let masses = vec![2.0];
        let rg = radius_of_gyration(&pos, &masses);
        assert!(rg.abs() < 1e-12, "Rg of single atom should be 0, got {rg}");
    }

    // --- RMSD clustering tests ---

    fn make_frame(positions: Vec<[f64; 3]>, t: f64) -> TrajectoryFrame {
        TrajectoryFrame {
            positions,
            velocities: None,
            box_lengths: [10.0, 10.0, 10.0],
            time: t,
        }
    }

    #[test]
    fn test_rmsd_cluster_identical_frames_one_cluster() {
        let pos = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let frames = vec![
            make_frame(pos.clone(), 0.0),
            make_frame(pos.clone(), 1.0),
            make_frame(pos.clone(), 2.0),
        ];
        let (assignments, reps) = rmsd_cluster(&frames, 0.5);
        assert_eq!(reps.len(), 1, "Identical frames should form 1 cluster");
        assert!(
            assignments.iter().all(|&a| a == 0),
            "All should be in cluster 0"
        );
    }

    #[test]
    fn test_rmsd_cluster_separate_frames_multiple_clusters() {
        let frames = vec![
            make_frame(vec![[0.0, 0.0, 0.0]], 0.0),
            make_frame(vec![[10.0, 0.0, 0.0]], 1.0), // far away
            make_frame(vec![[0.1, 0.0, 0.0]], 2.0),  // near first
        ];
        let (assignments, reps) = rmsd_cluster(&frames, 0.5);
        assert_eq!(reps.len(), 2, "Should have 2 clusters");
        assert_eq!(assignments[0], 0);
        assert_eq!(assignments[1], 1);
        assert_eq!(assignments[2], 0);
    }

    #[test]
    fn test_cluster_populations_correct() {
        let assignments = vec![0, 1, 0, 0, 1, 2];
        let pop = cluster_populations(&assignments, 3);
        assert_eq!(pop, vec![3, 2, 1]);
    }

    // --- PCA tests ---

    #[test]
    fn test_trajectory_covariance_single_frame_gives_zero() {
        let frames = vec![make_frame(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], 0.0)];
        let cov = trajectory_covariance(&frames);
        // Single frame: all deviations zero → covariance should be zero
        for row in &cov {
            for &v in row {
                assert!(
                    v.abs() < 1e-12,
                    "Single-frame covariance should be 0, got {v}"
                );
            }
        }
    }

    #[test]
    fn test_trajectory_covariance_size() {
        let frames = vec![
            make_frame(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], 0.0),
            make_frame(vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], 1.0),
        ];
        let cov = trajectory_covariance(&frames);
        let n_atoms = 2;
        let dim = 3 * n_atoms;
        assert_eq!(cov.len(), dim);
        assert_eq!(cov[0].len(), dim);
    }

    #[test]
    fn test_frame_to_vector_length() {
        let frame = make_frame(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], 0.0);
        let v = frame_to_vector(&frame);
        assert_eq!(v.len(), 6);
        assert!((v[0] - 1.0).abs() < 1e-12);
        assert!((v[3] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_project_onto_pcs_length() {
        let frames = vec![
            make_frame(vec![[0.0, 0.0, 0.0]], 0.0),
            make_frame(vec![[1.0, 0.0, 0.0]], 1.0),
        ];
        let pc = vec![vec![1.0, 0.0, 0.0]]; // one PC along x
        let proj = project_onto_pcs(&frames, &pc);
        assert_eq!(proj.len(), 2);
        assert_eq!(proj[0].len(), 1);
    }

    // --- DCCM tests ---

    #[test]
    fn test_dccm_diagonal_is_one() {
        let frames: Vec<TrajectoryFrame> = (0..10)
            .map(|t| make_frame(vec![[t as f64, 0.0, 0.0], [0.0, t as f64, 0.0]], t as f64))
            .collect();
        let dccm = dynamic_cross_correlation(&frames);
        assert_eq!(dccm.len(), 2);
        assert!((dccm[0][0] - 1.0).abs() < 1e-6, "DCCM diagonal should be 1");
        assert!((dccm[1][1] - 1.0).abs() < 1e-6, "DCCM diagonal should be 1");
    }

    #[test]
    fn test_dccm_range() {
        let frames: Vec<TrajectoryFrame> = (0..10)
            .map(|t| {
                let x = t as f64 * 0.1;
                make_frame(vec![[x, 0.0, 0.0], [-x, 0.0, 0.0]], x)
            })
            .collect();
        let dccm = dynamic_cross_correlation(&frames);
        for (i, row) in dccm.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (-1.0 - 1e-10..=1.0 + 1e-10).contains(&val),
                    "DCCM[{i}][{j}] = {val} out of range"
                );
            }
        }
    }

    #[test]
    fn test_dccm_anticorrelated_atoms() {
        // Two atoms moving in opposite directions → DCCM[0][1] ≈ -1
        let frames: Vec<TrajectoryFrame> = (0..20)
            .map(|t| {
                let x = t as f64 * 0.1;
                make_frame(vec![[x, 0.0, 0.0], [-x, 0.0, 0.0]], x)
            })
            .collect();
        let dccm = dynamic_cross_correlation(&frames);
        assert!(
            dccm[0][1] < -0.9,
            "Anticorrelated atoms should have DCCM ≈ -1, got {}",
            dccm[0][1]
        );
    }

    // --- Contact frequency map tests ---

    #[test]
    fn test_contact_frequency_self_is_always_one() {
        let frames: Vec<TrajectoryFrame> = (0..5)
            .map(|t| make_frame(vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]], t as f64))
            .collect();
        let cf = contact_frequency_map(&frames, 1.0);
        assert!(
            (cf[0][0] - 1.0).abs() < 1e-12,
            "Self-contact should always be 1"
        );
        assert!(
            (cf[1][1] - 1.0).abs() < 1e-12,
            "Self-contact should always be 1"
        );
    }

    #[test]
    fn test_contact_frequency_distant_always_zero() {
        let frames: Vec<TrajectoryFrame> = (0..5)
            .map(|t| make_frame(vec![[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]], t as f64))
            .collect();
        let cf = contact_frequency_map(&frames, 1.0);
        assert!(
            cf[0][1] < 1e-12,
            "Very distant atoms should have contact frequency 0"
        );
    }

    #[test]
    fn test_contact_frequency_intermittent() {
        // Frame 0: close; frames 1-4: far
        let frames: Vec<TrajectoryFrame> = (0..5)
            .map(|t| {
                let x = if t == 0 { 0.5 } else { 100.0 };
                make_frame(vec![[0.0, 0.0, 0.0], [x, 0.0, 0.0]], t as f64)
            })
            .collect();
        let cf = contact_frequency_map(&frames, 1.0);
        // Contact in 1 out of 5 frames → frequency = 0.2
        assert!(
            (cf[0][1] - 0.2).abs() < 1e-10,
            "Contact frequency should be 0.2, got {}",
            cf[0][1]
        );
    }

    // --- Representative structure tests ---

    #[test]
    fn test_representative_structure_single_frame() {
        let frames = vec![make_frame(vec![[1.0, 0.0, 0.0]], 0.0)];
        let idx = representative_structure(&frames);
        assert_eq!(idx, Some(0));
    }

    #[test]
    fn test_representative_structure_middle_frame() {
        // Mean is (1,0,0); frame at (1,0,0) is closest
        let frames = vec![
            make_frame(vec![[0.0, 0.0, 0.0]], 0.0),
            make_frame(vec![[1.0, 0.0, 0.0]], 1.0),
            make_frame(vec![[2.0, 0.0, 0.0]], 2.0),
        ];
        // Mean = (1,0,0) → closest is frame 1
        let idx = representative_structure(&frames);
        assert_eq!(idx, Some(1));
    }

    // --- RMSF tests ---

    #[test]
    fn test_rmsf_constant_trajectory() {
        let frames: Vec<TrajectoryFrame> = (0..10)
            .map(|t| make_frame(vec![[1.0, 2.0, 3.0]], t as f64))
            .collect();
        let f = rmsf(&frames);
        assert_eq!(f.len(), 1);
        assert!(
            f[0].abs() < 1e-12,
            "RMSF of constant position should be 0, got {}",
            f[0]
        );
    }

    #[test]
    fn test_rmsf_oscillating_atom() {
        // Atom oscillates ±A along x → RMSF = A
        let amp = 2.0;
        let frames: Vec<TrajectoryFrame> = (0..100)
            .map(|t| {
                let x = amp * ((t as f64 * std::f64::consts::PI / 50.0).cos());
                make_frame(vec![[x, 0.0, 0.0]], t as f64)
            })
            .collect();
        let f = rmsf(&frames);
        // RMSF for cos is amplitude / sqrt(2) ≈ 1.414 for amp=2
        assert!(f[0] > 0.5, "RMSF of oscillating atom should be > 0");
        assert!(f[0] < amp * 1.1, "RMSF should be less than amplitude");
    }

    #[test]
    fn test_rmsf_length_matches_n_atoms() {
        let frames: Vec<TrajectoryFrame> = (0..5)
            .map(|t| {
                make_frame(
                    vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
                    t as f64,
                )
            })
            .collect();
        let f = rmsf(&frames);
        assert_eq!(f.len(), 3);
    }
}
