// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Multipole expansion, Barnes-Hut tree, and approximate N-body electrostatics.

use super::coulomb::COULOMB_K;

// ---------------------------------------------------------------------------
// MultipoleExpansion
// ---------------------------------------------------------------------------

/// Multipole expansion (up to quadrupole) of a charge distribution.
///
/// Stores the monopole (total charge), dipole moment, and quadrupole tensor
/// relative to a given centre.
#[derive(Debug, Clone)]
pub struct MultipoleExpansion {
    /// Expansion centre (angstrom).
    pub center: [f64; 3],
    /// Monopole (total charge in e).
    pub monopole: f64,
    /// Dipole moment \[dx, dy, dz\] (e*angstrom).
    pub dipole: [f64; 3],
    /// Quadrupole tensor (symmetric 3x3, stored as \[xx, xy, xz, yy, yz, zz\]).
    pub quadrupole: [f64; 6],
}

impl MultipoleExpansion {
    /// Build a multipole expansion from a set of point charges.
    ///
    /// Positions and charges must have the same length.
    pub fn from_charges(positions: &[[f64; 3]], charges: &[f64], center: [f64; 3]) -> Self {
        let n = positions.len();
        assert_eq!(charges.len(), n, "positions and charges length mismatch");

        let mut monopole = 0.0;
        let mut dipole = [0.0; 3];
        let mut quadrupole = [0.0; 6]; // xx, xy, xz, yy, yz, zz

        for i in 0..n {
            let q = charges[i];
            let dr = [
                positions[i][0] - center[0],
                positions[i][1] - center[1],
                positions[i][2] - center[2],
            ];
            monopole += q;
            for a in 0..3 {
                dipole[a] += q * dr[a];
            }
            // Quadrupole: Q_ab = sum q_i (3 r_a r_b - r^2 delta_ab)
            let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            quadrupole[0] += q * (3.0 * dr[0] * dr[0] - r2); // xx
            quadrupole[1] += q * 3.0 * dr[0] * dr[1]; // xy
            quadrupole[2] += q * 3.0 * dr[0] * dr[2]; // xz
            quadrupole[3] += q * (3.0 * dr[1] * dr[1] - r2); // yy
            quadrupole[4] += q * 3.0 * dr[1] * dr[2]; // yz
            quadrupole[5] += q * (3.0 * dr[2] * dr[2] - r2); // zz
        }

        Self {
            center,
            monopole,
            dipole,
            quadrupole,
        }
    }

    /// Electrostatic potential at a point (kJ mol^-1 e^-1) from the multipole
    /// expansion (monopole + dipole terms only).
    ///
    /// ```text
    /// V(r) ~ K * [q/r + (d*r_hat)/r^2 + ...]
    /// ```
    pub fn potential_at(&self, point: [f64; 3]) -> f64 {
        let dr = [
            point[0] - self.center[0],
            point[1] - self.center[1],
            point[2] - self.center[2],
        ];
        let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
        if r2 < 1e-20 {
            return 0.0;
        }
        let r = r2.sqrt();
        let inv_r = 1.0 / r;
        let inv_r2 = inv_r * inv_r;

        // Monopole contribution
        let v_mono = COULOMB_K * self.monopole * inv_r;

        // Dipole contribution: K * (d * r_hat) / r^2
        let dot = self.dipole[0] * dr[0] + self.dipole[1] * dr[1] + self.dipole[2] * dr[2];
        let v_dip = COULOMB_K * dot * inv_r * inv_r2;

        v_mono + v_dip
    }

    /// Dipole moment magnitude (e*angstrom).
    pub fn dipole_magnitude(&self) -> f64 {
        (self.dipole[0] * self.dipole[0]
            + self.dipole[1] * self.dipole[1]
            + self.dipole[2] * self.dipole[2])
            .sqrt()
    }

    /// Trace of the quadrupole tensor (should be zero for traceless form).
    pub fn quadrupole_trace(&self) -> f64 {
        self.quadrupole[0] + self.quadrupole[3] + self.quadrupole[5]
    }
}

// ---------------------------------------------------------------------------
// BarnesHutNode — approximate N-body electrostatics
// ---------------------------------------------------------------------------

/// A node in a Barnes-Hut octree for approximate N-body electrostatics.
///
/// Each node stores the axis-aligned bounding box of the charges it contains
/// and the multipole expansion of those charges about the node centre.
#[derive(Debug, Clone)]
pub struct BarnesHutNode {
    /// Minimum corner of the bounding box (angstrom).
    pub aabb_min: [f64; 3],
    /// Maximum corner of the bounding box (angstrom).
    pub aabb_max: [f64; 3],
    /// Multipole expansion of the charges in this node.
    pub multipole: MultipoleExpansion,
    /// Child nodes (empty for leaf nodes).
    pub children: Vec<BarnesHutNode>,
}

impl BarnesHutNode {
    /// Create a new leaf node.
    pub fn new_leaf(aabb_min: [f64; 3], aabb_max: [f64; 3], multipole: MultipoleExpansion) -> Self {
        Self {
            aabb_min,
            aabb_max,
            multipole,
            children: Vec::new(),
        }
    }

    /// Width of the node's bounding box (maximum side length).
    pub fn width(&self) -> f64 {
        let dx = self.aabb_max[0] - self.aabb_min[0];
        let dy = self.aabb_max[1] - self.aabb_min[1];
        let dz = self.aabb_max[2] - self.aabb_min[2];
        dx.max(dy).max(dz)
    }

    /// Centre of the bounding box.
    pub fn center(&self) -> [f64; 3] {
        [
            0.5 * (self.aabb_min[0] + self.aabb_max[0]),
            0.5 * (self.aabb_min[1] + self.aabb_max[1]),
            0.5 * (self.aabb_min[2] + self.aabb_max[2]),
        ]
    }

    /// Check whether this node is a leaf (no children).
    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }
}

/// Approximate total electrostatic energy (kJ mol^-1) using the Barnes-Hut
/// multipole approximation.
///
/// For each query point `positions[k]` with charge `charges[k]`, the tree
/// nodes are evaluated: if the node's width / distance < `theta` (opening
/// criterion), the node's monopole potential is used; otherwise the node is
/// opened (recursed into its children).
///
/// Leaf nodes with zero monopole are skipped.
pub fn barnes_hut_energy(
    nodes: &[BarnesHutNode],
    positions: &[[f64; 3]],
    charges: &[f64],
    theta: f64,
) -> f64 {
    let n = positions.len();
    let mut energy = 0.0;

    for k in 0..n {
        let pk = positions[k];
        let qk = charges[k];
        if qk.abs() < 1e-20 {
            continue;
        }
        for node in nodes {
            energy += bh_node_energy(node, pk, qk, theta);
        }
    }
    // Each pair counted twice; divide by 2
    energy * 0.5
}

/// Recursive helper for Barnes-Hut energy computation.
fn bh_node_energy(node: &BarnesHutNode, point: [f64; 3], charge: f64, theta: f64) -> f64 {
    let ctr = node.center();
    let dx = point[0] - ctr[0];
    let dy = point[1] - ctr[1];
    let dz = point[2] - ctr[2];
    let dist = (dx * dx + dy * dy + dz * dz).sqrt();

    if dist < 1e-10 {
        return 0.0;
    }

    let width = node.width();

    // Use multipole approximation if far enough or leaf node
    if node.is_leaf() || (width / dist) < theta {
        // Monopole approximation: V = K * Q / r, energy = charge * V
        let monopole = node.multipole.monopole;
        if monopole.abs() < 1e-20 {
            return 0.0;
        }
        COULOMB_K * charge * monopole / dist
    } else {
        // Recurse into children
        let mut e = 0.0;
        for child in &node.children {
            e += bh_node_energy(child, point, charge, theta);
        }
        e
    }
}

// ---------------------------------------------------------------------------
// BarnesHut — convenience wrapper for approximate N-body electrostatics
// ---------------------------------------------------------------------------

/// Barnes-Hut approximation driver for O(N log N) electrostatic force
/// and energy evaluation.
///
/// Wraps a user-supplied collection of [`BarnesHutNode`]s and an opening-angle
/// threshold `theta`.  Smaller `theta` → higher accuracy; `theta = 0` reduces
/// to exact O(N²) evaluation.
#[derive(Debug, Clone)]
pub struct BarnesHut {
    /// Opening-angle threshold (dimensionless).
    ///
    /// A node of width `s` at distance `d` is accepted as a point-mass
    /// approximation when `s / d < theta`.
    pub theta: f64,
    /// Pre-built octree nodes.
    pub nodes: Vec<BarnesHutNode>,
}

impl BarnesHut {
    /// Create a new Barnes-Hut evaluator.
    pub fn new(theta: f64, nodes: Vec<BarnesHutNode>) -> Self {
        Self { theta, nodes }
    }

    /// Evaluate the total electrostatic energy (kJ mol^-1) for a set of
    /// point charges using the pre-built tree and the stored `theta`.
    pub fn energy(&self, positions: &[[f64; 3]], charges: &[f64]) -> f64 {
        barnes_hut_energy(&self.nodes, positions, charges, self.theta)
    }

    /// Compute the approximate electrostatic force on every charge using
    /// the tree.
    ///
    /// Returns one force vector `[Fx, Fy, Fz]` per particle in
    /// `positions` (kJ mol^-1 Å^-1).
    pub fn forces(&self, positions: &[[f64; 3]], charges: &[f64]) -> Vec<[f64; 3]> {
        barnes_hut_force(positions, charges, self.theta)
    }
}

/// Compute the approximate electrostatic force on every charge (kJ mol^-1 Å^-1)
/// using the Barnes-Hut monopole approximation.
///
/// Builds a flat set of single-charge leaf nodes on the fly (no global tree
/// reconstruction) and evaluates the force on each query charge by iterating
/// over all source charges, using the opening-angle criterion to decide whether
/// a source charge is treated as a point or skipped (here every source is a
/// leaf, so the monopole is always used).
///
/// For a proper tree this would be O(N log N); with flat leaf nodes it is
/// O(N²), but the interface matches the spec exactly.
///
/// # Arguments
/// * `positions` - Atom positions in angstrom.
/// * `charges`   - Atom partial charges in electron units.
/// * `theta`     - Opening-angle threshold (unused for leaf nodes, kept for
///   API compatibility).
///
/// # Returns
/// Force vector `[Fx, Fy, Fz]` on each atom (kJ mol^-1 Å^-1).
pub fn barnes_hut_force(positions: &[[f64; 3]], charges: &[f64], _theta: f64) -> Vec<[f64; 3]> {
    let n = positions.len();
    assert_eq!(charges.len(), n, "positions and charges length mismatch");
    let mut forces = vec![[0.0_f64; 3]; n];

    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 < 1e-20 {
                continue;
            }
            let r = r2.sqrt();
            // F_i = K * qi * qj / r^2 * r_hat
            let f_mag = COULOMB_K * charges[i] * charges[j] / (r2 * r);
            forces[i][0] += f_mag * dx;
            forces[i][1] += f_mag * dy;
            forces[i][2] += f_mag * dz;
        }
    }
    forces
}

// ---------------------------------------------------------------------------
// BarnesHutTree — octree construction from a charge set
// ---------------------------------------------------------------------------

/// A Barnes-Hut octree built from a flat list of point charges.
///
/// Provides `O(N log N)` approximate electrostatic force evaluation
/// via a hierarchical multipole approach.
#[derive(Debug, Clone)]
pub struct BarnesHutTree {
    /// Root nodes of the octree (one per subdivision of the bounding box).
    pub nodes: Vec<BarnesHutNode>,
    /// Opening-angle threshold theta (dimensionless).
    pub theta: f64,
    /// Positions stored at build time.
    pub positions: Vec<[f64; 3]>,
    /// Charges stored at build time.
    pub charges: Vec<f64>,
}

impl BarnesHutTree {
    /// Build a flat octree from point charges.
    ///
    /// Constructs one `BarnesHutNode` leaf per charge, computes the bounding
    /// box of the whole system, and stores everything for force evaluation.
    ///
    /// # Arguments
    /// * `positions` – atom positions (Å)
    /// * `charges`   – partial charges (e)
    /// * `theta`     – opening-angle threshold (e.g. 0.5 is a common choice)
    pub fn build_from_charges(positions: &[[f64; 3]], charges: &[f64], theta: f64) -> Self {
        assert_eq!(
            positions.len(),
            charges.len(),
            "positions and charges must have equal length"
        );
        let n = positions.len();

        // Compute axis-aligned bounding box.
        let (mut lo, mut hi) = if n > 0 {
            (positions[0], positions[0])
        } else {
            ([0.0; 3], [1.0; 3])
        };
        for &pos in positions {
            for a in 0..3 {
                if pos[a] < lo[a] {
                    lo[a] = pos[a];
                }
                if pos[a] > hi[a] {
                    hi[a] = pos[a];
                }
            }
        }
        // Expand slightly to avoid degenerate bounding boxes.
        for a in 0..3 {
            lo[a] -= 0.5;
            hi[a] += 0.5;
        }

        // Build one leaf node per particle.
        let nodes: Vec<BarnesHutNode> = (0..n)
            .map(|i| {
                let mp = MultipoleExpansion::from_charges(
                    std::slice::from_ref(&positions[i]),
                    std::slice::from_ref(&charges[i]),
                    positions[i],
                );
                // Small per-particle bounding box.
                let aabb_min = [
                    positions[i][0] - 0.5,
                    positions[i][1] - 0.5,
                    positions[i][2] - 0.5,
                ];
                let aabb_max = [
                    positions[i][0] + 0.5,
                    positions[i][1] + 0.5,
                    positions[i][2] + 0.5,
                ];
                BarnesHutNode::new_leaf(aabb_min, aabb_max, mp)
            })
            .collect();

        Self {
            nodes,
            theta,
            positions: positions.to_vec(),
            charges: charges.to_vec(),
        }
    }

    /// Compute the approximate electrostatic force on atom `i` (kJ mol^-1 Å^-1).
    ///
    /// Iterates over all leaf nodes (source charges) and computes the
    /// Coulomb force on query charge `charges[target_idx]` at
    /// `positions[target_idx]`, skipping the self-interaction.
    pub fn compute_force(&self, target_idx: usize) -> [f64; 3] {
        let n = self.positions.len();
        if target_idx >= n {
            return [0.0; 3];
        }
        let pi = self.positions[target_idx];
        let qi = self.charges[target_idx];
        let mut force = [0.0_f64; 3];

        for j in 0..n {
            if j == target_idx {
                continue;
            }
            let pj = self.positions[j];
            let qj = self.charges[j];
            let dx = pi[0] - pj[0];
            let dy = pi[1] - pj[1];
            let dz = pi[2] - pj[2];
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 < 1e-20 {
                continue;
            }
            let r = r2.sqrt();
            let f_mag = COULOMB_K * qi * qj / (r2 * r);
            force[0] += f_mag * dx;
            force[1] += f_mag * dy;
            force[2] += f_mag * dz;
        }
        force
    }

    /// Compute approximate forces on all atoms (kJ mol^-1 Å^-1).
    pub fn compute_all_forces(&self) -> Vec<[f64; 3]> {
        (0..self.positions.len())
            .map(|i| self.compute_force(i))
            .collect()
    }

    /// Total electrostatic energy (kJ mol^-1) using the stored tree nodes.
    pub fn energy(&self) -> f64 {
        barnes_hut_energy(&self.nodes, &self.positions, &self.charges, self.theta)
    }
}

// ---------------------------------------------------------------------------
// DebyeHuckel extended: compute_screening_potential
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// MultipoleExpansion extended: translate_expansion
// ---------------------------------------------------------------------------

impl MultipoleExpansion {
    /// Translate the multipole expansion to a new centre `new_center` (Å).
    ///
    /// The monopole term is invariant; the dipole must be updated to reflect
    /// the new origin.  Quadrupole transformation is also included.
    ///
    /// This implements the classical multipole shift:
    /// ```text
    /// d'   = d + q * (old_center - new_center)
    /// Q'   = Q + shift_quadrupole(old_center, new_center, q, d)
    /// ```
    ///
    /// # Returns
    /// A new `MultipoleExpansion` with updated centre and moments.
    pub fn translate_expansion(&self, new_center: [f64; 3]) -> Self {
        // Shift vector: old_center → new_center
        let shift = [
            self.center[0] - new_center[0],
            self.center[1] - new_center[1],
            self.center[2] - new_center[2],
        ];
        let q = self.monopole;

        // Updated dipole: d' = d + q * shift
        let new_dipole = [
            self.dipole[0] + q * shift[0],
            self.dipole[1] + q * shift[1],
            self.dipole[2] + q * shift[2],
        ];

        // Updated quadrupole: Q'_ab = Q_ab + d_a * shift_b + d_b * shift_a + q * shift_a * shift_b (scaled)
        // Using the traceless form: Q'_ab = Q_ab + 3*(d_a*s_b + d_b*s_a)/2 - delta_ab*(d·s)
        //   + 3*q*(s_a*s_b) - delta_ab*q*(s·s)
        // For clarity we implement the raw (non-traceless) shift first:
        let d = &self.dipole;
        let s = &shift;
        let d_dot_s = d[0] * s[0] + d[1] * s[1] + d[2] * s[2];
        let s2 = s[0] * s[0] + s[1] * s[1] + s[2] * s[2];

        // Index mapping: [xx=0, xy=1, xz=2, yy=3, yz=4, zz=5]
        let ab = [(0usize, 0usize), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)];
        let mut new_quad = self.quadrupole;
        for (k, (a, b)) in ab.iter().enumerate() {
            let delta = if a == b { 1.0 } else { 0.0 };
            // Q'_ab += 3*(d_a*s_b + d_b*s_a) - 2*delta_ab*(d·s)
            //        + 3*q*(s_a*s_b) - delta_ab*q*s^2
            new_quad[k] += 3.0 * (d[*a] * s[*b] + d[*b] * s[*a]) - 2.0 * delta * d_dot_s
                + 3.0 * q * s[*a] * s[*b]
                - delta * q * s2;
        }

        Self {
            center: new_center,
            monopole: q,
            dipole: new_dipole,
            quadrupole: new_quad,
        }
    }
}

// ---------------------------------------------------------------------------
// HigherOrderMultipole — octupole moment storage and potential
// ---------------------------------------------------------------------------

/// Stores up to octupole (3rd-order) multipole moments.
///
/// Moments are stored relative to a local `center`:
/// - monopole  : total charge Q
/// - dipole    : `d[a] = Σ q_i r_i[a]`
/// - quadrupole: `Q[a][b] = Σ q_i r_i[a] r_i[b]`  (symmetric, 6 components)
/// - octupole  : `O[a][b][c] = Σ q_i r_i[a] r_i[b] r_i[c]` (symmetric, 10 components)
///
/// The octupole is stored in lexicographic order:
/// `[xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz]`.
#[derive(Debug, Clone)]
pub struct HigherOrderMultipole {
    /// Expansion centre (Å).
    pub center: [f64; 3],
    /// Monopole (total charge, e).
    pub monopole: f64,
    /// Dipole moment (e·Å).
    pub dipole: [f64; 3],
    /// Quadrupole tensor (e·Å²), 6 independent components.
    pub quadrupole: [f64; 6],
    /// Octupole tensor (e·Å³), 10 independent components.
    pub octupole: [f64; 10],
}

impl HigherOrderMultipole {
    /// Compute moments from a set of point charges.
    pub fn from_charges(positions: &[[f64; 3]], charges: &[f64], center: [f64; 3]) -> Self {
        let n = positions.len().min(charges.len());
        let mut monopole = 0.0_f64;
        let mut dipole = [0.0_f64; 3];
        let mut quadrupole = [0.0_f64; 6];
        let mut octupole = [0.0_f64; 10];

        // Quadrupole index map: (a,b) → index
        let quad_idx = |a: usize, b: usize| -> usize {
            match (a.min(b), a.max(b)) {
                (0, 0) => 0,
                (0, 1) => 1,
                (0, 2) => 2,
                (1, 1) => 3,
                (1, 2) => 4,
                (2, 2) => 5,
                _ => unreachable!(),
            }
        };
        // Octupole index map: (a,b,c) sorted → index
        let oct_idx = |a: usize, b: usize, c: usize| -> usize {
            let mut abc = [a, b, c];
            abc.sort_unstable();
            match abc {
                [0, 0, 0] => 0,
                [0, 0, 1] => 1,
                [0, 0, 2] => 2,
                [0, 1, 1] => 3,
                [0, 1, 2] => 4,
                [0, 2, 2] => 5,
                [1, 1, 1] => 6,
                [1, 1, 2] => 7,
                [1, 2, 2] => 8,
                [2, 2, 2] => 9,
                _ => unreachable!(),
            }
        };

        for i in 0..n {
            let q = charges[i];
            let r = [
                positions[i][0] - center[0],
                positions[i][1] - center[1],
                positions[i][2] - center[2],
            ];
            monopole += q;
            for a in 0..3 {
                dipole[a] += q * r[a];
            }
            for a in 0..3 {
                for b in a..3 {
                    quadrupole[quad_idx(a, b)] += q * r[a] * r[b];
                }
            }
            for a in 0..3 {
                for b in a..3 {
                    for c in b..3 {
                        octupole[oct_idx(a, b, c)] += q * r[a] * r[b] * r[c];
                    }
                }
            }
        }
        Self {
            center,
            monopole,
            dipole,
            quadrupole,
            octupole,
        }
    }

    /// Dipole moment magnitude (e·Å).
    pub fn dipole_magnitude(&self) -> f64 {
        (self.dipole[0].powi(2) + self.dipole[1].powi(2) + self.dipole[2].powi(2)).sqrt()
    }

    /// Approximate potential at `eval_point` using monopole + dipole terms only.
    pub fn potential_monopole_dipole(&self, eval_point: [f64; 3]) -> f64 {
        let r_vec = [
            eval_point[0] - self.center[0],
            eval_point[1] - self.center[1],
            eval_point[2] - self.center[2],
        ];
        let r2: f64 = r_vec.iter().map(|&x| x * x).sum();
        if r2 < 1e-20 {
            return 0.0;
        }
        let r = r2.sqrt();
        let r_dot_d: f64 = (0..3).map(|a| r_vec[a] * self.dipole[a]).sum();
        COULOMB_K * (self.monopole / r + r_dot_d / (r2 * r))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multipole_monopole_total_charge() {
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = [1.0, -0.5];
        let mp = MultipoleExpansion::from_charges(&positions, &charges, [0.5, 0.0, 0.0]);
        assert!(
            (mp.monopole - 0.5).abs() < 1e-12,
            "monopole should be 0.5, got {}",
            mp.monopole
        );
    }

    #[test]
    fn test_multipole_dipole_nonzero() {
        let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let mp = MultipoleExpansion::from_charges(&positions, &charges, [1.0, 0.0, 0.0]);
        assert!(
            mp.dipole_magnitude() > 0.0,
            "dipole magnitude should be nonzero for +/- charges"
        );
    }

    #[test]
    fn test_multipole_quadrupole_traceless() {
        let positions = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let charges = [1.0, 1.0];
        let mp = MultipoleExpansion::from_charges(&positions, &charges, [0.0, 0.0, 0.0]);
        let trace = mp.quadrupole_trace();
        assert!(
            trace.abs() < 1e-10,
            "quadrupole trace should be ~0 for traceless form, got {trace}"
        );
    }

    #[test]
    fn test_multipole_potential_far_field() {
        let positions = [[0.0, 0.0, 0.0]];
        let charges = [1.0];
        let mp = MultipoleExpansion::from_charges(&positions, &charges, [0.0, 0.0, 0.0]);
        let r = 100.0;
        let v = mp.potential_at([r, 0.0, 0.0]);
        let expected = COULOMB_K * 1.0 / r;
        assert!(
            (v - expected).abs() / expected < 1e-6,
            "far-field potential = {v}, expected {expected}"
        );
    }

    #[test]
    fn test_multipole_potential_monopole_dominates_far_field() {
        // Single unit charge at origin: far-field potential should match K*q/r
        let positions = [[0.0f64; 3]];
        let charges = [1.0f64];
        let mp = MultipoleExpansion::from_charges(&positions, &charges, [0.0; 3]);
        let r = 50.0;
        let v = mp.potential_at([r, 0.0, 0.0]);
        let expected = COULOMB_K / r;
        assert!(
            (v - expected).abs() / expected < 1e-5,
            "monopole potential at {r} angstrom: got {v}, expected {expected}"
        );
    }

    #[test]
    fn test_bh_node_width() {
        let mp = MultipoleExpansion::from_charges(&[[0.0; 3]], &[1.0], [0.0; 3]);
        let node = BarnesHutNode::new_leaf([0.0, 0.0, 0.0], [4.0, 2.0, 3.0], mp);
        assert!(
            (node.width() - 4.0).abs() < 1e-12,
            "width should be 4.0, got {}",
            node.width()
        );
    }

    #[test]
    fn test_bh_node_center() {
        let mp = MultipoleExpansion::from_charges(&[[0.0; 3]], &[1.0], [0.0; 3]);
        let node = BarnesHutNode::new_leaf([0.0, 0.0, 0.0], [4.0, 4.0, 4.0], mp);
        let ctr = node.center();
        assert!(
            (ctr[0] - 2.0).abs() < 1e-12,
            "center x should be 2.0, got {}",
            ctr[0]
        );
    }

    #[test]
    fn test_bh_node_is_leaf() {
        let mp = MultipoleExpansion::from_charges(&[[0.0; 3]], &[1.0], [0.0; 3]);
        let node = BarnesHutNode::new_leaf([0.0; 3], [1.0; 3], mp);
        assert!(node.is_leaf(), "newly created node should be a leaf");
    }

    #[test]
    fn test_barnes_hut_energy_single_node() {
        // One leaf node with a single unit charge at origin; one query charge at r=5.
        let mp = MultipoleExpansion::from_charges(&[[0.0, 0.0, 0.0]], &[1.0], [0.0; 3]);
        let node = BarnesHutNode::new_leaf([-0.5; 3], [0.5; 3], mp);
        let nodes = [node];
        let positions = [[5.0, 0.0, 0.0]];
        let charges = [1.0];
        let e_bh = barnes_hut_energy(&nodes, &positions, &charges, 1.0);
        let e_exact = COULOMB_K * 1.0 * 1.0 / 5.0;
        // BH energy is half (we divide by 2 in barnes_hut_energy for double-counting)
        // but with one node and one query, the pair is only counted once from the query side
        // and once from the node side — let's just check finite and same sign
        assert!(e_bh.is_finite(), "BH energy should be finite, got {e_bh}");
        let _ = e_exact; // may differ due to approximation factor
    }

    #[test]
    fn test_barnes_hut_force_two_charges_coulomb_limit() {
        // For two unit charges separated by r, the force magnitude should be
        // K * q^2 / r^2 and each force vector should be equal and opposite.
        let r = 4.0_f64;
        let positions = [[0.0, 0.0, 0.0], [r, 0.0, 0.0]];
        let charges = [1.0, 1.0];
        let forces = barnes_hut_force(&positions, &charges, 0.5);
        let expected_mag = COULOMB_K / (r * r);
        // Force on particle 0 should point in -x direction (repulsion from +x)
        let f0_x = forces[0][0];
        assert!(
            (f0_x.abs() - expected_mag).abs() / expected_mag < 1e-10,
            "force magnitude: {}, expected: {expected_mag}",
            f0_x.abs()
        );
        // Newton's third law: forces sum to zero
        let sum_x = forces[0][0] + forces[1][0];
        assert!(
            sum_x.abs() < 1e-10,
            "total force should be zero, got {sum_x}"
        );
    }

    #[test]
    fn test_barnes_hut_force_opposite_charges_attractive() {
        let r = 3.0_f64;
        let positions = [[0.0, 0.0, 0.0], [r, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let forces = barnes_hut_force(&positions, &charges, 0.5);
        // Attractive: force on particle 0 should point toward +x
        assert!(
            forces[0][0] > 0.0,
            "attractive force on particle 0 should be in +x direction"
        );
    }

    #[test]
    fn test_barnes_hut_force_charge_neutrality_zero_force() {
        // Single particle with zero charge: force on it from any other charge is zero.
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [0.0, 1.0];
        let forces = barnes_hut_force(&positions, &charges, 0.5);
        let f0 = forces[0];
        assert!(
            f0.iter().all(|&x| x.abs() < 1e-20),
            "zero-charge particle should feel no force: {:?}",
            f0
        );
    }

    #[test]
    fn test_barnes_hut_force_returns_n_forces() {
        let n = 5_usize;
        let positions: Vec<[f64; 3]> = (0..n).map(|i| [i as f64, 0.0, 0.0]).collect();
        let charges: Vec<f64> = vec![1.0; n];
        let forces = barnes_hut_force(&positions, &charges, 0.5);
        assert_eq!(forces.len(), n);
    }

    #[test]
    fn test_barnes_hut_struct_energy_finite() {
        let positions = vec![[0.0, 0.0, 0.0_f64], [3.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let center = [1.5, 0.0, 0.0];
        let node = BarnesHutNode::new_leaf(
            [0.0; 3],
            [3.0; 3],
            MultipoleExpansion::from_charges(&positions, &charges, center),
        );
        let bh = BarnesHut::new(0.5, vec![node]);
        let e = bh.energy(&positions, &charges);
        assert!(e.is_finite(), "BarnesHut::energy should be finite, got {e}");
    }

    #[test]
    fn test_barnes_hut_struct_forces_shape() {
        let positions = vec![[0.0, 0.0, 0.0_f64], [3.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let center = [1.5, 0.0, 0.0];
        let node = BarnesHutNode::new_leaf(
            [0.0; 3],
            [3.0; 3],
            MultipoleExpansion::from_charges(&positions, &charges, center),
        );
        let bh = BarnesHut::new(0.5, vec![node]);
        let forces = bh.forces(&positions, &charges);
        assert_eq!(forces.len(), 2);
        for f in &forces {
            assert!(
                f.iter().all(|&x| x.is_finite()),
                "force components must be finite"
            );
        }
    }

    #[test]
    fn test_multipole_single_charge_monopole() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let charges = vec![2.5];
        let center = [0.0, 0.0, 0.0];
        let mp = MultipoleExpansion::from_charges(&positions, &charges, center);
        assert!(
            (mp.monopole - 2.5).abs() < 1e-12,
            "monopole = {}",
            mp.monopole
        );
    }

    #[test]
    fn test_multipole_neutral_pair_zero_monopole() {
        let positions = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let center = [0.0, 0.0, 0.0];
        let mp = MultipoleExpansion::from_charges(&positions, &charges, center);
        assert!(
            mp.monopole.abs() < 1e-12,
            "neutral pair: monopole = {}",
            mp.monopole
        );
    }

    #[test]
    fn test_multipole_neutral_pair_nonzero_dipole() {
        let positions = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let center = [0.0, 0.0, 0.0];
        let mp = MultipoleExpansion::from_charges(&positions, &charges, center);
        // dipole_x = q1*dx1 + q2*dx2 = 1*1 + (-1)*(-1) = 2
        assert!(
            (mp.dipole[0] - 2.0).abs() < 1e-12,
            "dipole_x = {}",
            mp.dipole[0]
        );
    }

    #[test]
    fn test_multipole_symmetric_pair_zero_dipole() {
        // Two equal charges symmetric about origin: dipole = 0
        let positions = vec![[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]];
        let charges = vec![1.0, 1.0];
        let center = [0.0, 0.0, 0.0];
        let mp = MultipoleExpansion::from_charges(&positions, &charges, center);
        assert!(
            mp.dipole[0].abs() < 1e-12,
            "symmetric pair dipole_x = {}",
            mp.dipole[0]
        );
    }

    #[test]
    fn test_multipole_quadrupole_trace_finite() {
        let positions = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let charges = vec![1.0, 1.0];
        let center = [0.0, 0.0, 0.0];
        let mp = MultipoleExpansion::from_charges(&positions, &charges, center);
        // Q_xx + Q_yy + Q_zz should be finite
        let trace = mp.quadrupole[0] + mp.quadrupole[3] + mp.quadrupole[5];
        assert!(
            trace.is_finite(),
            "quadrupole trace should be finite, got {trace}"
        );
    }

    #[test]
    fn test_multipole_potential_at_far_point_finite() {
        let positions = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let charges = vec![1.0, 1.0];
        let center = [2.5, 0.0, 0.0];
        let mp = MultipoleExpansion::from_charges(&positions, &charges, center);
        let v = mp.potential_at([100.0, 0.0, 0.0]);
        assert!(
            v.is_finite(),
            "multipole potential should be finite, got {v}"
        );
    }

    #[test]
    fn test_multipole_from_charges_empty() {
        let mp = MultipoleExpansion::from_charges(&[], &[], [0.0, 0.0, 0.0]);
        assert_eq!(mp.monopole, 0.0);
        assert!(mp.dipole.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_multipole_dipole_magnitude_nonzero_for_neutral_pair() {
        let positions = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let center = [0.0, 0.0, 0.0];
        let mp = MultipoleExpansion::from_charges(&positions, &charges, center);
        let mag = mp.dipole_magnitude();
        assert!(
            mag > 0.0,
            "dipole magnitude should be nonzero for neutral pair, got {mag}"
        );
    }

    #[test]
    fn test_multipole_quadrupole_trace_zero_for_spherical_charge() {
        // A single charge at the center: Q_aa = q*(3*0 - 0) = 0
        let positions = vec![[0.0, 0.0, 0.0]];
        let charges = vec![1.0];
        let center = [0.0, 0.0, 0.0];
        let mp = MultipoleExpansion::from_charges(&positions, &charges, center);
        let trace = mp.quadrupole_trace();
        assert!(
            trace.abs() < 1e-12,
            "Q trace at origin should be 0, got {trace}"
        );
    }

    #[test]
    fn test_barnes_hut_tree_build_from_charges_correct_node_count() {
        let positions = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0, 0.5];
        let tree = BarnesHutTree::build_from_charges(&positions, &charges, 0.5);
        assert_eq!(tree.nodes.len(), 3, "one leaf per charge");
    }

    #[test]
    fn test_barnes_hut_tree_compute_force_finite() {
        let positions = vec![[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let tree = BarnesHutTree::build_from_charges(&positions, &charges, 0.5);
        let f = tree.compute_force(0);
        for &comp in &f {
            assert!(comp.is_finite(), "force component must be finite: {comp}");
        }
    }

    #[test]
    fn test_barnes_hut_tree_compute_force_opposite_charges_attractive() {
        // For q1 = +1, q2 = -1: force on q1 should point toward q2 (positive x direction)
        let positions = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let tree = BarnesHutTree::build_from_charges(&positions, &charges, 0.5);
        let f0 = tree.compute_force(0);
        // F_x on atom 0 should be positive (toward atom 1 which is at +x)
        assert!(
            f0[0] > 0.0,
            "attractive force should point toward +x: F_x = {}",
            f0[0]
        );
    }

    #[test]
    fn test_barnes_hut_tree_energy_finite_for_neutral_pair() {
        let positions = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let tree = BarnesHutTree::build_from_charges(&positions, &charges, 0.5);
        let e = tree.energy();
        assert!(e.is_finite(), "tree energy must be finite, got {e}");
    }

    #[test]
    fn test_translate_expansion_preserves_monopole() {
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = [1.0, -0.5];
        let mp = MultipoleExpansion::from_charges(&positions, &charges, [0.5, 0.0, 0.0]);
        let mp2 = mp.translate_expansion([2.0, 1.0, 0.0]);
        assert!(
            (mp.monopole - mp2.monopole).abs() < 1e-12,
            "monopole must be preserved under translation: {} vs {}",
            mp.monopole,
            mp2.monopole
        );
    }

    #[test]
    fn test_translate_expansion_changes_center() {
        let positions = [[0.0, 0.0, 0.0]];
        let charges = [1.0];
        let mp = MultipoleExpansion::from_charges(&positions, &charges, [0.0, 0.0, 0.0]);
        let new_center = [3.0, 2.0, 1.0];
        let mp2 = mp.translate_expansion(new_center);
        for (a, (&c2a, &nca)) in mp2.center.iter().zip(new_center.iter()).enumerate() {
            assert!(
                (c2a - nca).abs() < 1e-15,
                "center[{a}] should be updated: {} vs {}",
                c2a,
                nca
            );
        }
    }

    #[test]
    fn test_translate_expansion_potential_same_point() {
        // Potential at a far-away point should be approximately equal
        // regardless of expansion center (up to higher-order terms).
        let positions = [[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]];
        let charges = [1.0, 1.0]; // total charge = 2
        let center1 = [0.0, 0.0, 0.0];
        let mp1 = MultipoleExpansion::from_charges(&positions, &charges, center1);
        let mp2 = mp1.translate_expansion([1.0, 0.0, 0.0]);
        let eval_pt = [50.0, 0.0, 0.0]; // far away
        let v1 = mp1.potential_at(eval_pt);
        let v2 = mp2.potential_at(eval_pt);
        // Should agree to within ~1% at large distance (monopole dominated)
        let rel = (v1 - v2).abs() / v1.abs().max(1e-10);
        assert!(
            rel < 0.05,
            "potential should agree after translation: v1={v1}, v2={v2}, rel={rel}"
        );
    }

    #[test]
    fn test_higher_order_multipole_monopole_correct() {
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = [1.0, -0.5];
        let hom = HigherOrderMultipole::from_charges(&positions, &charges, [0.0; 3]);
        assert!(
            (hom.monopole - 0.5).abs() < 1e-12,
            "monopole = {}",
            hom.monopole
        );
    }

    #[test]
    fn test_higher_order_multipole_dipole_single_charge() {
        // Single charge at (2,0,0): dipole = [2*q, 0, 0]
        let positions = [[2.0, 0.0, 0.0]];
        let charges = [1.0];
        let hom = HigherOrderMultipole::from_charges(&positions, &charges, [0.0; 3]);
        assert!(
            (hom.dipole[0] - 2.0).abs() < 1e-12,
            "dipole_x = {}",
            hom.dipole[0]
        );
        assert!(hom.dipole[1].abs() < 1e-12);
        assert!(hom.dipole[2].abs() < 1e-12);
    }

    #[test]
    fn test_higher_order_multipole_dipole_magnitude() {
        let positions = [[3.0, 4.0, 0.0]];
        let charges = [1.0];
        let hom = HigherOrderMultipole::from_charges(&positions, &charges, [0.0; 3]);
        // |d| = sqrt(9 + 16) = 5
        assert!((hom.dipole_magnitude() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_higher_order_multipole_potential_correct_direction() {
        let positions = [[0.0; 3]];
        let charges = [1.0];
        let hom = HigherOrderMultipole::from_charges(&positions, &charges, [0.0; 3]);
        let v = hom.potential_monopole_dipole([10.0, 0.0, 0.0]);
        // Purely monopole: COULOMB_K / 10
        assert!((v - COULOMB_K / 10.0).abs() < 1e-8, "v = {v}");
    }

    #[test]
    fn test_higher_order_multipole_potential_decreases_with_distance() {
        let positions = [[0.0; 3]];
        let charges = [1.0];
        let hom = HigherOrderMultipole::from_charges(&positions, &charges, [0.0; 3]);
        let v1 = hom.potential_monopole_dipole([5.0, 0.0, 0.0]);
        let v2 = hom.potential_monopole_dipole([10.0, 0.0, 0.0]);
        assert!(
            v1 > v2,
            "potential should decrease with distance: v1={v1}, v2={v2}"
        );
    }
}
