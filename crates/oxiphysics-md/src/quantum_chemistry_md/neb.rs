// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Nudged elastic band (NEB) and climbing-image NEB on full 3N-dimensional
//! configurations.

use super::dist;

/// Norm below which a configuration-space vector is treated as zero.
const TINY_NORM: f64 = 1e-30;

/// Nudged elastic band (NEB) method for reaction path optimization.
///
/// Finds the minimum-energy path (MEP) between two stable configurations
/// by optimizing a chain of images connected by harmonic spring forces.
/// The climbing-image variant (CI-NEB) converges to the exact saddle point.
///
/// Every image is a point in the 3N-dimensional configuration space; tangents,
/// projections and spring forces are all evaluated in that space.
#[derive(Debug, Clone)]
pub struct NebMethod {
    /// Chain of images: `images[i]` is the i-th image's atomic coordinates.
    pub images: Vec<Vec<[f64; 3]>>,
    /// Spring constant connecting adjacent images \[eV/Å²\].
    pub spring_k: f64,
    /// Potential energy at each image \[eV\].
    ///
    /// Used for the energy-weighted (Henkelman-Jonsson) tangent and to select
    /// the climbing image.  When all energies are equal (the initial state),
    /// the tangent falls back to the bisection of the two neighbour directions.
    pub energies: Vec<f64>,
}

impl NebMethod {
    /// Create a NEB chain by linearly interpolating between two endpoint images.
    ///
    /// # Arguments
    /// * `start` – Reactant geometry.
    /// * `end` – Product geometry.
    /// * `n_images` – Total number of images including endpoints.
    /// * `spring_k` – Spring constant in eV/Å².
    pub fn new(start: Vec<[f64; 3]>, end: Vec<[f64; 3]>, n_images: usize, spring_k: f64) -> Self {
        assert!(n_images >= 2, "Need at least 2 images");
        let mut images = Vec::with_capacity(n_images);
        for img in 0..n_images {
            let t = img as f64 / (n_images - 1) as f64;
            let interp: Vec<[f64; 3]> = start
                .iter()
                .zip(end.iter())
                .map(|(s, e)| {
                    [
                        s[0] + t * (e[0] - s[0]),
                        s[1] + t * (e[1] - s[1]),
                        s[2] + t * (e[2] - s[2]),
                    ]
                })
                .collect();
            images.push(interp);
        }
        let energies = vec![0.0; n_images];
        Self {
            images,
            spring_k,
            energies,
        }
    }

    /// Compute the NEB force on each interior image from the supplied true
    /// (potential) atomic forces.
    ///
    /// For each interior image `i` the NEB force is
    ///
    /// F_i = F_i^⊥(true) + F_i^∥(spring),
    ///
    /// where the true force is projected onto the hyperplane perpendicular to
    /// the 3N-dimensional unit tangent τ̂_i and the spring force acts along τ̂_i:
    ///
    /// F^⊥ = F_true − (F_true·τ̂)τ̂,   F^∥ = k(|R_{i+1}−R_i| − |R_i−R_{i−1}|) τ̂,
    ///
    /// with the dot product and the distances taken over all 3N coordinates.
    ///
    /// `true_forces[i][a]` is the physical force (−∇V) on atom `a` of image `i`,
    /// which the caller must obtain from the potential-energy surface — this
    /// method evaluates no PES of its own. Endpoint images are fixed and receive
    /// zero NEB force; missing entries in `true_forces` are treated as zero.
    pub fn neb_force(&self, true_forces: &[Vec<[f64; 3]>]) -> Vec<Vec<[f64; 3]>> {
        let n = self.images.len();
        let n_atoms = self.images.first().map(|im| im.len()).unwrap_or(0);
        let mut forces = vec![vec![[0.0f64; 3]; n_atoms]; n];
        if n < 3 {
            return forces;
        }
        for (i, force_row) in forces.iter_mut().enumerate().take(n - 1).skip(1) {
            let tau = self.local_tangent(i);
            let d_next = configuration_distance(&self.images[i + 1], &self.images[i]);
            let d_prev = configuration_distance(&self.images[i], &self.images[i - 1]);
            let f_spring = self.spring_k * (d_next - d_prev);
            let ft = true_force_row(true_forces, i, force_row.len());
            let f_dot_tau = dot_3n(&ft, &tau);
            for ((fa, ft_a), tau_a) in force_row.iter_mut().zip(ft.iter()).zip(tau.iter()) {
                for k in 0..3 {
                    fa[k] = ft_a[k] - f_dot_tau * tau_a[k] + f_spring * tau_a[k];
                }
            }
        }
        forces
    }

    /// Compute the climbing-image NEB force on the highest-energy interior image.
    ///
    /// The climbing image feels the full true force with its tangential
    /// component inverted (and no spring force), so it is driven uphill to the
    /// saddle point:
    ///
    /// F^CI = F_true − 2 (F_true·τ̂) τ̂,
    ///
    /// where F_true·τ̂ is the 3N-dimensional dot product.
    ///
    /// Returns the per-atom climbing force for the highest-energy image.
    /// `true_forces` supplies the physical atomic forces (the caller evaluates
    /// the PES). Returns an empty vector when there are fewer than three images.
    pub fn climbing_image_neb(&self, true_forces: &[Vec<[f64; 3]>]) -> Vec<[f64; 3]> {
        let n = self.images.len();
        if n < 3 {
            return Vec::new();
        }
        let ci = self.highest_energy_image();
        let tau = self.local_tangent(ci);
        let ft = true_force_row(true_forces, ci, self.images[ci].len());
        let f_dot_tau = dot_3n(&ft, &tau);
        ft.iter()
            .zip(tau.iter())
            .map(|(f, t)| {
                [
                    f[0] - 2.0 * f_dot_tau * t[0],
                    f[1] - 2.0 * f_dot_tau * t[1],
                    f[2] - 2.0 * f_dot_tau * t[2],
                ]
            })
            .collect()
    }

    /// Reaction coordinate (arc length) along the NEB path in Angstrom.
    pub fn reaction_coordinate(&self) -> Vec<f64> {
        let n = self.images.len();
        let mut s = vec![0.0f64; n];
        for i in 1..n {
            s[i] = s[i - 1] + configuration_distance(&self.images[i], &self.images[i - 1]);
        }
        s
    }

    /// Activation energy estimated as the maximum energy along the path.
    ///
    /// Returns E_max − E_reactant in eV.
    pub fn activation_energy(&self) -> f64 {
        let e_max = self
            .energies
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let e_min = self.energies.first().cloned().unwrap_or(0.0);
        e_max - e_min
    }

    // -- private helpers --

    /// Energy-weighted unit tangent (Henkelman & Jónsson, J. Chem. Phys. 113,
    /// 9978 (2000)) of interior image `i`, as a 3N-dimensional vector.
    ///
    /// * E_{i+1} > E_i > E_{i−1}: τ = R_{i+1} − R_i.
    /// * E_{i+1} < E_i < E_{i−1}: τ = R_i − R_{i−1}.
    /// * Extremum: τ = ΔE_max τ⁺ + ΔE_min τ⁻ (if E_{i+1} > E_{i−1}) or
    ///   ΔE_min τ⁺ + ΔE_max τ⁻ (otherwise).
    /// * Flat energies (ΔE_max = 0, e.g. before any energy is known): the
    ///   bisection τ⁺/|τ⁺| + τ⁻/|τ⁻|.
    ///
    /// The result is normalised over all 3N components (zero if degenerate).
    fn local_tangent(&self, i: usize) -> Vec<[f64; 3]> {
        let cur = &self.images[i];
        let tau_plus = difference(&self.images[i + 1], cur);
        let tau_minus = difference(cur, &self.images[i - 1]);
        let e_prev = self.energy_at(i - 1);
        let e_cur = self.energy_at(i);
        let e_next = self.energy_at(i + 1);

        let mut tau = if e_next > e_cur && e_cur > e_prev {
            tau_plus
        } else if e_next < e_cur && e_cur < e_prev {
            tau_minus
        } else {
            let de_next = (e_next - e_cur).abs();
            let de_prev = (e_cur - e_prev).abs();
            let de_max = de_next.max(de_prev);
            let de_min = de_next.min(de_prev);
            if de_max > 0.0 {
                let (w_plus, w_minus) = if e_next > e_prev {
                    (de_max, de_min)
                } else {
                    (de_min, de_max)
                };
                combine(&tau_plus, w_plus, &tau_minus, w_minus)
            } else {
                let n_plus = norm_3n(&tau_plus);
                let n_minus = norm_3n(&tau_minus);
                let w_plus = if n_plus > TINY_NORM {
                    1.0 / n_plus
                } else {
                    0.0
                };
                let w_minus = if n_minus > TINY_NORM {
                    1.0 / n_minus
                } else {
                    0.0
                };
                combine(&tau_plus, w_plus, &tau_minus, w_minus)
            }
        };

        let norm = norm_3n(&tau);
        if norm > TINY_NORM {
            for t in &mut tau {
                for c in t.iter_mut() {
                    *c /= norm;
                }
            }
        } else {
            for t in &mut tau {
                *t = [0.0; 3];
            }
        }
        tau
    }

    /// Energy of image `i`, or 0 when `energies` is shorter than the chain.
    fn energy_at(&self, i: usize) -> f64 {
        self.energies.get(i).copied().unwrap_or(0.0)
    }

    /// Index of the highest-energy interior image (1 if all are equal).
    fn highest_energy_image(&self) -> usize {
        let n = self.images.len();
        if n < 3 {
            return 0;
        }
        (1..(n - 1))
            .max_by(|&a, &b| {
                self.energy_at(a)
                    .partial_cmp(&self.energy_at(b))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(1)
    }
}

/// Per-atom true force row of image `i`, zero-filled to `n_atoms` atoms.
fn true_force_row(true_forces: &[Vec<[f64; 3]>], i: usize, n_atoms: usize) -> Vec<[f64; 3]> {
    (0..n_atoms)
        .map(|a| {
            true_forces
                .get(i)
                .and_then(|im| im.get(a))
                .copied()
                .unwrap_or([0.0; 3])
        })
        .collect()
}

/// Euclidean distance between two configurations in 3N-space.
fn configuration_distance(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(pa, pb)| dist(pa, pb).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Component-wise difference `a − b` of two configurations.
fn difference(a: &[[f64; 3]], b: &[[f64; 3]]) -> Vec<[f64; 3]> {
    a.iter()
        .zip(b.iter())
        .map(|(pa, pb)| [pa[0] - pb[0], pa[1] - pb[1], pa[2] - pb[2]])
        .collect()
}

/// Linear combination `wa·a + wb·b` of two configuration vectors.
fn combine(a: &[[f64; 3]], wa: f64, b: &[[f64; 3]], wb: f64) -> Vec<[f64; 3]> {
    a.iter()
        .zip(b.iter())
        .map(|(va, vb)| {
            [
                wa * va[0] + wb * vb[0],
                wa * va[1] + wb * vb[1],
                wa * va[2] + wb * vb[2],
            ]
        })
        .collect()
}

/// 3N-dimensional dot product.
fn dot_3n(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(va, vb)| va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2])
        .sum()
}

/// 3N-dimensional Euclidean norm.
fn norm_3n(a: &[[f64; 3]]) -> f64 {
    dot_3n(a, a).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_neb_interpolation_endpoints() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[4.0f64, 0.0, 0.0]];
        let neb = NebMethod::new(start, end, 5, 1.0);
        assert_eq!(neb.images[0][0][0], 0.0);
        assert!((neb.images[4][0][0] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_neb_image_count() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[1.0f64, 0.0, 0.0]];
        let neb = NebMethod::new(start, end, 7, 2.0);
        assert_eq!(neb.images.len(), 7);
    }

    #[test]
    fn test_neb_reaction_coordinate_monotone() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[3.0f64, 0.0, 0.0]];
        let neb = NebMethod::new(start, end, 4, 1.0);
        let s = neb.reaction_coordinate();
        for i in 1..s.len() {
            assert!(s[i] >= s[i - 1]);
        }
        assert!((s[3] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_neb_activation_energy_flat_path() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[1.0f64, 0.0, 0.0]];
        let neb = NebMethod::new(start, end, 3, 1.0);
        assert_eq!(neb.activation_energy(), 0.0);
    }

    #[test]
    fn test_neb_activation_energy_with_barrier() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[2.0f64, 0.0, 0.0]];
        let mut neb = NebMethod::new(start, end, 3, 1.0);
        neb.energies = vec![0.0, 0.5, 0.1];
        assert!((neb.activation_energy() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_neb_force_shape() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[3.0f64, 0.0, 0.0]];
        let neb = NebMethod::new(start, end, 5, 1.0);
        let true_forces = vec![vec![[0.0f64; 3]]; 5];
        let f = neb.neb_force(&true_forces);
        assert_eq!(f.len(), 5);
    }

    #[test]
    fn test_neb_force_removes_tangential_component() {
        // 3 equally spaced images on the x-axis: tangent at the middle image is
        // +x. A true force [3, 5, 0] must lose its tangential (x) part and keep
        // the perpendicular (y) part; equal spacing => zero spring force.
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[2.0f64, 0.0, 0.0]];
        let neb = NebMethod::new(start, end, 3, 1.0);
        let mut true_forces = vec![vec![[0.0f64; 3]]; 3];
        true_forces[1][0] = [3.0, 5.0, 0.0];
        let f = neb.neb_force(&true_forces);
        assert!(f[1][0][0].abs() < 1e-9);
        assert!((f[1][0][1] - 5.0).abs() < 1e-9);
        assert_eq!(f[0][0], [0.0, 0.0, 0.0]);
        assert_eq!(f[2][0], [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_neb_climbing_image_inverts_tangential_force() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[3.0f64, 0.0, 0.0]];
        let mut neb = NebMethod::new(start, end, 5, 1.0);
        neb.energies = vec![0.0, 0.3, 0.8, 0.4, 0.0];
        let mut true_forces = vec![vec![[0.0f64; 3]]; 5];
        true_forces[2][0] = [2.0, 1.0, 0.0];
        let f_ci = neb.climbing_image_neb(&true_forces);
        assert_eq!(f_ci.len(), 1);
        assert!((f_ci[0][0] + 2.0).abs() < 1e-9);
        assert!((f_ci[0][1] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_neb_tangent_opposite_moving_atoms_is_nonzero_3n_vector() {
        // Two atoms moving apart along ±x: the per-atom displacements cancel
        // when summed into a single 3-vector, but the 3N tangent is
        // (−1, 0, 0 | +1, 0, 0)/√2.
        let start = vec![[-1.0f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let end = vec![[-2.0f64, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let neb = NebMethod::new(start, end, 3, 1.0);
        let tau = neb.local_tangent(1);
        let s = 1.0 / 2.0_f64.sqrt();
        assert!((tau[0][0] + s).abs() < 1e-12);
        assert!((tau[1][0] - s).abs() < 1e-12);
        assert!((norm_3n(&tau) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_neb_force_projects_with_full_3n_dot_product() {
        let start = vec![[-1.0f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let end = vec![[-2.0f64, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let neb = NebMethod::new(start, end, 3, 1.0);
        // Stretching force along the path plus a perpendicular y-component.
        let mut true_forces = vec![vec![[0.0f64; 3]; 2]; 3];
        true_forces[1] = vec![[-1.0, 0.5, 0.0], [1.0, 0.0, 0.0]];
        let f = neb.neb_force(&true_forces);
        // F·τ = √2 → F⊥ = F − √2 τ = [[0, 0.5, 0], [0, 0, 0]]; equal spacing
        // gives no spring force.
        assert!(f[1][0][0].abs() < 1e-12);
        assert!((f[1][0][1] - 0.5).abs() < 1e-12);
        assert!(f[1][1][0].abs() < 1e-12);
        // The projected force is orthogonal to the 3N tangent.
        let tau = neb.local_tangent(1);
        assert!(dot_3n(&f[1], &tau).abs() < 1e-12);
    }

    #[test]
    fn test_neb_spring_force_independent_of_atom_count() {
        // Two atoms translating together with unequal image spacing: the 3N
        // spring force has magnitude k·(d_next − d_prev), i.e. k per atom here,
        // not k·√N per atom.
        let k = 2.0;
        let mut neb = NebMethod::new(vec![[0.0; 3]; 2], vec![[0.0; 3]; 2], 3, k);
        neb.images[0] = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        neb.images[1] = vec![[1.0, 0.0, 0.0], [6.0, 0.0, 0.0]];
        neb.images[2] = vec![[3.0, 0.0, 0.0], [8.0, 0.0, 0.0]];
        let f = neb.neb_force(&vec![vec![[0.0f64; 3]; 2]; 3]);
        for atom in 0..2 {
            assert!((f[1][atom][0] - k).abs() < 1e-12);
            assert!(f[1][atom][1].abs() < 1e-12);
        }
        let total = norm_3n(&f[1]);
        let expected = k * (8.0_f64.sqrt() - 2.0_f64.sqrt());
        assert!((total - expected).abs() < 1e-12);
    }

    #[test]
    fn test_neb_tangent_follows_uphill_neighbour() {
        let mut neb = NebMethod::new(vec![[0.0; 3]], vec![[0.0; 3]], 3, 1.0);
        neb.images[0] = vec![[0.0, 0.0, 0.0]];
        neb.images[1] = vec![[1.0, 0.0, 0.0]];
        neb.images[2] = vec![[1.0, 1.0, 0.0]];
        neb.energies = vec![0.0, 1.0, 2.0];
        let tau = neb.local_tangent(1);
        assert!((tau[0][1] - 1.0).abs() < 1e-12 && tau[0][0].abs() < 1e-12);
        neb.energies = vec![2.0, 1.0, 0.0];
        let tau = neb.local_tangent(1);
        assert!((tau[0][0] - 1.0).abs() < 1e-12 && tau[0][1].abs() < 1e-12);
    }

    #[test]
    fn test_neb_climbing_image_opposite_moving_atoms() {
        let start = vec![[-1.0f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let end = vec![[-3.0f64, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let mut neb = NebMethod::new(start, end, 5, 1.0);
        neb.energies = vec![0.0, 0.2, 0.9, 0.3, 0.0];
        let mut true_forces = vec![vec![[0.0f64; 3]; 2]; 5];
        // Pure compression force along the path (pushes atoms together).
        true_forces[2] = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let f_ci = neb.climbing_image_neb(&true_forces);
        // The tangential component is inverted: the atoms are pushed apart.
        assert!((f_ci[0][0] + 1.0).abs() < 1e-12);
        assert!((f_ci[1][0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_neb_two_images() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[1.0f64, 0.0, 0.0]];
        let neb = NebMethod::new(start, end, 2, 1.0);
        assert_eq!(neb.images.len(), 2);
        assert_eq!(neb.images[0][0][0], 0.0);
        assert!((neb.images[1][0][0] - 1.0).abs() < 1e-10);
        assert!(neb.climbing_image_neb(&[]).is_empty());
    }

    #[test]
    fn test_neb_spring_k_stored() {
        let start = vec![[0.0f64, 0.0, 0.0]];
        let end = vec![[1.0f64, 0.0, 0.0]];
        let neb = NebMethod::new(start, end, 3, 5.0);
        assert_eq!(neb.spring_k, 5.0);
    }
}
