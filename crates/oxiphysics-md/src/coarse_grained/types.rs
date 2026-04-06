//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
/// Elastic Network Model node.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct EnmNode {
    /// Node position (nm).
    pub pos: [f64; 3],
    /// Node name (optional).
    pub name: String,
}
impl EnmNode {
    /// Create a new ENM node.
    #[allow(dead_code)]
    pub fn new(pos: [f64; 3], name: impl Into<String>) -> Self {
        Self {
            pos,
            name: name.into(),
        }
    }
}
/// Builders for simple polymer topologies.
pub struct PolymerChain;
impl PolymerChain {
    /// Build a linear polymer of `n_beads` beads placed along the x-axis.
    ///
    /// Returns a `CgSystem` with N-1 harmonic bonds.
    pub fn linear_chain(n_beads: usize, bead_mass: f64, r0: f64, kb: f64) -> CgSystem {
        let cutoff = 2.0;
        let mut sys = CgSystem::new(cutoff);
        for idx in 0..n_beads {
            let x = idx as f64 * r0;
            let bead = CgBead::new(format!("B{}", idx), 0, bead_mass, 0.0, [x, 0.0, 0.0]);
            sys.add_bead(bead);
        }
        for idx in 0..(n_beads - 1) {
            sys.bonds.push(CgBond::new(idx, idx + 1, r0, kb));
        }
        sys
    }
    /// Build a ring polymer of `n_beads` beads placed on a circle.
    ///
    /// Returns a `CgSystem` with N harmonic bonds (including the wrap-around bond).
    pub fn ring_polymer(n_beads: usize, bead_mass: f64, r0: f64, kb: f64) -> CgSystem {
        let cutoff = 2.0;
        let mut sys = CgSystem::new(cutoff);
        let angle_step = 2.0 * std::f64::consts::PI / n_beads as f64;
        let radius = if n_beads > 1 {
            r0 / (2.0 * (angle_step / 2.0).sin())
        } else {
            r0
        };
        for idx in 0..n_beads {
            let theta = idx as f64 * angle_step;
            let pos = [radius * theta.cos(), radius * theta.sin(), 0.0];
            let bead = CgBead::new(format!("B{}", idx), 0, bead_mass, 0.0, pos);
            sys.add_bead(bead);
        }
        for idx in 0..n_beads {
            let next = (idx + 1) % n_beads;
            sys.bonds.push(CgBond::new(idx, next, r0, kb));
        }
        sys
    }
}
/// Mapping from atomistic to coarse-grained beads.
///
/// Each CG bead is a geometric centroid of a subset of atomistic atoms.
pub struct BeadMapping {
    /// For each CG bead: list of atomistic atom indices that belong to it.
    pub atom_groups: Vec<Vec<usize>>,
    /// Number of CG beads.
    pub n_beads: usize,
}
impl BeadMapping {
    /// Create a bead mapping with the given atom groups.
    pub fn new(atom_groups: Vec<Vec<usize>>) -> Self {
        let n_beads = atom_groups.len();
        Self {
            atom_groups,
            n_beads,
        }
    }
    /// Compute CG bead positions as the centroid of each atom group.
    pub fn map_positions(&self, aa_positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        self.atom_groups
            .iter()
            .map(|group| {
                if group.is_empty() {
                    return [0.0; 3];
                }
                let n = group.len() as f64;
                let mut centroid = [0.0; 3];
                for &idx in group {
                    for k in 0..3 {
                        centroid[k] += aa_positions[idx][k];
                    }
                }
                for k in 0..3 {
                    centroid[k] /= n;
                }
                centroid
            })
            .collect()
    }
    /// Compute CG bead velocities as mass-weighted average of atom velocities.
    ///
    /// `masses[i]` is the mass of atomistic atom `i`.
    pub fn map_velocities(&self, aa_velocities: &[[f64; 3]], masses: &[f64]) -> Vec<[f64; 3]> {
        self.atom_groups
            .iter()
            .map(|group| {
                if group.is_empty() {
                    return [0.0; 3];
                }
                let total_mass: f64 = group.iter().map(|&i| masses[i]).sum();
                if total_mass < 1e-30 {
                    return [0.0; 3];
                }
                let mut vel = [0.0; 3];
                for &idx in group {
                    for k in 0..3 {
                        vel[k] += masses[idx] * aa_velocities[idx][k];
                    }
                }
                for k in 0..3 {
                    vel[k] /= total_mass;
                }
                vel
            })
            .collect()
    }
    /// Compute total mass of each CG bead.
    pub fn bead_masses(&self, aa_masses: &[f64]) -> Vec<f64> {
        self.atom_groups
            .iter()
            .map(|group| group.iter().map(|&i| aa_masses[i]).sum())
            .collect()
    }
}
/// A complete coarse-grained system with beads, bonds, angles, dihedrals,
/// and non-bonded interactions.
pub struct CgSystem {
    /// All beads in the system.
    pub beads: Vec<CgBead>,
    /// Harmonic bonds.
    pub bonds: Vec<CgBond>,
    /// Harmonic angles.
    pub angles: Vec<CgAngle>,
    /// Periodic dihedrals.
    pub dihedrals: Vec<CgDihedral>,
    /// Non-bonded cutoff in nm.
    pub cutoff: f64,
}
impl CgSystem {
    /// Create an empty CG system with the given non-bonded cutoff.
    pub fn new(cutoff: f64) -> Self {
        Self {
            beads: Vec::new(),
            bonds: Vec::new(),
            angles: Vec::new(),
            dihedrals: Vec::new(),
            cutoff,
        }
    }
    /// Add a bead and return its index.
    pub fn add_bead(&mut self, bead: CgBead) -> usize {
        let idx = self.beads.len();
        self.beads.push(bead);
        idx
    }
    /// Total bonded potential energy (bonds + angles + dihedrals).
    pub fn compute_bonded_energy(&self) -> f64 {
        let e_bond: f64 = self.bonds.iter().map(|b| b.energy(&self.beads)).sum();
        let e_angle: f64 = self.angles.iter().map(|a| a.energy(&self.beads)).sum();
        let e_dih: f64 = self.dihedrals.iter().map(|d| d.energy(&self.beads)).sum();
        e_bond + e_angle + e_dih
    }
    /// Total non-bonded LJ energy (all pairs within cutoff, excluding 1-2 bonded pairs).
    pub fn compute_nonbonded_energy(&self) -> f64 {
        let n = self.beads.len();
        let mut excluded: std::collections::HashSet<(usize, usize)> =
            std::collections::HashSet::new();
        for b in &self.bonds {
            let (a, c) = if b.i < b.j { (b.i, b.j) } else { (b.j, b.i) };
            excluded.insert((a, c));
        }
        let mut energy = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                if excluded.contains(&(i, j)) {
                    continue;
                }
                let rij = sub(self.beads[j].pos, self.beads[i].pos);
                let r = length(rij);
                if r > self.cutoff {
                    continue;
                }
                let (sigma, epsilon) =
                    MartiniLj::get_lj_params(self.beads[i].bead_type, self.beads[j].bead_type);
                energy += MartiniLj::energy(r, sigma, epsilon);
            }
        }
        energy
    }
    /// Compute all forces and accumulate into `bead.force`.
    ///
    /// Forces are reset to zero before accumulation.
    pub fn compute_forces(&mut self) {
        for bead in &mut self.beads {
            bead.force = [0.0; 3];
        }
        for bond in &self.bonds {
            let (fi, fj) = bond.forces(&self.beads);
            self.beads[bond.i].force = add(self.beads[bond.i].force, fi);
            self.beads[bond.j].force = add(self.beads[bond.j].force, fj);
        }
        for angle in &self.angles {
            let (fi, fj, fk) = angle.forces(&self.beads);
            self.beads[angle.i].force = add(self.beads[angle.i].force, fi);
            self.beads[angle.j].force = add(self.beads[angle.j].force, fj);
            self.beads[angle.k].force = add(self.beads[angle.k].force, fk);
        }
        let n = self.beads.len();
        let mut excluded: std::collections::HashSet<(usize, usize)> =
            std::collections::HashSet::new();
        for b in &self.bonds {
            let (a, c) = if b.i < b.j { (b.i, b.j) } else { (b.j, b.i) };
            excluded.insert((a, c));
        }
        let mut nb_forces: Vec<[f64; 3]> = vec![[0.0; 3]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                if excluded.contains(&(i, j)) {
                    continue;
                }
                let rij = sub(self.beads[j].pos, self.beads[i].pos);
                let r = length(rij);
                if r < 1e-10 || r > self.cutoff {
                    continue;
                }
                let (sigma, epsilon) =
                    MartiniLj::get_lj_params(self.beads[i].bead_type, self.beads[j].bead_type);
                let fmag = MartiniLj::force_magnitude(r, sigma, epsilon);
                let fvec = scale(normalize(rij), fmag);
                nb_forces[i] = add(nb_forces[i], fvec);
                nb_forces[j] = sub(nb_forces[j], fvec);
            }
        }
        for i in 0..n {
            self.beads[i].force = add(self.beads[i].force, nb_forces[i]);
        }
    }
    /// Velocity Verlet integration step.
    ///
    /// 1. v(t + dt/2) = v(t) + 0.5 * F(t)/m * dt
    /// 2. r(t + dt)   = r(t) + v(t + dt/2) * dt
    /// 3. Recompute F(t + dt)
    /// 4. v(t + dt)   = v(t + dt/2) + 0.5 * F(t+dt)/m * dt
    pub fn step_verlet(&mut self, dt: f64) {
        let n = self.beads.len();
        for i in 0..n {
            let inv_mass = 1.0 / self.beads[i].mass;
            for a in 0..3 {
                self.beads[i].vel[a] += 0.5 * self.beads[i].force[a] * inv_mass * dt;
                self.beads[i].pos[a] += self.beads[i].vel[a] * dt;
            }
        }
        self.compute_forces();
        for i in 0..n {
            let inv_mass = 1.0 / self.beads[i].mass;
            for a in 0..3 {
                self.beads[i].vel[a] += 0.5 * self.beads[i].force[a] * inv_mass * dt;
            }
        }
    }
}
/// Periodic dihedral: V = kd * (1 + cos(n*phi - phi0)).
#[derive(Debug, Clone)]
pub struct CgDihedral {
    /// Indices of the four beads.
    pub i: usize,
    /// Indices of the four beads.
    pub j: usize,
    /// Indices of the four beads.
    pub k: usize,
    /// Indices of the four beads.
    pub l: usize,
    /// Phase angle in radians.
    pub phi0: f64,
    /// Dihedral stiffness in kJ/mol.
    pub kd: f64,
    /// Multiplicity (periodicity).
    pub multiplicity: u32,
}
impl CgDihedral {
    /// Create a new periodic dihedral.
    pub fn new(
        i: usize,
        j: usize,
        k: usize,
        l: usize,
        phi0: f64,
        kd: f64,
        multiplicity: u32,
    ) -> Self {
        Self {
            i,
            j,
            k,
            l,
            phi0,
            kd,
            multiplicity,
        }
    }
    /// Compute the dihedral angle phi (in radians, range -pi..pi).
    pub fn dihedral_angle(&self, beads: &[CgBead]) -> f64 {
        let b1 = sub(beads[self.j].pos, beads[self.i].pos);
        let b2 = sub(beads[self.k].pos, beads[self.j].pos);
        let b3 = sub(beads[self.l].pos, beads[self.k].pos);
        let n1 = cross(b1, b2);
        let n2 = cross(b2, b3);
        let len_n1 = length(n1);
        let len_n2 = length(n2);
        if len_n1 < 1e-15 || len_n2 < 1e-15 {
            return 0.0;
        }
        let cos_phi = dot(n1, n2) / (len_n1 * len_n2);
        let cos_phi = cos_phi.clamp(-1.0, 1.0);
        let m = cross(n1, b2);
        let sign = if dot(m, n2) < 0.0 { -1.0 } else { 1.0 };
        sign * cos_phi.acos()
    }
    /// Potential energy of this dihedral.
    pub fn energy(&self, beads: &[CgBead]) -> f64 {
        let phi = self.dihedral_angle(beads);
        self.kd * (1.0 + (self.multiplicity as f64 * phi - self.phi0).cos())
    }
}
/// Elastic Network Model for protein dynamics and normal mode analysis.
///
/// All pairs within `cutoff` nm are connected by harmonic springs.
pub struct ElasticNetworkModel {
    /// Reference positions (nm).
    pub positions: Vec<[f64; 3]>,
    /// Distance cutoff for spring connections (nm).
    pub cutoff: f64,
    /// Spring constant (kJ/mol/nm^2).
    pub spring_constant: f64,
}
impl ElasticNetworkModel {
    /// Build an ENM from reference positions.
    pub fn new(positions: Vec<[f64; 3]>, cutoff: f64, k: f64) -> Self {
        Self {
            positions,
            cutoff,
            spring_constant: k,
        }
    }
    /// Return all pairs (i, j) with i < j whose reference distance is within cutoff.
    pub fn contact_map(&self) -> Vec<(usize, usize)> {
        let n = self.positions.len();
        let mut contacts = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let r = length(sub(self.positions[j], self.positions[i]));
                if r <= self.cutoff {
                    contacts.push((i, j));
                }
            }
        }
        contacts
    }
    /// Potential energy given current positions.
    ///
    /// V = 0.5 * k * sum_{(i,j) in contacts} (r_ij - r0_ij)^2
    pub fn potential_energy(&self, current_pos: &[[f64; 3]]) -> f64 {
        let contacts = self.contact_map();
        contacts
            .iter()
            .map(|&(i, j)| {
                let r0 = length(sub(self.positions[j], self.positions[i]));
                let r = length(sub(current_pos[j], current_pos[i]));
                let dr = r - r0;
                0.5 * self.spring_constant * dr * dr
            })
            .sum()
    }
    /// Forces on all beads given current positions.
    pub fn forces(&self, current_pos: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let n = current_pos.len();
        let mut forces = vec![[0.0f64; 3]; n];
        let contacts = self.contact_map();
        for (i, j) in contacts {
            let r0 = length(sub(self.positions[j], self.positions[i]));
            let rij = sub(current_pos[j], current_pos[i]);
            let r = length(rij);
            if r < 1e-15 {
                continue;
            }
            let dr = r - r0;
            let mag = self.spring_constant * dr / r;
            let fij = scale(rij, mag);
            forces[i] = add(forces[i], fij);
            forces[j] = sub(forces[j], fij);
        }
        forces
    }
    /// Build the 3N×3N Hessian matrix (mass-unweighted).
    ///
    /// The Hessian is used for normal mode analysis.
    pub fn hessian(&self) -> Vec<Vec<f64>> {
        let n = self.positions.len();
        let dim = 3 * n;
        let mut h = vec![vec![0.0f64; dim]; dim];
        let contacts = self.contact_map();
        for (i, j) in contacts {
            let rij = sub(self.positions[j], self.positions[i]);
            let r2 = dot(rij, rij);
            if r2 < 1e-15 {
                continue;
            }
            for a in 0..3 {
                for b in 0..3 {
                    let val = self.spring_constant * rij[a] * rij[b] / r2;
                    h[3 * i + a][3 * j + b] -= val;
                    h[3 * j + b][3 * i + a] -= val;
                    h[3 * i + a][3 * i + b] += val;
                    h[3 * j + a][3 * j + b] += val;
                }
            }
        }
        h
    }
}
/// A single coarse-grained water bead (MARTINI W bead: 4 water molecules).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CgWaterBead {
    /// Position (nm).
    pub pos: [f64; 3],
    /// Velocity (nm/ps).
    pub vel: [f64; 3],
    /// Force (kJ/mol/nm).
    pub force: [f64; 3],
}
impl CgWaterBead {
    /// Create a new CG water bead at a given position.
    #[allow(dead_code)]
    pub fn new(pos: [f64; 3]) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            force: [0.0; 3],
        }
    }
    /// MARTINI W bead mass (4 × 18.015 = 72.06 Da).
    #[allow(dead_code)]
    pub fn mass() -> f64 {
        72.06
    }
}
/// MARTINI 2.x non-bonded interaction levels (sigma in nm, epsilon in kJ/mol).
///
/// Each row encodes a bead-type pair interaction following the MARTINI
/// interaction matrix (C-level approximation).
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct MartiniInteraction {
    /// Lennard-Jones sigma parameter (nm).
    pub sigma: f64,
    /// Lennard-Jones epsilon parameter (kJ/mol).
    pub epsilon: f64,
}
impl MartiniInteraction {
    /// Construct from sigma and epsilon.
    #[allow(dead_code)]
    pub fn new(sigma: f64, epsilon: f64) -> Self {
        Self { sigma, epsilon }
    }
    /// Self-interaction for a MARTINI C-level bead (sigma = 0.47 nm, epsilon = 5.6 kJ/mol).
    #[allow(dead_code)]
    pub fn c_level() -> Self {
        Self {
            sigma: 0.47,
            epsilon: 5.6,
        }
    }
    /// Self-interaction for a MARTINI D-level bead (weak: epsilon = 1.5 kJ/mol).
    #[allow(dead_code)]
    pub fn d_level() -> Self {
        Self {
            sigma: 0.47,
            epsilon: 1.5,
        }
    }
    /// LJ 12-6 energy between two beads at distance r.
    #[allow(dead_code)]
    pub fn energy(&self, r: f64) -> f64 {
        MartiniLj::energy(r, self.sigma, self.epsilon)
    }
    /// LJ 12-6 force magnitude (negative = attractive) at distance r.
    #[allow(dead_code)]
    pub fn force_magnitude(&self, r: f64) -> f64 {
        MartiniLj::force_magnitude(r, self.sigma, self.epsilon)
    }
    /// Minimum energy distance r_min = 2^(1/6) * sigma.
    #[allow(dead_code)]
    pub fn r_min(&self) -> f64 {
        2.0_f64.powf(1.0 / 6.0) * self.sigma
    }
}
/// Cosine-series Fourier CG dihedral potential.
///
/// V(φ) = sum_{n=0}^{N} k_n * (1 + cos(n*φ - δ_n))
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CgFourierDihedral {
    /// Bead index a.
    pub a: usize,
    /// Bead index b.
    pub b: usize,
    /// Bead index c.
    pub c: usize,
    /// Bead index d.
    pub d: usize,
    /// Force constants k_n (kJ/mol).
    pub k: Vec<f64>,
    /// Phase offsets δ_n (radians).
    pub delta: Vec<f64>,
}
impl CgFourierDihedral {
    /// Create a new Fourier dihedral.
    #[allow(dead_code)]
    pub fn new(a: usize, b: usize, c: usize, d: usize, k: Vec<f64>, delta: Vec<f64>) -> Self {
        Self {
            a,
            b,
            c,
            d,
            k,
            delta,
        }
    }
    /// Compute the dihedral angle for four bead positions.
    #[allow(dead_code)]
    pub fn dihedral_angle(beads: &[CgBead], a: usize, b: usize, c: usize, d: usize) -> f64 {
        let b1 = sub(beads[b].pos, beads[a].pos);
        let b2 = sub(beads[c].pos, beads[b].pos);
        let b3 = sub(beads[d].pos, beads[c].pos);
        let n1 = cross(b1, b2);
        let n2 = cross(b2, b3);
        let m1 = cross(n1, b2);
        let n1_n = length(n1);
        let n2_n = length(n2);
        let m1_n = length(m1);
        if n1_n < 1e-15 || n2_n < 1e-15 || m1_n < 1e-15 {
            return 0.0;
        }
        let x = (dot(n1, n2) / (n1_n * n2_n)).clamp(-1.0, 1.0);
        let y = dot(m1, n2) / (m1_n * n2_n);
        x.acos().copysign(y)
    }
    /// Potential energy.
    #[allow(dead_code)]
    pub fn energy(&self, beads: &[CgBead]) -> f64 {
        let phi = Self::dihedral_angle(beads, self.a, self.b, self.c, self.d);
        self.k
            .iter()
            .zip(self.delta.iter())
            .enumerate()
            .map(|(n, (k, delta))| k * (1.0 + (n as f64 * phi - delta).cos()))
            .sum()
    }
}
/// Harmonic angle between three beads: V = 0.5 * ka * (theta - theta0)^2.
#[derive(Debug, Clone)]
pub struct CgAngle {
    /// Index of first bead.
    pub i: usize,
    /// Index of central bead.
    pub j: usize,
    /// Index of third bead.
    pub k: usize,
    /// Equilibrium angle in radians.
    pub theta0: f64,
    /// Angle stiffness in kJ/mol/rad^2.
    pub ka: f64,
}
impl CgAngle {
    /// Create a new harmonic angle.
    pub fn new(i: usize, j: usize, k: usize, theta0: f64, ka: f64) -> Self {
        Self {
            i,
            j,
            k,
            theta0,
            ka,
        }
    }
    /// Potential energy of this angle.
    pub fn energy(&self, beads: &[CgBead]) -> f64 {
        let theta = self.angle(beads);
        let dtheta = theta - self.theta0;
        0.5 * self.ka * dtheta * dtheta
    }
    /// Compute the angle theta at bead j.
    fn angle(&self, beads: &[CgBead]) -> f64 {
        let rji = sub(beads[self.i].pos, beads[self.j].pos);
        let rjk = sub(beads[self.k].pos, beads[self.j].pos);
        let cos_theta = dot(rji, rjk) / (length(rji) * length(rjk) + 1e-15);
        cos_theta.clamp(-1.0, 1.0).acos()
    }
    /// Forces on beads i, j, k. Returns (fi, fj, fk).
    pub fn forces(&self, beads: &[CgBead]) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let rji = sub(beads[self.i].pos, beads[self.j].pos);
        let rjk = sub(beads[self.k].pos, beads[self.j].pos);
        let len_ji = length(rji);
        let len_jk = length(rjk);
        if len_ji < 1e-15 || len_jk < 1e-15 {
            return ([0.0; 3], [0.0; 3], [0.0; 3]);
        }
        let u_ji = scale(rji, 1.0 / len_ji);
        let u_jk = scale(rjk, 1.0 / len_jk);
        let cos_theta = dot(u_ji, u_jk).clamp(-1.0, 1.0);
        let theta = cos_theta.acos();
        let sin_theta = theta.sin();
        if sin_theta.abs() < 1e-10 {
            return ([0.0; 3], [0.0; 3], [0.0; 3]);
        }
        let dtheta = theta - self.theta0;
        let prefactor = -self.ka * dtheta / sin_theta;
        let fi: [f64; 3] = scale(sub(scale(u_ji, cos_theta), u_jk), prefactor / len_ji);
        let fk: [f64; 3] = scale(sub(scale(u_jk, cos_theta), u_ji), prefactor / len_jk);
        let fj: [f64; 3] = [-(fi[0] + fk[0]), -(fi[1] + fk[1]), -(fi[2] + fk[2])];
        (fi, fj, fk)
    }
}
/// A Go-like model for protein folding simulation.
///
/// Native contacts are defined from a reference structure; non-native
/// contacts use a purely repulsive potential.
pub struct GoModel {
    /// Reference (native) positions for defining contacts.
    pub native_positions: Vec<[f64; 3]>,
    /// Cutoff distance for defining native contacts (nm).
    pub cutoff: f64,
    /// Attractive well depth ε for native contacts (kJ/mol).
    pub epsilon: f64,
    /// Repulsion coefficient for non-native contacts.
    pub sigma_rep: f64,
    /// Go-model bond stiffness (kJ/mol/nm^2).
    pub k_bond: f64,
    /// Equilibrium bond length (nm).
    pub r0_bond: f64,
}
impl GoModel {
    /// Build a Go model from a reference structure.
    pub fn new(native_positions: Vec<[f64; 3]>, cutoff: f64, epsilon: f64) -> Self {
        Self {
            native_positions,
            cutoff,
            epsilon,
            sigma_rep: 0.4,
            k_bond: 100.0,
            r0_bond: 0.38,
        }
    }
    /// Identify native contact pairs (i, j) with i < j and d_native <= cutoff.
    pub fn native_contacts(&self) -> Vec<(usize, usize, f64)> {
        let n = self.native_positions.len();
        let mut contacts = Vec::new();
        for i in 0..n {
            for j in (i + 2)..n {
                let r = length(sub(self.native_positions[j], self.native_positions[i]));
                if r <= self.cutoff {
                    contacts.push((i, j, r));
                }
            }
        }
        contacts
    }
    /// Potential energy for a native contact pair using a 12-10 Go potential.
    ///
    /// V(r) = ε \[ 5*(r0/r)^12 - 6*(r0/r)^10 \]
    pub fn native_contact_energy(epsilon: f64, r0: f64, r: f64) -> f64 {
        if r < 1e-15 {
            return f64::INFINITY;
        }
        let x = r0 / r;
        epsilon * (5.0 * x.powi(12) - 6.0 * x.powi(10))
    }
    /// Potential energy for a non-native contact (purely repulsive).
    pub fn repulsive_energy(sigma: f64, r: f64) -> f64 {
        if r < 1e-15 {
            return f64::INFINITY;
        }
        let x = sigma / r;
        x.powi(12)
    }
    /// Total Go-model potential energy given current positions.
    pub fn potential_energy(&self, positions: &[[f64; 3]]) -> f64 {
        let n = positions.len();
        let contacts = self.native_contacts();
        let native_set: std::collections::HashSet<(usize, usize)> =
            contacts.iter().map(|&(i, j, _)| (i, j)).collect();
        let mut energy = 0.0;
        for i in 0..(n - 1) {
            let r = length(sub(positions[i + 1], positions[i]));
            let dr = r - self.r0_bond;
            energy += 0.5 * self.k_bond * dr * dr;
        }
        for &(i, j, r0) in &contacts {
            let r = length(sub(positions[j], positions[i]));
            energy += Self::native_contact_energy(self.epsilon, r0, r);
        }
        for i in 0..n {
            for j in (i + 2)..n {
                if !native_set.contains(&(i, j)) {
                    let r = length(sub(positions[j], positions[i]));
                    energy += Self::repulsive_energy(self.sigma_rep, r);
                }
            }
        }
        energy
    }
    /// Q-value: fraction of native contacts formed (within 1.2*r0).
    pub fn q_value(&self, positions: &[[f64; 3]]) -> f64 {
        let contacts = self.native_contacts();
        if contacts.is_empty() {
            return 0.0;
        }
        let formed = contacts
            .iter()
            .filter(|&&(i, j, r0)| {
                let r = length(sub(positions[j], positions[i]));
                r <= 1.2 * r0
            })
            .count();
        formed as f64 / contacts.len() as f64
    }
}
/// Shifted MARTINI LJ with smooth cutoff at r_cut.
///
/// V_shifted(r) = V_LJ(r) - V_LJ(r_cut) - (r - r_cut) * dV/dr|_{r_cut}
#[allow(dead_code)]
pub struct ShiftedMartiniLj {
    /// Base interaction parameters.
    pub params: MartiniInteraction,
    /// Cutoff distance (nm).
    pub r_cut: f64,
    /// V(r_cut) for shifting.
    pub(super) v_cut: f64,
    /// dV/dr at r_cut for shifting.
    pub(super) dv_cut: f64,
}
impl ShiftedMartiniLj {
    /// Construct a shifted LJ potential.
    #[allow(dead_code)]
    pub fn new(params: MartiniInteraction, r_cut: f64) -> Self {
        let v_cut = MartiniLj::energy(r_cut, params.sigma, params.epsilon);
        let dv_cut = -MartiniLj::force_magnitude(r_cut, params.sigma, params.epsilon);
        Self {
            params,
            r_cut,
            v_cut,
            dv_cut,
        }
    }
    /// Shifted potential energy at r.  Returns 0 for r >= r_cut.
    #[allow(dead_code)]
    pub fn energy(&self, r: f64) -> f64 {
        if r >= self.r_cut {
            return 0.0;
        }
        let v = self.params.energy(r);
        v - self.v_cut - (r - self.r_cut) * self.dv_cut
    }
}
/// Harmonic bond between two beads: V = 0.5 * kb * (r - r0)^2.
#[derive(Debug, Clone)]
pub struct CgBond {
    /// Index of first bead.
    pub i: usize,
    /// Index of second bead.
    pub j: usize,
    /// Equilibrium distance in nm.
    pub r0: f64,
    /// Bond stiffness in kJ/mol/nm^2.
    pub kb: f64,
}
impl CgBond {
    /// Create a new harmonic bond.
    pub fn new(i: usize, j: usize, r0: f64, kb: f64) -> Self {
        Self { i, j, r0, kb }
    }
    /// Potential energy of this bond.
    pub fn energy(&self, beads: &[CgBead]) -> f64 {
        let rij = sub(beads[self.j].pos, beads[self.i].pos);
        let r = length(rij);
        let dr = r - self.r0;
        0.5 * self.kb * dr * dr
    }
    /// Forces on beads i and j (returns (fi, fj)).
    pub fn forces(&self, beads: &[CgBead]) -> ([f64; 3], [f64; 3]) {
        let rij = sub(beads[self.j].pos, beads[self.i].pos);
        let r = length(rij);
        if r < 1e-15 {
            return ([0.0; 3], [0.0; 3]);
        }
        let dr = r - self.r0;
        let mag = self.kb * dr / r;
        let fi = scale(rij, mag);
        let fj = scale(rij, -mag);
        (fi, fj)
    }
}
/// MARTINI-inspired Lennard-Jones 12-6 non-bonded interactions.
pub struct MartiniLj;
impl MartiniLj {
    /// LJ 12-6 energy: V(r) = 4*epsilon*\[(sigma/r)^12 - (sigma/r)^6\].
    pub fn energy(r: f64, sigma: f64, epsilon: f64) -> f64 {
        if r < 1e-10 {
            return f64::INFINITY;
        }
        let sr6 = (sigma / r).powi(6);
        4.0 * epsilon * (sr6 * sr6 - sr6)
    }
    /// Magnitude of the LJ force (positive = repulsive).
    /// F(r) = 24*epsilon/r * \[2*(sigma/r)^12 - (sigma/r)^6\].
    pub fn force_magnitude(r: f64, sigma: f64, epsilon: f64) -> f64 {
        if r < 1e-10 {
            return f64::INFINITY;
        }
        let sr6 = (sigma / r).powi(6);
        24.0 * epsilon / r * (2.0 * sr6 * sr6 - sr6)
    }
    /// Return (sigma, epsilon) for a pair of bead types.
    ///
    /// Uses a 5x5 MARTINI-like epsilon table with sigma = 0.47 nm for all pairs.
    /// Bead types: 0=polar(P), 1=nonpolar(N), 2=apolar(C), 3=charged_pos(Qd), 4=charged_neg(Qa).
    pub fn get_lj_params(type_i: u8, type_j: u8) -> (f64, f64) {
        const EPS: [[f64; 5]; 5] = [
            [5.6, 5.0, 2.7, 5.6, 5.6],
            [5.0, 4.0, 3.1, 4.5, 4.5],
            [2.7, 3.1, 3.5, 2.0, 2.0],
            [5.6, 4.5, 2.0, 5.6, 2.0],
            [5.6, 4.5, 2.0, 2.0, 5.6],
        ];
        let i = (type_i as usize).min(4);
        let j = (type_j as usize).min(4);
        (0.47, EPS[i][j])
    }
}
/// Simplified Langevin thermostat parameters.
pub struct LangevinThermostat {
    /// Target temperature (K).
    pub temperature: f64,
    /// Friction coefficient γ (ps^-1).
    pub gamma: f64,
    /// Boltzmann constant (kJ/mol/K).
    pub kb: f64,
}
impl LangevinThermostat {
    /// Create a Langevin thermostat.
    pub fn new(temperature: f64, gamma: f64) -> Self {
        Self {
            temperature,
            gamma,
            kb: 8.314e-3,
        }
    }
    /// Standard deviation of random force for one degree of freedom.
    ///
    /// σ = sqrt(2 * m * γ * k_B * T / dt)
    pub fn random_force_sigma(&self, mass: f64, dt: f64) -> f64 {
        (2.0 * mass * self.gamma * self.kb * self.temperature / dt).sqrt()
    }
    /// Velocity scaling factor per step for the friction term.
    ///
    /// scale = exp(-γ * dt)  (Ornstein-Uhlenbeck integrator)
    pub fn velocity_scale(&self, dt: f64) -> f64 {
        (-self.gamma * dt).exp()
    }
}
/// Extended Elastic Network Model with named nodes (Cα residue-level ENM).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct EnmNodeModel {
    /// Nodes (typically Cα atoms).
    pub nodes: Vec<EnmNode>,
    /// Spring constant (kJ/mol/nm²).
    pub spring_k: f64,
    /// Cutoff distance for spring formation (nm).
    pub cutoff: f64,
}
impl EnmNodeModel {
    /// Create a new named-node ENM.
    #[allow(dead_code)]
    pub fn new(nodes: Vec<EnmNode>, spring_k: f64, cutoff: f64) -> Self {
        Self {
            nodes,
            spring_k,
            cutoff,
        }
    }
    /// Count the number of elastic contacts (springs).
    #[allow(dead_code)]
    pub fn n_contacts(&self) -> usize {
        let n = self.nodes.len();
        let mut count = 0;
        for i in 0..n {
            for j in (i + 1)..n {
                let d = length(sub(self.nodes[i].pos, self.nodes[j].pos));
                if d <= self.cutoff {
                    count += 1;
                }
            }
        }
        count
    }
    /// Compute the displaced energy relative to a reference configuration.
    ///
    /// `ref_nodes` should be the native (equilibrium) configuration.
    #[allow(dead_code)]
    pub fn displaced_energy(&self, ref_nodes: &[EnmNode]) -> f64 {
        let n = self.nodes.len().min(ref_nodes.len());
        let mut e = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r0 = length(sub(ref_nodes[i].pos, ref_nodes[j].pos));
                if r0 > self.cutoff {
                    continue;
                }
                let r = length(sub(self.nodes[i].pos, self.nodes[j].pos));
                let d = r - r0;
                e += 0.5 * self.spring_k * d * d;
            }
        }
        e
    }
    /// Compute the connectivity (number of springs per node).
    #[allow(dead_code)]
    pub fn connectivity(&self) -> Vec<usize> {
        let n = self.nodes.len();
        let mut conn = vec![0usize; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let d = length(sub(self.nodes[i].pos, self.nodes[j].pos));
                if d <= self.cutoff {
                    conn[i] += 1;
                    conn[j] += 1;
                }
            }
        }
        conn
    }
}
/// Cosine-based harmonic angle: V = 0.5 * ka * (cos(theta) - cos(theta0))^2.
#[derive(Debug, Clone)]
pub struct CgCosineAngle {
    /// Index of first bead.
    pub i: usize,
    /// Index of central bead.
    pub j: usize,
    /// Index of third bead.
    pub k: usize,
    /// Equilibrium angle in radians.
    pub theta0: f64,
    /// Force constant (kJ/mol).
    pub ka: f64,
}
impl CgCosineAngle {
    /// Create a new cosine angle.
    pub fn new(i: usize, j: usize, k: usize, theta0: f64, ka: f64) -> Self {
        Self {
            i,
            j,
            k,
            theta0,
            ka,
        }
    }
    /// Potential energy.
    pub fn energy(&self, beads: &[CgBead]) -> f64 {
        let rji = sub(beads[self.i].pos, beads[self.j].pos);
        let rjk = sub(beads[self.k].pos, beads[self.j].pos);
        let len_ji = length(rji);
        let len_jk = length(rjk);
        if len_ji < 1e-15 || len_jk < 1e-15 {
            return 0.0;
        }
        let cos_theta = (dot(rji, rjk) / (len_ji * len_jk)).clamp(-1.0, 1.0);
        let cos_theta0 = self.theta0.cos();
        let d = cos_theta - cos_theta0;
        0.5 * self.ka * d * d
    }
}
/// Tabulated CG pair potential (discretised V(r) and F(r)).
pub struct CgPairPotential {
    /// Minimum r value.
    pub r_min: f64,
    /// Maximum r value.
    pub r_max: f64,
    /// Tabulated V(r) values.
    pub energy_table: Vec<f64>,
    /// Tabulated F(r) = -dV/dr values.
    pub force_table: Vec<f64>,
}
impl CgPairPotential {
    /// Build a tabulated potential from a function V(r).
    pub fn from_function(
        r_min: f64,
        r_max: f64,
        n_points: usize,
        v_fn: impl Fn(f64) -> f64,
    ) -> Self {
        assert!(n_points >= 2);
        let dr = (r_max - r_min) / (n_points - 1) as f64;
        let rs: Vec<f64> = (0..n_points).map(|i| r_min + i as f64 * dr).collect();
        let energy_table: Vec<f64> = rs.iter().map(|&r| v_fn(r)).collect();
        let mut force_table = vec![0.0; n_points];
        for i in 1..(n_points - 1) {
            force_table[i] = -(energy_table[i + 1] - energy_table[i - 1]) / (2.0 * dr);
        }
        force_table[0] = force_table[1];
        force_table[n_points - 1] = force_table[n_points - 2];
        Self {
            r_min,
            r_max,
            energy_table,
            force_table,
        }
    }
    /// Interpolate the energy V(r).
    pub fn energy(&self, r: f64) -> f64 {
        self.interpolate(&self.energy_table, r)
    }
    /// Interpolate the force F(r) = -dV/dr.
    pub fn force(&self, r: f64) -> f64 {
        self.interpolate(&self.force_table, r)
    }
    fn interpolate(&self, table: &[f64], r: f64) -> f64 {
        let n = table.len();
        if r <= self.r_min {
            return table[0];
        }
        if r >= self.r_max {
            return table[n - 1];
        }
        let dr = (self.r_max - self.r_min) / (n - 1) as f64;
        let idx = ((r - self.r_min) / dr) as usize;
        let idx = idx.min(n - 2);
        let t = (r - self.r_min - idx as f64 * dr) / dr;
        table[idx] * (1.0 - t) + table[idx + 1] * t
    }
    /// Number of tabulation points.
    pub fn n_points(&self) -> usize {
        self.energy_table.len()
    }
}
/// A coarse-grained bead (MARTINI-like, representing ~4 heavy atoms).
#[derive(Debug, Clone)]
pub struct CgBead {
    /// Bead name (e.g. "BB", "SC1").
    pub name: String,
    /// Interaction type: 0=polar, 1=nonpolar, 2=apolar, 3=charged_pos, 4=charged_neg.
    pub bead_type: u8,
    /// Mass in Da (~72 Da typical for 4-to-1 mapping).
    pub mass: f64,
    /// Partial charge.
    pub charge: f64,
    /// Position in nm.
    pub pos: [f64; 3],
    /// Velocity in nm/ps.
    pub vel: [f64; 3],
    /// Accumulated force in kJ/mol/nm.
    pub force: [f64; 3],
}
impl CgBead {
    /// Create a new bead at the given position with zero velocity and force.
    pub fn new(
        name: impl Into<String>,
        bead_type: u8,
        mass: f64,
        charge: f64,
        pos: [f64; 3],
    ) -> Self {
        Self {
            name: name.into(),
            bead_type,
            mass,
            charge,
            pos,
            vel: [0.0; 3],
            force: [0.0; 3],
        }
    }
}
/// Full MARTINI bead type parameter set.
#[derive(Debug, Clone)]
pub struct MartiniBeadType {
    /// Bead type name (e.g. "P5", "C1", "Qa").
    pub name: &'static str,
    /// Bead type index used in the epsilon table.
    pub type_idx: u8,
    /// Default mass (Da).
    pub mass: f64,
}
impl MartiniBeadType {
    /// List of standard MARTINI 2.2 bead types.
    pub fn standard_types() -> Vec<MartiniBeadType> {
        vec![
            MartiniBeadType {
                name: "P5",
                type_idx: 0,
                mass: 72.0,
            },
            MartiniBeadType {
                name: "N0",
                type_idx: 1,
                mass: 72.0,
            },
            MartiniBeadType {
                name: "C1",
                type_idx: 2,
                mass: 72.0,
            },
            MartiniBeadType {
                name: "Qd",
                type_idx: 3,
                mass: 72.0,
            },
            MartiniBeadType {
                name: "Qa",
                type_idx: 4,
                mass: 72.0,
            },
        ]
    }
}
