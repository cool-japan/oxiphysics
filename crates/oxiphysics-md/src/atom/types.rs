//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::math::Vec3;

/// Predicate-based atom selection (MDAnalysis-style).
#[allow(dead_code)]
pub struct AtomSelection;
impl AtomSelection {
    /// Select all atoms with `element == elem`.
    #[allow(dead_code)]
    pub fn by_element(atoms: &[Atom], elem: Element) -> AtomGroup {
        let indices = atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.element == elem)
            .map(|(i, _)| i)
            .collect();
        AtomGroup::new(format!("element_{elem:?}"), indices)
    }
    /// Select all atoms with `residue_id == rid`.
    #[allow(dead_code)]
    pub fn by_residue(atoms: &[Atom], rid: u32) -> AtomGroup {
        let indices = atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.residue_id == rid)
            .map(|(i, _)| i)
            .collect();
        AtomGroup::new(format!("resid_{rid}"), indices)
    }
    /// Select all atoms whose atom_type contains the substring `needle`.
    #[allow(dead_code)]
    pub fn by_atom_type(atoms: &[Atom], needle: &str) -> AtomGroup {
        let indices = atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.atom_type.contains(needle))
            .map(|(i, _)| i)
            .collect();
        AtomGroup::new(format!("type_{needle}"), indices)
    }
    /// Select atoms with charge greater than `min_charge`.
    #[allow(dead_code)]
    pub fn by_charge_gt(atoms: &[Atom], min_charge: f64) -> AtomGroup {
        let indices = atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.charge > min_charge)
            .map(|(i, _)| i)
            .collect();
        AtomGroup::new("charge_gt".to_string(), indices)
    }
    /// Select atoms within `radius` nm of a point.
    #[allow(dead_code)]
    pub fn within_sphere(atoms: &[Atom], center: [f64; 3], radius: f64) -> AtomGroup {
        let r2 = radius * radius;
        let indices = atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| {
                let dx = a.position[0] - center[0];
                let dy = a.position[1] - center[1];
                let dz = a.position[2] - center[2];
                dx * dx + dy * dy + dz * dz <= r2
            })
            .map(|(i, _)| i)
            .collect();
        AtomGroup::new("sphere_sel", indices)
    }
    /// Select all heavy atoms (non-hydrogen).
    #[allow(dead_code)]
    pub fn heavy_atoms(atoms: &[Atom]) -> AtomGroup {
        let indices = atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.element != Element::H)
            .map(|(i, _)| i)
            .collect();
        AtomGroup::new("heavy", indices)
    }
}
/// A chain is a linear sequence of residues (e.g., a protein chain).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Chain {
    /// Chain identifier (e.g., 'A').
    pub id: char,
    /// Residue indices (into a topology's residue list).
    pub residue_indices: Vec<usize>,
}
impl Chain {
    /// Create a new empty chain with the given identifier.
    #[allow(dead_code)]
    pub fn new(id: char) -> Self {
        Self {
            id,
            residue_indices: Vec::new(),
        }
    }
    /// Add a residue index.
    #[allow(dead_code)]
    pub fn add_residue(&mut self, idx: usize) {
        self.residue_indices.push(idx);
    }
    /// Number of residues in this chain.
    #[allow(dead_code)]
    pub fn n_residues(&self) -> usize {
        self.residue_indices.len()
    }
}
/// Common chemical elements used in MD simulations.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Element {
    /// Hydrogen
    H,
    /// Carbon
    C,
    /// Nitrogen
    N,
    /// Oxygen
    O,
    /// Fluorine
    F,
    /// Phosphorus
    P,
    /// Sulfur
    S,
    /// Chlorine
    Cl,
    /// Bromine
    Br,
    /// Iodine
    I,
    /// Sodium
    Na,
    /// Potassium
    K,
    /// Calcium
    Ca,
    /// Magnesium
    Mg,
    /// Iron
    Fe,
    /// Zinc
    Zn,
}
impl Element {
    /// Atomic mass in atomic mass units (u).
    #[allow(dead_code)]
    pub fn mass(self) -> f64 {
        match self {
            Element::H => 1.008,
            Element::C => 12.011,
            Element::N => 14.007,
            Element::O => 15.999,
            Element::F => 18.998,
            Element::P => 30.974,
            Element::S => 32.060,
            Element::Cl => 35.450,
            Element::Br => 79.904,
            Element::I => 126.904,
            Element::Na => 22.990,
            Element::K => 39.098,
            Element::Ca => 40.078,
            Element::Mg => 24.305,
            Element::Fe => 55.845,
            Element::Zn => 65.380,
        }
    }
    /// Van der Waals radius in Angstrom (Å).
    #[allow(dead_code)]
    pub fn van_der_waals_radius(self) -> f64 {
        match self {
            Element::H => 1.20,
            Element::C => 1.70,
            Element::N => 1.55,
            Element::O => 1.52,
            Element::F => 1.47,
            Element::P => 1.80,
            Element::S => 1.80,
            Element::Cl => 1.75,
            Element::Br => 1.85,
            Element::I => 1.98,
            Element::Na => 2.27,
            Element::K => 2.75,
            Element::Ca => 2.31,
            Element::Mg => 1.73,
            Element::Fe => 2.04,
            Element::Zn => 2.10,
        }
    }
    /// Covalent radius in Angstrom (Å).
    #[allow(dead_code)]
    pub fn covalent_radius(self) -> f64 {
        match self {
            Element::H => 0.31,
            Element::C => 0.76,
            Element::N => 0.71,
            Element::O => 0.66,
            Element::F => 0.57,
            Element::P => 1.07,
            Element::S => 1.05,
            Element::Cl => 1.02,
            Element::Br => 1.20,
            Element::I => 1.39,
            Element::Na => 1.66,
            Element::K => 2.03,
            Element::Ca => 1.76,
            Element::Mg => 1.41,
            Element::Fe => 1.32,
            Element::Zn => 1.22,
        }
    }
    /// Pauling electronegativity (dimensionless).
    #[allow(dead_code)]
    pub fn electronegativity(self) -> f64 {
        match self {
            Element::H => 2.20,
            Element::C => 2.55,
            Element::N => 3.04,
            Element::O => 3.44,
            Element::F => 3.98,
            Element::P => 2.19,
            Element::S => 2.58,
            Element::Cl => 3.16,
            Element::Br => 2.96,
            Element::I => 2.66,
            Element::Na => 0.93,
            Element::K => 0.82,
            Element::Ca => 1.00,
            Element::Mg => 1.31,
            Element::Fe => 1.83,
            Element::Zn => 1.65,
        }
    }
    /// Atomic number (proton count).
    #[allow(dead_code)]
    pub fn atomic_number(self) -> u32 {
        match self {
            Element::H => 1,
            Element::C => 6,
            Element::N => 7,
            Element::O => 8,
            Element::F => 9,
            Element::P => 15,
            Element::S => 16,
            Element::Cl => 17,
            Element::Br => 35,
            Element::I => 53,
            Element::Na => 11,
            Element::K => 19,
            Element::Ca => 20,
            Element::Mg => 12,
            Element::Fe => 26,
            Element::Zn => 30,
        }
    }
}
/// A collection of atoms stored in Structure-of-Arrays (SoA) layout.
///
/// This layout is cache-friendly for MD simulations where
/// we typically iterate over one property at a time (e.g., all positions).
#[derive(Debug, Clone)]
pub struct AtomSet {
    /// Positions of all atoms.
    pub positions: Vec<Vec3>,
    /// Velocities of all atoms.
    pub velocities: Vec<Vec3>,
    /// Forces acting on each atom.
    pub forces: Vec<Vec3>,
    /// Mass of each atom.
    pub masses: Vec<f64>,
    /// Charge of each atom.
    pub charges: Vec<f64>,
    /// Type index of each atom (used to look up potential parameters).
    pub types: Vec<u32>,
}
impl AtomSet {
    /// Create a new empty `AtomSet`.
    pub fn new() -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            forces: Vec::new(),
            masses: Vec::new(),
            charges: Vec::new(),
            types: Vec::new(),
        }
    }
    /// Create an `AtomSet` with pre-allocated capacity.
    pub fn with_capacity(cap: usize) -> Self {
        Self {
            positions: Vec::with_capacity(cap),
            velocities: Vec::with_capacity(cap),
            forces: Vec::with_capacity(cap),
            masses: Vec::with_capacity(cap),
            charges: Vec::with_capacity(cap),
            types: Vec::with_capacity(cap),
        }
    }
    /// Add an atom with the given properties.
    pub fn add_atom(
        &mut self,
        position: Vec3,
        velocity: Vec3,
        mass: f64,
        charge: f64,
        atom_type: u32,
    ) {
        self.positions.push(position);
        self.velocities.push(velocity);
        self.forces.push(Vec3::zeros());
        self.masses.push(mass);
        self.charges.push(charge);
        self.types.push(atom_type);
    }
    /// Remove the atom at index `idx`, swapping with the last element.
    ///
    /// O(1) removal; does not preserve ordering.
    #[allow(dead_code)]
    pub fn remove_atom(&mut self, idx: usize) {
        let last = self.positions.len() - 1;
        self.positions.swap(idx, last);
        self.velocities.swap(idx, last);
        self.forces.swap(idx, last);
        self.masses.swap(idx, last);
        self.charges.swap(idx, last);
        self.types.swap(idx, last);
        self.positions.pop();
        self.velocities.pop();
        self.forces.pop();
        self.masses.pop();
        self.charges.pop();
        self.types.pop();
    }
    /// Return the number of atoms.
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    /// Return `true` if there are no atoms.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
    /// Zero out all forces.
    pub fn clear_forces(&mut self) {
        for f in &mut self.forces {
            *f = Vec3::zeros();
        }
    }
    /// Clear all atoms, resetting to empty.
    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.positions.clear();
        self.velocities.clear();
        self.forces.clear();
        self.masses.clear();
        self.charges.clear();
        self.types.clear();
    }
    /// Compute the total kinetic energy: sum of 0.5 * m * v^2.
    pub fn kinetic_energy(&self) -> f64 {
        self.masses
            .iter()
            .zip(self.velocities.iter())
            .map(|(m, v)| 0.5 * m * v.norm_squared())
            .sum()
    }
    /// Compute the instantaneous temperature from the equipartition theorem.
    ///
    /// T = 2 * KE / (N_dof * k_B), where N_dof = 3*N for N atoms.
    /// Returns 0.0 if there are no atoms.
    pub fn temperature(&self, boltzmann_k: f64) -> f64 {
        let n = self.len();
        if n == 0 {
            return 0.0;
        }
        let ndof = 3 * n;
        2.0 * self.kinetic_energy() / (ndof as f64 * boltzmann_k)
    }
    /// Compute the center of mass position.
    ///
    /// Returns the zero vector if total mass is zero or there are no atoms.
    pub fn center_of_mass(&self) -> Vec3 {
        let mut total_mass = 0.0;
        let mut com = Vec3::zeros();
        for (m, r) in self.masses.iter().zip(self.positions.iter()) {
            com += r * *m;
            total_mass += m;
        }
        if total_mass > 0.0 {
            com / total_mass
        } else {
            Vec3::zeros()
        }
    }
    /// Compute the total momentum: sum of m * v.
    pub fn total_momentum(&self) -> Vec3 {
        let mut p = Vec3::zeros();
        for (m, v) in self.masses.iter().zip(self.velocities.iter()) {
            p += v * *m;
        }
        p
    }
    /// Remove the center-of-mass velocity from all atoms.
    pub fn remove_com_velocity(&mut self) {
        let p = self.total_momentum();
        let total_mass: f64 = self.masses.iter().sum();
        if total_mass > 0.0 {
            let v_com = p / total_mass;
            for v in &mut self.velocities {
                *v -= v_com;
            }
        }
    }
    /// Compute the radius of gyration about the center of mass.
    ///
    /// Rg = sqrt( sum_i m_i |r_i - r_com|^2 / M_total )
    ///
    /// Returns 0.0 for empty or zero-mass systems.
    #[allow(dead_code)]
    pub fn radius_of_gyration(&self) -> f64 {
        let total_mass: f64 = self.masses.iter().sum();
        if total_mass <= 0.0 || self.positions.is_empty() {
            return 0.0;
        }
        let com = self.center_of_mass();
        let sum_mr2: f64 = self
            .masses
            .iter()
            .zip(self.positions.iter())
            .map(|(m, r)| {
                let dr = *r - com;
                m * dr.norm_squared()
            })
            .sum();
        (sum_mr2 / total_mass).sqrt()
    }
    /// Apply periodic boundary conditions: wrap all positions into \[0, L).
    #[allow(dead_code)]
    pub fn apply_pbc(&mut self, box_lengths: [f64; 3]) {
        for pos in &mut self.positions {
            for k in 0..3 {
                let l = box_lengths[k];
                if l > 0.0 {
                    pos[k] -= (pos[k] / l).floor() * l;
                }
            }
        }
    }
    /// Compute the total charge of the system.
    #[allow(dead_code)]
    pub fn total_charge(&self) -> f64 {
        self.charges.iter().sum()
    }
    /// Compute the total mass of the system.
    #[allow(dead_code)]
    pub fn total_mass(&self) -> f64 {
        self.masses.iter().sum()
    }
}
impl AtomSet {
    /// Compute a coarse-grained contact map: returns an n×n boolean matrix
    /// where entry \[i\]\[j\] is `true` if atoms i and j are within `cutoff`.
    #[allow(dead_code)]
    pub fn contact_map(&self, cutoff: f64) -> Vec<Vec<bool>> {
        let n = self.len();
        let mut map = vec![vec![false; n]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let p = &self.positions[i];
                let q = &self.positions[j];
                let dx = p.x - q.x;
                let dy = p.y - q.y;
                let dz = p.z - q.z;
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 <= cutoff * cutoff {
                    map[i][j] = true;
                    map[j][i] = true;
                }
            }
        }
        map
    }
    /// Scale all velocities to match a target temperature using the equipartition theorem.
    ///
    /// T_current must be > 0. Does nothing if T_current ≈ 0.
    #[allow(dead_code)]
    pub fn rescale_velocities(&mut self, target_temperature: f64, boltzmann_k: f64) {
        let t_current = self.temperature(boltzmann_k);
        if t_current < 1e-10 {
            return;
        }
        let scale = (target_temperature / t_current).sqrt();
        for v in &mut self.velocities {
            v.x *= scale;
            v.y *= scale;
            v.z *= scale;
        }
    }
    /// Compute the root-mean-square velocity (nm/ps).
    #[allow(dead_code)]
    pub fn rms_velocity(&self) -> f64 {
        if self.is_empty() {
            return 0.0;
        }
        let sum_v2: f64 = self.velocities.iter().map(|v| v.norm_squared()).sum();
        (sum_v2 / self.len() as f64).sqrt()
    }
    /// Translate all atoms by a displacement vector.
    #[allow(dead_code)]
    pub fn translate(&mut self, displacement: [f64; 3]) {
        for pos in &mut self.positions {
            pos.x += displacement[0];
            pos.y += displacement[1];
            pos.z += displacement[2];
        }
    }
    /// Return the indices of atoms within a spherical cutoff of a reference point.
    #[allow(dead_code)]
    pub fn atoms_in_sphere(&self, center: [f64; 3], radius: f64) -> Vec<usize> {
        let r2 = radius * radius;
        self.positions
            .iter()
            .enumerate()
            .filter(|(_, p)| {
                let dx = p.x - center[0];
                let dy = p.y - center[1];
                let dz = p.z - center[2];
                dx * dx + dy * dy + dz * dz <= r2
            })
            .map(|(i, _)| i)
            .collect()
    }
}
impl AtomSet {
    /// Return the minimum and maximum coordinates across all dimensions.
    ///
    /// Returns `([xmin, ymin, zmin], [xmax, ymax, zmax])`.
    /// Returns zeros if empty.
    #[allow(dead_code)]
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        if self.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut lo = [f64::MAX; 3];
        let mut hi = [f64::MIN; 3];
        for p in &self.positions {
            let arr = [p.x, p.y, p.z];
            for k in 0..3 {
                if arr[k] < lo[k] {
                    lo[k] = arr[k];
                }
                if arr[k] > hi[k] {
                    hi[k] = arr[k];
                }
            }
        }
        (lo, hi)
    }
    /// Return the maximum pairwise distance among all atoms.
    #[allow(dead_code)]
    pub fn max_pairwise_distance(&self) -> f64 {
        let n = self.len();
        let mut max_d2 = 0.0_f64;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = self.positions[i].x - self.positions[j].x;
                let dy = self.positions[i].y - self.positions[j].y;
                let dz = self.positions[i].z - self.positions[j].z;
                let d2 = dx * dx + dy * dy + dz * dz;
                if d2 > max_d2 {
                    max_d2 = d2;
                }
            }
        }
        max_d2.sqrt()
    }
    /// Compute the moment of inertia tensor (3×3, row-major) about the center of mass.
    ///
    /// Returns a flat array \[I_xx, I_xy, I_xz, I_yx, I_yy, I_yz, I_zx, I_zy, I_zz\].
    #[allow(dead_code)]
    pub fn inertia_tensor(&self) -> [f64; 9] {
        let com = self.center_of_mass();
        let mut tensor = [0.0_f64; 9];
        for (m, p) in self.masses.iter().zip(self.positions.iter()) {
            let r = [p.x - com.x, p.y - com.y, p.z - com.z];
            let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
            for a in 0..3 {
                for b in 0..3 {
                    let delta = if a == b { 1.0 } else { 0.0 };
                    tensor[a * 3 + b] += m * (r2 * delta - r[a] * r[b]);
                }
            }
        }
        tensor
    }
    /// Rotate all positions by a rotation matrix R (row-major 3×3) about the origin.
    #[allow(dead_code)]
    pub fn rotate(&mut self, r: [f64; 9]) {
        for pos in &mut self.positions {
            let x = r[0] * pos.x + r[1] * pos.y + r[2] * pos.z;
            let y = r[3] * pos.x + r[4] * pos.y + r[5] * pos.z;
            let z = r[6] * pos.x + r[7] * pos.y + r[8] * pos.z;
            pos.x = x;
            pos.y = y;
            pos.z = z;
        }
    }
    /// Compute the root-mean-square displacement (RMSD) of positions from a reference set.
    ///
    /// Returns 0.0 if atoms are empty or lengths differ.
    #[allow(dead_code)]
    pub fn rmsd(&self, reference: &[Vec3]) -> f64 {
        let n = self.len();
        if n == 0 || n != reference.len() {
            return 0.0;
        }
        let sum_d2: f64 = self
            .positions
            .iter()
            .zip(reference.iter())
            .map(|(p, r)| {
                let dx = p.x - r.x;
                let dy = p.y - r.y;
                let dz = p.z - r.z;
                dx * dx + dy * dy + dz * dz
            })
            .sum();
        (sum_d2 / n as f64).sqrt()
    }
    /// Scale all positions relative to the origin by a scalar factor.
    #[allow(dead_code)]
    pub fn scale_positions(&mut self, factor: f64) {
        for pos in &mut self.positions {
            pos.x *= factor;
            pos.y *= factor;
            pos.z *= factor;
        }
    }
    /// Compute the dipole moment (sum of q_i * r_i).
    #[allow(dead_code)]
    pub fn dipole_moment(&self) -> Vec3 {
        let mut d = Vec3::zeros();
        for (q, p) in self.charges.iter().zip(self.positions.iter()) {
            d += *p * *q;
        }
        d
    }
    /// Return the index of the atom with maximum speed.
    #[allow(dead_code)]
    pub fn fastest_atom(&self) -> Option<usize> {
        self.velocities
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                a.norm_squared()
                    .partial_cmp(&b.norm_squared())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
    }
    /// Return the index of the heaviest atom.
    #[allow(dead_code)]
    pub fn heaviest_atom(&self) -> Option<usize> {
        self.masses
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
    }
    /// Return indices of all atoms with charge != 0.
    #[allow(dead_code)]
    pub fn charged_atoms(&self) -> Vec<usize> {
        self.charges
            .iter()
            .enumerate()
            .filter(|(_, q)| q.abs() > 1e-15)
            .map(|(i, _)| i)
            .collect()
    }
    /// Compute the total force magnitude (norm of force sum).
    #[allow(dead_code)]
    pub fn total_force_norm(&self) -> f64 {
        let mut fx = 0.0;
        let mut fy = 0.0;
        let mut fz = 0.0;
        for f in &self.forces {
            fx += f.x;
            fy += f.y;
            fz += f.z;
        }
        (fx * fx + fy * fy + fz * fz).sqrt()
    }
    /// Compute the maximum force magnitude among all atoms.
    #[allow(dead_code)]
    pub fn max_force(&self) -> f64 {
        self.forces
            .iter()
            .map(|f| f.norm_squared())
            .fold(0.0_f64, f64::max)
            .sqrt()
    }
    /// Mirror all positions about the XY plane (flip Z coordinates).
    #[allow(dead_code)]
    pub fn mirror_z(&mut self) {
        for pos in &mut self.positions {
            pos.z = -pos.z;
        }
    }
}
/// An atom group is a named set of atom indices, analogous to MDAnalysis AtomGroup.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AtomGroup {
    /// Group name.
    pub name: String,
    /// Atom indices in the group.
    pub indices: Vec<usize>,
}
impl AtomGroup {
    /// Create a new named group.
    #[allow(dead_code)]
    pub fn new(name: impl Into<String>, indices: Vec<usize>) -> Self {
        Self {
            name: name.into(),
            indices,
        }
    }
    /// Number of atoms in the group.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.indices.len()
    }
    /// Returns true if the group is empty.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }
    /// Return the union of two groups (no duplicate indices).
    #[allow(dead_code)]
    pub fn union(&self, other: &AtomGroup) -> AtomGroup {
        let mut idx: Vec<usize> = self.indices.clone();
        for &i in &other.indices {
            if !idx.contains(&i) {
                idx.push(i);
            }
        }
        idx.sort_unstable();
        AtomGroup::new(format!("{}_union_{}", self.name, other.name), idx)
    }
    /// Return the intersection of two groups.
    #[allow(dead_code)]
    pub fn intersection(&self, other: &AtomGroup) -> AtomGroup {
        let idx: Vec<usize> = self
            .indices
            .iter()
            .filter(|&&i| other.indices.contains(&i))
            .copied()
            .collect();
        AtomGroup::new(format!("{}_inter_{}", self.name, other.name), idx)
    }
    /// Return the difference (self minus other).
    #[allow(dead_code)]
    pub fn difference(&self, other: &AtomGroup) -> AtomGroup {
        let idx: Vec<usize> = self
            .indices
            .iter()
            .filter(|&&i| !other.indices.contains(&i))
            .copied()
            .collect();
        AtomGroup::new(format!("{}_diff_{}", self.name, other.name), idx)
    }
    /// Extract the positions from an AtomSet for this group.
    #[allow(dead_code)]
    pub fn positions<'a>(&self, atoms: &'a AtomSet) -> Vec<&'a Vec3> {
        self.indices.iter().map(|&i| &atoms.positions[i]).collect()
    }
    /// Compute the centroid position (geometric center) of this group.
    #[allow(dead_code)]
    pub fn centroid(&self, atoms: &AtomSet) -> Vec3 {
        if self.indices.is_empty() {
            return Vec3::zeros();
        }
        let mut sum = Vec3::zeros();
        for &i in &self.indices {
            sum += atoms.positions[i];
        }
        sum / (self.indices.len() as f64)
    }
    /// Compute the center of mass of this group.
    #[allow(dead_code)]
    pub fn center_of_mass(&self, atoms: &AtomSet) -> Vec3 {
        let mut total_mass = 0.0;
        let mut com = Vec3::zeros();
        for &i in &self.indices {
            let m = atoms.masses[i];
            com += atoms.positions[i] * m;
            total_mass += m;
        }
        if total_mass > 0.0 {
            com / total_mass
        } else {
            Vec3::zeros()
        }
    }
    /// Compute the total mass of the group.
    #[allow(dead_code)]
    pub fn total_mass(&self, atoms: &AtomSet) -> f64 {
        self.indices.iter().map(|&i| atoms.masses[i]).sum()
    }
}
/// A bond entry in the graph: connects atom `i` and atom `j` with a given bond order.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BondEntry {
    /// Index of first atom.
    pub atom_i: usize,
    /// Index of second atom.
    pub atom_j: usize,
    /// Bond order (1.0 = single, 2.0 = double, etc.).
    pub bond_order: f64,
}
/// A single atom with all physical properties.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Atom {
    /// Cartesian position in nm.
    pub position: [f64; 3],
    /// Velocity in nm/ps.
    pub velocity: [f64; 3],
    /// Force in kJ/mol/nm.
    pub force: [f64; 3],
    /// Mass in atomic mass units.
    pub mass: f64,
    /// Partial charge in elementary charge units.
    pub charge: f64,
    /// Chemical element.
    pub element: Element,
    /// Atom type string (e.g. "CA", "NH1").
    pub atom_type: String,
    /// Residue identifier.
    pub residue_id: u32,
}
impl Atom {
    /// Create a new atom with given element, position, and residue id.
    /// Velocity and force are zeroed; mass and charge taken from element defaults.
    #[allow(dead_code)]
    pub fn new(element: Element, position: [f64; 3], residue_id: u32) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            mass: element.mass(),
            charge: 0.0,
            element,
            atom_type: format!("{element:?}"),
            residue_id,
        }
    }
    /// Kinetic energy of this single atom: 0.5 * m * |v|^2.
    #[allow(dead_code)]
    pub fn kinetic_energy(&self) -> f64 {
        let v = &self.velocity;
        0.5 * self.mass * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    }
}
/// A residue groups atoms belonging to the same chemical unit (e.g., amino acid).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Residue {
    /// Unique residue identifier.
    pub id: u32,
    /// Three-letter residue name (e.g., "ALA", "GLY").
    pub name: String,
    /// Chain identifier (e.g., 'A').
    pub chain_id: char,
    /// Atom indices belonging to this residue.
    pub atom_indices: Vec<usize>,
}
impl Residue {
    /// Create a new residue with the given id, name, and chain.
    #[allow(dead_code)]
    pub fn new(id: u32, name: impl Into<String>, chain_id: char) -> Self {
        Self {
            id,
            name: name.into(),
            chain_id,
            atom_indices: Vec::new(),
        }
    }
    /// Add an atom index to this residue.
    #[allow(dead_code)]
    pub fn add_atom(&mut self, idx: usize) {
        self.atom_indices.push(idx);
    }
    /// Number of atoms in this residue.
    #[allow(dead_code)]
    pub fn n_atoms(&self) -> usize {
        self.atom_indices.len()
    }
}
/// Adjacency-list bond graph for a set of atoms.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct BondGraph {
    /// All bonds stored as a flat list.
    pub bonds: Vec<BondEntry>,
    /// Number of atoms the graph covers.
    pub n_atoms: usize,
}
impl BondGraph {
    /// Create an empty bond graph for `n_atoms` atoms.
    #[allow(dead_code)]
    pub fn new(n_atoms: usize) -> Self {
        Self {
            bonds: Vec::new(),
            n_atoms,
        }
    }
    /// Add a bond between atoms `i` and `j` with given order.
    #[allow(dead_code)]
    pub fn add_bond(&mut self, atom_i: usize, atom_j: usize, bond_order: f64) {
        self.bonds.push(BondEntry {
            atom_i,
            atom_j,
            bond_order,
        });
    }
    /// Return all bonds involving atom `idx`.
    #[allow(dead_code)]
    pub fn bonds_of(&self, idx: usize) -> Vec<&BondEntry> {
        self.bonds
            .iter()
            .filter(|b| b.atom_i == idx || b.atom_j == idx)
            .collect()
    }
    /// Return the number of bonds in the graph.
    #[allow(dead_code)]
    pub fn num_bonds(&self) -> usize {
        self.bonds.len()
    }
    /// Return the bond order between atoms `i` and `j`, or 0.0 if none.
    #[allow(dead_code)]
    pub fn bond_order(&self, i: usize, j: usize) -> f64 {
        self.bonds
            .iter()
            .find(|b| (b.atom_i == i && b.atom_j == j) || (b.atom_i == j && b.atom_j == i))
            .map(|b| b.bond_order)
            .unwrap_or(0.0)
    }
    /// Check whether atoms `i` and `j` are bonded.
    #[allow(dead_code)]
    pub fn are_bonded(&self, i: usize, j: usize) -> bool {
        self.bond_order(i, j) > 0.0
    }
    /// Clear all bonds.
    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.bonds.clear();
    }
}
/// Molecular topology: atoms → residues → chains.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct Topology {
    /// All atoms.
    pub atoms: Vec<Atom>,
    /// All residues.
    pub residues: Vec<Residue>,
    /// All chains.
    pub chains: Vec<Chain>,
    /// Bond graph.
    pub bonds: BondGraph,
}
impl Topology {
    /// Create an empty topology.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self::default()
    }
    /// Add an atom to the topology, returning its index.
    #[allow(dead_code)]
    pub fn add_atom(&mut self, atom: Atom) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(atom);
        idx
    }
    /// Add a residue and return its index.
    #[allow(dead_code)]
    pub fn add_residue(&mut self, residue: Residue) -> usize {
        let idx = self.residues.len();
        self.residues.push(residue);
        idx
    }
    /// Add a chain and return its index.
    #[allow(dead_code)]
    pub fn add_chain(&mut self, chain: Chain) -> usize {
        let idx = self.chains.len();
        self.chains.push(chain);
        idx
    }
    /// Return atoms belonging to a given residue id.
    #[allow(dead_code)]
    pub fn atoms_in_residue(&self, residue_id: u32) -> Vec<&Atom> {
        self.atoms
            .iter()
            .filter(|a| a.residue_id == residue_id)
            .collect()
    }
    /// Return residues belonging to a given chain id.
    #[allow(dead_code)]
    pub fn residues_in_chain(&self, chain_id: char) -> Vec<&Residue> {
        self.residues
            .iter()
            .filter(|r| r.chain_id == chain_id)
            .collect()
    }
    /// Total number of atoms.
    #[allow(dead_code)]
    pub fn n_atoms(&self) -> usize {
        self.atoms.len()
    }
    /// Total number of residues.
    #[allow(dead_code)]
    pub fn n_residues(&self) -> usize {
        self.residues.len()
    }
    /// Total number of chains.
    #[allow(dead_code)]
    pub fn n_chains(&self) -> usize {
        self.chains.len()
    }
}
/// Exclusion policy for non-bonded pairs.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExclusionPolicy {
    /// Exclude 1-2 bonded pairs only.
    OneTwoOnly,
    /// Exclude 1-2 and 1-3 pairs.
    OneTwoThree,
    /// Exclude 1-2, 1-3, and 1-4 pairs.
    OneTwoThreeFour,
}
