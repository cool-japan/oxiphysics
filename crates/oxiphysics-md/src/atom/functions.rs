//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{Atom, BondGraph, Element, ExclusionPolicy};

/// Return the atomic (covalent) radius in Angstrom for an element.
#[allow(dead_code)]
pub fn atomic_radius(element: Element) -> f64 {
    element.covalent_radius()
}
/// Return the van der Waals radius in Angstrom for an element.
#[allow(dead_code)]
pub fn van_der_waals_radius(element: Element) -> f64 {
    element.van_der_waals_radius()
}
/// Check whether two atoms are likely bonded based on the sum of their covalent radii.
///
/// Uses a tolerance factor of 1.3 to account for slightly stretched bonds.
#[allow(dead_code)]
pub fn are_likely_bonded(a: &Atom, b: &Atom) -> bool {
    let dx = a.position[0] - b.position[0];
    let dy = a.position[1] - b.position[1];
    let dz = a.position[2] - b.position[2];
    let dist_nm = (dx * dx + dy * dy + dz * dz).sqrt();
    let r_sum_nm = (a.element.covalent_radius() + b.element.covalent_radius()) * 0.1;
    dist_nm <= r_sum_nm * 1.3
}
/// Compute the distance between two atoms in nm.
#[allow(dead_code)]
pub fn atom_distance(a: &Atom, b: &Atom) -> f64 {
    let dx = a.position[0] - b.position[0];
    let dy = a.position[1] - b.position[1];
    let dz = a.position[2] - b.position[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
/// Compute the bond angle at atom `b` formed by atoms `a`-`b`-`c` (in radians).
///
/// Returns 0.0 if any distance is zero.
#[allow(dead_code)]
pub fn bond_angle(a: &Atom, b: &Atom, c: &Atom) -> f64 {
    let u = [
        a.position[0] - b.position[0],
        a.position[1] - b.position[1],
        a.position[2] - b.position[2],
    ];
    let v = [
        c.position[0] - b.position[0],
        c.position[1] - b.position[1],
        c.position[2] - b.position[2],
    ];
    let u_norm = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
    let v_norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if u_norm < 1e-15 || v_norm < 1e-15 {
        return 0.0;
    }
    let cos_theta = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (u_norm * v_norm);
    cos_theta.clamp(-1.0, 1.0).acos()
}
/// Compute the dihedral angle for four atoms a-b-c-d (in radians).
///
/// Uses the Praxedes-Hückel sign convention.
/// Returns a value in (-π, π].
#[allow(dead_code)]
pub fn dihedral_angle(a: &Atom, b: &Atom, c: &Atom, d: &Atom) -> f64 {
    let b1 = [
        b.position[0] - a.position[0],
        b.position[1] - a.position[1],
        b.position[2] - a.position[2],
    ];
    let b2 = [
        c.position[0] - b.position[0],
        c.position[1] - b.position[1],
        c.position[2] - b.position[2],
    ];
    let b3 = [
        d.position[0] - c.position[0],
        d.position[1] - c.position[1],
        d.position[2] - c.position[2],
    ];
    let n1 = cross3(b1, b2);
    let n2 = cross3(b2, b3);
    let m1 = cross3(n1, b2);
    let n1_n = norm3(n1);
    let n2_n = norm3(n2);
    let m1_n = norm3(m1);
    if n1_n < 1e-15 || n2_n < 1e-15 || m1_n < 1e-15 {
        return 0.0;
    }
    let x = dot3(n1, n2) / (n1_n * n2_n);
    let y = dot3(m1, n2) / (m1_n * n2_n);
    x.clamp(-1.0, 1.0).acos().copysign(y)
}
/// Sort atoms by their element atomic number (ascending).
///
/// Returns a vector of indices into the original atom slice in sorted order.
#[allow(dead_code)]
pub fn sort_atoms_by_atomic_number(atoms: &[Atom]) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..atoms.len()).collect();
    indices.sort_by_key(|&i| atoms[i].element.atomic_number());
    indices
}
/// Sort atoms by residue ID.
#[allow(dead_code)]
pub fn sort_atoms_by_residue(atoms: &[Atom]) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..atoms.len()).collect();
    indices.sort_by_key(|&i| atoms[i].residue_id);
    indices
}
/// Compute the geometric centroid of a set of atom positions (in nm).
#[allow(dead_code)]
pub fn geometric_centroid(atoms: &[Atom]) -> [f64; 3] {
    if atoms.is_empty() {
        return [0.0; 3];
    }
    let n = atoms.len() as f64;
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    for a in atoms {
        cx += a.position[0];
        cy += a.position[1];
        cz += a.position[2];
    }
    [cx / n, cy / n, cz / n]
}
/// Check if two atoms are within van der Waals contact distance.
///
/// Contact is detected when the interatomic distance is less than the sum
/// of their van der Waals radii (converted from Å to nm).
#[allow(dead_code)]
pub fn in_vdw_contact(a: &Atom, b: &Atom) -> bool {
    let dist = atom_distance(a, b);
    let r_sum_nm = (a.element.van_der_waals_radius() + b.element.van_der_waals_radius()) * 0.1;
    dist < r_sum_nm
}
/// Return all pairs of atom indices that are within a cutoff distance (nm).
#[allow(dead_code)]
pub fn neighbors_within(atoms: &[Atom], cutoff_nm: f64) -> Vec<(usize, usize)> {
    let n = atoms.len();
    let mut pairs = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if atom_distance(&atoms[i], &atoms[j]) <= cutoff_nm {
                pairs.push((i, j));
            }
        }
    }
    pairs
}
/// Count atoms of each element in a slice and return (element, count) pairs.
#[allow(dead_code)]
pub fn element_composition(atoms: &[Atom]) -> std::collections::HashMap<u32, usize> {
    let mut map = std::collections::HashMap::new();
    for a in atoms {
        *map.entry(a.element.atomic_number()).or_insert(0) += 1;
    }
    map
}
/// Compute the coordination number of atom `idx` within `cutoff` nm.
///
/// Returns the count of other atoms whose distance to atom `idx` is <= cutoff.
#[allow(dead_code)]
pub fn coordination_number(atoms: &[Atom], idx: usize, cutoff_nm: f64) -> usize {
    let ref_pos = atoms[idx].position;
    atoms
        .iter()
        .enumerate()
        .filter(|(i, a)| {
            if *i == idx {
                return false;
            }
            let dx = a.position[0] - ref_pos[0];
            let dy = a.position[1] - ref_pos[1];
            let dz = a.position[2] - ref_pos[2];
            dx * dx + dy * dy + dz * dz <= cutoff_nm * cutoff_nm
        })
        .count()
}
/// Compute the coordination number for every atom in the slice.
#[allow(dead_code)]
pub fn all_coordination_numbers(atoms: &[Atom], cutoff_nm: f64) -> Vec<usize> {
    (0..atoms.len())
        .map(|i| coordination_number(atoms, i, cutoff_nm))
        .collect()
}
/// Build a list of non-bonded atom pairs, respecting an exclusion policy.
///
/// Returns pairs (i, j) with i < j that are NOT excluded by the bond graph.
#[allow(dead_code)]
pub fn non_bonded_pairs(
    n_atoms: usize,
    bonds: &BondGraph,
    policy: ExclusionPolicy,
) -> Vec<(usize, usize)> {
    let mut excluded: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    for b in &bonds.bonds {
        let (i, j) = (b.atom_i.min(b.atom_j), b.atom_i.max(b.atom_j));
        excluded.insert((i, j));
        if policy == ExclusionPolicy::OneTwoThree || policy == ExclusionPolicy::OneTwoThreeFour {
            for b2 in &bonds.bonds {
                if b2.atom_i == b.atom_i || b2.atom_j == b.atom_i {
                    let k = if b2.atom_i == b.atom_i {
                        b2.atom_j
                    } else {
                        b2.atom_i
                    };
                    if k != b.atom_j {
                        let p = (k.min(b.atom_j), k.max(b.atom_j));
                        excluded.insert(p);
                    }
                }
                if b2.atom_i == b.atom_j || b2.atom_j == b.atom_j {
                    let k = if b2.atom_i == b.atom_j {
                        b2.atom_j
                    } else {
                        b2.atom_i
                    };
                    if k != b.atom_i {
                        let p = (k.min(b.atom_i), k.max(b.atom_i));
                        excluded.insert(p);
                    }
                }
            }
        }
    }
    let mut pairs = Vec::new();
    for i in 0..n_atoms {
        for j in (i + 1)..n_atoms {
            if !excluded.contains(&(i, j)) {
                pairs.push((i, j));
            }
        }
    }
    pairs
}
/// Count bonded neighbors for each atom in the bond graph.
#[allow(dead_code)]
pub fn bond_degree(bonds: &BondGraph, n_atoms: usize) -> Vec<usize> {
    let mut degree = vec![0usize; n_atoms];
    for b in &bonds.bonds {
        if b.atom_i < n_atoms {
            degree[b.atom_i] += 1;
        }
        if b.atom_j < n_atoms {
            degree[b.atom_j] += 1;
        }
    }
    degree
}
/// Compute bond angle (radians) at position `b` given positions `a`, `b`, `c`.
#[allow(dead_code)]
pub fn bond_angle_pos(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let u = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let v = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
    let un = norm3(u);
    let vn = norm3(v);
    if un < 1e-15 || vn < 1e-15 {
        return 0.0;
    }
    let cos = (dot3(u, v) / (un * vn)).clamp(-1.0, 1.0);
    cos.acos()
}
/// Compute the dihedral angle (radians) for four positions a-b-c-d.
#[allow(dead_code)]
pub fn dihedral_angle_pos(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let b1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let b2 = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
    let b3 = [d[0] - c[0], d[1] - c[1], d[2] - c[2]];
    let n1 = cross3(b1, b2);
    let n2 = cross3(b2, b3);
    let m1 = cross3(n1, b2);
    let n1n = norm3(n1);
    let n2n = norm3(n2);
    let m1n = norm3(m1);
    if n1n < 1e-15 || n2n < 1e-15 || m1n < 1e-15 {
        return 0.0;
    }
    let x = (dot3(n1, n2) / (n1n * n2n)).clamp(-1.0, 1.0);
    let y = dot3(m1, n2) / (m1n * n2n);
    x.acos().copysign(y)
}
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn norm3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    use oxiphysics_core::Vec3;
    #[test]
    fn test_element_mass_hydrogen() {
        let m = Element::H.mass();
        assert!((m - 1.008).abs() < 1e-6, "H mass should be ~1.008, got {m}");
    }
    #[test]
    fn test_element_mass_carbon() {
        let m = Element::C.mass();
        assert!(
            (m - 12.011).abs() < 1e-6,
            "C mass should be ~12.011, got {m}"
        );
    }
    #[test]
    fn test_element_vdw_oxygen() {
        let r = Element::O.van_der_waals_radius();
        assert!(
            (r - 1.52).abs() < 1e-9,
            "O vdW radius should be 1.52 Å, got {r}"
        );
    }
    #[test]
    fn test_element_covalent_radius_nitrogen() {
        let r = Element::N.covalent_radius();
        assert!(
            r > 0.0 && r < 2.0,
            "N covalent radius should be reasonable, got {r}"
        );
    }
    #[test]
    fn test_element_electronegativity_fluorine() {
        let en_f = Element::F.electronegativity();
        let en_k = Element::K.electronegativity();
        assert!(en_f > en_k, "F should be more electronegative than K");
    }
    #[test]
    fn test_element_atomic_number() {
        assert_eq!(Element::H.atomic_number(), 1);
        assert_eq!(Element::C.atomic_number(), 6);
        assert_eq!(Element::O.atomic_number(), 8);
        assert_eq!(Element::Fe.atomic_number(), 26);
    }
    #[test]
    fn test_atom_new_defaults() {
        let atom = Atom::new(Element::C, [1.0, 2.0, 3.0], 1);
        assert_eq!(atom.element, Element::C);
        assert_eq!(atom.residue_id, 1);
        assert!((atom.mass - 12.011).abs() < 1e-6);
        assert!(atom.velocity.iter().all(|&v| v == 0.0));
    }
    #[test]
    fn test_atom_kinetic_energy() {
        let mut atom = Atom::new(Element::H, [0.0; 3], 0);
        atom.velocity = [1.0, 0.0, 0.0];
        atom.mass = 2.0;
        assert!((atom.kinetic_energy() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_atom_set_add_and_len() {
        let mut atoms = AtomSet::new();
        assert!(atoms.is_empty());
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        assert_eq!(atoms.len(), 1);
    }
    #[test]
    fn test_atom_set_remove_atom() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(0.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 2.0, 0.0, 1);
        atoms.add_atom(Vec3::new(2.0, 0.0, 0.0), Vec3::zeros(), 3.0, 0.0, 2);
        atoms.remove_atom(0);
        assert_eq!(atoms.len(), 2);
    }
    #[test]
    fn test_atom_set_clear() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::zeros(), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.clear();
        assert!(atoms.is_empty());
    }
    #[test]
    fn test_kinetic_energy() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0), 2.0, 0.0, 0);
        let ke = atoms.kinetic_energy();
        assert!((ke - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_temperature_equipartition() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0), 1.0, 0.0, 0);
        let temp = atoms.temperature(1.0);
        assert!((temp - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_center_of_mass() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(0.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::new(2.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        let com = atoms.center_of_mass();
        assert!((com.x - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_total_momentum() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0), 2.0, 0.0, 0);
        atoms.add_atom(Vec3::zeros(), Vec3::new(-1.0, 0.0, 0.0), 2.0, 0.0, 0);
        let p = atoms.total_momentum();
        assert!(p.norm() < 1e-12);
    }
    #[test]
    fn test_radius_of_gyration_two_atoms() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(-1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        let rg = atoms.radius_of_gyration();
        assert!((rg - 1.0).abs() < 1e-12, "Rg should be 1.0, got {rg}");
    }
    #[test]
    fn test_radius_of_gyration_empty() {
        let atoms = AtomSet::new();
        assert_eq!(atoms.radius_of_gyration(), 0.0);
    }
    #[test]
    fn test_apply_pbc_wraps_positions() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(12.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.apply_pbc([10.0, 10.0, 10.0]);
        let x = atoms.positions[0].x;
        assert!(
            (x - 2.0).abs() < 1e-10,
            "PBC: x=12 in box=10 should wrap to 2, got {x}"
        );
    }
    #[test]
    fn test_apply_pbc_already_inside() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(5.0, 3.0, 1.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.apply_pbc([10.0, 10.0, 10.0]);
        let pos = atoms.positions[0];
        assert!((pos.x - 5.0).abs() < 1e-12);
        assert!((pos.y - 3.0).abs() < 1e-12);
        assert!((pos.z - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_bond_graph_add_and_query() {
        let mut bg = BondGraph::new(3);
        bg.add_bond(0, 1, 1.0);
        bg.add_bond(1, 2, 2.0);
        assert_eq!(bg.num_bonds(), 2);
        assert!(bg.are_bonded(0, 1));
        assert!(!bg.are_bonded(0, 2));
        assert!((bg.bond_order(1, 2) - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_bond_graph_bonds_of() {
        let mut bg = BondGraph::new(4);
        bg.add_bond(0, 1, 1.0);
        bg.add_bond(0, 2, 1.5);
        bg.add_bond(1, 3, 1.0);
        let bonds_of_0 = bg.bonds_of(0);
        assert_eq!(bonds_of_0.len(), 2, "Atom 0 should have 2 bonds");
    }
    #[test]
    fn test_bond_graph_clear() {
        let mut bg = BondGraph::new(3);
        bg.add_bond(0, 1, 1.0);
        bg.add_bond(1, 2, 1.0);
        bg.clear();
        assert_eq!(bg.num_bonds(), 0);
    }
    #[test]
    fn test_atom_distance() {
        let a = Atom::new(Element::C, [0.0, 0.0, 0.0], 0);
        let b = Atom::new(Element::C, [0.3, 0.4, 0.0], 0);
        let d = atom_distance(&a, &b);
        assert!((d - 0.5).abs() < 1e-12, "expected 0.5 nm, got {d}");
    }
    #[test]
    fn test_bond_angle_linear() {
        let a = Atom::new(Element::C, [-0.1, 0.0, 0.0], 0);
        let b = Atom::new(Element::C, [0.0, 0.0, 0.0], 0);
        let c = Atom::new(Element::C, [0.1, 0.0, 0.0], 0);
        let angle = bond_angle(&a, &b, &c);
        assert!(
            (angle - std::f64::consts::PI).abs() < 1e-10,
            "expected π, got {angle}"
        );
    }
    #[test]
    fn test_bond_angle_right_angle() {
        let a = Atom::new(Element::C, [1.0, 0.0, 0.0], 0);
        let b = Atom::new(Element::C, [0.0, 0.0, 0.0], 0);
        let c = Atom::new(Element::C, [0.0, 1.0, 0.0], 0);
        let angle = bond_angle(&a, &b, &c);
        assert!(
            (angle - std::f64::consts::PI / 2.0).abs() < 1e-10,
            "expected π/2, got {angle}"
        );
    }
    #[test]
    fn test_dihedral_angle_zero() {
        let a = Atom::new(Element::C, [0.0, 1.0, 0.0], 0);
        let b = Atom::new(Element::C, [0.0, 0.0, 0.0], 0);
        let c = Atom::new(Element::C, [1.0, 0.0, 0.0], 0);
        let d = Atom::new(Element::C, [1.0, 1.0, 0.0], 0);
        let phi = dihedral_angle(&a, &b, &c, &d);
        assert!(phi.abs() < 1e-10, "expected 0 dihedral, got {phi}");
    }
    #[test]
    fn test_dihedral_angle_90deg() {
        let a = Atom::new(Element::C, [0.0, 0.0, 1.0], 0);
        let b = Atom::new(Element::C, [0.0, 0.0, 0.0], 0);
        let c = Atom::new(Element::C, [1.0, 0.0, 0.0], 0);
        let d = Atom::new(Element::C, [1.0, 1.0, 0.0], 0);
        let phi = dihedral_angle(&a, &b, &c, &d);
        assert!(
            (phi.abs() - std::f64::consts::PI / 2.0).abs() < 1e-8,
            "expected ±π/2, got {phi}"
        );
    }
    #[test]
    fn test_sort_atoms_by_atomic_number() {
        let atoms = vec![
            Atom::new(Element::Fe, [0.0; 3], 0),
            Atom::new(Element::H, [0.0; 3], 0),
            Atom::new(Element::C, [0.0; 3], 0),
        ];
        let order = sort_atoms_by_atomic_number(&atoms);
        assert_eq!(order[0], 1, "H should be first (Z=1)");
        assert_eq!(order[1], 2, "C should be second (Z=6)");
        assert_eq!(order[2], 0, "Fe should be last (Z=26)");
    }
    #[test]
    fn test_sort_atoms_by_residue() {
        let atoms = vec![
            Atom::new(Element::C, [0.0; 3], 5),
            Atom::new(Element::N, [0.0; 3], 1),
            Atom::new(Element::O, [0.0; 3], 3),
        ];
        let order = sort_atoms_by_residue(&atoms);
        assert_eq!(order[0], 1, "residue 1 first");
        assert_eq!(order[1], 2, "residue 3 second");
        assert_eq!(order[2], 0, "residue 5 last");
    }
    #[test]
    fn test_geometric_centroid() {
        let atoms = vec![
            Atom::new(Element::C, [-1.0, 0.0, 0.0], 0),
            Atom::new(Element::C, [1.0, 0.0, 0.0], 0),
        ];
        let cen = geometric_centroid(&atoms);
        assert!(
            cen[0].abs() < 1e-12,
            "centroid x should be 0, got {}",
            cen[0]
        );
    }
    #[test]
    fn test_in_vdw_contact() {
        let a = Atom::new(Element::H, [0.0, 0.0, 0.0], 0);
        let b = Atom::new(Element::H, [0.18, 0.0, 0.0], 0);
        assert!(
            in_vdw_contact(&a, &b),
            "H atoms at 0.18 nm should be in vdW contact (sum=0.24)"
        );
    }
    #[test]
    fn test_not_in_vdw_contact() {
        let a = Atom::new(Element::H, [0.0, 0.0, 0.0], 0);
        let b = Atom::new(Element::H, [0.5, 0.0, 0.0], 0);
        assert!(
            !in_vdw_contact(&a, &b),
            "H atoms at 0.5 nm should NOT be in vdW contact"
        );
    }
    #[test]
    fn test_are_likely_bonded() {
        let c_atom = Atom::new(Element::C, [0.0, 0.0, 0.0], 0);
        let h_atom = Atom::new(Element::H, [0.11, 0.0, 0.0], 0);
        assert!(
            are_likely_bonded(&c_atom, &h_atom),
            "C-H at 0.11 nm should be bonded"
        );
    }
    #[test]
    fn test_neighbors_within() {
        let atoms = vec![
            Atom::new(Element::C, [0.0, 0.0, 0.0], 0),
            Atom::new(Element::C, [0.2, 0.0, 0.0], 0),
            Atom::new(Element::C, [0.6, 0.0, 0.0], 0),
        ];
        let pairs = neighbors_within(&atoms, 0.3);
        assert_eq!(pairs.len(), 1, "expected 1 pair, got {}", pairs.len());
        assert_eq!(pairs[0], (0, 1));
    }
    #[test]
    fn test_element_composition() {
        let atoms = vec![
            Atom::new(Element::C, [0.0; 3], 0),
            Atom::new(Element::C, [0.0; 3], 0),
            Atom::new(Element::H, [0.0; 3], 0),
        ];
        let comp = element_composition(&atoms);
        assert_eq!(*comp.get(&6).unwrap(), 2, "2 carbons");
        assert_eq!(*comp.get(&1).unwrap(), 1, "1 hydrogen");
    }
    #[test]
    fn test_contact_map_basic() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(0.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::new(0.2, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        let map = atoms.contact_map(0.3);
        assert!(
            map[0][1],
            "atoms 0 and 1 (0.2 nm) should be in contact at cutoff 0.3"
        );
        assert!(
            !map[0][2],
            "atoms 0 and 2 (1.0 nm) should NOT be in contact at cutoff 0.3"
        );
        assert_eq!(map[0][1], map[1][0]);
    }
    #[test]
    fn test_rescale_velocities() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0), 1.0, 0.0, 0);
        let kb = 1.0;
        let initial_temp = atoms.temperature(kb);
        let target = initial_temp * 4.0;
        atoms.rescale_velocities(target, kb);
        let new_temp = atoms.temperature(kb);
        assert!(
            (new_temp - target).abs() < 1e-10,
            "Temperature after rescaling should be {target}, got {new_temp}"
        );
    }
    #[test]
    fn test_rms_velocity() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::new(3.0, 4.0, 0.0), 1.0, 0.0, 0);
        let rms = atoms.rms_velocity();
        assert!(
            (rms - 5.0).abs() < 1e-12,
            "RMS velocity should be 5.0, got {rms}"
        );
    }
    #[test]
    fn test_translate_atomset() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 2.0, 3.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.translate([0.5, -1.0, 2.0]);
        let p = atoms.positions[0];
        assert!((p.x - 1.5).abs() < 1e-12);
        assert!((p.y - 1.0).abs() < 1e-12);
        assert!((p.z - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_atoms_in_sphere() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(0.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::new(0.5, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::new(2.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        let inside = atoms.atoms_in_sphere([0.0, 0.0, 0.0], 1.0);
        assert_eq!(
            inside.len(),
            2,
            "atoms 0 and 1 should be in sphere of radius 1"
        );
    }
    #[test]
    fn test_residue_add_atoms() {
        let mut res = Residue::new(1, "ALA", 'A');
        res.add_atom(0);
        res.add_atom(1);
        assert_eq!(res.n_atoms(), 2);
        assert_eq!(res.id, 1);
    }
    #[test]
    fn test_chain_add_residues() {
        let mut chain = Chain::new('B');
        chain.add_residue(0);
        chain.add_residue(1);
        assert_eq!(chain.n_residues(), 2);
        assert_eq!(chain.id, 'B');
    }
    #[test]
    fn test_topology_add_and_query() {
        let mut top = Topology::new();
        let a0 = Atom::new(Element::C, [0.0; 3], 1);
        let a1 = Atom::new(Element::N, [0.1; 3], 1);
        top.add_atom(a0);
        top.add_atom(a1);
        let r0 = Residue::new(1, "GLY", 'A');
        top.add_residue(r0);
        assert_eq!(top.n_atoms(), 2);
        assert_eq!(top.n_residues(), 1);
    }
    #[test]
    fn test_topology_atoms_in_residue() {
        let mut top = Topology::new();
        top.add_atom(Atom::new(Element::C, [0.0; 3], 5));
        top.add_atom(Atom::new(Element::N, [0.1; 3], 5));
        top.add_atom(Atom::new(Element::O, [0.2; 3], 7));
        let in_5 = top.atoms_in_residue(5);
        assert_eq!(in_5.len(), 2);
    }
    #[test]
    fn test_topology_chains_and_residues() {
        let mut top = Topology::new();
        let mut c = Chain::new('A');
        c.add_residue(0);
        top.add_chain(c);
        let r = Residue::new(1, "ALA", 'A');
        top.add_residue(r);
        let chains = top.residues_in_chain('A');
        assert_eq!(chains.len(), 1);
    }
    #[test]
    fn test_atom_group_len_and_empty() {
        let g = AtomGroup::new("empty", vec![]);
        assert!(g.is_empty());
        let g2 = AtomGroup::new("two", vec![0, 1]);
        assert_eq!(g2.len(), 2);
    }
    #[test]
    fn test_atom_group_union() {
        let g1 = AtomGroup::new("g1", vec![0, 1, 2]);
        let g2 = AtomGroup::new("g2", vec![2, 3, 4]);
        let u = g1.union(&g2);
        assert_eq!(u.len(), 5);
        assert!(u.indices.contains(&0));
        assert!(u.indices.contains(&4));
    }
    #[test]
    fn test_atom_group_intersection() {
        let g1 = AtomGroup::new("g1", vec![0, 1, 2, 3]);
        let g2 = AtomGroup::new("g2", vec![2, 3, 4]);
        let inter = g1.intersection(&g2);
        assert_eq!(inter.len(), 2);
        assert!(inter.indices.contains(&2));
        assert!(inter.indices.contains(&3));
    }
    #[test]
    fn test_atom_group_difference() {
        let g1 = AtomGroup::new("g1", vec![0, 1, 2, 3]);
        let g2 = AtomGroup::new("g2", vec![2, 3]);
        let diff = g1.difference(&g2);
        assert_eq!(diff.len(), 2);
        assert!(diff.indices.contains(&0));
        assert!(diff.indices.contains(&1));
    }
    #[test]
    fn test_atom_group_centroid() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(0.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::new(2.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::new(10.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        let g = AtomGroup::new("g", vec![0, 1]);
        let c = g.centroid(&atoms);
        assert!((c.x - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_atom_group_com() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(0.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::new(4.0, 0.0, 0.0), Vec3::zeros(), 3.0, 0.0, 0);
        let g = AtomGroup::new("g", vec![0, 1]);
        let com = g.center_of_mass(&atoms);
        assert!((com.x - 3.0).abs() < 1e-12, "weighted COM: {}", com.x);
    }
    #[test]
    fn test_atom_group_total_mass() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::zeros(), Vec3::zeros(), 4.0, 0.0, 0);
        let g = AtomGroup::new("g", vec![0, 1]);
        assert!((g.total_mass(&atoms) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_selection_by_element() {
        let atoms = vec![
            Atom::new(Element::C, [0.0; 3], 1),
            Atom::new(Element::N, [0.0; 3], 1),
            Atom::new(Element::C, [0.0; 3], 2),
        ];
        let g = AtomSelection::by_element(&atoms, Element::C);
        assert_eq!(g.len(), 2);
    }
    #[test]
    fn test_selection_by_residue() {
        let atoms = vec![
            Atom::new(Element::C, [0.0; 3], 1),
            Atom::new(Element::N, [0.0; 3], 2),
            Atom::new(Element::O, [0.0; 3], 1),
        ];
        let g = AtomSelection::by_residue(&atoms, 1);
        assert_eq!(g.len(), 2);
    }
    #[test]
    fn test_selection_heavy_atoms() {
        let atoms = vec![
            Atom::new(Element::H, [0.0; 3], 1),
            Atom::new(Element::C, [0.0; 3], 1),
            Atom::new(Element::H, [0.0; 3], 1),
            Atom::new(Element::N, [0.0; 3], 1),
        ];
        let g = AtomSelection::heavy_atoms(&atoms);
        assert_eq!(g.len(), 2);
    }
    #[test]
    fn test_selection_within_sphere() {
        let atoms = vec![
            Atom::new(Element::C, [0.0, 0.0, 0.0], 1),
            Atom::new(Element::C, [0.5, 0.0, 0.0], 1),
            Atom::new(Element::C, [5.0, 0.0, 0.0], 1),
        ];
        let g = AtomSelection::within_sphere(&atoms, [0.0, 0.0, 0.0], 1.0);
        assert_eq!(g.len(), 2);
    }
    #[test]
    fn test_selection_by_charge_gt() {
        let mut a0 = Atom::new(Element::Na, [0.0; 3], 1);
        a0.charge = 1.0;
        let mut a1 = Atom::new(Element::Cl, [0.0; 3], 1);
        a1.charge = -1.0;
        let a2 = Atom::new(Element::C, [0.0; 3], 1);
        let atoms = vec![a0, a1, a2];
        let g = AtomSelection::by_charge_gt(&atoms, 0.0);
        assert_eq!(g.len(), 1);
    }
    #[test]
    fn test_coordination_number_basic() {
        let atoms = vec![
            Atom::new(Element::C, [0.0, 0.0, 0.0], 1),
            Atom::new(Element::C, [0.15, 0.0, 0.0], 1),
            Atom::new(Element::C, [0.30, 0.0, 0.0], 1),
            Atom::new(Element::C, [5.0, 0.0, 0.0], 1),
        ];
        let cn = coordination_number(&atoms, 0, 0.35);
        assert_eq!(cn, 2, "atom 0 should have 2 neighbors within 0.35 nm");
    }
    #[test]
    fn test_coordination_number_isolated() {
        let atoms = vec![
            Atom::new(Element::C, [0.0, 0.0, 0.0], 1),
            Atom::new(Element::C, [10.0, 0.0, 0.0], 1),
        ];
        let cn = coordination_number(&atoms, 0, 0.5);
        assert_eq!(cn, 0);
    }
    #[test]
    fn test_all_coordination_numbers_len() {
        let atoms = vec![
            Atom::new(Element::C, [0.0; 3], 1),
            Atom::new(Element::C, [0.1; 3], 1),
            Atom::new(Element::C, [0.2; 3], 1),
        ];
        let cns = all_coordination_numbers(&atoms, 0.3);
        assert_eq!(cns.len(), 3);
    }
    #[test]
    fn test_non_bonded_pairs_no_bonds() {
        let bg = BondGraph::new(4);
        let pairs = non_bonded_pairs(4, &bg, ExclusionPolicy::OneTwoOnly);
        assert_eq!(pairs.len(), 6);
    }
    #[test]
    fn test_non_bonded_pairs_excludes_12() {
        let mut bg = BondGraph::new(4);
        bg.add_bond(0, 1, 1.0);
        let pairs = non_bonded_pairs(4, &bg, ExclusionPolicy::OneTwoOnly);
        assert_eq!(pairs.len(), 5);
        assert!(!pairs.contains(&(0, 1)));
    }
    #[test]
    fn test_bond_degree_basic() {
        let mut bg = BondGraph::new(3);
        bg.add_bond(0, 1, 1.0);
        bg.add_bond(0, 2, 1.0);
        let deg = bond_degree(&bg, 3);
        assert_eq!(deg[0], 2);
        assert_eq!(deg[1], 1);
        assert_eq!(deg[2], 1);
    }
    #[test]
    fn test_bounding_box() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, -1.0, 2.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::new(-1.0, 3.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        let (lo, hi) = atoms.bounding_box();
        assert!((lo[0] - (-1.0)).abs() < 1e-12);
        assert!((hi[0] - 1.0).abs() < 1e-12);
        assert!((lo[1] - (-1.0)).abs() < 1e-12);
        assert!((hi[1] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_max_pairwise_distance() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(0.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::new(3.0, 4.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        let d = atoms.max_pairwise_distance();
        assert!((d - 5.0).abs() < 1e-12, "expected 5.0, got {d}");
    }
    #[test]
    fn test_inertia_tensor_symmetry() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::new(-1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        let inertia = atoms.inertia_tensor();
        assert!((inertia[1] - inertia[3]).abs() < 1e-12);
        assert!(inertia[0].abs() < 1e-12, "I_xx = {} expected 0", inertia[0]);
    }
    #[test]
    fn test_rmsd_identical() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 2.0, 3.0), Vec3::zeros(), 1.0, 0.0, 0);
        let ref_pos = atoms.positions.clone();
        let rmsd = atoms.rmsd(&ref_pos);
        assert!(rmsd.abs() < 1e-12, "RMSD to self should be 0, got {rmsd}");
    }
    #[test]
    fn test_rmsd_displaced() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(0.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        let ref_pos = vec![Vec3::new(1.0, 0.0, 0.0)];
        let rmsd = atoms.rmsd(&ref_pos);
        assert!((rmsd - 1.0).abs() < 1e-12, "RMSD should be 1.0, got {rmsd}");
    }
    #[test]
    fn test_scale_positions() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 2.0, 3.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.scale_positions(2.0);
        let p = atoms.positions[0];
        assert!((p.x - 2.0).abs() < 1e-12);
        assert!((p.y - 4.0).abs() < 1e-12);
        assert!((p.z - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_dipole_moment() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 1.0, 0);
        atoms.add_atom(Vec3::new(-1.0, 0.0, 0.0), Vec3::zeros(), 1.0, -1.0, 0);
        let d = atoms.dipole_moment();
        assert!((d.x - 2.0).abs() < 1e-12, "dipole x = {}", d.x);
    }
    #[test]
    fn test_fastest_atom() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::zeros(), Vec3::new(5.0, 0.0, 0.0), 1.0, 0.0, 0);
        let idx = atoms.fastest_atom().unwrap();
        assert_eq!(idx, 1);
    }
    #[test]
    fn test_heaviest_atom() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::zeros(), Vec3::zeros(), 100.0, 0.0, 0);
        let idx = atoms.heaviest_atom().unwrap();
        assert_eq!(idx, 1);
    }
    #[test]
    fn test_charged_atoms() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::zeros(), 1.0, 1.0, 0);
        atoms.add_atom(Vec3::zeros(), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::zeros(), Vec3::zeros(), 1.0, -1.0, 0);
        let charged = atoms.charged_atoms();
        assert_eq!(charged.len(), 2);
    }
    #[test]
    fn test_mirror_z() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 2.0, 3.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.mirror_z();
        assert!((atoms.positions[0].z - (-3.0)).abs() < 1e-12);
    }
    #[test]
    fn test_rotate_identity() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 2.0, 3.0), Vec3::zeros(), 1.0, 0.0, 0);
        let rot = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        atoms.rotate(rot);
        let p = atoms.positions[0];
        assert!((p.x - 1.0).abs() < 1e-12);
        assert!((p.y - 2.0).abs() < 1e-12);
        assert!((p.z - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_bond_angle_pos_right_angle() {
        let angle = bond_angle_pos([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((angle - std::f64::consts::FRAC_PI_2).abs() < 1e-10);
    }
    #[test]
    fn test_dihedral_angle_pos_zero() {
        let phi = dihedral_angle_pos(
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
        );
        assert!(phi.abs() < 1e-10, "expected 0 dihedral, got {phi}");
    }
    #[test]
    fn test_total_force_norm_zero_forces() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::zeros(), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::zeros(), Vec3::zeros(), 1.0, 0.0, 0);
        assert!(atoms.total_force_norm() < 1e-15);
    }
    #[test]
    fn test_max_force_zero() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::zeros(), 1.0, 0.0, 0);
        assert!(atoms.max_force() < 1e-15);
    }
    #[test]
    fn test_remove_com_velocity_zeroes_momentum() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::zeros(), Vec3::new(3.0, 0.0, 0.0), 1.0, 0.0, 0);
        atoms.remove_com_velocity();
        let p = atoms.total_momentum();
        assert!(p.norm() < 1e-12, "momentum after removal: {}", p.norm());
    }
    #[test]
    fn test_bounding_box_empty() {
        let atoms = AtomSet::new();
        let (lo, hi) = atoms.bounding_box();
        for k in 0..3 {
            assert_eq!(lo[k], 0.0);
            assert_eq!(hi[k], 0.0);
        }
    }
}
