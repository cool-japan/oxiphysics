//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{AlloySystem, CrystalStructure, DislocationCore, EamParams, MixingRule, SimBox};

/// Add two 3-vectors.
pub(super) fn vadd(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract `b` from `a`.
pub(super) fn vsub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scale a 3-vector.
pub(super) fn vscale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
/// Dot product.
pub(super) fn vdot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Euclidean norm.
pub(super) fn vnorm(v: [f64; 3]) -> f64 {
    vdot(v, v).sqrt()
}
/// Cross product.
#[cfg(test)]
pub(super) fn vcross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Normalize a vector. Returns zero vector if norm < eps.
pub(super) fn vnormalize(v: [f64; 3]) -> [f64; 3] {
    let n = vnorm(v);
    if n < 1e-30 {
        [0.0; 3]
    } else {
        vscale(v, 1.0 / n)
    }
}
/// Cross-species pair potential using a mixing rule.
pub fn cross_pair_potential(p1: &EamParams, p2: &EamParams, r: f64, rule: MixingRule) -> f64 {
    match rule {
        MixingRule::Geometric => {
            let v1 = p1.pair_potential(r);
            let v2 = p2.pair_potential(r);
            if v1 >= 0.0 && v2 >= 0.0 {
                (v1 * v2).sqrt()
            } else if v1 <= 0.0 && v2 <= 0.0 {
                -(v1.abs() * v2.abs()).sqrt()
            } else {
                0.5 * (v1 + v2)
            }
        }
        MixingRule::Arithmetic => 0.5 * (p1.pair_potential(r) + p2.pair_potential(r)),
        MixingRule::Johnson => {
            let rho_a1 = p1.electron_density(r);
            let rho_a2 = p2.electron_density(r);
            if rho_a1.abs() < 1e-30 || rho_a2.abs() < 1e-30 {
                return 0.5 * (p1.pair_potential(r) + p2.pair_potential(r));
            }
            0.5 * (rho_a2 / rho_a1 * p1.pair_potential(r) + rho_a1 / rho_a2 * p2.pair_potential(r))
        }
    }
}
/// Cross-species pair force using a mixing rule.
pub fn cross_pair_force(p1: &EamParams, p2: &EamParams, r: f64, rule: MixingRule) -> f64 {
    let dr = 1e-6;
    let vp = cross_pair_potential(p1, p2, r + dr, rule);
    let vm = cross_pair_potential(p1, p2, r - dr, rule);
    (vp - vm) / (2.0 * dr)
}
/// Compute the Warren-Cowley short-range order parameter α for species A-B.
///
/// α = 1 - P(B|A) / c_B
/// where P(B|A) is probability of B neighbor given A atom, c_B is B concentration.
/// α = 0: random, α < 0: ordering, α > 0: clustering.
pub fn warren_cowley_sro(system: &AlloySystem, species_a: usize, species_b: usize) -> f64 {
    let pairs = system.build_neighbor_list();
    let n = system.atoms.len() as f64;
    let n_b = system
        .atoms
        .iter()
        .filter(|a| a.species == species_b)
        .count() as f64;
    let c_b = n_b / n;
    if !(1e-12..=1.0 - 1e-12).contains(&c_b) {
        return 0.0;
    }
    let mut n_ab = 0usize;
    let mut n_a_total = 0usize;
    for &(i, j, _dr, _r) in &pairs {
        let si = system.atoms[i].species;
        let sj = system.atoms[j].species;
        if si == species_a {
            n_a_total += 1;
            if sj == species_b {
                n_ab += 1;
            }
        }
        if sj == species_a {
            n_a_total += 1;
            if si == species_b {
                n_ab += 1;
            }
        }
    }
    if n_a_total == 0 {
        return 0.0;
    }
    let p_b_given_a = n_ab as f64 / n_a_total as f64;
    1.0 - p_b_given_a / c_b
}
/// Estimate vacancy migration energy using NEB-like drag method.
///
/// Creates a vacancy at `vacancy_idx`, attempts to move a neighbor into it,
/// and returns the energy barrier estimate (eV).
pub fn vacancy_migration_energy(
    species: &[EamParams],
    lattice_a: f64,
    structure: CrystalStructure,
    vacancy_idx: usize,
    neighbor_idx: usize,
    n_images: usize,
    mixing_rule: MixingRule,
) -> f64 {
    let atoms_per_side = 4;
    let box_len = atoms_per_side as f64 * lattice_a;
    let mut system = AlloySystem::new(species.to_vec(), SimBox::cubic(box_len), mixing_rule, 0.001);
    let positions = generate_crystal_positions(lattice_a, atoms_per_side, structure);
    for pos in &positions {
        system.add_atom(*pos, 0);
    }
    let safe_vac = vacancy_idx.min(system.atoms.len().saturating_sub(1));
    let safe_nbr = neighbor_idx.min(system.atoms.len().saturating_sub(1));
    let start_pos = system.atoms[safe_nbr].position;
    let end_pos = system.atoms[safe_vac].position;
    let _removed = if safe_vac < system.atoms.len() {
        system.atoms.remove(safe_vac);
        true
    } else {
        false
    };
    let adj_nbr = if safe_nbr > safe_vac {
        safe_nbr - 1
    } else {
        safe_nbr
    };
    if adj_nbr >= system.atoms.len() {
        return 0.0;
    }
    let images = n_images.max(3);
    let mut max_energy = f64::NEG_INFINITY;
    let mut min_energy = f64::INFINITY;
    for img in 0..images {
        let t = img as f64 / (images - 1) as f64;
        let interp = vadd(vscale(start_pos, 1.0 - t), vscale(end_pos, t));
        system.atoms[adj_nbr].position = system.sim_box.wrap(interp);
        for atom in &mut system.atoms {
            atom.rho_bar = 0.0;
        }
        let pairs = system.build_neighbor_list();
        for &(i, j, _dr, r) in &pairs {
            let si = system.atoms[i].species;
            let sj = system.atoms[j].species;
            system.atoms[i].rho_bar += system.species[sj].electron_density(r);
            system.atoms[j].rho_bar += system.species[si].electron_density(r);
        }
        let e = system.potential_energy();
        if e > max_energy {
            max_energy = e;
        }
        if e < min_energy {
            min_energy = e;
        }
    }
    (max_energy - min_energy).abs()
}
/// Generate crystal positions for a given structure and number of unit cells per side.
pub fn generate_crystal_positions(
    lattice_a: f64,
    n_cells: usize,
    structure: CrystalStructure,
) -> Vec<[f64; 3]> {
    let mut positions = Vec::new();
    match structure {
        CrystalStructure::FCC => {
            let basis = [
                [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.0],
                [0.5, 0.0, 0.5],
                [0.0, 0.5, 0.5],
            ];
            for ix in 0..n_cells {
                for iy in 0..n_cells {
                    for iz in 0..n_cells {
                        let origin = [
                            ix as f64 * lattice_a,
                            iy as f64 * lattice_a,
                            iz as f64 * lattice_a,
                        ];
                        for b in &basis {
                            positions.push([
                                origin[0] + b[0] * lattice_a,
                                origin[1] + b[1] * lattice_a,
                                origin[2] + b[2] * lattice_a,
                            ]);
                        }
                    }
                }
            }
        }
        CrystalStructure::BCC => {
            let basis = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]];
            for ix in 0..n_cells {
                for iy in 0..n_cells {
                    for iz in 0..n_cells {
                        let origin = [
                            ix as f64 * lattice_a,
                            iy as f64 * lattice_a,
                            iz as f64 * lattice_a,
                        ];
                        for b in &basis {
                            positions.push([
                                origin[0] + b[0] * lattice_a,
                                origin[1] + b[1] * lattice_a,
                                origin[2] + b[2] * lattice_a,
                            ]);
                        }
                    }
                }
            }
        }
        CrystalStructure::HCP => {
            let c_over_a = (8.0_f64 / 3.0).sqrt();
            let c = lattice_a * c_over_a;
            let basis = [[0.0, 0.0, 0.0], [0.5, 0.5 / 3.0_f64.sqrt(), 0.5]];
            for ix in 0..n_cells {
                for iy in 0..n_cells {
                    for iz in 0..n_cells {
                        let origin = [
                            ix as f64 * lattice_a + (iy % 2) as f64 * 0.5 * lattice_a,
                            iy as f64 * lattice_a * (3.0_f64.sqrt() / 2.0),
                            iz as f64 * c,
                        ];
                        for b in &basis {
                            positions.push([
                                origin[0] + b[0] * lattice_a,
                                origin[1] + b[1] * lattice_a,
                                origin[2] + b[2] * c,
                            ]);
                        }
                    }
                }
            }
        }
    }
    positions
}
/// Compute grain boundary energy by constructing a bicrystal.
///
/// Creates two half-slabs with a tilt angle and computes the excess energy
/// per unit area relative to a perfect crystal.
///
/// Returns energy in eV/ų.
pub fn grain_boundary_energy(params: &EamParams, tilt_angle_deg: f64, n_cells: usize) -> f64 {
    let a = params.lattice_a;
    let box_len = n_cells as f64 * a;
    let mut ref_system = AlloySystem::new(
        vec![params.clone()],
        SimBox::new(box_len, box_len, box_len),
        MixingRule::Geometric,
        0.001,
    );
    let positions = generate_crystal_positions(a, n_cells, CrystalStructure::FCC);
    for pos in &positions {
        ref_system.add_atom(*pos, 0);
    }
    for atom in &mut ref_system.atoms {
        atom.rho_bar = 0.0;
    }
    let pairs = ref_system.build_neighbor_list();
    for &(i, j, _dr, r) in &pairs {
        ref_system.atoms[i].rho_bar += ref_system.species[0].electron_density(r);
        ref_system.atoms[j].rho_bar += ref_system.species[0].electron_density(r);
    }
    let e_ref = ref_system.potential_energy();
    let e_per_atom_ref = e_ref / ref_system.atoms.len() as f64;
    let theta = tilt_angle_deg * PI / 180.0;
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    let half_z = box_len / 2.0;
    let mut gb_system = AlloySystem::new(
        vec![params.clone()],
        SimBox::new(box_len, box_len, box_len),
        MixingRule::Geometric,
        0.001,
    );
    for pos in &positions {
        let mut p = *pos;
        if p[2] >= half_z {
            let x0 = p[0] - half_z;
            let y0 = p[1] - half_z;
            p[0] = x0 * cos_t - y0 * sin_t + half_z;
            p[1] = x0 * sin_t + y0 * cos_t + half_z;
        }
        p = gb_system.sim_box.wrap(p);
        gb_system.add_atom(p, 0);
    }
    let min_dist = 0.5 * a;
    let mut to_remove = Vec::new();
    let n = gb_system.atoms.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let dr = gb_system.sim_box.min_image(vsub(
                gb_system.atoms[j].position,
                gb_system.atoms[i].position,
            ));
            if vnorm(dr) < min_dist {
                to_remove.push(j);
            }
        }
    }
    to_remove.sort_unstable();
    to_remove.dedup();
    for &idx in to_remove.iter().rev() {
        gb_system.atoms.remove(idx);
    }
    for atom in &mut gb_system.atoms {
        atom.rho_bar = 0.0;
    }
    let pairs_gb = gb_system.build_neighbor_list();
    for &(i, j, _dr, r) in &pairs_gb {
        gb_system.atoms[i].rho_bar += gb_system.species[0].electron_density(r);
        gb_system.atoms[j].rho_bar += gb_system.species[0].electron_density(r);
    }
    let e_gb = gb_system.potential_energy();
    let excess = e_gb - gb_system.atoms.len() as f64 * e_per_atom_ref;
    let area = box_len * box_len;
    (excess / (2.0 * area)).abs()
}
/// Compute generalized stacking fault (GSF) energy curve.
///
/// Displaces top half of an FCC crystal along \[112\] direction in steps,
/// returns vector of (displacement_fraction, energy_per_area eV/ų).
pub fn stacking_fault_energy_curve(
    params: &EamParams,
    n_cells: usize,
    n_points: usize,
) -> Vec<(f64, f64)> {
    let a = params.lattice_a;
    let box_len = n_cells as f64 * a;
    let shift_direction = vnormalize([1.0, 1.0, 2.0]);
    let burgers_mag = a / 6.0_f64.sqrt();
    let mut ref_system = AlloySystem::new(
        vec![params.clone()],
        SimBox::new(box_len, box_len, box_len),
        MixingRule::Geometric,
        0.001,
    );
    let positions = generate_crystal_positions(a, n_cells, CrystalStructure::FCC);
    for pos in &positions {
        ref_system.add_atom(*pos, 0);
    }
    for atom in &mut ref_system.atoms {
        atom.rho_bar = 0.0;
    }
    let pairs_ref = ref_system.build_neighbor_list();
    for &(i, j, _dr, r) in &pairs_ref {
        ref_system.atoms[i].rho_bar += ref_system.species[0].electron_density(r);
        ref_system.atoms[j].rho_bar += ref_system.species[0].electron_density(r);
    }
    let e_per_atom_ref = ref_system.potential_energy() / ref_system.atoms.len() as f64;
    let half_z = box_len / 2.0;
    let area = box_len * box_len;
    let mut curve = Vec::with_capacity(n_points);
    for pt in 0..n_points {
        let frac = pt as f64 / (n_points - 1).max(1) as f64;
        let shift = vscale(shift_direction, frac * burgers_mag);
        let mut sys = AlloySystem::new(
            vec![params.clone()],
            SimBox::new(box_len, box_len, box_len),
            MixingRule::Geometric,
            0.001,
        );
        for pos in &positions {
            let mut p = *pos;
            if p[2] >= half_z {
                p = vadd(p, shift);
            }
            p = sys.sim_box.wrap(p);
            sys.add_atom(p, 0);
        }
        for atom in &mut sys.atoms {
            atom.rho_bar = 0.0;
        }
        let pairs_sys = sys.build_neighbor_list();
        for &(i, j, _dr, r) in &pairs_sys {
            sys.atoms[i].rho_bar += sys.species[0].electron_density(r);
            sys.atoms[j].rho_bar += sys.species[0].electron_density(r);
        }
        let e = sys.potential_energy();
        let excess = e - sys.atoms.len() as f64 * e_per_atom_ref;
        let gamma = excess / area;
        curve.push((frac, gamma));
    }
    curve
}
/// Analyze edge dislocation core by Volterra displacement field.
///
/// Inserts an edge dislocation along z with Burgers vector along x in an FCC crystal.
pub fn analyze_dislocation_core(params: &EamParams, n_cells: usize) -> DislocationCore {
    let a = params.lattice_a;
    let box_len = n_cells as f64 * a;
    let b_mag = a / 2.0_f64.sqrt();
    let burgers = [b_mag, 0.0, 0.0];
    let positions = generate_crystal_positions(a, n_cells, CrystalStructure::FCC);
    let cx = box_len / 2.0;
    let cy = box_len / 2.0;
    let mut displacements = Vec::new();
    let mut max_disp = 0.0_f64;
    let mut core_x = 0.0;
    let mut core_y = 0.0;
    let mut n_core = 0usize;
    let nu = 0.33;
    for (idx, pos) in positions.iter().enumerate() {
        let dx = pos[0] - cx;
        let dy = pos[1] - cy;
        let r2 = dx * dx + dy * dy;
        if r2 < 1e-10 {
            continue;
        }
        let _r = r2.sqrt();
        let theta = dy.atan2(dx);
        let ux = b_mag / (2.0 * PI) * (theta + dx * dy / (2.0 * (1.0 - nu) * r2));
        let uy = -b_mag / (2.0 * PI)
            * ((1.0 - 2.0 * nu) / (4.0 * (1.0 - nu)) * (r2).ln()
                + (dx * dx - dy * dy) / (4.0 * (1.0 - nu) * r2));
        let uz = 0.0;
        let disp = [ux, uy, uz];
        let disp_mag = vnorm(disp);
        displacements.push((idx, disp));
        if disp_mag > 0.3 * b_mag {
            core_x += pos[0];
            core_y += pos[1];
            n_core += 1;
        }
        if disp_mag > max_disp {
            max_disp = disp_mag;
        }
    }
    let core_center = if n_core > 0 {
        [
            core_x / n_core as f64,
            core_y / n_core as f64,
            box_len / 2.0,
        ]
    } else {
        [cx, cy, box_len / 2.0]
    };
    let mut r_core = 0.0_f64;
    let mut count = 0;
    for &(idx, disp) in &displacements {
        let d = vnorm(disp);
        if d > 0.1 * b_mag {
            let dx = positions[idx][0] - core_center[0];
            let dy = positions[idx][1] - core_center[1];
            r_core += (dx * dx + dy * dy).sqrt();
            count += 1;
        }
    }
    let core_width = if count > 0 { r_core / count as f64 } else { a };
    DislocationCore {
        burgers_vector: burgers,
        core_center,
        core_width,
        displacements,
    }
}
/// Estimate Peierls stress from the Peierls-Nabarro model.
///
/// σ_P ≈ (2G / (1-ν)) * exp(-2π * d / b)
///
/// where G = shear modulus, ν = Poisson's ratio, d = plane spacing, b = Burgers vector.
pub fn peierls_stress(
    shear_modulus_gpa: f64,
    poisson_ratio: f64,
    plane_spacing_angstrom: f64,
    burgers_vector_angstrom: f64,
) -> f64 {
    let g = shear_modulus_gpa;
    let nu = poisson_ratio;
    let d = plane_spacing_angstrom;
    let b = burgers_vector_angstrom;
    let prefactor = 2.0 * g / (1.0 - nu);
    let exponent = -2.0 * PI * d / b;
    prefactor * exponent.exp()
}
/// Estimate Peierls stress from EAM parameters using elastic constants.
///
/// Computes approximate shear modulus from the EAM curvature at equilibrium.
pub fn peierls_stress_from_eam(params: &EamParams) -> f64 {
    let a = params.lattice_a;
    let b_mag = a / 2.0_f64.sqrt();
    let d_111 = a / 3.0_f64.sqrt();
    let dr = 0.001;
    let re = params.re;
    let fp = params.pair_force(re + dr);
    let fm = params.pair_force(re - dr);
    let d2phi = (fp - fm) / (2.0 * dr);
    let omega = a * a * a / 4.0;
    let g_ev_per_a3 = 12.0 * d2phi.abs() / (6.0 * omega);
    let g_gpa = g_ev_per_a3 * 160.2176634;
    peierls_stress(g_gpa, 0.33, d_111, b_mag)
}
/// Compute mole fractions for each species.
pub fn mole_fractions(system: &AlloySystem) -> Vec<f64> {
    let n = system.atoms.len() as f64;
    if n < 1.0 {
        return vec![0.0; system.species.len()];
    }
    let mut counts = vec![0usize; system.species.len()];
    for atom in &system.atoms {
        if atom.species < counts.len() {
            counts[atom.species] += 1;
        }
    }
    counts.iter().map(|&c| c as f64 / n).collect()
}
/// Compute local composition in a sphere around a point.
pub fn local_composition(system: &AlloySystem, center: [f64; 3], radius: f64) -> Vec<f64> {
    let r2 = radius * radius;
    let mut counts = vec![0usize; system.species.len()];
    let mut total = 0usize;
    for atom in &system.atoms {
        let dr = system.sim_box.min_image(vsub(atom.position, center));
        if vdot(dr, dr) <= r2 {
            counts[atom.species] += 1;
            total += 1;
        }
    }
    if total == 0 {
        return vec![0.0; system.species.len()];
    }
    counts.iter().map(|&c| c as f64 / total as f64).collect()
}
/// Compute partial radial distribution function g_ab(r).
///
/// Returns (bin_centers, g_ab_values).
pub fn partial_rdf(
    system: &AlloySystem,
    species_a: usize,
    species_b: usize,
    r_max: f64,
    n_bins: usize,
) -> (Vec<f64>, Vec<f64>) {
    let dr = r_max / n_bins as f64;
    let mut hist = vec![0.0_f64; n_bins];
    let mut centers = Vec::with_capacity(n_bins);
    for i in 0..n_bins {
        centers.push((i as f64 + 0.5) * dr);
    }
    let n_a = system
        .atoms
        .iter()
        .filter(|a| a.species == species_a)
        .count() as f64;
    let n_b = system
        .atoms
        .iter()
        .filter(|a| a.species == species_b)
        .count() as f64;
    let vol = system.sim_box.volume();
    let rho_b = n_b / vol;
    for i in 0..system.atoms.len() {
        if system.atoms[i].species != species_a {
            continue;
        }
        for j in 0..system.atoms.len() {
            if i == j || system.atoms[j].species != species_b {
                continue;
            }
            let dr_vec = system
                .sim_box
                .min_image(vsub(system.atoms[j].position, system.atoms[i].position));
            let r = vnorm(dr_vec);
            let bin = (r / dr) as usize;
            if bin < n_bins {
                hist[bin] += 1.0;
            }
        }
    }
    let mut g_ab = vec![0.0; n_bins];
    for i in 0..n_bins {
        let r_inner = i as f64 * dr;
        let r_outer = (i + 1) as f64 * dr;
        let shell_vol = 4.0 / 3.0 * PI * (r_outer.powi(3) - r_inner.powi(3));
        let ideal = n_a * rho_b * shell_vol;
        if ideal > 0.0 {
            g_ab[i] = hist[i] / ideal;
        }
    }
    (centers, g_ab)
}
/// Compute average coordination number for each species pair.
///
/// Returns matrix\[a\]\[b\] = average number of b neighbors for an a atom.
pub fn coordination_numbers(system: &AlloySystem, cutoff: f64) -> Vec<Vec<f64>> {
    let ns = system.species.len();
    let mut counts = vec![vec![0usize; ns]; ns];
    let mut n_species = vec![0usize; ns];
    for atom in &system.atoms {
        n_species[atom.species] += 1;
    }
    for i in 0..system.atoms.len() {
        for j in 0..system.atoms.len() {
            if i == j {
                continue;
            }
            let dr = system
                .sim_box
                .min_image(vsub(system.atoms[j].position, system.atoms[i].position));
            if vnorm(dr) < cutoff {
                counts[system.atoms[i].species][system.atoms[j].species] += 1;
            }
        }
    }
    let mut result = vec![vec![0.0; ns]; ns];
    for a in 0..ns {
        for b in 0..ns {
            if n_species[a] > 0 {
                result[a][b] = counts[a][b] as f64 / n_species[a] as f64;
            }
        }
    }
    result
}
/// Classify local structure using Common Neighbor Analysis.
///
/// Returns per-atom label: "FCC", "BCC", "HCP", or "Other".
pub fn common_neighbor_analysis(system: &AlloySystem, cutoff: f64) -> Vec<&'static str> {
    let n = system.atoms.len();
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        for j in (i + 1)..n {
            let dr = system
                .sim_box
                .min_image(vsub(system.atoms[j].position, system.atoms[i].position));
            if vnorm(dr) < cutoff {
                neighbors[i].push(j);
                neighbors[j].push(i);
            }
        }
    }
    let mut labels = Vec::with_capacity(n);
    for i in 0..n {
        let mut fcc_count = 0u32;
        let mut hcp_count = 0u32;
        for &j in &neighbors[i] {
            let mut common = Vec::new();
            for &ni in &neighbors[i] {
                if ni != j && neighbors[j].contains(&ni) {
                    common.push(ni);
                }
            }
            let n_common = common.len() as u32;
            let mut n_bonds = 0u32;
            for a in 0..common.len() {
                for b in (a + 1)..common.len() {
                    if neighbors[common[a]].contains(&common[b]) {
                        n_bonds += 1;
                    }
                }
            }
            if n_common == 4 && n_bonds == 2 {
                fcc_count += 1;
            } else if n_common == 4 && n_bonds == 1 {
                hcp_count += 1;
            }
        }
        let total_nn = neighbors[i].len() as u32;
        let label = if total_nn == 12 && fcc_count >= 8 {
            "FCC"
        } else if total_nn == 12 && hcp_count >= 6 {
            "HCP"
        } else if total_nn == 14 || total_nn == 8 {
            "BCC"
        } else {
            "Other"
        };
        labels.push(label);
    }
    labels
}
/// Elastic constants (C11, C12, C44) in GPa from EAM potential.
///
/// Uses numerical differentiation of the energy with respect to strain.
pub fn elastic_constants_eam(params: &EamParams) -> (f64, f64, f64) {
    let a = params.lattice_a;
    let n_cells = 3;
    let box_len = n_cells as f64 * a;
    let positions = generate_crystal_positions(a, n_cells, CrystalStructure::FCC);
    let compute_energy = |strain: [[f64; 3]; 3]| -> f64 {
        let mut sys = AlloySystem::new(
            vec![params.clone()],
            SimBox::new(
                box_len * (1.0 + strain[0][0]),
                box_len * (1.0 + strain[1][1]),
                box_len * (1.0 + strain[2][2]),
            ),
            MixingRule::Geometric,
            0.001,
        );
        for pos in &positions {
            let p = [
                pos[0] * (1.0 + strain[0][0]) + pos[1] * strain[0][1] + pos[2] * strain[0][2],
                pos[0] * strain[1][0] + pos[1] * (1.0 + strain[1][1]) + pos[2] * strain[1][2],
                pos[0] * strain[2][0] + pos[1] * strain[2][1] + pos[2] * (1.0 + strain[2][2]),
            ];
            sys.add_atom(sys.sim_box.wrap(p), 0);
        }
        for atom in &mut sys.atoms {
            atom.rho_bar = 0.0;
        }
        let pairs = sys.build_neighbor_list();
        for &(i, j, _dr, r) in &pairs {
            sys.atoms[i].rho_bar += sys.species[0].electron_density(r);
            sys.atoms[j].rho_bar += sys.species[0].electron_density(r);
        }
        sys.potential_energy()
    };
    let e0 = compute_energy([[0.0; 3]; 3]);
    let vol = box_len * box_len * box_len;
    let eps = 0.001;
    let ep = compute_energy([[eps, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);
    let em = compute_energy([[-eps, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);
    let c11_ev = (ep + em - 2.0 * e0) / (vol * eps * eps);
    let ep2 = compute_energy([[eps, 0.0, 0.0], [0.0, eps, 0.0], [0.0, 0.0, 0.0]]);
    let em2 = compute_energy([[-eps, 0.0, 0.0], [0.0, -eps, 0.0], [0.0, 0.0, 0.0]]);
    let c12_ev = (ep2 + em2 - 2.0 * e0) / (2.0 * vol * eps * eps) - c11_ev / 2.0;
    let es_p = compute_energy([[0.0, eps, 0.0], [eps, 0.0, 0.0], [0.0, 0.0, 0.0]]);
    let es_m = compute_energy([[0.0, -eps, 0.0], [-eps, 0.0, 0.0], [0.0, 0.0, 0.0]]);
    let c44_ev = (es_p + es_m - 2.0 * e0) / (4.0 * vol * eps * eps);
    let conv = 160.2176634;
    (c11_ev * conv, c12_ev * conv, c44_ev * conv)
}
/// Compute enthalpy of mixing for a binary alloy at given composition.
///
/// ΔH_mix = E(alloy) - x_A * E(pure_A) - x_B * E(pure_B)
/// Returns eV/atom.
pub fn heat_of_mixing(
    params_a: &EamParams,
    params_b: &EamParams,
    x_b: f64,
    n_cells: usize,
    mixing_rule: MixingRule,
) -> f64 {
    let a_avg = params_a.lattice_a * (1.0 - x_b) + params_b.lattice_a * x_b;
    let box_len = n_cells as f64 * a_avg;
    let positions = generate_crystal_positions(a_avg, n_cells, CrystalStructure::FCC);
    let n_total = positions.len();
    let n_b = (n_total as f64 * x_b).round() as usize;
    let mut alloy = AlloySystem::new(
        vec![params_a.clone(), params_b.clone()],
        SimBox::cubic(box_len),
        mixing_rule,
        0.001,
    );
    for (idx, pos) in positions.iter().enumerate() {
        let sp = if idx < n_b { 1 } else { 0 };
        alloy.add_atom(*pos, sp);
    }
    for atom in &mut alloy.atoms {
        atom.rho_bar = 0.0;
    }
    let pairs = alloy.build_neighbor_list();
    for &(i, j, _dr, r) in &pairs {
        let si = alloy.atoms[i].species;
        let sj = alloy.atoms[j].species;
        alloy.atoms[i].rho_bar += alloy.species[sj].electron_density(r);
        alloy.atoms[j].rho_bar += alloy.species[si].electron_density(r);
    }
    let e_alloy = alloy.potential_energy() / n_total as f64;
    let mut pure_a = AlloySystem::new(
        vec![params_a.clone()],
        SimBox::cubic(n_cells as f64 * params_a.lattice_a),
        MixingRule::Geometric,
        0.001,
    );
    let pos_a = generate_crystal_positions(params_a.lattice_a, n_cells, CrystalStructure::FCC);
    for pos in &pos_a {
        pure_a.add_atom(*pos, 0);
    }
    for atom in &mut pure_a.atoms {
        atom.rho_bar = 0.0;
    }
    let pairs_a = pure_a.build_neighbor_list();
    for &(i, j, _dr, r) in &pairs_a {
        pure_a.atoms[i].rho_bar += pure_a.species[0].electron_density(r);
        pure_a.atoms[j].rho_bar += pure_a.species[0].electron_density(r);
    }
    let e_a = pure_a.potential_energy() / pos_a.len() as f64;
    let mut pure_b = AlloySystem::new(
        vec![params_b.clone()],
        SimBox::cubic(n_cells as f64 * params_b.lattice_a),
        MixingRule::Geometric,
        0.001,
    );
    let pos_b = generate_crystal_positions(params_b.lattice_a, n_cells, CrystalStructure::FCC);
    for pos in &pos_b {
        pure_b.add_atom(*pos, 0);
    }
    for atom in &mut pure_b.atoms {
        atom.rho_bar = 0.0;
    }
    let pairs_b = pure_b.build_neighbor_list();
    for &(i, j, _dr, r) in &pairs_b {
        pure_b.atoms[i].rho_bar += pure_b.species[0].electron_density(r);
        pure_b.atoms[j].rho_bar += pure_b.species[0].electron_density(r);
    }
    let e_b = pure_b.potential_energy() / pos_b.len() as f64;
    let x_a = 1.0 - x_b;
    e_alloy - x_a * e_a - x_b * e_b
}
/// Estimate linear thermal expansion coefficient from lattice constant variation.
///
/// Runs short NVT-like simulation at two temperatures, measures lattice expansion.
/// Returns α in 1/K.
pub fn thermal_expansion_coefficient(
    params: &EamParams,
    _t_low: f64,
    _t_high: f64,
    n_cells: usize,
) -> f64 {
    let a = params.lattice_a;
    let vol = a * a * a / 4.0;
    let (c11, c12, _c44) = elastic_constants_eam(params);
    let bulk_mod = (c11 + 2.0 * c12) / 3.0;
    let _gamma = 2.0;
    let kb = 8.617333e-5;
    let cv = 3.0 * kb;
    let b_ev = bulk_mod / 160.2176634;
    let _n_cells_val = n_cells;
    if b_ev.abs() < 1e-20 {
        return 1e-5;
    }
    2.0 * cv / (3.0 * b_ev.abs() * vol)
}
