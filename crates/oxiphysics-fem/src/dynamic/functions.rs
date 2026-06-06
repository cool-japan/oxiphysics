//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{ModalSuperpositionResult, ResponseSpectrumInput};
use crate::mesh::TetrahedralMesh;
use crate::sparse::CsrMatrix;

/// Compute the 12x12 consistent mass matrix for a single linear tetrahedron.
///
/// For a 4-node linear tetrahedron with density `rho` and volume `V`, the
/// consistent mass matrix in 3D is:
///
/// `M_e = rho * V / 20 * [[2I, I, I, I\], [I, 2I, I, I], [I, I, 2I, I], [I, I, I, 2I]]`
///
/// where `I` is the 3x3 identity matrix.
pub fn element_mass_matrix(nodes: &[[f64; 3]; 4], density: f64) -> [[f64; 12]; 12] {
    let x10 = [
        nodes[1][0] - nodes[0][0],
        nodes[1][1] - nodes[0][1],
        nodes[1][2] - nodes[0][2],
    ];
    let x20 = [
        nodes[2][0] - nodes[0][0],
        nodes[2][1] - nodes[0][1],
        nodes[2][2] - nodes[0][2],
    ];
    let x30 = [
        nodes[3][0] - nodes[0][0],
        nodes[3][1] - nodes[0][1],
        nodes[3][2] - nodes[0][2],
    ];
    let det = x10[0] * (x20[1] * x30[2] - x20[2] * x30[1])
        - x10[1] * (x20[0] * x30[2] - x20[2] * x30[0])
        + x10[2] * (x20[0] * x30[1] - x20[1] * x30[0]);
    let volume = det.abs() / 6.0;
    let factor = density * volume / 20.0;
    let mut me = [[0.0f64; 12]; 12];
    for i in 0..4 {
        for j in 0..4 {
            let scale = if i == j { 2.0 * factor } else { factor };
            for d in 0..3 {
                me[i * 3 + d][j * 3 + d] = scale;
            }
        }
    }
    me
}
/// Assemble the global consistent mass matrix from a tetrahedral mesh.
pub fn assemble_mass_matrix(mesh: &TetrahedralMesh, density: f64) -> CsrMatrix {
    let ndof = mesh.num_nodes() * 3;
    let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
    for elem in &mesh.elements {
        let nodes = [
            [
                mesh.nodes[elem[0]].x,
                mesh.nodes[elem[0]].y,
                mesh.nodes[elem[0]].z,
            ],
            [
                mesh.nodes[elem[1]].x,
                mesh.nodes[elem[1]].y,
                mesh.nodes[elem[1]].z,
            ],
            [
                mesh.nodes[elem[2]].x,
                mesh.nodes[elem[2]].y,
                mesh.nodes[elem[2]].z,
            ],
            [
                mesh.nodes[elem[3]].x,
                mesh.nodes[elem[3]].y,
                mesh.nodes[elem[3]].z,
            ],
        ];
        let me = element_mass_matrix(&nodes, density);
        for local_i in 0..4 {
            for local_j in 0..4 {
                for di in 0..3 {
                    for dj in 0..3 {
                        let global_i = elem[local_i] * 3 + di;
                        let global_j = elem[local_j] * 3 + dj;
                        let val = me[local_i * 3 + di][local_j * 3 + dj];
                        if val.abs() > 1e-60 {
                            triplets.push((global_i, global_j, val));
                        }
                    }
                }
            }
        }
    }
    CsrMatrix::from_triplets(ndof, ndof, &triplets)
}
/// Compute a lumped (diagonal) mass matrix from a consistent mass matrix
/// using row-sum lumping.
///
/// Returns the diagonal as a `Vec`f64`.
pub fn lump_mass_row_sum(mass: &CsrMatrix) -> Vec<f64> {
    let n = mass.nrows;
    let mut diag = vec![0.0; n];
    for (row, diag_row) in diag.iter_mut().enumerate().take(n) {
        let start = mass.row_ptr[row];
        let end = mass.row_ptr[row + 1];
        for idx in start..end {
            *diag_row += mass.values[idx];
        }
    }
    diag
}
/// Compute a lumped mass matrix using diagonal scaling (HRZ lumping).
///
/// The HRZ method scales the diagonal entries so that the total mass
/// is preserved: `m_lumped\[i\] = m_diag\[i\] * (total_mass / sum_of_diag)`.
pub fn lump_mass_hrz(mass: &CsrMatrix) -> Vec<f64> {
    let n = mass.nrows;
    let mut diag = vec![0.0; n];
    for (row, diag_row) in diag.iter_mut().enumerate().take(n) {
        let start = mass.row_ptr[row];
        let end = mass.row_ptr[row + 1];
        for idx in start..end {
            if mass.col_indices[idx] == row {
                *diag_row = mass.values[idx];
            }
        }
    }
    let total_mass: f64 = (0..n)
        .map(|row| {
            let start = mass.row_ptr[row];
            let end = mass.row_ptr[row + 1];
            mass.values[start..end].iter().sum::<f64>()
        })
        .sum::<f64>()
        / 3.0;
    let sum_diag: f64 = diag.iter().sum();
    if sum_diag.abs() < 1e-60 {
        return diag;
    }
    let scale = total_mass * 3.0 / sum_diag;
    for d in &mut diag {
        *d *= scale;
    }
    diag
}
/// Compute the inverse of a diagonal mass vector, setting zero masses
/// to zero inverse (static DOFs).
pub fn invert_diagonal_mass(m_diag: &[f64]) -> Vec<f64> {
    m_diag
        .iter()
        .map(|&m| if m.abs() > 1e-60 { 1.0 / m } else { 0.0 })
        .collect()
}
/// Compute the kinetic energy: `T = 0.5 * v^T * M * v`.
///
/// Uses a diagonal mass for efficiency.
pub fn kinetic_energy_diagonal(v: &[f64], m_diag: &[f64]) -> f64 {
    let mut t = 0.0;
    for (vi, &mi) in v.iter().zip(m_diag.iter()) {
        t += mi * vi * vi;
    }
    0.5 * t
}
/// Compute the kinetic energy using a full (CSR) mass matrix.
pub fn kinetic_energy_full(v: &[f64], mass: &CsrMatrix) -> f64 {
    let mv = mass.mul_vec(v);
    let mut t = 0.0;
    for (vi, mvi) in v.iter().zip(mv.iter()) {
        t += vi * mvi;
    }
    0.5 * t
}
/// Compute the strain (potential) energy: `U = 0.5 * u^T * K * u`.
pub fn strain_energy(u: &[f64], stiffness: &CsrMatrix) -> f64 {
    let ku = stiffness.mul_vec(u);
    let mut e = 0.0;
    for (ui, kui) in u.iter().zip(ku.iter()) {
        e += ui * kui;
    }
    0.5 * e
}
/// Compute the total mechanical energy (kinetic + potential).
pub fn total_energy(u: &[f64], v: &[f64], m_diag: &[f64], stiffness: &CsrMatrix) -> f64 {
    kinetic_energy_diagonal(v, m_diag) + strain_energy(u, stiffness)
}
/// Assemble the Rayleigh damping matrix: `C = alpha * M + beta * K`.
///
/// For typical structural applications:
/// - `alpha` (mass proportional) damps low-frequency modes
/// - `beta` (stiffness proportional) damps high-frequency modes
pub fn rayleigh_damping_matrix(
    mass: &CsrMatrix,
    stiffness: &CsrMatrix,
    alpha: f64,
    beta: f64,
) -> CsrMatrix {
    let n = mass.nrows;
    let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
    for row in 0..mass.nrows {
        let start = mass.row_ptr[row];
        let end = mass.row_ptr[row + 1];
        for idx in start..end {
            let val = alpha * mass.values[idx];
            if val.abs() > 1e-60 {
                triplets.push((row, mass.col_indices[idx], val));
            }
        }
    }
    for row in 0..stiffness.nrows {
        let start = stiffness.row_ptr[row];
        let end = stiffness.row_ptr[row + 1];
        for idx in start..end {
            let val = beta * stiffness.values[idx];
            if val.abs() > 1e-60 {
                triplets.push((row, stiffness.col_indices[idx], val));
            }
        }
    }
    CsrMatrix::from_triplets(n, n, &triplets)
}
/// Compute Rayleigh coefficients alpha and beta to achieve a given damping ratio `zeta`
/// at a single frequency `omega_i` (single-mode targeting).
///
/// From `zeta = alpha/(2*omega) + beta*omega/2`:
/// - If only targeting one mode, there are infinite solutions.
/// - This function returns the stiffness-only solution: `alpha = 0, beta = 2*zeta/omega`.
pub fn rayleigh_from_modal(omega: f64, zeta: f64) -> (f64, f64) {
    if omega < 1e-15 {
        return (0.0, 0.0);
    }
    let beta = 2.0 * zeta / omega;
    let alpha = 0.0;
    (alpha, beta)
}
/// Compute Rayleigh coefficients from two target frequencies and damping ratios.
///
/// Given `(omega1, zeta1)` and `(omega2, zeta2)`, solve the 2×2 system:
/// ```text
/// [1/(2*omega1)  omega1/2] [alpha]   [zeta1]
/// [1/(2*omega2)  omega2/2] [beta ] = [zeta2]
/// ```
pub fn rayleigh_from_two_modes(omega1: f64, zeta1: f64, omega2: f64, zeta2: f64) -> (f64, f64) {
    let a11 = 0.5 / omega1;
    let a12 = 0.5 * omega1;
    let a21 = 0.5 / omega2;
    let a22 = 0.5 * omega2;
    let det = a11 * a22 - a12 * a21;
    if det.abs() < 1e-30 {
        return (0.0, 0.0);
    }
    let alpha = (a22 * zeta1 - a12 * zeta2) / det;
    let beta = (-a21 * zeta1 + a11 * zeta2) / det;
    (alpha, beta)
}
/// Compute the effective seismic force vector for base excitation.
///
/// For ground acceleration `a_g` in direction `d` (0=x, 1=y, 2=z),
/// the equivalent force is: `F_eff = -M * r * a_g`
///
/// where `r` is the influence vector (1 for DOFs in the excitation direction, 0 otherwise).
///
/// `n_dof_per_node` = 3 (3D) or 2 (2D).
pub fn seismic_base_excitation_force(
    mass: &CsrMatrix,
    ground_acceleration: &[f64],
    n_dof_per_node: usize,
) -> Vec<f64> {
    let n = mass.nrows;
    let ndim = ground_acceleration.len().min(n_dof_per_node);
    let mut r = vec![0.0_f64; n];
    for (i, r_i) in r.iter_mut().enumerate().take(n) {
        let dof_dir = i % n_dof_per_node;
        if dof_dir < ndim {
            *r_i = ground_acceleration[dof_dir];
        }
    }
    let mr = mass.mul_vec(&r);
    mr.iter().map(|&v| -v).collect()
}
/// Compute the contact impulse for a DOF that has penetrated its constraint.
///
/// Uses the impulse-based contact model:
/// `J = -(1 + e) * v_n / (1/m_a + 1/m_b)`
///
/// where `e` is the coefficient of restitution, `v_n` is the normal velocity,
/// and `m_a`, `m_b` are the masses of the two bodies (here, body B is a rigid wall with infinite mass).
///
/// Returns the impulse to apply to the DOF. Returns 0 if gap > 0 (no contact).
pub fn contact_impulse(
    gap: f64,
    velocities: &[f64],
    m_diag: &[f64],
    contact_dof: usize,
    restitution: f64,
) -> f64 {
    if gap >= 0.0 {
        return 0.0;
    }
    if contact_dof >= velocities.len() {
        return 0.0;
    }
    let v_n = velocities[contact_dof];
    if v_n >= 0.0 {
        return 0.0;
    }
    let m = m_diag[contact_dof];
    if m < 1e-60 {
        return 0.0;
    }
    -(1.0 + restitution) * v_n * m
}
/// Apply contact impulses to a velocity vector for a set of contact DOFs.
///
/// For each contact DOF, if penetration exists, applies the impulse.
pub fn apply_contact_impulses(
    velocities: &mut [f64],
    gaps: &[f64],
    m_diag: &[f64],
    contact_dofs: &[usize],
    restitution: f64,
) {
    for &dof in contact_dofs {
        if dof >= gaps.len() {
            continue;
        }
        let j = contact_impulse(gaps[dof], velocities, m_diag, dof, restitution);
        if dof < velocities.len() && dof < m_diag.len() && m_diag[dof] > 1e-60 {
            velocities[dof] += j / m_diag[dof];
        }
    }
}
/// Estimate the critical time step for the central difference method.
///
/// For a linear tetrahedron: `dt_crit = L_min / c` where
/// `c = sqrt(E/rho)` is the wave speed and `L_min` is the
/// minimum element edge length.
pub fn critical_time_step(min_edge_length: f64, youngs_modulus: f64, density: f64) -> f64 {
    let c = (youngs_modulus / density).sqrt();
    if c < 1e-60 {
        return f64::MAX;
    }
    min_edge_length / c
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::assemble_stiffness;
    use crate::constitutive::LinearElasticMaterial;
    use crate::dynamic::*;
    #[test]
    fn test_element_mass_matrix_symmetry() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let density = 7800.0;
        let me = element_mass_matrix(&nodes, density);
        for (i, row) in me.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let diff = (val - me[j][i]).abs();
                let scale = val.abs().max(me[j][i].abs()).max(1e-30);
                assert!(
                    diff / scale < 1e-12,
                    "Mass matrix not symmetric at ({i},{j}): {} vs {}",
                    val,
                    me[j][i]
                );
            }
        }
    }
    #[test]
    fn test_element_mass_matrix_trace() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let density = 7800.0;
        let volume = 1.0 / 6.0;
        let total_mass = density * volume;
        let me = element_mass_matrix(&nodes, density);
        let trace: f64 = (0..12).map(|i| me[i][i]).sum();
        let expected_trace = 1.2 * total_mass;
        assert!(
            (trace - expected_trace).abs() / expected_trace < 1e-10,
            "trace = {trace}, expected {expected_trace}"
        );
    }
    fn make_bar_mesh() -> TetrahedralMesh {
        use oxiphysics_core::math::Vec3;
        let nodes = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        TetrahedralMesh::from_nodes_and_elements(nodes, vec![[0, 1, 2, 3]])
    }
    #[test]
    fn test_newmark_free_vibration() {
        let mesh = make_bar_mesh();
        let density = 7800.0;
        let material = LinearElasticMaterial::new(200.0e9, 0.3);
        let mass = assemble_mass_matrix(&mesh, density);
        let stiffness = assemble_stiffness(&mesh, &material);
        let ndof = mesh.num_nodes() * 3;
        let mut state = NewmarkState::new(ndof);
        let fixed_dofs: Vec<usize> = vec![0, 1, 2];
        let mut force = vec![0.0; ndof];
        force[3] = 1.0e6;
        state.step(&mass, &stiffness, &force, 1e-5, &fixed_dofs);
        let zero_force = vec![0.0; ndof];
        let initial_u3 = state.u[3];
        let mut changed_sign = false;
        for _ in 0..50 {
            state.step(&mass, &stiffness, &zero_force, 1e-5, &fixed_dofs);
            if state.u[3] * initial_u3 < 0.0 {
                changed_sign = true;
                break;
            }
        }
        assert!(
            changed_sign,
            "Expected free vibration to change sign within 50 steps, initial_u3={initial_u3:.3e}"
        );
    }
    #[test]
    fn test_newmark_static_limit() {
        let k_spring = 1000.0_f64;
        let mass_val = 1.0_f64;
        let f0 = 10.0_f64;
        let stiffness = CsrMatrix::from_triplets(
            2,
            2,
            &[
                (0, 0, k_spring),
                (0, 1, -k_spring),
                (1, 0, -k_spring),
                (1, 1, k_spring),
            ],
        );
        let mass = CsrMatrix::from_triplets(2, 2, &[(0, 0, mass_val), (1, 1, mass_val)]);
        let force = vec![0.0, f0];
        let fixed_dofs = vec![0usize];
        let dt = 0.01;
        let mut state = NewmarkState::new(2);
        for _ in 0..5000 {
            state.step(&mass, &stiffness, &force, dt, &fixed_dofs);
        }
        let u_static = f0 / k_spring;
        let u_dynamic = state.u[1];
        assert!(
            u_dynamic > 0.0,
            "Displacement should be positive, got {u_dynamic:.4e}"
        );
        assert!(
            u_dynamic < 2.5 * u_static,
            "Displacement should stay within 2x static solution for undamped system, \
             got {u_dynamic:.4e}, static={u_static:.4e}"
        );
    }
    /// HHT-alpha with alpha=0 should behave like standard Newmark.
    #[test]
    fn test_hht_alpha_zero_is_newmark() {
        let k_spring = 1000.0_f64;
        let mass_val = 1.0_f64;
        let f0 = 10.0_f64;
        let stiffness = CsrMatrix::from_triplets(
            2,
            2,
            &[
                (0, 0, k_spring),
                (0, 1, -k_spring),
                (1, 0, -k_spring),
                (1, 1, k_spring),
            ],
        );
        let mass = CsrMatrix::from_triplets(2, 2, &[(0, 0, mass_val), (1, 1, mass_val)]);
        let force = vec![0.0, f0];
        let fixed_dofs = vec![0usize];
        let dt = 0.01;
        let mut hht = HhtAlphaState::new(2, 0.0);
        for _ in 0..100 {
            let f_old = force.clone();
            hht.step(&mass, &stiffness, &force, &f_old, dt, &fixed_dofs);
        }
        let u_static = f0 / k_spring;
        assert!(hht.u[1] > 0.0, "HHT u should be positive, got {}", hht.u[1]);
        assert!(
            hht.u[1] < 2.5 * u_static,
            "HHT u should be bounded, got {}",
            hht.u[1]
        );
    }
    /// HHT-alpha with negative alpha provides dissipation.
    #[test]
    fn test_hht_alpha_parameters() {
        let state = HhtAlphaState::new(1, -0.1);
        assert!((state.gamma - 0.6).abs() < 1e-12, "gamma = {}", state.gamma);
        assert!((state.beta - 0.3025).abs() < 1e-12, "beta = {}", state.beta);
    }
    /// Generalized-alpha with rho_inf=1 should have specific parameter values.
    #[test]
    fn test_gen_alpha_no_dissipation() {
        let state = GeneralizedAlphaState::new(1, 1.0);
        assert!(
            (state.alpha_m - 0.5).abs() < 1e-12,
            "alpha_m = {}",
            state.alpha_m
        );
        assert!(
            (state.alpha_f - 0.5).abs() < 1e-12,
            "alpha_f = {}",
            state.alpha_f
        );
        assert!((state.gamma - 0.5).abs() < 1e-12, "gamma = {}", state.gamma);
    }
    /// Generalized-alpha with rho_inf=0 (max dissipation).
    #[test]
    fn test_gen_alpha_max_dissipation() {
        let state = GeneralizedAlphaState::new(1, 0.0);
        assert!(
            (state.alpha_m - (-1.0)).abs() < 1e-12,
            "alpha_m = {}",
            state.alpha_m
        );
        assert!((state.alpha_f).abs() < 1e-12, "alpha_f = {}", state.alpha_f);
        assert!((state.gamma - 1.5).abs() < 1e-12, "gamma = {}", state.gamma);
    }
    /// A single-DOF spring-mass system under the central difference method
    /// should oscillate with the analytical frequency.
    #[test]
    fn test_central_difference_spring_mass() {
        let k_val = 100.0_f64;
        let m_val = 1.0_f64;
        let stiffness = CsrMatrix::from_triplets(
            2,
            2,
            &[(0, 0, k_val), (0, 1, -k_val), (1, 0, -k_val), (1, 1, k_val)],
        );
        let inv_m = vec![0.0, 1.0 / m_val];
        let mut state = CentralDifferenceState::new(2, inv_m);
        state.u[1] = 0.1;
        state.u_prev[1] = 0.1;
        let force = vec![0.0; 2];
        let fixed_dofs = vec![0usize];
        let dt = 0.001;
        let mut changed_sign = false;
        for _ in 0..1000 {
            state.step(&stiffness, &force, dt, &fixed_dofs);
            if state.u[1] < 0.0 {
                changed_sign = true;
                break;
            }
        }
        assert!(
            changed_sign,
            "Expected oscillation to cross zero, u[1] = {}",
            state.u[1]
        );
    }
    /// Row-sum lumping of the consistent mass matrix should preserve total mass.
    #[test]
    fn test_lump_mass_row_sum_preserves_total() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let density = 7800.0;
        let me = element_mass_matrix(&nodes, density);
        let mut triplets = Vec::new();
        for (i, row) in me.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                if val.abs() > 1e-60 {
                    triplets.push((i, j, val));
                }
            }
        }
        let mass = CsrMatrix::from_triplets(12, 12, &triplets);
        let lumped = lump_mass_row_sum(&mass);
        let total_lumped: f64 = lumped.iter().sum();
        let total_consistent: f64 = me.iter().flat_map(|r| r.iter()).sum();
        assert!(
            (total_lumped - total_consistent).abs() / total_consistent.abs() < 1e-10,
            "lumped total = {total_lumped}, consistent total = {total_consistent}"
        );
    }
    /// All lumped masses should be positive.
    #[test]
    fn test_lump_mass_all_positive() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let density = 7800.0;
        let me = element_mass_matrix(&nodes, density);
        let mut triplets = Vec::new();
        for (i, row) in me.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                if val.abs() > 1e-60 {
                    triplets.push((i, j, val));
                }
            }
        }
        let mass = CsrMatrix::from_triplets(12, 12, &triplets);
        let lumped = lump_mass_row_sum(&mass);
        for (i, &m) in lumped.iter().enumerate() {
            assert!(m > 0.0, "lumped mass[{i}] = {m}, should be positive");
        }
    }
    /// Inverse diagonal mass should be finite and positive for positive masses.
    #[test]
    fn test_invert_diagonal_mass() {
        let m = vec![1.0, 2.0, 0.5, 0.0];
        let inv = invert_diagonal_mass(&m);
        assert!((inv[0] - 1.0).abs() < 1e-12);
        assert!((inv[1] - 0.5).abs() < 1e-12);
        assert!((inv[2] - 2.0).abs() < 1e-12);
        assert!((inv[3]).abs() < 1e-30);
    }
    /// Kinetic energy of unit mass at unit velocity = 0.5.
    #[test]
    fn test_kinetic_energy_diagonal() {
        let v = vec![1.0, 0.0, 0.0];
        let m = vec![1.0, 1.0, 1.0];
        let ke = kinetic_energy_diagonal(&v, &m);
        assert!((ke - 0.5).abs() < 1e-12, "KE = {ke}, expected 0.5");
    }
    /// Strain energy of a unit displacement under unit stiffness = 0.5.
    #[test]
    fn test_strain_energy() {
        let stiffness = CsrMatrix::from_triplets(1, 1, &[(0, 0, 1.0)]);
        let u = vec![1.0];
        let se = strain_energy(&u, &stiffness);
        assert!((se - 0.5).abs() < 1e-12, "SE = {se}, expected 0.5");
    }
    /// Energy monitor should track drift correctly.
    #[test]
    fn test_energy_monitor_drift() {
        let mut mon = EnergyMonitor::new();
        assert!(mon.is_empty());
        mon.record(0.0, 50.0, 50.0);
        mon.record(0.1, 48.0, 52.0);
        mon.record(0.2, 51.0, 50.0);
        assert_eq!(mon.len(), 3);
        let drift = mon.max_relative_drift();
        assert!(
            (drift - 0.01).abs() < 1e-10,
            "drift = {drift}, expected 0.01"
        );
    }
    /// Critical time step for steel: E=200 GPa, rho=7800, L=0.01 m.
    #[test]
    fn test_critical_time_step() {
        let e = 200.0e9;
        let rho = 7800.0;
        let l_min = 0.01;
        let dt_crit = critical_time_step(l_min, e, rho);
        let c = (e / rho).sqrt();
        let expected = l_min / c;
        assert!(
            (dt_crit - expected).abs() < 1e-15,
            "dt_crit = {dt_crit}, expected {expected}"
        );
        assert!(dt_crit < 1e-5, "dt_crit should be small, got {dt_crit}");
    }
    /// HRZ lumping should produce positive diagonal entries.
    #[test]
    fn test_hrz_lumping_positive() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let density = 7800.0;
        let me = element_mass_matrix(&nodes, density);
        let mut triplets = Vec::new();
        for (i, row) in me.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                if val.abs() > 1e-60 {
                    triplets.push((i, j, val));
                }
            }
        }
        let mass = CsrMatrix::from_triplets(12, 12, &triplets);
        let lumped = lump_mass_hrz(&mass);
        for (i, &m) in lumped.iter().enumerate() {
            assert!(m > 0.0, "HRZ lumped mass[{i}] = {m}, should be positive");
        }
    }
    /// Total energy should be constant for free vibration with Newmark.
    #[test]
    fn test_total_energy_function() {
        let u = vec![1.0, 0.0];
        let v = vec![0.0, 1.0];
        let m = vec![1.0, 1.0];
        let stiffness = CsrMatrix::from_triplets(2, 2, &[(0, 0, 1.0), (1, 1, 1.0)]);
        let te = total_energy(&u, &v, &m, &stiffness);
        assert!(
            (te - 1.0).abs() < 1e-12,
            "total energy = {te}, expected 1.0"
        );
    }
    /// Rayleigh damping matrix: C = alpha * M + beta * K.
    #[test]
    fn test_rayleigh_damping_construction() {
        let mass = CsrMatrix::from_triplets(2, 2, &[(0, 0, 2.0), (1, 1, 3.0)]);
        let stiff = CsrMatrix::from_triplets(2, 2, &[(0, 0, 100.0), (1, 1, 200.0)]);
        let c = rayleigh_damping_matrix(&mass, &stiff, 0.01, 0.001);
        let cv = c.mul_vec(&[1.0, 0.0]);
        assert!((cv[0] - 0.12).abs() < 1e-12, "C[0,0] = {}", cv[0]);
        let cv2 = c.mul_vec(&[0.0, 1.0]);
        assert!((cv2[1] - 0.23).abs() < 1e-12, "C[1,1] = {}", cv2[1]);
    }
    #[test]
    fn test_rayleigh_damping_zero_coeffs() {
        let mass = CsrMatrix::from_triplets(2, 2, &[(0, 0, 1.0), (1, 1, 1.0)]);
        let stiff = CsrMatrix::from_triplets(2, 2, &[(0, 0, 1.0), (1, 1, 1.0)]);
        let c = rayleigh_damping_matrix(&mass, &stiff, 0.0, 0.0);
        let cv = c.mul_vec(&[1.0, 1.0]);
        assert!(
            cv.iter().all(|&v| v.abs() < 1e-12),
            "Zero-coeff damping should give zero C"
        );
    }
    #[test]
    fn test_seismic_excitation_force_nonzero() {
        let mass = CsrMatrix::from_triplets(2, 2, &[(0, 0, 1.0), (1, 1, 2.0)]);
        let acc = vec![9.8, 0.0, 0.0];
        let f = seismic_base_excitation_force(&mass, &acc, 3);
        let max_f = f.iter().map(|&v| v.abs()).fold(0.0_f64, f64::max);
        assert!(max_f > 0.0, "Seismic force should be nonzero");
    }
    #[test]
    fn test_seismic_excitation_zero_acceleration() {
        let mass = CsrMatrix::from_triplets(2, 2, &[(0, 0, 1.0), (1, 1, 2.0)]);
        let acc = vec![0.0, 0.0, 0.0];
        let f = seismic_base_excitation_force(&mass, &acc, 3);
        assert!(
            f.iter().all(|&v| v.abs() < 1e-15),
            "Zero acceleration → zero force"
        );
    }
    #[test]
    fn test_contact_impulse_positive() {
        let v = vec![-1.0, 0.0];
        let m_diag = vec![1.0, 1.0];
        let contact_dof = 0;
        let gap = -0.01;
        let imp = contact_impulse(gap, &v, &m_diag, contact_dof, 0.8);
        assert!(
            imp > 0.0,
            "Impact impulse should be positive (restoring): {imp}"
        );
    }
    #[test]
    fn test_contact_impulse_no_penetration() {
        let v = vec![1.0, 0.0];
        let m_diag = vec![1.0, 1.0];
        let imp = contact_impulse(0.01, &v, &m_diag, 0, 0.8);
        assert!(
            imp.abs() < 1e-14,
            "Separating contact → zero impulse: {imp}"
        );
    }
    #[test]
    fn test_modal_damping_diagonal() {
        let omega = 10.0_f64;
        let zeta = 0.05;
        let (alpha, beta) = rayleigh_from_modal(omega, zeta);
        let zeta_check = alpha / (2.0 * omega) + beta * omega / 2.0;
        assert!(
            (zeta_check - zeta).abs() < 1e-10,
            "zeta check = {zeta_check}"
        );
    }
}
/// Modal superposition solution for M*u'' + C*u' + K*u = F(t).
///
/// Transforms to modal coordinates: q_r'' + 2*zeta_r*omega_r*q_r' + omega_r^2*q_r = F_r(t)/M_r
/// using the Newmark-beta method in each decoupled modal equation.
///
/// # Arguments
/// * `omega_r` – natural frequency of mode r (rad/s)
/// * `zeta_r` – modal damping ratio for mode r
/// * `modal_mass` – generalised mass M_r = phi^T M phi
/// * `modal_force_history` – time history of generalised force F_r(t), one value per step
/// * `dt` – time step (s)
///
/// Returns the time history of modal displacement q_r(t).
pub fn modal_newmark_1dof(
    omega_r: f64,
    zeta_r: f64,
    modal_mass: f64,
    modal_force_history: &[f64],
    dt: f64,
) -> Vec<f64> {
    let n = modal_force_history.len();
    if n == 0 {
        return Vec::new();
    }
    let beta = 0.25_f64;
    let gamma = 0.5_f64;
    let k_modal = omega_r * omega_r * modal_mass;
    let c_modal = 2.0 * zeta_r * omega_r * modal_mass;
    let inv_beta_dt2 = 1.0 / (beta * dt * dt);
    let k_eff = k_modal + c_modal * gamma / (beta * dt) + modal_mass * inv_beta_dt2;
    let mut q = vec![0.0_f64; n];
    let mut dq = vec![0.0_f64; n];
    let mut ddq = vec![0.0_f64; n];
    if modal_mass.abs() > f64::EPSILON {
        ddq[0] = (modal_force_history[0] - c_modal * dq[0] - k_modal * q[0]) / modal_mass;
    }
    let inv_beta_dt = 1.0 / (beta * dt);
    let gamma_over_beta = gamma / beta;
    for i in 1..n {
        let a1 = inv_beta_dt2;
        let a2 = inv_beta_dt;
        let a3 = 0.5 / beta - 1.0;
        let a4 = gamma_over_beta * inv_beta_dt;
        let a5 = gamma_over_beta - 1.0;
        let a6 = 0.5 * dt * (gamma_over_beta - 2.0);
        let f_eff = modal_force_history[i]
            + modal_mass * (a1 * q[i - 1] + a2 * dq[i - 1] + a3 * ddq[i - 1])
            + c_modal * (a4 * q[i - 1] + a5 * dq[i - 1] + a6 * ddq[i - 1]);
        if k_eff.abs() < f64::EPSILON {
            q[i] = 0.0;
        } else {
            q[i] = f_eff / k_eff;
        }
        ddq[i] = a1 * (q[i] - q[i - 1]) - a2 * dq[i - 1] - a3 * ddq[i - 1];
        dq[i] = dq[i - 1] + dt * ((1.0 - gamma) * ddq[i - 1] + gamma * ddq[i]);
    }
    q
}
/// Compute the response spectrum (pseudo-spectral acceleration Sa vs. period T)
/// for a given earthquake acceleration record.
///
/// For each natural period T in the range, an SDOF oscillator is subjected to
/// the ground acceleration record. The maximum absolute acceleration is the Sa.
///
/// # Arguments
/// * `ground_acc` – time history of ground acceleration (m/s²)
/// * `dt` – sampling interval (s)
/// * `periods` – slice of natural periods (s) at which to compute Sa
/// * `zeta` – viscous damping ratio (typically 0.05)
///
/// Returns a vector of (T, Sa) pairs.
pub fn response_spectrum(
    ground_acc: &[f64],
    dt: f64,
    periods: &[f64],
    zeta: f64,
) -> Vec<(f64, f64)> {
    let n = ground_acc.len();
    periods
        .iter()
        .map(|&t_period| {
            if t_period < 1e-6 {
                let sa = ground_acc.iter().map(|a| a.abs()).fold(0.0_f64, f64::max);
                return (t_period, sa);
            }
            let omega = 2.0 * std::f64::consts::PI / t_period;
            let modal_force: Vec<f64> = ground_acc.iter().map(|&ag| -ag).collect();
            let q = modal_newmark_1dof(omega, zeta, 1.0, &modal_force, dt);
            let q_max = q.iter().take(n).map(|v| v.abs()).fold(0.0_f64, f64::max);
            let sa = omega * omega * q_max;
            (t_period, sa)
        })
        .collect()
}
/// Generate a half-sine pulse ground acceleration record.
///
/// The pulse has the form: a(t) = a_max * sin(π * t / t_d)  for 0 ≤ t ≤ t_d
///                                0                          otherwise
///
/// # Arguments
/// * `a_max` – peak acceleration (m/s²)
/// * `t_d` – pulse duration (s)
/// * `dt` – time step (s)
/// * `t_total` – total record duration (s)
///
/// Returns a vector of acceleration values.
pub fn half_sine_pulse(a_max: f64, t_d: f64, dt: f64, t_total: f64) -> Vec<f64> {
    let n = (t_total / dt.max(f64::EPSILON)) as usize + 1;
    (0..n)
        .map(|i| {
            let t = i as f64 * dt;
            if t <= t_d {
                a_max * (std::f64::consts::PI * t / t_d.max(f64::EPSILON)).sin()
            } else {
                0.0
            }
        })
        .collect()
}
/// Generate a symmetric triangular pulse ground acceleration.
///
/// The pulse rises linearly from 0 to a_max over [0, t_d/2] and falls
/// linearly back to 0 over [t_d/2, t_d].
///
/// # Arguments
/// * `a_max` – peak acceleration (m/s²)
/// * `t_d` – pulse duration (s)
/// * `dt` – time step (s)
/// * `t_total` – total record duration (s)
pub fn triangular_pulse(a_max: f64, t_d: f64, dt: f64, t_total: f64) -> Vec<f64> {
    let n = (t_total / dt.max(f64::EPSILON)) as usize + 1;
    let t_half = t_d / 2.0;
    (0..n)
        .map(|i| {
            let t = i as f64 * dt;
            if t <= 0.0 || t >= t_d {
                0.0
            } else if t <= t_half {
                a_max * t / t_half.max(f64::EPSILON)
            } else {
                a_max * (t_d - t) / t_half.max(f64::EPSILON)
            }
        })
        .collect()
}
/// Generate a rectangular (boxcar) ground acceleration pulse.
///
/// The acceleration is `a_max` for 0 ≤ t ≤ t_d and 0 otherwise.
pub fn rectangular_pulse(a_max: f64, t_d: f64, dt: f64, t_total: f64) -> Vec<f64> {
    let n = (t_total / dt.max(f64::EPSILON)) as usize + 1;
    (0..n)
        .map(|i| {
            let t = i as f64 * dt;
            if t <= t_d { a_max } else { 0.0 }
        })
        .collect()
}
/// Dynamic stress intensity factor K_I^dyn using the Freund formula.
///
/// The dynamic SIF is related to the static SIF via a universal function
/// of the crack velocity:
///
/// K_I^dyn(v) = k(v) * K_I^static
///
/// where k(v) ≈ (1 - v / c_R) / sqrt(1 - v / c_d)  (Freund, 1990)
///
/// # Arguments
/// * `k1_static` – static mode-I SIF (Pa√m)
/// * `crack_velocity` – current crack velocity (m/s)
/// * `c_r` – Rayleigh wave speed (m/s)
/// * `c_d` – dilatational (P-wave) speed (m/s)
pub fn dynamic_sif_freund(k1_static: f64, crack_velocity: f64, c_r: f64, c_d: f64) -> f64 {
    let v = crack_velocity.clamp(0.0, c_r * 0.9999);
    let cr = c_r.max(f64::EPSILON);
    let cd = c_d.max(f64::EPSILON);
    let k_v = (1.0 - v / cr) / (1.0 - v / cd).abs().sqrt();
    k1_static * k_v
}
/// Compute the dynamic energy release rate G_dyn from the dynamic SIF.
///
/// For plane strain: G_dyn = (1 - ν²) / E * K_I^dyn²
pub fn dynamic_energy_release_rate(k1_dyn: f64, young_modulus: f64, poisson: f64) -> f64 {
    let e_prime = young_modulus / (1.0 - poisson * poisson).max(f64::EPSILON);
    k1_dyn * k1_dyn / e_prime
}
/// Compute the modal participation factor Γ_r for mode r.
///
/// Γ_r = phi_r^T * M * r_influence / (phi_r^T * M * phi_r)
///
/// where phi_r is the mode shape vector and r_influence is the influence
/// vector (1 for DOFs in excitation direction, 0 otherwise).
///
/// # Arguments
/// * `mode_shape` – eigenvector of mode r (length n_dof)
/// * `m_diag` – diagonal mass vector (length n_dof)
/// * `influence_vec` – seismic influence vector (length n_dof)
pub fn modal_participation_factor(
    mode_shape: &[f64],
    m_diag: &[f64],
    influence_vec: &[f64],
) -> f64 {
    let n = mode_shape.len().min(m_diag.len()).min(influence_vec.len());
    let numerator: f64 = (0..n)
        .map(|i| mode_shape[i] * m_diag[i] * influence_vec[i])
        .sum();
    let denominator: f64 = (0..n)
        .map(|i| mode_shape[i] * m_diag[i] * mode_shape[i])
        .sum();
    if denominator.abs() < f64::EPSILON {
        return 0.0;
    }
    numerator / denominator
}
/// Compute the effective modal mass for mode r.
///
/// M_eff_r = (phi_r^T * M * r)² / (phi_r^T * M * phi_r)
///
/// The sum of all effective modal masses equals the total mass.
pub fn effective_modal_mass(mode_shape: &[f64], m_diag: &[f64], influence_vec: &[f64]) -> f64 {
    let n = mode_shape.len().min(m_diag.len()).min(influence_vec.len());
    let numerator: f64 = (0..n)
        .map(|i| mode_shape[i] * m_diag[i] * influence_vec[i])
        .sum();
    let denominator: f64 = (0..n)
        .map(|i| mode_shape[i] * m_diag[i] * mode_shape[i])
        .sum();
    if denominator.abs() < f64::EPSILON {
        return 0.0;
    }
    numerator * numerator / denominator
}
/// Combine peak modal responses using the SRSS rule.
///
/// u_SRSS = sqrt( sum_r (u_r²) )
///
/// where u_r are the peak responses of individual modes.
pub fn srss_combination(modal_peaks: &[f64]) -> f64 {
    modal_peaks.iter().map(|v| v * v).sum::<f64>().sqrt()
}
/// Combine peak modal responses using the CQC rule.
///
/// u_CQC = sqrt( sum_r sum_s rho_rs * u_r * u_s )
///
/// The cross-correlation coefficient (Der Kiureghian, 1981):
/// rho_rs = 8 * zeta^2 * (1 + beta_rs) * beta_rs^(3/2)
///          / ((1 - beta_rs²)² + 4*zeta²*beta_rs*(1+beta_rs)²)
///
/// where beta_rs = omega_r / omega_s.
///
/// # Arguments
/// * `modal_peaks` – peak responses for each mode
/// * `frequencies` – natural frequencies (rad/s) for each mode
/// * `zeta` – uniform damping ratio
pub fn cqc_combination(modal_peaks: &[f64], frequencies: &[f64], zeta: f64) -> f64 {
    let n = modal_peaks.len().min(frequencies.len());
    if n == 0 {
        return 0.0;
    }
    let mut sum = 0.0;
    for r in 0..n {
        for s in 0..n {
            let beta = if frequencies[s].abs() < f64::EPSILON {
                1.0
            } else {
                frequencies[r] / frequencies[s]
            };
            let zz = zeta * zeta;
            let num = 8.0 * zz * (1.0 + beta) * beta.powf(1.5);
            let den = (1.0 - beta * beta).powi(2) + 4.0 * zz * beta * (1.0 + beta).powi(2);
            let rho = if den.abs() < f64::EPSILON {
                if r == s { 1.0 } else { 0.0 }
            } else {
                num / den
            };
            sum += rho * modal_peaks[r] * modal_peaks[s];
        }
    }
    sum.sqrt()
}
/// Compute the damped natural frequency ω_d = ω_n * sqrt(1 - ζ²).
pub fn damped_frequency(omega_n: f64, zeta: f64) -> f64 {
    omega_n * (1.0 - zeta * zeta).max(0.0).sqrt()
}
/// Estimate the damping ratio from logarithmic decrement.
///
/// δ = ln(u_1 / u_{1+n}) / n
/// ζ = δ / sqrt(4π² + δ²)
///
/// # Arguments
/// * `u1` – amplitude at cycle 1
/// * `u1n` – amplitude at cycle 1+n (n cycles later)
/// * `n_cycles` – number of cycles between measurements
pub fn damping_ratio_from_log_decrement(u1: f64, u1n: f64, n_cycles: usize) -> f64 {
    if u1n.abs() < f64::EPSILON || n_cycles == 0 {
        return 0.0;
    }
    let delta = (u1 / u1n).abs().ln() / n_cycles as f64;
    let pi2 = 4.0 * std::f64::consts::PI * std::f64::consts::PI;
    delta / (pi2 + delta * delta).sqrt()
}
#[cfg(test)]
mod tests_dynamic_expanded {
    use super::*;
    use crate::dynamic::*;
    #[test]
    fn test_modal_newmark_zero_force() {
        let force = vec![0.0; 50];
        let q = modal_newmark_1dof(10.0, 0.0, 1.0, &force, 0.001);
        for (i, &v) in q.iter().enumerate() {
            assert!(v.abs() < 1e-12, "q[{i}] = {v}");
        }
    }
    #[test]
    fn test_modal_newmark_response_finite() {
        let omega_n = 2.0 * std::f64::consts::PI * 5.0;
        let omega_f = 2.0 * std::f64::consts::PI * 2.0;
        let force: Vec<f64> = (0..500)
            .map(|i| (omega_f * i as f64 * 0.001).sin())
            .collect();
        let q = modal_newmark_1dof(omega_n, 0.05, 1.0, &force, 0.001);
        assert_eq!(q.len(), 500);
        let max_q = q.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        assert!(max_q.is_finite(), "max q = {max_q}");
        assert!(max_q > 0.0, "should have nonzero response to nonzero force");
    }
    #[test]
    fn test_half_sine_pulse_peak() {
        let a_max = 9.81;
        let t_d = 0.1;
        let dt = 0.001;
        let acc = half_sine_pulse(a_max, t_d, dt, 0.5);
        let peak = acc.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        assert!((peak - a_max).abs() / a_max < 0.01, "peak = {peak}");
    }
    #[test]
    fn test_half_sine_pulse_zero_after_duration() {
        let a_max = 9.81;
        let t_d = 0.1;
        let dt = 0.001;
        let acc = half_sine_pulse(a_max, t_d, dt, 0.5);
        let n_td = (t_d / dt) as usize + 2;
        for (i, &a) in acc.iter().enumerate().skip(n_td) {
            assert!(
                a.abs() < 1e-10,
                "acc[{i}] = {} should be zero after pulse",
                a
            );
        }
    }
    #[test]
    fn test_triangular_pulse_peak_at_midpoint() {
        let a_max = 5.0;
        let t_d = 0.2;
        let dt = 0.001;
        let acc = triangular_pulse(a_max, t_d, dt, 0.5);
        let peak_idx = acc
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        let t_peak = peak_idx as f64 * dt;
        assert!(
            (t_peak - t_d / 2.0).abs() < dt * 2.0,
            "peak should be at t_d/2 = {}, got {}",
            t_d / 2.0,
            t_peak
        );
    }
    #[test]
    fn test_rectangular_pulse_constant() {
        let a_max = 3.0;
        let t_d = 0.1;
        let dt = 0.01;
        let acc = rectangular_pulse(a_max, t_d, dt, 0.3);
        let n_td = (t_d / dt) as usize;
        for (i, &a) in acc.iter().enumerate().take(n_td + 1) {
            assert!((a - a_max).abs() < 1e-10, "acc[{i}] = {}", a);
        }
    }
    #[test]
    fn test_triangular_pulse_zero_at_boundaries() {
        let acc = triangular_pulse(5.0, 0.2, 0.001, 0.5);
        assert!(acc[0].abs() < 1e-10, "acc[0] = {}", acc[0]);
    }
    #[test]
    fn test_response_spectrum_rigid_system() {
        let acc = vec![1.0, -2.0, 3.0, -1.5];
        let rs = response_spectrum(&acc, 0.01, &[1e-9], 0.05);
        assert_eq!(rs.len(), 1);
        let (_, sa) = rs[0];
        assert!((sa - 3.0).abs() < 1e-10, "rigid Sa = {sa}");
    }
    #[test]
    fn test_response_spectrum_positive_sa() {
        let acc: Vec<f64> = (0..200).map(|i| (i as f64 * 0.05).sin()).collect();
        let periods = vec![0.1, 0.5, 1.0];
        let rs = response_spectrum(&acc, 0.005, &periods, 0.05);
        assert_eq!(rs.len(), 3);
        for &(_t, sa) in &rs {
            assert!(sa >= 0.0, "Sa should be non-negative");
            assert!(sa.is_finite(), "Sa should be finite");
        }
    }
    #[test]
    fn test_dynamic_sif_zero_velocity() {
        let k1_dyn = dynamic_sif_freund(1e6, 0.0, 3000.0, 5000.0);
        assert!((k1_dyn - 1e6).abs() / 1e6 < 1e-6, "K_dyn(v=0) = {k1_dyn}");
    }
    #[test]
    fn test_dynamic_sif_decreases_with_velocity() {
        let k1_dyn_slow = dynamic_sif_freund(1e6, 100.0, 3000.0, 5000.0);
        let k1_dyn_fast = dynamic_sif_freund(1e6, 2000.0, 3000.0, 5000.0);
        assert!(
            k1_dyn_fast < k1_dyn_slow,
            "dynamic SIF should decrease with crack velocity"
        );
    }
    #[test]
    fn test_dynamic_energy_release_rate_positive() {
        let g = dynamic_energy_release_rate(30e6, 200e9, 0.3);
        assert!(g > 0.0 && g.is_finite(), "G_dyn = {g}");
    }
    #[test]
    fn test_modal_participation_factor_normalized() {
        let phi = vec![1.0, 1.0, 1.0];
        let m = vec![1.0, 1.0, 1.0];
        let r = vec![1.0, 1.0, 1.0];
        let gamma = modal_participation_factor(&phi, &m, &r);
        assert!((gamma - 1.0).abs() < 1e-12, "Γ = {gamma}");
    }
    #[test]
    fn test_effective_modal_mass_sum_equals_total() {
        let phi = vec![1.0];
        let m = vec![5.0];
        let r = vec![1.0];
        let m_eff = effective_modal_mass(&phi, &m, &r);
        assert!((m_eff - 5.0).abs() < 1e-10, "M_eff = {m_eff}");
    }
    #[test]
    fn test_srss_combination_pythagoras() {
        let peaks = vec![3.0, 4.0];
        let srss = srss_combination(&peaks);
        assert!((srss - 5.0).abs() < 1e-10, "SRSS = {srss}");
    }
    #[test]
    fn test_cqc_combination_single_mode_equals_peak() {
        let peaks = vec![5.0];
        let freqs = vec![10.0];
        let cqc = cqc_combination(&peaks, &freqs, 0.05);
        assert!((cqc - 5.0).abs() < 1e-10, "CQC single mode = {cqc}");
    }
    #[test]
    fn test_cqc_combination_well_separated_modes() {
        let peaks = vec![3.0, 4.0];
        let freqs = vec![1.0, 100.0];
        let cqc = cqc_combination(&peaks, &freqs, 0.05);
        let srss = srss_combination(&peaks);
        assert!(
            (cqc - srss).abs() / srss < 0.05,
            "CQC = {cqc}, SRSS = {srss}"
        );
    }
    #[test]
    fn test_damped_frequency_underdamped() {
        let omega_n = 10.0;
        let zeta = 0.1;
        let wd = damped_frequency(omega_n, zeta);
        let expected = omega_n * (1.0 - 0.01_f64).sqrt();
        assert!((wd - expected).abs() < 1e-12, "wd = {wd}");
    }
    #[test]
    fn test_damped_frequency_undamped_equals_natural() {
        let wd = damped_frequency(10.0, 0.0);
        assert!((wd - 10.0).abs() < 1e-12, "wd(zeta=0) = {wd}");
    }
    #[test]
    fn test_log_decrement_round_trip() {
        let zeta = 0.05;
        let omega_n = 10.0;
        let wd = damped_frequency(omega_n, zeta);
        let n_cyc = 3usize;
        let t_n = 2.0 * std::f64::consts::PI / wd * n_cyc as f64;
        let u1 = 1.0;
        let u1n = (-zeta * omega_n * t_n).exp();
        let zeta_recovered = damping_ratio_from_log_decrement(u1, u1n, n_cyc);
        assert!(
            (zeta_recovered - zeta).abs() < 1e-5,
            "zeta = {zeta_recovered}"
        );
    }
    #[test]
    fn test_newmark_linear_accel_critical_dt() {
        let omega_max = 100.0;
        let dt_crit = NewmarkLinearAccelerationParams::critical_dt(omega_max);
        let expected = 3.0_f64.sqrt() / 100.0;
        assert!((dt_crit - expected).abs() < 1e-12, "dt_crit = {dt_crit}");
    }
    #[test]
    fn test_rayleigh_from_two_modes() {
        let omega1 = 5.0;
        let zeta1 = 0.03;
        let omega2 = 50.0;
        let zeta2 = 0.05;
        let (alpha, beta) = rayleigh_from_two_modes(omega1, zeta1, omega2, zeta2);
        let zeta_check1 = alpha / (2.0 * omega1) + beta * omega1 / 2.0;
        let zeta_check2 = alpha / (2.0 * omega2) + beta * omega2 / 2.0;
        assert!(
            (zeta_check1 - zeta1).abs() < 1e-10,
            "zeta1 check = {zeta_check1}"
        );
        assert!(
            (zeta_check2 - zeta2).abs() < 1e-10,
            "zeta2 check = {zeta_check2}"
        );
    }
    #[test]
    fn test_generalized_alpha_step_output_length() {
        let n = 4;
        let mut solver = DynamicSolver::new(n);
        let triplets: Vec<(usize, usize, f64)> = (0..n).map(|i| (i, i, 10.0)).collect();
        let mass = CsrMatrix::from_triplets(n, n, &triplets);
        let stiff = CsrMatrix::from_triplets(n, n, &triplets);
        let f = vec![1.0_f64; n];
        solver.generalized_alpha(&mass, &stiff, &f, 0.01, 0.9);
        assert_eq!(solver.u.len(), n);
        assert_eq!(solver.v.len(), n);
        assert_eq!(solver.a.len(), n);
    }
    #[test]
    fn test_generalized_alpha_zero_force_stays_zero() {
        let n = 3;
        let mut solver = DynamicSolver::new(n);
        let triplets: Vec<(usize, usize, f64)> = (0..n).map(|i| (i, i, 5.0)).collect();
        let mass = CsrMatrix::from_triplets(n, n, &triplets);
        let stiff = CsrMatrix::from_triplets(n, n, &triplets);
        let f = vec![0.0_f64; n];
        solver.generalized_alpha(&mass, &stiff, &f, 0.01, 0.8);
        for i in 0..n {
            assert!(solver.u[i].abs() < 1e-14, "u[{i}] = {}", solver.u[i]);
        }
    }
    #[test]
    fn test_generalized_alpha_rho_inf_zero_maximum_dissipation() {
        let n = 2;
        let mut solver = DynamicSolver::new(n);
        let m_tri = vec![(0, 0, 2.0), (1, 1, 2.0)];
        let k_tri = vec![(0, 0, 4.0), (0, 1, -1.0), (1, 0, -1.0), (1, 1, 4.0)];
        let mass = CsrMatrix::from_triplets(n, n, &m_tri);
        let stiff = CsrMatrix::from_triplets(n, n, &k_tri);
        let f = vec![1.0_f64, 0.0];
        solver.generalized_alpha(&mass, &stiff, &f, 0.01, 0.0);
        for &ui in &solver.u {
            assert!(ui.is_finite(), "u = {ui}");
        }
    }
    #[test]
    fn test_rayleigh_damping_diagonal_positive() {
        let n = 3;
        let triplets: Vec<(usize, usize, f64)> = (0..n).map(|i| (i, i, 10.0)).collect();
        let mass = CsrMatrix::from_triplets(n, n, &triplets);
        let stiff = CsrMatrix::from_triplets(n, n, &triplets);
        let (alpha, beta) = compute_rayleigh_coefficients(2.0, 0.05, 20.0, 0.05);
        let c = rayleigh_damping_matrix(&mass, &stiff, alpha, beta);
        for i in 0..n {
            let row_start = c.row_ptr[i];
            let row_end = c.row_ptr[i + 1];
            let diag = (row_start..row_end)
                .find(|&idx| c.col_indices[idx] == i)
                .map(|idx| c.values[idx])
                .unwrap_or(0.0);
            assert!(diag > 0.0, "C[{i},{i}] = {diag}");
        }
    }
    #[test]
    fn test_rayleigh_damping_scales_with_alpha() {
        let n = 2;
        let triplets: Vec<(usize, usize, f64)> = (0..n).map(|i| (i, i, 5.0)).collect();
        let mass = CsrMatrix::from_triplets(n, n, &triplets);
        let stiff = CsrMatrix::from_triplets(n, n, &triplets);
        let c1 = rayleigh_damping_matrix(&mass, &stiff, 1.0, 0.0);
        let c2 = rayleigh_damping_matrix(&mass, &stiff, 2.0, 0.0);
        let sum1: f64 = c1.values.iter().sum();
        let sum2: f64 = c2.values.iter().sum();
        assert!(
            (sum2 / sum1 - 2.0).abs() < 1e-10,
            "sum2/sum1 = {}",
            sum2 / sum1
        );
    }
    #[test]
    fn test_rayleigh_coefficients_round_trip() {
        let omega1 = 5.0;
        let omega2 = 50.0;
        let zeta = 0.04;
        let (alpha, beta) = compute_rayleigh_coefficients(omega1, zeta, omega2, zeta);
        let zeta_check = alpha / (2.0 * omega1) + beta * omega1 / 2.0;
        assert!(
            (zeta_check - zeta).abs() < 1e-10,
            "zeta check = {zeta_check}"
        );
    }
    #[test]
    fn test_modal_truncation_correction_zero_residual() {
        let omega_n = [10.0_f64, 20.0, 30.0];
        let phi: Vec<Vec<f64>> = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let f = vec![1.0_f64, 2.0, 3.0];
        let corr = modal_truncation_correction(&omega_n, &phi, &f, 3);
        for &c in &corr {
            assert!(c.abs() < 1e-10, "correction should be near-zero: {c}");
        }
    }
    #[test]
    fn test_modal_truncation_correction_nonzero_when_truncated() {
        let omega_n = [10.0_f64, 20.0, 300.0];
        let phi: Vec<Vec<f64>> = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let f = vec![0.0_f64, 0.0, 5.0];
        let corr = modal_truncation_correction(&omega_n, &phi, &f, 2);
        let norm: f64 = corr.iter().map(|&v| v * v).sum::<f64>().sqrt();
        assert!(
            norm > 0.0,
            "truncation correction should be non-zero: {norm}"
        );
    }
    #[test]
    fn test_modal_truncation_correction_finite() {
        let omega_n = [5.0_f64, 15.0, 50.0];
        let phi: Vec<Vec<f64>> = vec![
            vec![
                std::f64::consts::FRAC_1_SQRT_2,
                std::f64::consts::FRAC_1_SQRT_2,
                0.0,
            ],
            vec![
                -std::f64::consts::FRAC_1_SQRT_2,
                std::f64::consts::FRAC_1_SQRT_2,
                0.0,
            ],
            vec![0.0, 0.0, 1.0],
        ];
        let f = vec![1.0_f64, 0.5, 0.25];
        let corr = modal_truncation_correction(&omega_n, &phi, &f, 2);
        for &c in &corr {
            assert!(c.is_finite(), "correction must be finite: {c}");
        }
    }
}
/// Compute Rayleigh coefficients (alpha, beta) targeting two modes.
///
/// Solves the 2×2 system:
/// ```text
/// zeta1 = alpha/(2*omega1) + beta*omega1/2
/// zeta2 = alpha/(2*omega2) + beta*omega2/2
/// ```
pub fn compute_rayleigh_coefficients(
    omega1: f64,
    zeta1: f64,
    omega2: f64,
    zeta2: f64,
) -> (f64, f64) {
    let a11 = 0.5 / omega1;
    let a12 = 0.5 * omega1;
    let a21 = 0.5 / omega2;
    let a22 = 0.5 * omega2;
    let det = a11 * a22 - a12 * a21;
    if det.abs() < 1e-60 {
        return (0.0, 0.0);
    }
    let alpha = (a22 * zeta1 - a12 * zeta2) / det;
    let beta = (a11 * zeta2 - a21 * zeta1) / det;
    (alpha, beta)
}
/// Modal truncation correction (static correction for modes not included).
///
/// The residual flexibility accounts for the contribution of high-frequency
/// (truncated) modes via:
///   u_residual = Σ_{r=n_inc+1}^{N} φ_r * (φ_r^T F) / ω_r²
///
/// where the sum is over excluded modes.
///
/// # Arguments
/// * `omega_n`           – natural frequencies [rad/s], length N
/// * `phi`               – mode shapes, phi\[r\] is the r-th mode vector
/// * `force`             – applied force vector
/// * `n_modes_included`  – number of modes already used in the solution
///
/// Returns the static correction displacement vector.
pub fn modal_truncation_correction(
    omega_n: &[f64],
    phi: &[Vec<f64>],
    force: &[f64],
    n_modes_included: usize,
) -> Vec<f64> {
    let n_modes = omega_n.len().min(phi.len());
    let n_dof = force.len();
    let mut corr = vec![0.0_f64; n_dof];
    for r in n_modes_included..n_modes {
        let omega2 = omega_n[r] * omega_n[r];
        if omega2 < 1e-60 {
            continue;
        }
        let mode = &phi[r];
        if mode.len() != n_dof {
            continue;
        }
        let q_r: f64 = mode
            .iter()
            .zip(force.iter())
            .map(|(phi_i, f_i)| phi_i * f_i)
            .sum();
        let coeff = q_r / omega2;
        for i in 0..n_dof {
            corr[i] += coeff * mode[i];
        }
    }
    corr
}
/// Perform modal superposition time integration for a multi-DOF system.
///
/// Each mode is treated as an independent SDOF oscillator.  The modal
/// responses are combined via coordinate transformation `u = Φ q`.
///
/// # Arguments
/// * `omega_n`       – natural frequencies [rad/s]
/// * `zeta_modes`    – modal damping ratios (one per mode)
/// * `phi`           – mode shapes: `phi\[mode\]\[dof\]`
/// * `f_time`        – external force history: `f_time\[step\]\[dof\]`
/// * `dt`            – time step
/// * `u0`            – initial displacements
/// * `v0`            – initial velocities
pub fn modal_superposition(
    omega_n: &[f64],
    zeta_modes: &[f64],
    phi: &[Vec<f64>],
    f_time: &[Vec<f64>],
    dt: f64,
    u0: &[f64],
    v0: &[f64],
) -> ModalSuperpositionResult {
    let n_modes = omega_n.len().min(phi.len()).min(zeta_modes.len());
    let n_steps = f_time.len();
    let n_dof = if n_modes > 0 { phi[0].len() } else { 0 };
    let mut q = vec![0.0_f64; n_modes];
    let mut dq = vec![0.0_f64; n_modes];
    for r in 0..n_modes {
        let mode = &phi[r];
        q[r] = mode.iter().zip(u0.iter()).map(|(p, u)| p * u).sum();
        dq[r] = mode.iter().zip(v0.iter()).map(|(p, v)| p * v).sum();
    }
    let mut disp_history: Vec<Vec<f64>> = Vec::with_capacity(n_steps);
    let mut modal_history: Vec<Vec<f64>> = Vec::with_capacity(n_steps);
    for f in f_time {
        for r in 0..n_modes {
            let mode = &phi[r];
            let omega = omega_n[r];
            let zeta = zeta_modes[r];
            let omega_d = (omega * omega * (1.0 - zeta * zeta)).max(0.0).sqrt();
            let _ = omega_d;
            let f_modal: f64 = mode.iter().zip(f.iter()).map(|(p, fi)| p * fi).sum();
            let beta = 0.25_f64;
            let gamma = 0.5_f64;
            let a_old = -2.0 * zeta * omega * dq[r] - omega * omega * q[r] + f_modal;
            let k_eff =
                1.0 / (beta * dt * dt) + gamma / (beta * dt) * 2.0 * zeta * omega + omega * omega;
            let f_eff = f_modal
                + 1.0 / (beta * dt * dt) * q[r]
                + 1.0 / (beta * dt) * dq[r]
                + (1.0 / (2.0 * beta) - 1.0) * a_old
                + 2.0
                    * zeta
                    * omega
                    * (gamma / (beta * dt) * q[r]
                        + (gamma / beta - 1.0) * dq[r]
                        + dt * (gamma / (2.0 * beta) - 1.0) * a_old);
            let q_new = if k_eff.abs() < 1e-60 {
                q[r]
            } else {
                f_eff / k_eff
            };
            let dq_new = gamma / (beta * dt) * (q_new - q[r])
                - (gamma / beta - 1.0) * dq[r]
                - dt * (gamma / (2.0 * beta) - 1.0) * a_old;
            q[r] = q_new;
            dq[r] = dq_new;
        }
        let mut u_phys = vec![0.0_f64; n_dof];
        for r in 0..n_modes {
            for j in 0..n_dof {
                u_phys[j] += phi[r][j] * q[r];
            }
        }
        disp_history.push(u_phys);
        modal_history.push(q.clone());
    }
    ModalSuperpositionResult {
        disp: disp_history,
        modal_coords: modal_history,
    }
}
/// Perform a complete response spectrum analysis (SRSS combination).
///
/// Returns peak displacements and accelerations for each DOF using SRSS.
pub fn response_spectrum_analysis(input: &ResponseSpectrumInput) -> (Vec<f64>, Vec<f64>) {
    let n_modes = input.omega_n.len();
    let n_dof = if n_modes > 0 { input.phi[0].len() } else { 0 };
    let mut peak_u2 = vec![0.0_f64; n_dof];
    let mut peak_a2 = vec![0.0_f64; n_dof];
    for r in 0..n_modes
        .min(input.phi.len())
        .min(input.sd.len())
        .min(input.sa.len())
    {
        let gamma_r = if r < input.gamma.len() {
            input.gamma[r]
        } else {
            1.0
        };
        let sd_r = input.sd[r];
        let sa_r = input.sa[r];
        for j in 0..n_dof {
            let u_r = gamma_r * input.phi[r][j] * sd_r;
            let a_r = gamma_r * input.phi[r][j] * sa_r;
            peak_u2[j] += u_r * u_r;
            peak_a2[j] += a_r * a_r;
        }
    }
    let peak_u: Vec<f64> = peak_u2.iter().map(|v| v.sqrt()).collect();
    let peak_a: Vec<f64> = peak_a2.iter().map(|v| v.sqrt()).collect();
    (peak_u, peak_a)
}
/// Compute spectral displacement Sd from a ground acceleration record.
///
/// Uses Newmark-beta average acceleration to integrate each SDOF oscillator.
///
/// # Arguments
/// * `ag`    – ground acceleration time history \[m/s²\]
/// * `dt`    – time step \[s\]
/// * `omega` – natural frequency \[rad/s\]
/// * `zeta`  – damping ratio
///
/// Returns the peak absolute displacement Sd.
pub fn compute_sd(ag: &[f64], dt: f64, omega: f64, zeta: f64) -> f64 {
    let beta = 0.25_f64;
    let gamma_nb = 0.5_f64;
    let mut u = 0.0_f64;
    let mut v = 0.0_f64;
    let mut a = 0.0_f64;
    let mut sd = 0.0_f64;
    let k_eff =
        omega * omega + 2.0 * zeta * omega * gamma_nb / (beta * dt) + 1.0 / (beta * dt * dt);
    for &ag_i in ag {
        let f_eff = -ag_i
            + u / (beta * dt * dt)
            + v / (beta * dt)
            + (1.0 / (2.0 * beta) - 1.0) * a
            + 2.0
                * zeta
                * omega
                * (gamma_nb / (beta * dt) * u
                    + (gamma_nb / beta - 1.0) * v
                    + dt * (gamma_nb / (2.0 * beta) - 1.0) * a);
        let u_new = if k_eff.abs() < 1e-60 {
            u
        } else {
            f_eff / k_eff
        };
        let v_new = gamma_nb / (beta * dt) * (u_new - u)
            - (gamma_nb / beta - 1.0) * v
            - dt * (gamma_nb / (2.0 * beta) - 1.0) * a;
        let a_new =
            (u_new - u) / (beta * dt * dt) - v / (beta * dt) - (1.0 / (2.0 * beta) - 1.0) * a;
        u = u_new;
        v = v_new;
        a = a_new;
        if u.abs() > sd {
            sd = u.abs();
        }
    }
    sd
}
/// Compute spectral acceleration Sa from a ground acceleration record.
///
/// Returns the peak absolute acceleration Sa = ω² * Sd (pseudo-acceleration).
pub fn compute_sa(ag: &[f64], dt: f64, omega: f64, zeta: f64) -> f64 {
    let sd = compute_sd(ag, dt, omega, zeta);
    omega * omega * sd
}
/// Compute the mass-proportional damping coefficient α for a target modal damping.
///
/// For mass-proportional damping C = α M:
///   ζ_r = α / (2 ω_r)  →  α = 2 ζ_r ω_r
///
/// # Arguments
/// * `omega_r` – reference natural frequency [rad/s]
/// * `zeta_r`  – target damping ratio at `omega_r`
pub fn mass_proportional_damping_alpha(omega_r: f64, zeta_r: f64) -> f64 {
    2.0 * zeta_r * omega_r
}
/// Compute the stiffness-proportional damping coefficient β for a target modal damping.
///
/// For stiffness-proportional damping C = β K:
///   ζ_r = β ω_r / 2  →  β = 2 ζ_r / ω_r
///
/// # Arguments
/// * `omega_r` – reference natural frequency [rad/s]
/// * `zeta_r`  – target damping ratio at `omega_r`
pub fn stiffness_proportional_damping_beta(omega_r: f64, zeta_r: f64) -> f64 {
    2.0 * zeta_r / omega_r
}
/// Evaluate the effective damping ratio for a mode at frequency `omega`
/// given Rayleigh coefficients (α, β).
///
/// ζ(ω) = α / (2ω) + β ω / 2
pub fn rayleigh_modal_damping(omega: f64, alpha: f64, beta: f64) -> f64 {
    alpha / (2.0 * omega) + beta * omega / 2.0
}
/// Apply mass-proportional damping to a force vector (equivalent to C·v for diagonal M).
///
/// Returns the damping force vector d = α * M_diag * v.
pub fn mass_proportional_damping_force(m_diag: &[f64], v: &[f64], alpha: f64) -> Vec<f64> {
    m_diag
        .iter()
        .zip(v.iter())
        .map(|(m, vi)| alpha * m * vi)
        .collect()
}
/// Apply stiffness-proportional damping to a force vector: d = β * K_diag * v.
pub fn stiffness_proportional_damping_force(k_diag: &[f64], v: &[f64], beta: f64) -> Vec<f64> {
    k_diag
        .iter()
        .zip(v.iter())
        .map(|(k, vi)| beta * k * vi)
        .collect()
}
/// Cross-correlation coefficient ρ_ij for the CQC combination rule.
///
/// Uses the standard Wilson (1981) formula.
pub fn cqc_rho(omega_i: f64, omega_j: f64, zeta_i: f64, zeta_j: f64) -> f64 {
    let r = omega_j / omega_i;
    let zeta_sum = zeta_i + zeta_j;
    let num = 8.0 * (zeta_i * zeta_j).sqrt() * zeta_sum * r.powf(1.5);
    let den = (1.0 - r * r).powi(2)
        + 4.0 * zeta_i * zeta_j * r * (1.0 + r * r)
        + 4.0 * (zeta_i * zeta_i + zeta_j * zeta_j) * r * r;
    if den.abs() < 1e-30 { 1.0 } else { num / den }
}
/// CQC combination with individual modal damping ratios.
///
/// Returns the combined peak response from `modal_peaks` using the CQC rule.
pub fn cqc_combination_individual_zeta(modal_peaks: &[f64], omega: &[f64], zeta: &[f64]) -> f64 {
    let n = modal_peaks.len().min(omega.len()).min(zeta.len());
    let mut sum = 0.0_f64;
    for i in 0..n {
        for j in 0..n {
            let rho = cqc_rho(omega[i], omega[j], zeta[i], zeta[j]);
            sum += rho * modal_peaks[i] * modal_peaks[j];
        }
    }
    sum.max(0.0).sqrt()
}
/// Simulate free vibration decay for a damped SDOF system.
///
/// Returns displacement time history starting from u0=1, v0=0.
///
/// The exact analytical solution is:
///   u(t) = e^{-ζ ω_n t} [cos(ω_d t) + (ζ/√(1-ζ²)) sin(ω_d t)]
pub fn free_vibration_decay(omega_n: f64, zeta: f64, dt: f64, n_steps: usize) -> Vec<f64> {
    let omega_d = omega_n * (1.0 - zeta * zeta).max(0.0).sqrt();
    (0..n_steps)
        .map(|i| {
            let t = i as f64 * dt;
            let env = (-zeta * omega_n * t).exp();
            if omega_d > 1e-15 {
                env * (omega_d * t).cos()
                    + env * zeta / (1.0 - zeta * zeta).max(1e-30).sqrt() * (omega_d * t).sin()
            } else {
                env * (1.0 + omega_n * t)
            }
        })
        .collect()
}
/// Estimate damping ratio from successive peaks of a free vibration record.
///
/// Uses logarithmic decrement: δ = ln(u1/u2), ζ = δ / √(4π² + δ²).
pub fn damping_from_successive_peaks(u1: f64, u2: f64) -> f64 {
    if u2.abs() < 1e-60 || u1.abs() < 1e-60 {
        return 0.0;
    }
    let delta = (u1 / u2).abs().ln();
    delta / (4.0 * std::f64::consts::PI * std::f64::consts::PI + delta * delta).sqrt()
}
/// Compute the unit impulse response function h(t) for a damped SDOF.
///
/// h(t) = (1 / (m ω_d)) e^{-ζ ω_n t} sin(ω_d t)  for t ≥ 0
///
/// # Arguments
/// * `m`       – modal mass
/// * `omega_n` – natural frequency \[rad/s\]
/// * `zeta`    – damping ratio
/// * `dt`      – time step \[s\]
/// * `n`       – number of samples
pub fn impulse_response_function(m: f64, omega_n: f64, zeta: f64, dt: f64, n: usize) -> Vec<f64> {
    let omega_d = omega_n * (1.0 - zeta * zeta).max(0.0).sqrt();
    let denom = m * omega_d;
    (0..n)
        .map(|i| {
            let t = i as f64 * dt;
            if denom.abs() < 1e-60 {
                0.0
            } else {
                (-zeta * omega_n * t).exp() * (omega_d * t).sin() / denom
            }
        })
        .collect()
}
/// Compute the numerical energy ratio E_{n+1} / E_n for monitoring dissipation.
///
/// Returns a value ≤ 1 for stable (dissipative) schemes, ≈ 1 for conservative.
pub fn numerical_energy_ratio(
    u_new: &[f64],
    v_new: &[f64],
    u_old: &[f64],
    v_old: &[f64],
    m_diag: &[f64],
    k_diag: &[f64],
) -> f64 {
    let ke_new: f64 = v_new
        .iter()
        .zip(m_diag.iter())
        .map(|(v, m)| 0.5 * m * v * v)
        .sum();
    let pe_new: f64 = u_new
        .iter()
        .zip(k_diag.iter())
        .map(|(u, k)| 0.5 * k * u * u)
        .sum();
    let ke_old: f64 = v_old
        .iter()
        .zip(m_diag.iter())
        .map(|(v, m)| 0.5 * m * v * v)
        .sum();
    let pe_old: f64 = u_old
        .iter()
        .zip(k_diag.iter())
        .map(|(u, k)| 0.5 * k * u * u)
        .sum();
    let e_old = ke_old + pe_old;
    if e_old < 1e-60 {
        1.0
    } else {
        (ke_new + pe_new) / e_old
    }
}
