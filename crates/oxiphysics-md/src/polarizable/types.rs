//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// Thole short-range damping functions for dipole-dipole interactions.
pub struct TholeDamping;
impl TholeDamping {
    /// Thole screening function.
    pub fn thole_screening(u: f64, a: f64) -> f64 {
        let x = a * u * u * u;
        1.0 - (-x).exp() * (1.0 + x + 0.5 * x * x)
    }
    /// Dimensionless distance variable.
    pub fn u_ij(r: f64, alpha_i: f64, alpha_j: f64) -> f64 {
        let denom = (alpha_i * alpha_j).powf(1.0 / 6.0);
        if denom == 0.0 { 0.0 } else { r / denom }
    }
    /// Damped dipole interaction tensor.
    pub fn damped_dipole_tensor(
        r_vec: [f64; 3],
        alpha_i: f64,
        alpha_j: f64,
        a: f64,
    ) -> [[f64; 3]; 3] {
        let r = norm(r_vec);
        if r < 1e-12 {
            return [[0.0; 3]; 3];
        }
        let r2 = r * r;
        let r5 = r2 * r2 * r;
        let u = Self::u_ij(r, alpha_i, alpha_j);
        let s = Self::thole_screening(u, a);
        let mut t = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let delta_ij = if i == j { 1.0 } else { 0.0 };
                let bare = (3.0 * r_vec[i] * r_vec[j] - r2 * delta_ij) / r5;
                t[i][j] = s * bare;
            }
        }
        t
    }
    /// Linear Thole damping function (alternative to exponential).
    ///
    /// s_lin(u) = 4*a*u^3 - 3*a*u^4  for u < 1/a^(1/3)
    /// s_lin(u) = 1                    otherwise
    pub fn thole_screening_linear(u: f64, a: f64) -> f64 {
        let threshold = (1.0 / a).powf(1.0 / 3.0);
        if u >= threshold {
            1.0
        } else {
            let au3 = a * u * u * u;
            let au4 = au3 * u;
            4.0 * au3 - 3.0 * au4
        }
    }
}
impl TholeDamping {
    /// Thole-smeared (effective) charge at distance `r`.
    ///
    /// The Thole model replaces a point charge with a smeared charge
    /// distribution.  The effective charge at separation `r` is scaled
    /// by the Thole screening factor λ:
    ///
    /// ```text
    /// q_smeared(r) = q * λ(u)
    /// ```
    ///
    /// where `u = r / (α_i α_j)^{1/6}` and
    /// λ(u) = 1 − exp(−a·u³)(1 + a·u³ + ½(a·u³)²)
    ///
    /// # Arguments
    /// * `q`       – bare charge (e).
    /// * `r`       – intersite distance (Å).
    /// * `alpha_i` – polarizability of site i (Å³).
    /// * `alpha_j` – polarizability of site j (Å³).
    /// * `a`       – Thole damping parameter (dimensionless, typically 2.0–2.6).
    pub fn compute_smeared_charge(q: f64, r: f64, alpha_i: f64, alpha_j: f64, a: f64) -> f64 {
        let u = Self::u_ij(r, alpha_i, alpha_j);
        let lambda = Self::thole_screening(u, a);
        q * lambda
    }
    /// Smeared Thole charge gradient with respect to distance `r`.
    ///
    /// dq_smeared/dr = q * dλ/dr
    ///
    /// Using the chain rule and the analytic form of λ(u):
    /// dλ/du = a² u^4 exp(-a u³) / 2
    /// du/dr = 1 / (α_i α_j)^{1/6}
    pub fn compute_smeared_charge_gradient(
        q: f64,
        r: f64,
        alpha_i: f64,
        alpha_j: f64,
        a: f64,
    ) -> f64 {
        let denom = (alpha_i * alpha_j).powf(1.0 / 6.0);
        if denom < 1e-30 {
            return 0.0;
        }
        let u = r / denom;
        let au3 = a * u * u * u;
        let dl_du = 0.5 * a * a * u * u * u * u * (-au3).exp();
        let du_dr = 1.0 / denom;
        q * dl_du * du_dr
    }
}
/// Shell model atom for ionic systems.
///
/// The shell model splits each atom into a core and a massless shell
/// connected by a harmonic spring.  The shell carries the atomic
/// polarizability.
#[derive(Debug, Clone)]
pub struct ShellModelAtom {
    /// Core position.
    pub core_pos: [f64; 3],
    /// Shell position.
    pub shell_pos: [f64; 3],
    /// Core charge.
    pub core_charge: f64,
    /// Shell charge.
    pub shell_charge: f64,
    /// Core mass.
    pub core_mass: f64,
    /// Spring constant connecting core and shell.
    pub k_cs: f64,
    /// Polarizability (alpha = shell_charge^2 / k_cs).
    pub polarizability: f64,
}
/// Self-consistent point-dipole model.
#[derive(Debug, Clone)]
pub struct PointDipoleModel {
    /// Atom positions.
    pub positions: Vec<[f64; 3]>,
    /// Partial charges.
    pub charges: Vec<f64>,
    /// Isotropic polarizabilities.
    pub alphas: Vec<f64>,
    /// Induced dipole moments, solved by [`Self::solve_self_consistent`].
    pub dipoles: Vec<[f64; 3]>,
}
impl PointDipoleModel {
    /// Create a new `PointDipoleModel` with no atoms.
    pub fn new() -> Self {
        Self {
            positions: Vec::new(),
            charges: Vec::new(),
            alphas: Vec::new(),
            dipoles: Vec::new(),
        }
    }
    /// Add a site with position, charge, and polarizability.
    pub fn add_site(&mut self, pos: [f64; 3], charge: f64, alpha: f64) {
        self.positions.push(pos);
        self.charges.push(charge);
        self.alphas.push(alpha);
        self.dipoles.push([0.0; 3]);
    }
    /// Bare relay tensor between sites.
    fn relay_tensor(r_vec: [f64; 3]) -> [[f64; 3]; 3] {
        let r2 = norm2(r_vec);
        if r2 < 1e-24 {
            return [[0.0; 3]; 3];
        }
        let r = r2.sqrt();
        let r5 = r2 * r2 * r;
        let mut t = [[0.0f64; 3]; 3];
        for a in 0..3 {
            for b in 0..3 {
                let delta = if a == b { 1.0 } else { 0.0 };
                t[a][b] = (3.0 * r_vec[a] * r_vec[b] - r2 * delta) / r5;
            }
        }
        t
    }
    /// Solve for induced dipoles self-consistently by simple relaxation.
    ///
    /// Returns the number of iterations performed.
    pub fn solve_self_consistent(
        &mut self,
        e_ext: &[[f64; 3]],
        max_iter: usize,
        tol: f64,
    ) -> usize {
        let n = self.positions.len();
        for iter in 0..max_iter {
            let old_dipoles = self.dipoles.clone();
            let mut new_dipoles = vec![[0.0f64; 3]; n];
            for (i, (new_d, &e_ext_i)) in new_dipoles.iter_mut().zip(e_ext.iter()).enumerate() {
                let mut e_field = e_ext_i;
                for (j, &mu_j) in old_dipoles.iter().enumerate().take(n) {
                    if i == j {
                        continue;
                    }
                    let r_vec = sub(self.positions[i], self.positions[j]);
                    let t = Self::relay_tensor(r_vec);
                    for (ef_a, ta) in e_field.iter_mut().zip(t.iter()) {
                        *ef_a += ta[0] * mu_j[0] + ta[1] * mu_j[1] + ta[2] * mu_j[2];
                    }
                }
                *new_d = scale(self.alphas[i], e_field);
            }
            let mut rms = 0.0;
            for i in 0..n {
                let d = sub(new_dipoles[i], old_dipoles[i]);
                rms += norm2(d);
            }
            rms = (rms / (3 * n) as f64).sqrt();
            self.dipoles = new_dipoles;
            if rms < tol {
                return iter + 1;
            }
        }
        max_iter
    }
    /// Solve using SOR (Successive Over-Relaxation) for faster convergence.
    ///
    /// # Arguments
    /// * `omega` – relaxation parameter (1.0 = standard Gauss-Seidel, >1 = over-relaxation).
    pub fn solve_sor(
        &mut self,
        e_ext: &[[f64; 3]],
        max_iter: usize,
        tol: f64,
        omega: f64,
    ) -> usize {
        let n = self.positions.len();
        for iter in 0..max_iter {
            let mut max_change = 0.0_f64;
            for (i, &e_ext_i) in e_ext.iter().enumerate().take(n) {
                let mut e_field = e_ext_i;
                for j in 0..n {
                    if i == j {
                        continue;
                    }
                    let r_vec = sub(self.positions[i], self.positions[j]);
                    let t = Self::relay_tensor(r_vec);
                    let mu_j = self.dipoles[j];
                    for (ef_a, ta) in e_field.iter_mut().zip(t.iter()) {
                        *ef_a += ta[0] * mu_j[0] + ta[1] * mu_j[1] + ta[2] * mu_j[2];
                    }
                }
                let new_mu = scale(self.alphas[i], e_field);
                let old_mu = self.dipoles[i];
                let updated = add(scale(1.0 - omega, old_mu), scale(omega, new_mu));
                let change = norm(sub(updated, old_mu));
                if change > max_change {
                    max_change = change;
                }
                self.dipoles[i] = updated;
            }
            if max_change < tol {
                return iter + 1;
            }
        }
        max_iter
    }
    /// Polarization energy: `U_pol = -0.5 * sum_i mu_i . E_ext_i`.
    pub fn polarization_energy(&self, e_ext: &[[f64; 3]]) -> f64 {
        let mut u = 0.0;
        for (i, mu) in self.dipoles.iter().enumerate() {
            u -= 0.5 * dot(*mu, e_ext[i]);
        }
        u
    }
    /// Total induced dipole moment of the system.
    pub fn total_dipole(&self) -> [f64; 3] {
        let mut total = [0.0f64; 3];
        for mu in &self.dipoles {
            total = add(total, *mu);
        }
        total
    }
    /// Dipole-dipole interaction energy between all pairs.
    pub fn dipole_dipole_energy(&self) -> f64 {
        let n = self.dipoles.len();
        let mut energy = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r_vec = sub(self.positions[i], self.positions[j]);
                energy += InducibleDipoleForce::dipole_dipole_interaction(
                    self.dipoles[i],
                    self.dipoles[j],
                    r_vec,
                );
            }
        }
        energy
    }
}
/// Classical Drude oscillator force-field model.
#[derive(Debug, Clone)]
pub struct DrudeModel {
    /// All polarizable atoms in the system.
    pub atoms: Vec<PolarizableAtom>,
    /// Global Drude spring constant.
    pub k_drude: f64,
    /// Thole damping parameter (typically 2.6).
    pub thole_a: f64,
}
impl DrudeModel {
    /// Create a new `DrudeModel`.
    pub fn new(k_drude: f64, thole_a: f64) -> Self {
        Self {
            atoms: Vec::new(),
            k_drude,
            thole_a,
        }
    }
    /// Add a polarizable atom and automatically attach a Drude particle.
    pub fn add_atom(&mut self, pos: [f64; 3], charge: f64, mass: f64, alpha: f64) {
        let q_d = Self::drude_charge_from_alpha(alpha, self.k_drude);
        let drude = DrudeParticle {
            pos,
            charge: -q_d.abs(),
            mass: 0.4,
            k_spring: self.k_drude,
        };
        let atom = PolarizableAtom {
            pos,
            charge: charge + q_d.abs(),
            mass,
            polarizability: alpha,
            drude: Some(drude),
            dipole: [0.0; 3],
        };
        self.atoms.push(atom);
    }
    /// Compute the Drude charge from polarizability and spring constant:
    /// `q_D = sqrt(k * alpha)`.
    pub fn drude_charge_from_alpha(alpha: f64, k: f64) -> f64 {
        (k * alpha).sqrt()
    }
    /// Total spring (self-)energy over all Drude oscillators.
    pub fn spring_energy(&self) -> f64 {
        let mut energy = 0.0;
        for atom in &self.atoms {
            if let Some(ref drude) = atom.drude {
                let d = sub(drude.pos, atom.pos);
                energy += 0.5 * drude.k_spring * norm2(d);
            }
        }
        energy
    }
    /// Spring force on the *atom core* at index `atom_idx`.
    pub fn spring_force(&self, atom_idx: usize) -> [f64; 3] {
        let atom = &self.atoms[atom_idx];
        match &atom.drude {
            None => [0.0; 3],
            Some(drude) => {
                let d = sub(drude.pos, atom.pos);
                scale(drude.k_spring, d)
            }
        }
    }
    /// Displace the Drude particle for atom `atom_idx` by `displacement`.
    pub fn displace_drude(&mut self, atom_idx: usize, displacement: [f64; 3]) {
        if let Some(ref mut drude) = self.atoms[atom_idx].drude {
            drude.pos = add(drude.pos, displacement);
        }
    }
    /// Compute the induced dipole for each atom from Drude displacement.
    ///
    /// The induced dipole is mu = q_D * (r_D - r_core).
    pub fn compute_dipoles(&mut self) {
        for atom in &mut self.atoms {
            if let Some(ref drude) = atom.drude {
                let d = sub(drude.pos, atom.pos);
                atom.dipole = scale(drude.charge.abs(), d);
            }
        }
    }
    /// Total induced dipole moment of the system.
    pub fn total_dipole(&self) -> [f64; 3] {
        let mut total = [0.0f64; 3];
        for atom in &self.atoms {
            total = add(total, atom.dipole);
        }
        total
    }
    /// Total polarization energy = sum of spring energies.
    pub fn total_polarization_energy(&self) -> f64 {
        self.spring_energy()
    }
    /// Minimize Drude positions using steepest descent.
    ///
    /// Iterates to find equilibrium positions where spring forces balance
    /// external electric field forces.
    pub fn minimize_drude(
        &mut self,
        e_ext: &[[f64; 3]],
        max_iter: usize,
        step_size: f64,
        tol: f64,
    ) -> usize {
        for iter in 0..max_iter {
            let mut max_force = 0.0_f64;
            for (i, atom) in self.atoms.iter_mut().enumerate() {
                if let Some(ref mut drude) = atom.drude {
                    let d = sub(drude.pos, atom.pos);
                    let f_spring = scale(-drude.k_spring, d);
                    let f_field = scale(drude.charge.abs(), e_ext[i]);
                    let f_total = add(f_spring, f_field);
                    let f_mag = norm(f_total);
                    if f_mag > max_force {
                        max_force = f_mag;
                    }
                    drude.pos = add(drude.pos, scale(step_size, f_total));
                }
            }
            if max_force < tol {
                return iter + 1;
            }
        }
        max_iter
    }
}
/// An atom that can be polarized, optionally carrying a Drude oscillator.
#[derive(Debug, Clone)]
pub struct PolarizableAtom {
    /// Position of the atom core (Angstrom).
    pub pos: [f64; 3],
    /// Core (partial) charge in units of *e*.
    pub charge: f64,
    /// Mass of the atom core (Da).
    pub mass: f64,
    /// Isotropic polarizability alpha (Angstrom^3).
    pub polarizability: f64,
    /// Optional Drude oscillator attached to this atom.
    pub drude: Option<DrudeParticle>,
    /// Current induced dipole moment (e*Angstrom).
    pub dipole: [f64; 3],
}
/// Table of CHARMM Drude parameters by atom type.
#[derive(Debug, Default, Clone)]
pub struct CharmmDrudeParamTable {
    pub(super) entries: Vec<CharmmDrudeParams>,
}
impl CharmmDrudeParamTable {
    /// Create an empty table.
    pub fn new() -> Self {
        CharmmDrudeParamTable {
            entries: Vec::new(),
        }
    }
    /// Create a table pre-populated with common CHARMM atom types.
    pub fn default_charmm() -> Self {
        let mut t = Self::new();
        t.add(CharmmDrudeParams::cg331());
        t.add(CharmmDrudeParams::og311());
        t.add(CharmmDrudeParams::ng2s1());
        t
    }
    /// Add or replace an entry.
    pub fn add(&mut self, params: CharmmDrudeParams) {
        if let Some(e) = self
            .entries
            .iter_mut()
            .find(|e| e.atom_type == params.atom_type)
        {
            *e = params;
        } else {
            self.entries.push(params);
        }
    }
    /// Look up parameters by atom type.
    pub fn get(&self, atom_type: &str) -> Option<&CharmmDrudeParams> {
        self.entries.iter().find(|e| e.atom_type == atom_type)
    }
}
/// A polarizable site that tracks its own induced dipole moment.
#[derive(Debug, Clone)]
pub struct InducedDipole {
    /// Isotropic polarizability (Angstrom^3 or consistent units).
    pub polarizability: f64,
    /// Position of the site.
    pub position: [f64; 3],
    /// Current induced dipole moment (e·Angstrom or consistent units).
    pub dipole: [f64; 3],
}
impl InducedDipole {
    /// Create a new `InducedDipole` with zero dipole.
    pub fn new(polarizability: f64, position: [f64; 3]) -> Self {
        Self {
            polarizability,
            position,
            dipole: [0.0; 3],
        }
    }
    /// One SCF update step: p = alpha * E_ext.
    pub fn update(&mut self, e_ext: [f64; 3]) {
        self.dipole = scale(self.polarizability, e_ext);
    }
}
/// Simplified Drude oscillator model parameterised by a scalar
/// polarizability `alpha` and spring constant `k`.
#[derive(Debug, Clone)]
pub struct DruckerModel {
    /// Scalar polarizability (Angstrom^3 or consistent units).
    pub alpha: f64,
    /// Spring constant connecting core and Drude charge (kcal/mol/Angstrom^2).
    pub k: f64,
}
impl DruckerModel {
    /// Create a new `DruckerModel`.
    pub fn new(alpha: f64, k: f64) -> Self {
        Self { alpha, k }
    }
    /// Induced dipole moment: p = alpha * E.
    pub fn dipole_moment(&self, e_field: [f64; 3]) -> [f64; 3] {
        scale(self.alpha, e_field)
    }
    /// Spring force on the Drude charge: F = -k * (r_drude - r_core).
    pub fn drude_spring_force(&self, r_core: [f64; 3], r_drude: [f64; 3]) -> [f64; 3] {
        let d = sub(r_drude, r_core);
        scale(-self.k, d)
    }
}
/// Per-atom parameters for the charge equilibration (QEq) model.
///
/// The QEq model (Rappé & Goddard 1991) assigns partial charges by
/// minimising the total electrostatic energy subject to a charge
/// neutrality constraint, using per-element electronegativity χ and
/// hardness η.
#[derive(Debug, Clone)]
pub struct QeqAtom {
    /// Electronegativity χ (eV).
    pub chi: f64,
    /// Chemical hardness η (eV).
    pub eta: f64,
    /// Current partial charge (e).
    pub charge: f64,
}
impl QeqAtom {
    /// Create a new QEq atom.
    pub fn new(chi: f64, eta: f64) -> Self {
        QeqAtom {
            chi,
            eta,
            charge: 0.0,
        }
    }
}
/// Shell model for ionic systems (Dick-Overhauser model).
#[derive(Debug, Clone)]
pub struct ShellModel {
    /// All shell model atoms.
    pub atoms: Vec<ShellModelAtom>,
}
impl ShellModel {
    /// Create a new empty shell model.
    pub fn new() -> Self {
        Self { atoms: Vec::new() }
    }
    /// Add an atom with core and shell parameters.
    pub fn add_atom(
        &mut self,
        pos: [f64; 3],
        core_charge: f64,
        shell_charge: f64,
        core_mass: f64,
        k_cs: f64,
    ) {
        let polarizability = if k_cs.abs() > 1e-30 {
            shell_charge * shell_charge / k_cs
        } else {
            0.0
        };
        self.atoms.push(ShellModelAtom {
            core_pos: pos,
            shell_pos: pos,
            core_charge,
            shell_charge,
            core_mass,
            k_cs,
            polarizability,
        });
    }
    /// Total spring energy between cores and shells.
    pub fn spring_energy(&self) -> f64 {
        let mut energy = 0.0;
        for atom in &self.atoms {
            let d = sub(atom.shell_pos, atom.core_pos);
            energy += 0.5 * atom.k_cs * norm2(d);
        }
        energy
    }
    /// Core-shell displacement for atom `idx`.
    pub fn displacement(&self, idx: usize) -> [f64; 3] {
        sub(self.atoms[idx].shell_pos, self.atoms[idx].core_pos)
    }
    /// Induced dipole from core-shell displacement: mu = Y * d
    /// where Y is the shell charge.
    pub fn induced_dipole(&self, idx: usize) -> [f64; 3] {
        let d = self.displacement(idx);
        scale(self.atoms[idx].shell_charge, d)
    }
    /// Total polarization energy.
    pub fn total_polarization_energy(&self) -> f64 {
        self.spring_energy()
    }
    /// Relax shell positions to equilibrium given external field.
    pub fn relax_shells(
        &mut self,
        e_ext: &[[f64; 3]],
        max_iter: usize,
        step_size: f64,
        tol: f64,
    ) -> usize {
        for iter in 0..max_iter {
            let mut max_force = 0.0_f64;
            for (i, atom) in self.atoms.iter_mut().enumerate() {
                let d = sub(atom.shell_pos, atom.core_pos);
                let f_spring = scale(-atom.k_cs, d);
                let f_field = scale(atom.shell_charge, e_ext[i]);
                let f_total = add(f_spring, f_field);
                let f_mag = norm(f_total);
                if f_mag > max_force {
                    max_force = f_mag;
                }
                atom.shell_pos = add(atom.shell_pos, scale(step_size, f_total));
            }
            if max_force < tol {
                return iter + 1;
            }
        }
        max_iter
    }
}
/// Extended Thole damping with selectable model and per-pair parameters.
#[derive(Debug, Clone)]
pub struct TholeDampingExt {
    /// Smearing model.
    pub model: TholeModel,
    /// Damping parameter a (dimensionless, typically 1.1–2.3 for exponential).
    pub a: f64,
}
impl TholeDampingExt {
    /// Create an extended Thole damper.
    pub fn new(model: TholeModel, a: f64) -> Self {
        TholeDampingExt { model, a }
    }
    /// Reduced distance u = r / (α_i α_j)^{1/6}.
    pub fn reduced_distance(r: f64, alpha_i: f64, alpha_j: f64) -> f64 {
        let denom = (alpha_i * alpha_j).powf(1.0 / 6.0);
        if denom < 1e-30 {
            return f64::INFINITY;
        }
        r / denom
    }
    /// Damping function f(u) – scales the T2 dipole interaction tensor.
    pub fn damping_fn(&self, u: f64) -> f64 {
        let au3 = self.a * u * u * u;
        match self.model {
            TholeModel::Exponential => 1.0 - (-au3).exp(),
            TholeModel::Linear => {
                if au3 >= 1.0 {
                    1.0
                } else {
                    au3
                }
            }
            TholeModel::Gaussian => {
                let au2 = self.a * u * u;
                1.0 - (-au2).exp()
            }
        }
    }
    /// Damped T2 tensor element λ_{ij}(r).
    ///
    /// Returns the scalar damping factor to multiply the bare T2 tensor.
    pub fn t2_damping(&self, r: f64, alpha_i: f64, alpha_j: f64) -> f64 {
        let u = Self::reduced_distance(r, alpha_i, alpha_j);
        self.damping_fn(u)
    }
}
/// SWM4-NDP polarizable water model parameters.
///
/// Reference: Lamoureux *et al.* (2006).
#[derive(Debug, Clone)]
pub struct Swm4NdpParams {
    /// O-H bond length (Å).
    pub r_oh: f64,
    /// H-O-H angle (degrees).
    pub theta_hoh: f64,
    /// Partial charge on oxygen (e).
    pub q_o: f64,
    /// Partial charge on hydrogen (e).
    pub q_h: f64,
    /// Oxygen polarizability (Å³).
    pub alpha_o: f64,
    /// LJ ε for O-O (kcal/mol).
    pub eps_oo: f64,
    /// LJ σ for O-O (Å).
    pub sigma_oo: f64,
    /// Drude charge (e).
    pub q_drude: f64,
    /// Drude spring constant k_D (kcal mol^-1 Å^-2).
    pub k_drude: f64,
}
impl Swm4NdpParams {
    /// Standard SWM4-NDP parameters.
    pub fn new() -> Self {
        Swm4NdpParams {
            r_oh: 0.9572,
            theta_hoh: 104.52,
            q_o: -1.04844,
            q_h: 0.52422,
            alpha_o: 0.97825,
            eps_oo: -0.21094,
            sigma_oo: 3.18395,
            q_drude: -1.71636,
            k_drude: 1000.0,
        }
    }
    /// Effective charge on the Drude particle (same as q_drude by construction).
    pub fn drude_charge(&self) -> f64 {
        self.q_drude
    }
    /// Derived polarizability from α = q_D² / k_D (in Å³ assuming consistent units).
    pub fn derived_polarizability(&self) -> f64 {
        self.q_drude * self.q_drude / self.k_drude
    }
}
/// Extended-Lagrangian Drude integrator with dual thermostat.
///
/// The real atoms are coupled to a thermostat at temperature T_real,
/// while the Drude particles are coupled to a cold thermostat at T_drude
/// (typically 1 K) to keep them close to their SCF positions.
#[derive(Debug, Clone)]
pub struct DrudeExtLagrangian {
    /// Number of Drude pairs.
    pub n_pairs: usize,
    /// Core (real atom) positions (Å).
    pub core_pos: Vec<[f64; 3]>,
    /// Core velocities (Å ps^-1).
    pub core_vel: Vec<[f64; 3]>,
    /// Drude positions (Å).
    pub drude_pos: Vec<[f64; 3]>,
    /// Drude velocities (Å ps^-1).
    pub drude_vel: Vec<[f64; 3]>,
    /// Core masses (amu).
    pub core_mass: Vec<f64>,
    /// Drude masses (amu).
    pub drude_mass: Vec<f64>,
    /// Spring constant k_D (kcal mol^-1 Å^-2).
    pub k_spring: f64,
    /// Temperature for Drude thermostat (K).
    pub t_drude: f64,
    /// Langevin friction coefficient for Drude particles (ps^-1).
    pub gamma_drude: f64,
}
impl DrudeExtLagrangian {
    /// Create a new extended-Lagrangian Drude state.
    pub fn new(
        core_pos: Vec<[f64; 3]>,
        drude_pos: Vec<[f64; 3]>,
        core_mass: Vec<f64>,
        drude_mass: Vec<f64>,
        k_spring: f64,
        t_drude: f64,
        gamma_drude: f64,
    ) -> Self {
        let n = core_pos.len();
        DrudeExtLagrangian {
            n_pairs: n,
            core_pos,
            core_vel: vec![[0.0; 3]; n],
            drude_pos,
            drude_vel: vec![[0.0; 3]; n],
            core_mass,
            drude_mass,
            k_spring,
            t_drude,
            gamma_drude,
        }
    }
    /// Compute spring force on each Drude particle: F = -k_D (r_D - r_C).
    pub fn spring_forces(&self) -> Vec<[f64; 3]> {
        self.drude_pos
            .iter()
            .zip(self.core_pos.iter())
            .map(|(rd, rc)| {
                [
                    -self.k_spring * (rd[0] - rc[0]),
                    -self.k_spring * (rd[1] - rc[1]),
                    -self.k_spring * (rd[2] - rc[2]),
                ]
            })
            .collect()
    }
    /// Drude kinetic energy (sum over all Drude particles).
    pub fn drude_kinetic_energy(&self) -> f64 {
        self.drude_vel
            .iter()
            .zip(self.drude_mass.iter())
            .map(|(v, &m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum()
    }
    /// Apply Langevin thermostat kick to Drude velocities (Ornstein–Uhlenbeck step).
    ///
    /// Uses the analytical O-U propagator: v → v exp(-γ dt) + √(k_B T_D / m) √(1-exp(-2γ dt)) ξ.
    pub fn langevin_drude_kick(&mut self, dt: f64) {
        const KB_MD: f64 = 8.314_462_618e-3;
        let decay = (-self.gamma_drude * dt).exp();
        let noise_scale_sq = 1.0 - decay * decay;
        for i in 0..self.n_pairs {
            let sigma = (KB_MD * self.t_drude / self.drude_mass[i] * noise_scale_sq).sqrt();
            let mut rng = rand::rng();
            use rand::RngExt;
            for v in &mut self.drude_vel[i] {
                let xi: f64 = rng.random_range(-1.0_f64..1.0_f64);
                *v = decay * *v + sigma * xi;
            }
        }
    }
}
/// A simple container for a set of induced dipoles and their polarizabilities,
/// providing the classical polarization energy.
pub struct InducedDipoleModel {
    /// Isotropic polarizabilities α_i (Å³ or consistent units).
    pub polarizabilities: Vec<f64>,
    /// Induced dipole moments p_i (e·Å or consistent units).
    pub dipoles: Vec<[f64; 3]>,
}
impl InducedDipoleModel {
    /// Create a new `InducedDipoleModel`.
    pub fn new(polarizabilities: Vec<f64>, dipoles: Vec<[f64; 3]>) -> Self {
        Self {
            polarizabilities,
            dipoles,
        }
    }
    /// Classical polarization (self-induction) energy.
    ///
    /// For linear response (p = α E_ind):
    ///
    /// ```text
    /// U_pol = -½ Σ_i |p_i|² / α_i
    /// ```
    ///
    /// This is always ≤ 0 (polarization lowers the energy).
    pub fn compute_polarization_energy(&self) -> f64 {
        self.polarizabilities
            .iter()
            .zip(self.dipoles.iter())
            .filter(|(alpha, _)| **alpha > 1e-30)
            .map(|(alpha, p)| -0.5 * dot(*p, *p) / alpha)
            .sum()
    }
    /// Dipole-dipole interaction energy between two sites.
    ///
    /// ```text
    /// U_dd = (p_i · p_j) / r³ - 3 (p_i · r̂)(p_j · r̂) / r³
    /// ```
    ///
    /// # Arguments
    /// * `p_i`, `p_j` – dipole moments.
    /// * `r_vec`      – displacement vector from site j to site i.
    pub fn dipole_dipole_energy(p_i: [f64; 3], p_j: [f64; 3], r_vec: [f64; 3]) -> f64 {
        let r2 = norm2(r_vec);
        if r2 < 1e-24 {
            return 0.0;
        }
        let r = r2.sqrt();
        let r3 = r2 * r;
        let r_hat = scale(1.0 / r, r_vec);
        let pi_dot_pj = dot(p_i, p_j);
        let pi_dot_rhat = dot(p_i, r_hat);
        let pj_dot_rhat = dot(p_j, r_hat);
        (pi_dot_pj - 3.0 * pi_dot_rhat * pj_dot_rhat) / r3
    }
    /// Total pairwise dipole-dipole interaction energy among all sites.
    ///
    /// Uses the positions provided separately (not stored in this struct).
    pub fn total_dipole_dipole_energy(&self, positions: &[[f64; 3]]) -> f64 {
        let n = self.dipoles.len().min(positions.len());
        let mut energy = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r_vec = sub(positions[i], positions[j]);
                energy += Self::dipole_dipole_energy(self.dipoles[i], self.dipoles[j], r_vec);
            }
        }
        energy
    }
}
/// Fluctuating-charge (FQ) extended Lagrangian dynamics.
///
/// Each atom carries a charge degree of freedom q_i with a fictitious
/// mass m_q.  The equations of motion are integrated alongside atomic
/// positions to propagate charges adiabatically.
#[derive(Debug, Clone)]
pub struct FluctuatingCharge {
    /// Partial charges (e).
    pub charges: Vec<f64>,
    /// Charge velocities (e ps^-1).
    pub charge_velocities: Vec<f64>,
    /// Fictitious charge mass (amu Å² e^-2).
    pub charge_mass: f64,
    /// Per-atom electronegativity χ_i (eV).
    pub chi: Vec<f64>,
    /// Per-atom hardness η_i (eV).
    pub eta: Vec<f64>,
}
impl FluctuatingCharge {
    /// Create an FQ model for `n` atoms with the given parameters.
    pub fn new(n: usize, charge_mass: f64, chi: Vec<f64>, eta: Vec<f64>) -> Self {
        FluctuatingCharge {
            charges: vec![0.0; n],
            charge_velocities: vec![0.0; n],
            charge_mass,
            chi,
            eta,
        }
    }
    /// Compute the electronegativity force on each charge: -∂E/∂q_i.
    ///
    /// E = Σ χ_i q_i + Σ η_i q_i² + Σ_{i<j} q_i J_{ij} q_j
    pub fn electronegativity_force(&self, positions: &[[f64; 3]], shielding: f64) -> Vec<f64> {
        let n = self.charges.len();
        let mut f = vec![0.0f64; n];
        for i in 0..n {
            f[i] -= self.chi[i];
            f[i] -= 2.0 * self.eta[i] * self.charges[i];
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = positions[i][0] - positions[j][0];
                let dy = positions[i][1] - positions[j][1];
                let dz = positions[i][2] - positions[j][2];
                let r2 = dx * dx + dy * dy + dz * dz;
                let j_ij = 1.0 / (r2 + shielding * shielding).sqrt();
                f[i] -= j_ij * self.charges[j];
            }
        }
        f
    }
    /// Velocity Verlet half-step for charge velocities.
    pub fn velocity_half_step(&mut self, forces: &[f64], dt: f64) {
        for (cv, &f) in self.charge_velocities.iter_mut().zip(forces.iter()) {
            *cv += 0.5 * dt * f / self.charge_mass;
        }
    }
    /// Full position (charge) step.
    pub fn charge_step(&mut self, dt: f64) {
        for i in 0..self.charges.len() {
            self.charges[i] += dt * self.charge_velocities[i];
        }
    }
    /// Kinetic energy of the charge degrees of freedom.
    pub fn charge_kinetic_energy(&self) -> f64 {
        self.charge_velocities
            .iter()
            .map(|v| 0.5 * self.charge_mass * v * v)
            .sum()
    }
}
/// Utility functions for point-dipole electric fields and interaction energies.
pub struct InducibleDipoleForce;
impl InducibleDipoleForce {
    /// Electric field produced by a point dipole `mu` at displacement `r`.
    pub fn dipole_electric_field(mu: [f64; 3], r: [f64; 3]) -> [f64; 3] {
        let r2 = norm2(r);
        if r2 < 1e-24 {
            return [0.0; 3];
        }
        let r_norm = r2.sqrt();
        let r3 = r2 * r_norm;
        let mu_dot_r = dot(mu, r);
        let factor = 1.0 / r3;
        [
            factor * (3.0 * mu_dot_r * r[0] / r2 - mu[0]),
            factor * (3.0 * mu_dot_r * r[1] / r2 - mu[1]),
            factor * (3.0 * mu_dot_r * r[2] / r2 - mu[2]),
        ]
    }
    /// Interaction energy between a point charge `q` and dipole `mu`.
    pub fn charge_dipole_interaction(q: f64, mu: [f64; 3], r: [f64; 3]) -> f64 {
        let r2 = norm2(r);
        if r2 < 1e-24 {
            return 0.0;
        }
        -q * dot(mu, r) / (r2 * r2.sqrt())
    }
    /// Interaction energy between two point dipoles.
    pub fn dipole_dipole_interaction(mu_i: [f64; 3], mu_j: [f64; 3], r: [f64; 3]) -> f64 {
        let r2 = norm2(r);
        if r2 < 1e-24 {
            return 0.0;
        }
        let r_norm = r2.sqrt();
        let r3 = r2 * r_norm;
        let mi_dot_mj = dot(mu_i, mu_j);
        let mi_dot_r = dot(mu_i, r);
        let mj_dot_r = dot(mu_j, r);
        (mi_dot_mj - 3.0 * mi_dot_r * mj_dot_r / r2) / r3
    }
    /// Gradient of the dipole-dipole interaction energy with respect to
    /// the displacement vector r.
    ///
    /// Returns the force on dipole i due to dipole j.
    pub fn dipole_dipole_force(mu_i: [f64; 3], mu_j: [f64; 3], r: [f64; 3]) -> [f64; 3] {
        let r2 = norm2(r);
        if r2 < 1e-24 {
            return [0.0; 3];
        }
        let r_norm = r2.sqrt();
        let r5 = r2 * r2 * r_norm;
        let r7 = r5 * r2;
        let mi_dot_mj = dot(mu_i, mu_j);
        let mi_dot_r = dot(mu_i, r);
        let mj_dot_r = dot(mu_j, r);
        let mut f = [0.0f64; 3];
        for k in 0..3 {
            f[k] = -3.0 / r5 * (mi_dot_mj * r[k] + mu_i[k] * mj_dot_r + mu_j[k] * mi_dot_r)
                + 15.0 * mi_dot_r * mj_dot_r * r[k] / r7;
        }
        f
    }
}
/// Simplified charge-on-spring (COS) model.
///
/// Similar to the Drude model but with a simpler treatment where
/// the spring charge magnitude is a free parameter rather than
/// being derived from polarizability.
#[derive(Debug, Clone)]
pub struct ChargeOnSpring {
    /// Atom positions.
    pub positions: Vec<[f64; 3]>,
    /// Spring charge positions (auxiliary sites).
    pub spring_positions: Vec<[f64; 3]>,
    /// Magnitude of the spring charge.
    pub spring_charges: Vec<f64>,
    /// Spring constants.
    pub spring_constants: Vec<f64>,
    /// Damping coefficient for the charge displacement.
    pub damping: f64,
}
impl ChargeOnSpring {
    /// Create a new COS model with a given damping coefficient.
    pub fn new(damping: f64) -> Self {
        Self {
            positions: Vec::new(),
            spring_positions: Vec::new(),
            spring_charges: Vec::new(),
            spring_constants: Vec::new(),
            damping,
        }
    }
    /// Add an atom site.
    pub fn add_site(&mut self, pos: [f64; 3], q_spring: f64, k_spring: f64) {
        self.positions.push(pos);
        self.spring_positions.push(pos);
        self.spring_charges.push(q_spring);
        self.spring_constants.push(k_spring);
    }
    /// Compute the COS spring energy.
    pub fn spring_energy(&self) -> f64 {
        let mut energy = 0.0;
        for i in 0..self.positions.len() {
            let d = sub(self.spring_positions[i], self.positions[i]);
            energy += 0.5 * self.spring_constants[i] * norm2(d);
        }
        energy
    }
    /// Compute the effective polarizability from the COS parameters:
    /// alpha_eff = q^2 / k.
    pub fn effective_polarizability(&self, idx: usize) -> f64 {
        let q = self.spring_charges[idx];
        let k = self.spring_constants[idx];
        if k.abs() < 1e-30 {
            return 0.0;
        }
        q * q / k
    }
    /// Compute the induced dipole at site `idx`.
    pub fn induced_dipole(&self, idx: usize) -> [f64; 3] {
        let d = sub(self.spring_positions[idx], self.positions[idx]);
        scale(self.spring_charges[idx], d)
    }
    /// Update spring positions given external electric field.
    ///
    /// The equilibrium displacement is: d = q * E / k
    pub fn update_positions(&mut self, e_ext: &[[f64; 3]]) {
        for (i, sp) in self.spring_positions.iter_mut().enumerate() {
            let k = self.spring_constants[i];
            if k.abs() < 1e-30 {
                continue;
            }
            let q = self.spring_charges[i];
            let d = scale(q / k * self.damping, e_ext[i]);
            *sp = add(self.positions[i], d);
        }
    }
}
/// Utility methods related to Drude oscillator energetics.
pub struct Drude;
impl Drude {
    /// Drude self-energy correction.
    ///
    /// The harmonic spring energy stored in a single Drude pair when the Drude
    /// particle is displaced by (`dx`, `dy`, `dz`) relative to its parent core:
    ///
    /// ```text
    /// U_self = ½ k_D * (dx² + dy² + dz²)
    /// ```
    ///
    /// # Arguments
    /// * `dx`, `dy`, `dz` – displacement components (Å).
    /// * `k_spring`       – Drude spring constant (kcal mol⁻¹ Å⁻²).
    pub fn compute_self_energy(dx: f64, dy: f64, dz: f64, k_spring: f64) -> f64 {
        0.5 * k_spring * (dx * dx + dy * dy + dz * dz)
    }
    /// Total Drude self-energy for a collection of Drude pairs.
    ///
    /// Sums `compute_self_energy` over all `n` pairs.
    ///
    /// # Arguments
    /// * `core_pos`  – core particle positions (Å).
    /// * `drude_pos` – Drude particle positions (Å).
    /// * `k_spring`  – common spring constant (kcal mol⁻¹ Å⁻²).
    pub fn total_self_energy(core_pos: &[[f64; 3]], drude_pos: &[[f64; 3]], k_spring: f64) -> f64 {
        core_pos
            .iter()
            .zip(drude_pos.iter())
            .map(|(c, d)| {
                let dx = d[0] - c[0];
                let dy = d[1] - c[1];
                let dz = d[2] - c[2];
                Self::compute_self_energy(dx, dy, dz, k_spring)
            })
            .sum()
    }
}
/// Amoeba mutual-induction state for a collection of polarizable sites.
///
/// The Amoeba force field uses smeared Gaussian charge distributions and
/// iterative mutual polarisation.  This struct holds the current induced
/// dipoles and provides an iteration update step.
#[derive(Debug, Clone)]
pub struct AmoebaPolarization {
    /// Number of polarizable sites.
    pub n_sites: usize,
    /// Isotropic polarizability per site (Å³).
    pub alpha: Vec<f64>,
    /// Current induced dipoles (e·Å).
    pub mu_ind: Vec<[f64; 3]>,
    /// Thole damping parameter (Amoeba default ≈ 0.39).
    pub thole: f64,
}
impl AmoebaPolarization {
    /// Create a new Amoeba polarization state.
    pub fn new(alpha: Vec<f64>, thole: f64) -> Self {
        let n = alpha.len();
        AmoebaPolarization {
            n_sites: n,
            alpha,
            mu_ind: vec![[0.0; 3]; n],
            thole,
        }
    }
    /// Compute the electric field at site `i` due to induced dipoles at all other sites.
    ///
    /// Uses the bare dipole field (no Thole damping for simplicity).
    pub fn field_from_induced(&self, i: usize, positions: &[[f64; 3]]) -> [f64; 3] {
        let mut e = [0.0f64; 3];
        for j in 0..self.n_sites {
            if i == j {
                continue;
            }
            let r = [
                positions[i][0] - positions[j][0],
                positions[i][1] - positions[j][1],
                positions[i][2] - positions[j][2],
            ];
            let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
            let r1 = r2.sqrt();
            let r3 = r1 * r2;
            let r5 = r3 * r2;
            let mu = self.mu_ind[j];
            let mu_dot_r = mu[0] * r[0] + mu[1] * r[1] + mu[2] * r[2];
            for d in 0..3 {
                e[d] += 3.0 * mu_dot_r * r[d] / r5 - mu[d] / r3;
            }
        }
        e
    }
    /// Perform one Jacobi update of all induced dipoles.
    ///
    /// μ_i^{new} = α_i (E_perm_i + E_ind_i)
    /// where E_perm_i is the permanent field at site i.
    pub fn jacobi_update(&mut self, e_perm: &[[f64; 3]], positions: &[[f64; 3]]) {
        let mut mu_new = vec![[0.0f64; 3]; self.n_sites];
        for i in 0..self.n_sites {
            let e_ind = self.field_from_induced(i, positions);
            for d in 0..3 {
                mu_new[i][d] = self.alpha[i] * (e_perm[i][d] + e_ind[d]);
            }
        }
        self.mu_ind = mu_new;
    }
    /// Iterative SCF until convergence (max norm of Δμ < tol).
    pub fn scf_solve(
        &mut self,
        e_perm: &[[f64; 3]],
        positions: &[[f64; 3]],
        max_iter: usize,
        tol: f64,
    ) {
        for _ in 0..max_iter {
            let old = self.mu_ind.clone();
            self.jacobi_update(e_perm, positions);
            let mut max_delta = 0.0f64;
            for (mu_i, old_i) in self.mu_ind.iter().zip(old.iter()) {
                for d in 0..3 {
                    let delta = (mu_i[d] - old_i[d]).abs();
                    if delta > max_delta {
                        max_delta = delta;
                    }
                }
            }
            if max_delta < tol {
                break;
            }
        }
    }
    /// Induction energy = -1/2 Σ_i μ_i · E_perm_i.
    pub fn induction_energy(&self, e_perm: &[[f64; 3]]) -> f64 {
        let mut u = 0.0f64;
        for (mu_i, ep_i) in self.mu_ind.iter().zip(e_perm.iter()) {
            for d in 0..3 {
                u -= 0.5 * mu_i[d] * ep_i[d];
            }
        }
        u
    }
}
/// An auxiliary Drude charge tethered to its parent atom by a harmonic spring.
#[derive(Debug, Clone)]
pub struct DrudeParticle {
    /// Position of the Drude particle (Angstrom).
    pub pos: [f64; 3],
    /// Drude charge (usually -1 to -2 *e*).
    pub charge: f64,
    /// Mass of the Drude particle (~0.4 Da).
    pub mass: f64,
    /// Spring constant to the parent atom (kcal/mol/Angstrom^2).
    pub k_spring: f64,
}
/// Charge equilibration solver.
///
/// Minimises the energy E = Σ_i (χ_i q_i + η_i q_i²) + Σ_{i<j} q_i J_{ij} q_j
/// subject to Σ q_i = Q_total, where J_{ij} is the shielded Coulomb kernel.
#[derive(Debug, Clone)]
pub struct QeqSolver {
    /// Shielding parameter for off-diagonal Coulomb interactions (Å).
    pub shielding: f64,
    /// Maximum SCF iterations.
    pub max_iter: usize,
    /// Convergence threshold on charge change (e).
    pub tol: f64,
}
impl QeqSolver {
    /// Create a new QEq solver with default parameters.
    pub fn new() -> Self {
        QeqSolver {
            shielding: 0.5,
            max_iter: 200,
            tol: 1e-8,
        }
    }
    /// Shielded Coulomb interaction integral J_{ij} (in eV/e² when r in Å).
    ///
    /// Uses the shielded expression: J = 1 / sqrt(r² + λ²), with λ = shielding.
    pub fn coulomb_kernel(&self, r: f64) -> f64 {
        1.0 / (r * r + self.shielding * self.shielding).sqrt()
    }
    /// Build the electronegativity equalization matrix A and vector b.
    ///
    /// The system A q = b, with a Lagrange multiplier row/col for the
    /// charge-neutrality constraint, gives the equilibrium charges.
    ///
    /// # Arguments
    /// * `atoms` – slice of QEq atoms
    /// * `positions` – Cartesian positions (Å)
    /// * `q_total` – total charge (e)
    ///
    /// Returns (A, b) with size (n+1) × (n+1).
    pub fn build_system(
        &self,
        atoms: &[QeqAtom],
        positions: &[[f64; 3]],
        q_total: f64,
    ) -> (Vec<Vec<f64>>, Vec<f64>) {
        let n = atoms.len();
        let dim = n + 1;
        let mut a = vec![vec![0.0f64; dim]; dim];
        let mut b = vec![0.0f64; dim];
        for i in 0..n {
            a[i][i] = 2.0 * atoms[i].eta;
            b[i] = -atoms[i].chi;
            for j in (i + 1)..n {
                let dx = positions[i][0] - positions[j][0];
                let dy = positions[i][1] - positions[j][1];
                let dz = positions[i][2] - positions[j][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                let j_ij = self.coulomb_kernel(r);
                a[i][j] = j_ij;
                a[j][i] = j_ij;
            }
            a[i][n] = 1.0;
            a[n][i] = 1.0;
        }
        b[n] = q_total;
        (a, b)
    }
    /// Solve the linear system A x = b using Gaussian elimination (in-place).
    fn gaussian_solve(a: &mut [Vec<f64>], b: &mut [f64]) -> Option<Vec<f64>> {
        let n = b.len();
        for col in 0..n {
            let mut max_row = col;
            let mut max_val = a[col][col].abs();
            for (row, a_row) in a.iter().enumerate().take(n).skip(col + 1) {
                if a_row[col].abs() > max_val {
                    max_val = a_row[col].abs();
                    max_row = row;
                }
            }
            if max_val < 1e-14 {
                return None;
            }
            a.swap(col, max_row);
            b.swap(col, max_row);
            let pivot = a[col][col];
            for v in a[col].iter_mut().skip(col) {
                *v /= pivot;
            }
            b[col] /= pivot;
            for row in 0..n {
                if row == col {
                    continue;
                }
                let factor = a[row][col];
                let pivot_row_slice: Vec<f64> = a[col][col..n].to_vec();
                for (a_row_c, &pv) in a[row][col..n].iter_mut().zip(pivot_row_slice.iter()) {
                    *a_row_c -= factor * pv;
                }
                b[row] -= factor * b[col];
            }
        }
        Some(b.to_vec())
    }
    /// Solve for equilibrium charges.
    ///
    /// Returns None if the system is singular.
    pub fn solve(
        &self,
        atoms: &[QeqAtom],
        positions: &[[f64; 3]],
        q_total: f64,
    ) -> Option<Vec<f64>> {
        let (mut a, mut b) = self.build_system(atoms, positions, q_total);
        let sol = Self::gaussian_solve(&mut a, &mut b)?;
        let n = atoms.len();
        Some(sol[..n].to_vec())
    }
}
/// Thole smearing models for short-range dipole damping.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TholeModel {
    /// Exponential smearing: f(u) = 1 - exp(-a u³).
    Exponential,
    /// Linear smearing: f(u) = max(0, 1 - a u³).
    Linear,
    /// Tinker-style smearing (Gaussian): f(u) = 1 - exp(-a u²).
    Gaussian,
}
/// CHARMM Drude force field atom-type parameters.
///
/// Reference: Lamoureux & Roux (2003); MacKerell & Lopes (2015).
#[derive(Debug, Clone)]
pub struct CharmmDrudeParams {
    /// Atom type name (e.g. "CG331").
    pub atom_type: String,
    /// Drude charge q_D (e), always negative.
    pub drude_charge: f64,
    /// Drude spring constant k_D (kcal mol^-1 Å^-2).
    pub k_drude: f64,
    /// Atomic polarizability α (Å³).
    pub polarizability: f64,
    /// Thole damping parameter a (dimensionless).
    pub thole_a: f64,
}
impl CharmmDrudeParams {
    /// Create parameters for a common CHARMM Drude atom type.
    pub fn new(
        atom_type: impl Into<String>,
        drude_charge: f64,
        k_drude: f64,
        polarizability: f64,
        thole_a: f64,
    ) -> Self {
        CharmmDrudeParams {
            atom_type: atom_type.into(),
            drude_charge,
            k_drude,
            polarizability,
            thole_a,
        }
    }
    /// Return polarizability in SI units (C² s² kg^-1) from input in Å³.
    ///
    /// Conversion: 1 Å³ = 1e-30 m³; multiply by 4πε₀ = 1.1127e-10 C² s² kg^-1 m^-3.
    pub fn polarizability_si(&self) -> f64 {
        self.polarizability * 1e-30 * 1.112_650_056e-10
    }
    /// Derive Drude charge from spring constant and polarizability.
    ///
    /// α = q_D² / k_D  ⟹  q_D = sqrt(α k_D).
    /// The returned charge is in units consistent with input (kcal mol^-1 Å^-2, Å³, e).
    pub fn derived_drude_charge(polarizability: f64, k_drude: f64) -> f64 {
        (polarizability * k_drude).sqrt()
    }
    /// Typical CHARMM Drude parameters for a methyl carbon (CG331).
    pub fn cg331() -> Self {
        Self::new("CG331", -1.051, 1000.0, 1.105, 1.3)
    }
    /// Typical CHARMM Drude parameters for an oxygen (OG311).
    pub fn og311() -> Self {
        Self::new("OG311", -0.774, 500.0, 0.591, 1.3)
    }
    /// Typical CHARMM Drude parameters for a nitrogen (NG2S1).
    pub fn ng2s1() -> Self {
        Self::new("NG2S1", -0.843, 700.0, 1.013, 1.3)
    }
}
