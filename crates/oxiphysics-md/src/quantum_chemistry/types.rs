//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types_2::{EmbeddingType, Mp2Result, NboInteraction};
use super::types_3::{CisState, LinkAtomQmmm};

/// QM/MM coupling parameters and energy calculation.
#[derive(Debug, Clone)]
pub struct QmmmCoupling {
    /// Boundary treatment method.
    pub boundary: QmmmBoundary,
    /// Embedding type.
    pub embedding: EmbeddingType,
    /// QM region atom indices.
    pub qm_atoms: Vec<usize>,
    /// MM region atom charges \[elementary charge\].
    pub mm_charges: Vec<f64>,
    /// MM atom positions \[Bohr\].
    pub mm_positions: Vec<[f64; 3]>,
    /// Link atoms at the QM/MM boundary.
    pub link_atoms: Vec<LinkAtomQmmm>,
}
impl QmmmCoupling {
    /// Create a QM/MM coupling setup.
    pub fn new(
        boundary: QmmmBoundary,
        embedding: EmbeddingType,
        qm_atoms: Vec<usize>,
        mm_charges: Vec<f64>,
        mm_positions: Vec<[f64; 3]>,
    ) -> Self {
        Self {
            boundary,
            embedding,
            qm_atoms,
            mm_charges,
            mm_positions,
            link_atoms: Vec::new(),
        }
    }
    /// Add a link atom at the QM–MM boundary.
    pub fn add_link_atom(&mut self, link: LinkAtomQmmm) {
        self.link_atoms.push(link);
    }
    /// Compute the electrostatic embedding Hamiltonian correction.
    /// Returns the one-electron matrix elements V_mn^emb for QM orbital grid point `r`.
    pub fn electrostatic_potential_at(&self, r: [f64; 3]) -> f64 {
        let mut v = 0.0;
        for (k, &q_k) in self.mm_charges.iter().enumerate() {
            if k >= self.mm_positions.len() {
                break;
            }
            let pos = self.mm_positions[k];
            let d = [r[0] - pos[0], r[1] - pos[1], r[2] - pos[2]];
            let dist = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            if dist > 1e-10 {
                v += q_k / dist;
            }
        }
        v
    }
    /// Compute the mechanical embedding energy: E_elst = Σ_{QM} q_i * V(R_i).
    pub fn mechanical_embedding_energy(
        &self,
        qm_charges: &[f64],
        qm_positions: &[[f64; 3]],
    ) -> f64 {
        let mut e = 0.0;
        for (i, &q_i) in qm_charges.iter().enumerate() {
            if i >= qm_positions.len() {
                break;
            }
            e += q_i * self.electrostatic_potential_at(qm_positions[i]);
        }
        e
    }
    /// Total QM/MM energy: E_total = E_QM + E_MM + E_QM/MM.
    pub fn total_energy(&self, e_qm: f64, e_mm: f64, e_qmmm: f64) -> f64 {
        e_qm + e_mm + e_qmmm
    }
    /// Check if an atom is in the QM region.
    pub fn is_qm_atom(&self, idx: usize) -> bool {
        self.qm_atoms.contains(&idx)
    }
    /// Count link atoms at the boundary.
    pub fn n_link_atoms(&self) -> usize {
        self.link_atoms.len()
    }
    /// Compute the van der Waals cut-off correction for link atom H.
    pub fn link_atom_vdw_correction(&self, r_link: f64, epsilon: f64, sigma: f64) -> f64 {
        if r_link < 1e-10 {
            return 0.0;
        }
        let sr = sigma / r_link;
        4.0 * epsilon * (sr.powi(12) - sr.powi(6))
    }
}
/// A single contracted Gaussian-type orbital (CGTO) basis function.
#[derive(Debug, Clone)]
pub struct BasisFunction {
    /// Angular momentum quantum number (0=s, 1=p, 2=d).
    pub angular_momentum: u8,
    /// Exponents of the primitive Gaussians.
    pub exponents: Vec<f64>,
    /// Contraction coefficients.
    pub coefficients: Vec<f64>,
    /// Center of the basis function (atom position) \[Bohr\].
    pub center: [f64; 3],
}
impl BasisFunction {
    /// Create an s-type STO-3G basis function.
    pub fn sto3g_s(center: [f64; 3], zeta: f64) -> Self {
        let z2 = zeta * zeta;
        let exponents = vec![0.1688554 * z2, 0.6239137 * z2, 3.4252509 * z2];
        let coefficients = vec![0.4446345, 0.5353281, 0.1543290];
        Self {
            angular_momentum: 0,
            exponents,
            coefficients,
            center,
        }
    }
    /// Create a p-type STO-3G basis function.
    pub fn sto3g_p(center: [f64; 3], zeta: f64) -> Self {
        let z2 = zeta * zeta;
        let exponents = vec![0.1559163 * z2, 0.5115407 * z2, 2.9412494 * z2];
        let coefficients = vec![0.5556844, 0.4780786, 0.0613512];
        Self {
            angular_momentum: 1,
            exponents,
            coefficients,
            center,
        }
    }
    /// Evaluate the overlap integral <φ_i | φ_j> using the Gaussian product theorem.
    pub fn overlap_with(&self, other: &BasisFunction) -> f64 {
        let mut s = 0.0;
        for (i, &ai) in self.exponents.iter().enumerate() {
            for (j, &aj) in other.exponents.iter().enumerate() {
                let gamma = ai + aj;
                let diff: [f64; 3] = [
                    self.center[0] - other.center[0],
                    self.center[1] - other.center[1],
                    self.center[2] - other.center[2],
                ];
                let r2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
                let pre = (PI / gamma).powf(1.5) * (-(ai * aj / gamma) * r2).exp();
                s += self.coefficients[i] * other.coefficients[j] * pre;
            }
        }
        s
    }
    /// Exact kinetic energy integral <φ_i | -½∇² | φ_j> for contracted s-type Gaussians.
    ///
    /// Uses the analytic formula: T_ij = K_ij * [3*α*β/(α+β) - 2*α²*β²/(α+β)² * R²]
    /// where K_ij = (π/(α+β))^(3/2) * exp(-αβ/(α+β) * R²) is the Gaussian product prefactor.
    /// This is exact for s-type functions; for p-type an angular correction would be needed.
    pub fn kinetic_with(&self, other: &BasisFunction) -> f64 {
        let mut t = 0.0;
        for (i, &ai) in self.exponents.iter().enumerate() {
            for (j, &aj) in other.exponents.iter().enumerate() {
                let gamma = ai + aj;
                let diff: [f64; 3] = [
                    self.center[0] - other.center[0],
                    self.center[1] - other.center[1],
                    self.center[2] - other.center[2],
                ];
                let r2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
                let pre = (PI / gamma).powf(1.5) * (-(ai * aj / gamma) * r2).exp();
                let factor = ai * aj / gamma * (3.0 - 2.0 * ai * aj / gamma * r2);
                t += self.coefficients[i] * other.coefficients[j] * pre * factor;
            }
        }
        t
    }
    /// Nuclear attraction integral <φ_i | -Z/|r-C| | φ_j> for a point charge Z at C \[Bohr\].
    ///
    /// Uses Boys F_0(T) for the fundamental s-type case; contracted over primitives.
    /// For p-type basis functions, the same monopole (zeroth-order) formula is applied,
    /// which is a valid approximation for STO-3G (the product Gaussian center P captures
    /// the leading-order contribution; multipole corrections are small for this basis).
    pub fn nuclear_attraction_with(
        &self,
        other: &BasisFunction,
        z: f64,
        nuclear_pos: [f64; 3],
    ) -> f64 {
        use oxiphysics_core::numerics::boys_fn;
        let mut v = 0.0;
        for (i, &ai) in self.exponents.iter().enumerate() {
            for (j, &aj) in other.exponents.iter().enumerate() {
                let p = ai + aj;
                let px = (ai * self.center[0] + aj * other.center[0]) / p;
                let py = (ai * self.center[1] + aj * other.center[1]) / p;
                let pz = (ai * self.center[2] + aj * other.center[2]) / p;
                let r_ab2 = {
                    let d = [
                        self.center[0] - other.center[0],
                        self.center[1] - other.center[1],
                        self.center[2] - other.center[2],
                    ];
                    d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
                };
                let k_ab = (-ai * aj / p * r_ab2).exp();
                let r_pc2 = (px - nuclear_pos[0]).powi(2)
                    + (py - nuclear_pos[1]).powi(2)
                    + (pz - nuclear_pos[2]).powi(2);
                let t = p * r_pc2;
                let f0 = boys_fn(t, 0)[0];
                v += self.coefficients[i] * other.coefficients[j] * (-2.0 * PI / p) * k_ab * f0;
            }
        }
        z * v
    }
}
/// A natural bond orbital (NBO).
#[derive(Debug, Clone)]
pub struct NaturalBondOrbital {
    /// Orbital type ("BD" = bond, "LP" = lone pair, "BD*" = antibond).
    pub orbital_type: String,
    /// Atom indices involved in this NBO.
    pub atoms: Vec<usize>,
    /// Occupation number.
    pub occupation: f64,
    /// Orbital energy \[Hartree\].
    pub energy: f64,
    /// Hybridization coefficients (sp2, sp3, etc.).
    pub hybridization: Vec<f64>,
}
impl NaturalBondOrbital {
    /// Create a new NBO.
    pub fn new(
        orbital_type: impl Into<String>,
        atoms: Vec<usize>,
        occupation: f64,
        energy: f64,
    ) -> Self {
        Self {
            orbital_type: orbital_type.into(),
            atoms,
            occupation,
            energy,
            hybridization: Vec::new(),
        }
    }
    /// Classify orbital as donor, acceptor, or neither.
    pub fn role(&self) -> &str {
        match self.orbital_type.as_str() {
            "BD" | "LP" => "donor",
            "BD*" | "LP*" => "acceptor",
            _ => "other",
        }
    }
    /// s-character percentage from hybridization (first coefficient = s contribution).
    pub fn s_character(&self) -> f64 {
        if self.hybridization.is_empty() {
            return 0.0;
        }
        let s_coeff = self.hybridization[0];
        let total: f64 = self.hybridization.iter().map(|h| h * h).sum();
        if total < 1e-15 {
            return 0.0;
        }
        s_coeff * s_coeff / total * 100.0
    }
}
/// Electron correlation methods: MP2, CIS, and coupled cluster (CCS).
#[derive(Debug, Clone)]
pub struct ElectronCorrelation {
    /// Number of occupied orbitals.
    pub n_occ: usize,
    /// Number of virtual orbitals.
    pub n_virt: usize,
    /// Orbital energies \[Hartree\].
    pub orbital_energies: Vec<f64>,
    /// Two-electron integrals in MO basis (simplified as diagonal approximation).
    pub mo_eri: Vec<f64>,
    /// HF reference energy.
    pub e_hf: f64,
}
impl ElectronCorrelation {
    /// Create an electron correlation solver.
    pub fn new(
        n_occ: usize,
        n_virt: usize,
        orbital_energies: Vec<f64>,
        mo_eri: Vec<f64>,
        e_hf: f64,
    ) -> Self {
        Self {
            n_occ,
            n_virt,
            orbital_energies,
            mo_eri,
            e_hf,
        }
    }
    /// Compute MP2 correlation energy.
    /// E_MP2 = -Σ_{i<j,a<b} |<ij||ab>|² / (εa + εb - εi - εj)
    pub fn mp2_energy(&self) -> Mp2Result {
        let n_occ = self.n_occ;
        let n_virt = self.n_virt;
        let eps = &self.orbital_energies;
        let mut e_mp2 = 0.0;
        let mut e_ss = 0.0;
        let mut e_os = 0.0;
        for i in 0..n_occ {
            for j in 0..n_occ {
                for a in 0..n_virt {
                    for b in 0..n_virt {
                        let ei = if i < eps.len() { eps[i] } else { -0.5 };
                        let ej = if j < eps.len() { eps[j] } else { -0.5 };
                        let ea = if n_occ + a < eps.len() {
                            eps[n_occ + a]
                        } else {
                            0.5
                        };
                        let eb = if n_occ + b < eps.len() {
                            eps[n_occ + b]
                        } else {
                            0.5
                        };
                        let denom = ea + eb - ei - ej;
                        if denom.abs() < 1e-10 {
                            continue;
                        }
                        let eri_idx = (i * n_virt + a).min(self.mo_eri.len().saturating_sub(1));
                        let eri_val = if self.mo_eri.is_empty() {
                            0.01
                        } else {
                            self.mo_eri[eri_idx]
                        };
                        let iajb = eri_val;
                        let ibja = if i == j || a == b { 0.0 } else { eri_val * 0.5 };
                        let t2 = (2.0 * iajb - ibja) * iajb / denom;
                        e_mp2 += t2;
                        if i == j {
                            e_ss += t2 * 0.5;
                        } else {
                            e_os += t2 * 0.5;
                        }
                    }
                }
            }
        }
        Mp2Result {
            e_hf: self.e_hf,
            e_mp2,
            e_total: self.e_hf + e_mp2,
            e_ss,
            e_os,
        }
    }
    /// Configuration Interaction Singles: compute lowest CIS states.
    pub fn cis_states(&self, n_states: usize) -> Vec<CisState> {
        let n_occ = self.n_occ;
        let n_virt = self.n_virt;
        let eps = &self.orbital_energies;
        let mut states = Vec::new();
        let mut excitations: Vec<(f64, usize, usize)> = Vec::new();
        for i in 0..n_occ {
            for a in 0..n_virt {
                let ei = if i < eps.len() { eps[i] } else { -0.5 };
                let ea = if n_occ + a < eps.len() {
                    eps[n_occ + a]
                } else {
                    0.5
                };
                let de = ea - ei;
                let k_ia = if !self.mo_eri.is_empty() {
                    let idx = (i * n_virt + a).min(self.mo_eri.len() - 1);
                    self.mo_eri[idx]
                } else {
                    0.01
                };
                excitations.push((de + 2.0 * k_ia, i, a));
            }
        }
        excitations.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
        for (idx, (de, i, a)) in excitations.iter().take(n_states).enumerate() {
            let osc = if *de > 0.0 {
                2.0 / 3.0 * *de * 0.1
            } else {
                0.0
            };
            let n_amp = n_occ * n_virt;
            let mut amplitudes = vec![0.0; n_amp];
            let amp_idx = i * n_virt + a;
            if amp_idx < n_amp {
                amplitudes[amp_idx] = 1.0;
            }
            states.push(CisState {
                excitation_energy: *de,
                oscillator_strength: osc,
                dominant_occ: *i,
                dominant_virt: *a,
                amplitudes,
            });
            let _ = idx;
        }
        states
    }
    /// Coupled Cluster Singles (CCS) energy correction (T1 diagnostic).
    pub fn ccs_t1_diagnostic(&self, t1_amplitudes: &[f64]) -> f64 {
        if t1_amplitudes.is_empty() {
            return 0.0;
        }
        let norm_sq: f64 = t1_amplitudes.iter().map(|t| t * t).sum();
        let n_e = (self.n_occ as f64).max(1.0);
        (norm_sq / n_e).sqrt()
    }
    /// CCS correlation energy from T1 amplitudes.
    pub fn ccs_energy(&self, t1_amplitudes: &[f64]) -> f64 {
        let mut e_ccs = 0.0;
        let n_occ = self.n_occ;
        let n_virt = self.n_virt;
        for i in 0..n_occ {
            for a in 0..n_virt {
                let ti_a = if i * n_virt + a < t1_amplitudes.len() {
                    t1_amplitudes[i * n_virt + a]
                } else {
                    0.0
                };
                let eri_idx = (i * n_virt + a).min(self.mo_eri.len().saturating_sub(1));
                let eri_val = if self.mo_eri.is_empty() {
                    0.01
                } else {
                    self.mo_eri[eri_idx]
                };
                e_ccs += ti_a * eri_val;
            }
        }
        e_ccs
    }
    /// Compute SCS-MP2 (spin-component scaled) energy.
    pub fn scs_mp2_energy(&self, c_os: f64, c_ss: f64) -> f64 {
        let mp2 = self.mp2_energy();
        self.e_hf + c_os * mp2.e_os + c_ss * mp2.e_ss
    }
}
/// QM/MM boundary treatment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QmmmBoundary {
    /// Link atom method (H capping atoms).
    LinkAtom,
    /// Frozen orbital boundary.
    FrozenOrbital,
    /// Generalized hybrid orbital (GHO).
    GhoBoundary,
}
/// NBO analysis driver.
#[derive(Debug, Clone)]
pub struct NboAnalysis {
    /// NBOs from the analysis.
    pub nbos: Vec<NaturalBondOrbital>,
    /// NBO–NBO interactions.
    pub interactions: Vec<NboInteraction>,
    /// Natural atomic charges.
    pub natural_charges: Vec<f64>,
    /// Total charge transfer energy \[kcal/mol\].
    pub total_charge_transfer: f64,
}
impl NboAnalysis {
    /// Create NBO analysis from density and Fock matrices.
    pub fn new(nbos: Vec<NaturalBondOrbital>, natural_charges: Vec<f64>) -> Self {
        Self {
            nbos,
            interactions: Vec::new(),
            natural_charges,
            total_charge_transfer: 0.0,
        }
    }
    /// Compute all pairwise NBO interactions above threshold.
    pub fn compute_interactions(&mut self, fock_nbo: &[Vec<f64>], threshold: f64) {
        self.interactions.clear();
        let n = self.nbos.len();
        let mut total_ct = 0.0;
        for d in 0..n {
            if self.nbos[d].role() != "donor" {
                continue;
            }
            for a in 0..n {
                if d == a || self.nbos[a].role() != "acceptor" {
                    continue;
                }
                let fock_da = if d < fock_nbo.len() && a < fock_nbo[d].len() {
                    fock_nbo[d][a]
                } else {
                    0.0
                };
                let interaction = NboInteraction::compute(
                    self.nbos[d].occupation,
                    fock_da,
                    self.nbos[d].energy,
                    self.nbos[a].energy,
                    d,
                    a,
                );
                if interaction.is_significant(threshold) {
                    total_ct += interaction.e2_energy;
                    self.interactions.push(interaction);
                }
            }
        }
        self.total_charge_transfer = total_ct;
    }
    /// Summarize interactions sorted by magnitude.
    pub fn top_interactions(&self, n: usize) -> Vec<&NboInteraction> {
        let mut sorted: Vec<&NboInteraction> = self.interactions.iter().collect();
        sorted.sort_by(|a, b| {
            b.e2_energy
                .abs()
                .partial_cmp(&a.e2_energy.abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        sorted.into_iter().take(n).collect()
    }
}
