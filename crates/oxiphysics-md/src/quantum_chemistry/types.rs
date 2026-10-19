//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::integrals::{
    kinetic_integral, nuclear_attraction_integral, odd_double_factorial, overlap_integral,
    primitive_normalisation,
};
use super::types_2::{EmbeddingType, NboInteraction};
use super::types_3::LinkAtomQmmm;

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
/// Cartesian direction of a p-type basis function.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CartesianAxis {
    /// x direction.
    X,
    /// y direction.
    Y,
    /// z direction.
    Z,
}
impl CartesianAxis {
    /// The three axes in x, y, z order.
    pub const ALL: [CartesianAxis; 3] = [CartesianAxis::X, CartesianAxis::Y, CartesianAxis::Z];
    /// Index of the axis (0 = x, 1 = y, 2 = z).
    pub fn index(self) -> usize {
        match self {
            CartesianAxis::X => 0,
            CartesianAxis::Y => 1,
            CartesianAxis::Z => 2,
        }
    }
    /// Cartesian powers `(l, m, n)` of the p function `x`, `y` or `z` along this axis.
    pub fn p_powers(self) -> [u8; 3] {
        let mut powers = [0u8; 3];
        powers[self.index()] = 1;
        powers
    }
}
/// A single contracted Cartesian Gaussian-type orbital (CGTO) basis function
/// `φ(r) = Σ_k c_k (x−A_x)^l (y−A_y)^m (z−A_z)^n exp(−α_k |r−A|²)`.
#[derive(Debug, Clone)]
pub struct BasisFunction {
    /// Angular momentum quantum number (0=s, 1=p, 2=d); equals `l + m + n`.
    pub angular_momentum: u8,
    /// Cartesian powers `(l, m, n)` of the angular factor `x^l y^m z^n`
    /// (e.g. `[1, 0, 0]` for p_x).
    pub cartesian: [u8; 3],
    /// Exponents of the primitive Gaussians.
    pub exponents: Vec<f64>,
    /// Contraction coefficients for the *raw* primitives `x^l y^m z^n e^{-αr²}`
    /// (primitive normalisation and contraction normalisation included).
    pub coefficients: Vec<f64>,
    /// Center of the basis function (atom position) \[Bohr\].
    pub center: [f64; 3],
}
/// STO-3G least-squares expansion of a 1s Slater orbital with ζ = 1
/// (Hehre, Stewart & Pople 1969): primitive exponents.
const STO3G_1S_EXPONENTS: [f64; 3] = [0.109_817_510_4, 0.405_771_156_2, 2.227_660_584];
/// STO-3G 1s contraction coefficients (for normalised primitives).
const STO3G_1S_COEFFICIENTS: [f64; 3] = [0.444_634_542_2, 0.535_328_142_3, 0.154_328_967_3];
/// STO-3G expansion of the 2s/2p Slater orbitals with ζ = 1: shared primitive exponents.
const STO3G_2SP_EXPONENTS: [f64; 3] = [0.075_138_585_7, 0.231_031_333, 0.994_202_729];
/// STO-3G 2s contraction coefficients (for normalised primitives).
const STO3G_2S_COEFFICIENTS: [f64; 3] = [0.700_115_468_9, 0.399_512_826_1, -0.099_967_229_2];
/// STO-3G 2p contraction coefficients (for normalised primitives).
const STO3G_2P_COEFFICIENTS: [f64; 3] = [0.391_957_393_1, 0.607_683_718_6, 0.155_916_275];

/// Turn tabulated contraction coefficients (defined for *normalised* primitives)
/// into coefficients for the raw primitives `x^l y^m z^n e^{-αr²}` used by the
/// integral code, and renormalise the whole contraction to unit self-overlap.
///
/// Primitive normalisation:
/// N(α) = (2α/π)^{3/4} (4α)^{(l+m+n)/2} / √((2l−1)!!(2m−1)!!(2n−1)!!).
fn normalised_contraction(exponents: &[f64], coefficients: &[f64], cartesian: [u8; 3]) -> Vec<f64> {
    let mut c: Vec<f64> = exponents
        .iter()
        .zip(coefficients.iter())
        .map(|(&a, &ci)| ci * primitive_normalisation(a, cartesian))
        .collect();
    // Same-centre self-overlap: Σ c_i c_j (π/γ)^{3/2} Π_k (2l_k−1)!!/(2γ)^{l_k}.
    let mut self_overlap = 0.0;
    for (i, &ai) in exponents.iter().enumerate() {
        for (j, &aj) in exponents.iter().enumerate() {
            let gamma = ai + aj;
            let radial = (PI / gamma).powf(1.5);
            let angular: f64 = cartesian
                .iter()
                .map(|&l| odd_double_factorial(l) / (2.0 * gamma).powi(i32::from(l)))
                .product();
            self_overlap += c[i] * c[j] * radial * angular;
        }
    }
    if self_overlap > 0.0 {
        let scale = 1.0 / self_overlap.sqrt();
        for ci in &mut c {
            *ci *= scale;
        }
    }
    c
}

impl BasisFunction {
    /// Build a normalised contracted Cartesian Gaussian centred at `center`
    /// \[Bohr\] with angular factor `x^l y^m z^n` (`cartesian = [l, m, n]`).
    ///
    /// `coefficients` are the tabulated contraction weights of *normalised*
    /// primitives (the usual basis-set-library convention); the stored
    /// coefficients include the primitive normalisation and are rescaled so the
    /// contracted function has unit norm.  Extra entries of the longer slice
    /// are ignored.
    pub fn contracted(
        center: [f64; 3],
        cartesian: [u8; 3],
        exponents: &[f64],
        coefficients: &[f64],
    ) -> Self {
        let n = exponents.len().min(coefficients.len());
        let exponents = exponents[..n].to_vec();
        let coefficients = normalised_contraction(&exponents, &coefficients[..n], cartesian);
        Self {
            angular_momentum: cartesian[0] + cartesian[1] + cartesian[2],
            cartesian,
            exponents,
            coefficients,
            center,
        }
    }
    /// Create a 1s STO-3G basis function approximating a 1s Slater orbital of
    /// exponent `zeta` (Bohr⁻¹), e.g. ζ = 1.24 for molecular hydrogen.
    ///
    /// The primitive exponents are α_k = α_k(ζ=1)·ζ².  The stored `coefficients`
    /// already contain the primitive normalisation constants and the overall
    /// contraction normalisation, so that `overlap_with(self) == 1`.
    pub fn sto3g_s(center: [f64; 3], zeta: f64) -> Self {
        let z2 = zeta * zeta;
        let exponents: Vec<f64> = STO3G_1S_EXPONENTS.iter().map(|a| a * z2).collect();
        Self::contracted(center, [0, 0, 0], &exponents, &STO3G_1S_COEFFICIENTS)
    }
    /// Create a 2s STO-3G basis function approximating a 2s Slater orbital of
    /// exponent `zeta` (Bohr⁻¹).  Shares its primitive exponents with
    /// [`sto3g_p`](Self::sto3g_p) (the STO-3G "sp" shell).
    ///
    /// Coefficients are normalised as for [`sto3g_s`](Self::sto3g_s).
    pub fn sto3g_2s(center: [f64; 3], zeta: f64) -> Self {
        let z2 = zeta * zeta;
        let exponents: Vec<f64> = STO3G_2SP_EXPONENTS.iter().map(|a| a * z2).collect();
        Self::contracted(center, [0, 0, 0], &exponents, &STO3G_2S_COEFFICIENTS)
    }
    /// Create the `axis` Cartesian component (p_x, p_y or p_z) of a 2p STO-3G
    /// basis function approximating a 2p Slater orbital of exponent `zeta` (Bohr⁻¹).
    ///
    /// The stored `coefficients` include the p-type primitive normalisation
    /// (2α/π)^{3/4}·2√α and the contraction normalisation, so the component
    /// (e.g. x·Σ c_k e^{-α_k r²}) has unit norm.
    pub fn sto3g_p(center: [f64; 3], zeta: f64, axis: CartesianAxis) -> Self {
        let z2 = zeta * zeta;
        let exponents: Vec<f64> = STO3G_2SP_EXPONENTS.iter().map(|a| a * z2).collect();
        Self::contracted(center, axis.p_powers(), &exponents, &STO3G_2P_COEFFICIENTS)
    }
    /// Value φ(r) and gradient ∇φ(r) of the basis function at `r` \[Bohr\].
    pub fn value_and_gradient(&self, r: [f64; 3]) -> (f64, [f64; 3]) {
        let d = [
            r[0] - self.center[0],
            r[1] - self.center[1],
            r[2] - self.center[2],
        ];
        let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        // Radial sums Σ c e^{-αr²} and Σ c α e^{-αr²}.
        let mut radial = 0.0;
        let mut radial_alpha = 0.0;
        for (&a, &c) in self.exponents.iter().zip(self.coefficients.iter()) {
            let g = c * (-a * r2).exp();
            radial += g;
            radial_alpha += a * g;
        }
        let pow = |x: f64, n: u8| -> f64 { if n == 0 { 1.0 } else { x.powi(i32::from(n)) } };
        let ang = [
            pow(d[0], self.cartesian[0]),
            pow(d[1], self.cartesian[1]),
            pow(d[2], self.cartesian[2]),
        ];
        let angular = ang[0] * ang[1] * ang[2];
        let mut grad = [0.0_f64; 3];
        for k in 0..3 {
            let l = self.cartesian[k];
            // ∂/∂x_k [x_k^l] = l x_k^{l−1}
            let d_ang = if l == 0 {
                0.0
            } else {
                f64::from(l) * pow(d[k], l - 1)
            };
            let others: f64 = (0..3).filter(|&m| m != k).map(|m| ang[m]).product();
            grad[k] = d_ang * others * radial - 2.0 * d[k] * angular * radial_alpha;
        }
        (angular * radial, grad)
    }
    /// Overlap integral <φ_i | φ_j> (exact, McMurchie–Davidson).
    pub fn overlap_with(&self, other: &BasisFunction) -> f64 {
        overlap_integral(self, other)
    }
    /// Kinetic-energy integral <φ_i | -½∇² | φ_j> (exact for any angular momentum).
    pub fn kinetic_with(&self, other: &BasisFunction) -> f64 {
        kinetic_integral(self, other)
    }
    /// Nuclear attraction integral <φ_i | -Z/|r-C| | φ_j> for a point charge Z at C \[Bohr\].
    ///
    /// Exact for any angular momentum (Hermite–Coulomb integrals built on the
    /// Boys function).
    pub fn nuclear_attraction_with(
        &self,
        other: &BasisFunction,
        z: f64,
        nuclear_pos: [f64; 3],
    ) -> f64 {
        nuclear_attraction_integral(self, other, z, nuclear_pos)
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
