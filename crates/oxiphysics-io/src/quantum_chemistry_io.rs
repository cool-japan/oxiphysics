// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Quantum chemistry file format I/O.
//!
//! Implements readers and writers for the most common quantum chemistry
//! file formats: Gaussian (.com/.gjf/.log), ORCA (.inp/.out), along with
//! data structures for molecular orbitals, basis sets, normal modes,
//! electronic structure, excited states, and NBO analysis.

use std::collections::HashMap;
use std::fmt;

// ---------------------------------------------------------------------------
// Shared geometry types
// ---------------------------------------------------------------------------

/// A single atom with element symbol and Cartesian coordinates (Å).
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    /// Element symbol (e.g. "C", "H", "N").
    pub symbol: String,
    /// Cartesian coordinates in Ångströms: `[x, y, z]`.
    pub coords: [f64; 3],
    /// Optional atomic number (0 if unset).
    pub atomic_number: u8,
}

impl Atom {
    /// Create a new atom from symbol and coordinates.
    pub fn new(symbol: impl Into<String>, x: f64, y: f64, z: f64) -> Self {
        Self {
            symbol: symbol.into(),
            coords: [x, y, z],
            atomic_number: 0,
        }
    }

    /// Set the atomic number.
    pub fn with_atomic_number(mut self, z: u8) -> Self {
        self.atomic_number = z;
        self
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:<4} {:>12.6} {:>12.6} {:>12.6}",
            self.symbol, self.coords[0], self.coords[1], self.coords[2]
        )
    }
}

// ---------------------------------------------------------------------------
// GaussianInput — .com / .gjf format
// ---------------------------------------------------------------------------

/// Route section options for a Gaussian calculation.
#[derive(Debug, Clone)]
pub struct GaussianRouteSection {
    /// Method keyword (e.g., "B3LYP", "MP2", "CCSD(T)").
    pub method: String,
    /// Basis set keyword (e.g., "6-31G*", "cc-pVDZ").
    pub basis_set: String,
    /// Calculation type (e.g., "SP", "Opt", "Freq", "NMR").
    pub calc_type: String,
    /// Additional keywords as key=value pairs.
    pub keywords: HashMap<String, String>,
    /// Raw extra route tokens (e.g., "Geom=CheckPoint").
    pub extra_tokens: Vec<String>,
}

impl GaussianRouteSection {
    /// Create a new route section.
    pub fn new(
        method: impl Into<String>,
        basis_set: impl Into<String>,
        calc_type: impl Into<String>,
    ) -> Self {
        Self {
            method: method.into(),
            basis_set: basis_set.into(),
            calc_type: calc_type.into(),
            keywords: HashMap::new(),
            extra_tokens: Vec::new(),
        }
    }

    /// Add a keyword.
    pub fn add_keyword(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.keywords.insert(key.into(), value.into());
    }

    /// Render the route section as a Gaussian `#` line.
    pub fn render(&self) -> String {
        let mut parts = vec![format!(
            "# {}/{} {}",
            self.method, self.basis_set, self.calc_type
        )];
        for (k, v) in &self.keywords {
            if v.is_empty() {
                parts.push(k.clone());
            } else {
                parts.push(format!("{}={}", k, v));
            }
        }
        parts.extend(self.extra_tokens.iter().cloned());
        parts.join(" ")
    }
}

/// A complete Gaussian input file structure (.com / .gjf).
#[derive(Debug, Clone)]
pub struct GaussianInput {
    /// Link0 commands (e.g., `%chk=job.chk`).
    pub link0: Vec<String>,
    /// Route section.
    pub route: GaussianRouteSection,
    /// Job title / comment line.
    pub title: String,
    /// Molecular charge.
    pub charge: i32,
    /// Spin multiplicity (2S+1).
    pub multiplicity: u32,
    /// Cartesian atom list.
    pub atoms: Vec<Atom>,
    /// Optional Z-matrix string (used instead of Cartesian coordinates).
    pub z_matrix: Option<String>,
}

impl GaussianInput {
    /// Create a new Gaussian input with Cartesian coordinates.
    pub fn new(
        route: GaussianRouteSection,
        title: impl Into<String>,
        charge: i32,
        multiplicity: u32,
    ) -> Self {
        Self {
            link0: Vec::new(),
            route,
            title: title.into(),
            charge,
            multiplicity,
            atoms: Vec::new(),
            z_matrix: None,
        }
    }

    /// Add an atom to the geometry.
    pub fn add_atom(&mut self, atom: Atom) {
        self.atoms.push(atom);
    }

    /// Add a `%link0` directive (e.g., `%chk=job.chk`).
    pub fn add_link0(&mut self, directive: impl Into<String>) {
        self.link0.push(directive.into());
    }

    /// Render the full Gaussian input file as a string.
    pub fn render(&self) -> String {
        let mut out = String::new();
        for l0 in &self.link0 {
            out.push('%');
            out.push_str(l0);
            out.push('\n');
        }
        out.push_str(&self.route.render());
        out.push_str("\n\n");
        out.push_str(&self.title);
        out.push_str("\n\n");
        out.push_str(&format!("{} {}\n", self.charge, self.multiplicity));
        if let Some(zmat) = &self.z_matrix {
            out.push_str(zmat);
        } else {
            for atom in &self.atoms {
                out.push_str(&atom.to_string());
                out.push('\n');
            }
        }
        out.push('\n');
        out
    }

    /// Count the number of atoms.
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }
}

// ---------------------------------------------------------------------------
// GaussianOutput — parse .log files
// ---------------------------------------------------------------------------

/// A geometry optimization step extracted from a Gaussian log file.
#[derive(Debug, Clone)]
pub struct OptStep {
    /// Step number (1-based).
    pub step: usize,
    /// SCF energy at this geometry (Hartree).
    pub energy: f64,
    /// Maximum force component.
    pub max_force: f64,
    /// RMS force.
    pub rms_force: f64,
    /// Maximum displacement.
    pub max_displacement: f64,
    /// RMS displacement.
    pub rms_displacement: f64,
    /// Converged?
    pub converged: bool,
    /// Atomic coordinates at this step.
    pub atoms: Vec<Atom>,
}

impl OptStep {
    /// Create an optimization step record.
    pub fn new(step: usize, energy: f64) -> Self {
        Self {
            step,
            energy,
            max_force: 0.0,
            rms_force: 0.0,
            max_displacement: 0.0,
            rms_displacement: 0.0,
            converged: false,
            atoms: Vec::new(),
        }
    }
}

/// Parsed frequency data for a single normal mode.
#[derive(Debug, Clone)]
pub struct GaussianFrequency {
    /// Frequency in cm⁻¹.
    pub freq_cm: f64,
    /// IR intensity (km/mol).
    pub ir_intensity: f64,
    /// Raman activity (Å⁴/amu), if available.
    pub raman_activity: Option<f64>,
    /// Reduced mass (amu).
    pub reduced_mass: f64,
    /// Normal mode displacement vectors for each atom `[dx, dy, dz]`.
    pub displacements: Vec<[f64; 3]>,
}

impl GaussianFrequency {
    /// Create a frequency entry.
    pub fn new(freq_cm: f64, ir_intensity: f64, reduced_mass: f64) -> Self {
        Self {
            freq_cm,
            ir_intensity,
            raman_activity: None,
            reduced_mass,
            displacements: Vec::new(),
        }
    }

    /// Return `true` if this is an imaginary frequency (transition state).
    pub fn is_imaginary(&self) -> bool {
        self.freq_cm < 0.0
    }
}

/// NMR shielding tensor data for a single nucleus.
#[derive(Debug, Clone)]
pub struct NmrShielding {
    /// Atom index (0-based).
    pub atom_idx: usize,
    /// Isotropic chemical shielding (ppm).
    pub isotropic: f64,
    /// Anisotropy value (ppm).
    pub anisotropy: f64,
}

/// Parsed Gaussian output (.log) file data.
#[derive(Debug, Clone, Default)]
pub struct GaussianOutput {
    /// Final SCF energy (Hartree).
    pub scf_energy: Option<f64>,
    /// Whether the geometry optimization converged.
    pub opt_converged: bool,
    /// Optimization steps.
    pub opt_steps: Vec<OptStep>,
    /// Vibrational frequencies.
    pub frequencies: Vec<GaussianFrequency>,
    /// NMR shielding tensors.
    pub nmr_shieldings: Vec<NmrShielding>,
    /// Zero-point energy (Hartree).
    pub zero_point_energy: Option<f64>,
    /// Thermal correction to enthalpy (Hartree).
    pub thermal_correction_h: Option<f64>,
    /// Thermal correction to Gibbs free energy (Hartree).
    pub thermal_correction_g: Option<f64>,
    /// Final optimized geometry.
    pub final_geometry: Vec<Atom>,
    /// Number of imaginary frequencies (negative values).
    pub n_imaginary: usize,
}

impl GaussianOutput {
    /// Create an empty output record.
    pub fn new() -> Self {
        Self::default()
    }

    /// Parse a Gaussian log file content string and extract key data.
    ///
    /// This is a simplified parser that recognizes the most common output
    /// patterns. Returns a populated `GaussianOutput`.
    pub fn parse(content: &str) -> Self {
        let mut out = GaussianOutput::new();
        let mut current_opt_step: Option<OptStep> = None;

        for line in content.lines() {
            let trimmed = line.trim();

            // SCF energy
            if trimmed.starts_with("SCF Done:")
                && let Some(e) = parse_after_eq(trimmed)
            {
                out.scf_energy = Some(e);
            }

            // Optimization convergence
            if trimmed.contains("Optimization completed.") {
                out.opt_converged = true;
                if let Some(step) = current_opt_step.take() {
                    out.opt_steps.push(step);
                }
            }

            // Optimization step energy
            if trimmed.starts_with("Energy=")
                && let Some(e) = parse_value_after(trimmed, "Energy=")
            {
                let step_num = out.opt_steps.len() + 1;
                current_opt_step = Some(OptStep::new(step_num, e));
            }

            // Frequency
            if trimmed.starts_with("Frequencies --")
                && let Some(rest) = trimmed.strip_prefix("Frequencies --")
            {
                for tok in rest.split_whitespace() {
                    if let Ok(f) = tok.parse::<f64>() {
                        out.frequencies.push(GaussianFrequency::new(f, 0.0, 1.0));
                        if f < 0.0 {
                            out.n_imaginary += 1;
                        }
                    }
                }
            }

            // IR intensities
            if trimmed.starts_with("IR Inten    --")
                && let Some(rest) = trimmed.strip_prefix("IR Inten    --")
            {
                let intensities: Vec<f64> = rest
                    .split_whitespace()
                    .filter_map(|t| t.parse().ok())
                    .collect();
                let n_freqs = out.frequencies.len();
                let start = n_freqs.saturating_sub(intensities.len());
                for (i, &ir) in intensities.iter().enumerate() {
                    if start + i < n_freqs {
                        out.frequencies[start + i].ir_intensity = ir;
                    }
                }
            }

            // Zero-point energy
            if trimmed.starts_with("Zero-point correction=")
                && let Some(v) = parse_value_after(trimmed, "Zero-point correction=")
            {
                out.zero_point_energy = Some(v);
            }

            // Thermal correction to enthalpy
            if trimmed.starts_with("Thermal correction to Enthalpy=")
                && let Some(v) = parse_value_after(trimmed, "Thermal correction to Enthalpy=")
            {
                out.thermal_correction_h = Some(v);
            }

            // Thermal correction to Gibbs
            if trimmed.starts_with("Thermal correction to Gibbs Free Energy=")
                && let Some(v) =
                    parse_value_after(trimmed, "Thermal correction to Gibbs Free Energy=")
            {
                out.thermal_correction_g = Some(v);
            }
        }

        // Flush remaining opt step
        if let Some(step) = current_opt_step {
            out.opt_steps.push(step);
        }

        out
    }

    /// Return the lowest frequency (most negative if imaginary).
    pub fn lowest_frequency(&self) -> Option<f64> {
        self.frequencies.iter().map(|f| f.freq_cm).reduce(f64::min)
    }

    /// Return the number of frequencies stored.
    pub fn n_frequencies(&self) -> usize {
        self.frequencies.len()
    }
}

/// Parse the floating-point value after `=` in a string like `SCF Done: E(RB3LYP) = -78.123`.
fn parse_after_eq(s: &str) -> Option<f64> {
    let idx = s.rfind('=')?;
    s[idx + 1..].split_whitespace().next()?.parse().ok()
}

/// Parse a float value directly after a prefix.
fn parse_value_after(s: &str, prefix: &str) -> Option<f64> {
    s.strip_prefix(prefix)?
        .split_whitespace()
        .next()?
        .parse()
        .ok()
}

// ---------------------------------------------------------------------------
// OrcaInput — ORCA input format
// ---------------------------------------------------------------------------

/// Calculation settings block for ORCA.
#[derive(Debug, Clone)]
pub struct OrcaSettings {
    /// Maximum number of SCF iterations.
    pub max_scf_iter: usize,
    /// SCF energy convergence threshold (Eh).
    pub scf_conv: f64,
    /// Gradient convergence threshold.
    pub grad_conv: f64,
    /// Grid level (1–7 for DFT integration).
    pub grid: u8,
    /// Whether to use RIJCOSX (resolution-of-identity) approximation.
    pub rijcosx: bool,
    /// Auxiliary basis set for RI (e.g., "def2/J").
    pub aux_basis: Option<String>,
    /// Number of processors.
    pub n_procs: usize,
    /// Memory per process (MB).
    pub mem_per_proc: usize,
}

impl OrcaSettings {
    /// Create default ORCA settings.
    pub fn default_settings() -> Self {
        Self {
            max_scf_iter: 125,
            scf_conv: 1e-8,
            grad_conv: 1e-5,
            grid: 5,
            rijcosx: false,
            aux_basis: None,
            n_procs: 1,
            mem_per_proc: 2048,
        }
    }
}

/// A complete ORCA input file.
#[derive(Debug, Clone)]
pub struct OrcaInput {
    /// Method/keywords line (e.g., "! B3LYP def2-TZVP TightSCF Opt").
    pub keywords: Vec<String>,
    /// Charge of the molecule.
    pub charge: i32,
    /// Spin multiplicity.
    pub multiplicity: u32,
    /// Atom list (Cartesian coordinates, Å).
    pub atoms: Vec<Atom>,
    /// Calculation settings block.
    pub settings: OrcaSettings,
    /// Extra `%block ... end` sections as raw strings.
    pub extra_blocks: Vec<String>,
}

impl OrcaInput {
    /// Create a new ORCA input.
    pub fn new(charge: i32, multiplicity: u32) -> Self {
        Self {
            keywords: Vec::new(),
            charge,
            multiplicity,
            atoms: Vec::new(),
            settings: OrcaSettings::default_settings(),
            extra_blocks: Vec::new(),
        }
    }

    /// Add a keyword to the `!` line.
    pub fn add_keyword(&mut self, kw: impl Into<String>) {
        self.keywords.push(kw.into());
    }

    /// Add an atom.
    pub fn add_atom(&mut self, atom: Atom) {
        self.atoms.push(atom);
    }

    /// Render the full ORCA input file as a string.
    pub fn render(&self) -> String {
        let mut out = String::new();

        // Keywords line
        if !self.keywords.is_empty() {
            out.push_str("! ");
            out.push_str(&self.keywords.join(" "));
            out.push('\n');
        }

        // %pal block for parallel
        if self.settings.n_procs > 1 {
            out.push_str(&format!("%pal nprocs {} end\n", self.settings.n_procs));
        }

        // %maxcore block
        out.push_str(&format!("%maxcore {}\n", self.settings.mem_per_proc));

        // %scf block
        out.push_str(&format!(
            "%scf\n  MaxIter {}\n  ConvTol {:.2e}\nend\n",
            self.settings.max_scf_iter, self.settings.scf_conv
        ));

        // Extra blocks
        for block in &self.extra_blocks {
            out.push_str(block);
            out.push('\n');
        }

        // Coordinate block
        out.push_str(&format!("* xyz {} {}\n", self.charge, self.multiplicity));
        for atom in &self.atoms {
            out.push_str(&atom.to_string());
            out.push('\n');
        }
        out.push_str("*\n");
        out
    }
}

// ---------------------------------------------------------------------------
// OrcaOutput — parse ORCA output
// ---------------------------------------------------------------------------

/// Excited state from ORCA TD-DFT or EOM-CCSD output.
#[derive(Debug, Clone)]
pub struct OrcaExcitedState {
    /// State number (1-based).
    pub state: usize,
    /// Excitation energy in eV.
    pub energy_ev: f64,
    /// Oscillator strength (dimensionless).
    pub oscillator_strength: f64,
    /// Transition type label (e.g., "SINGLET-A").
    pub symmetry: String,
}

impl OrcaExcitedState {
    /// Create a new excited state record.
    pub fn new(
        state: usize,
        energy_ev: f64,
        oscillator_strength: f64,
        symmetry: impl Into<String>,
    ) -> Self {
        Self {
            state,
            energy_ev,
            oscillator_strength,
            symmetry: symmetry.into(),
        }
    }
}

/// Gradient vector for one atom.
#[derive(Debug, Clone)]
pub struct AtomGradient {
    /// Atom index (0-based).
    pub atom_idx: usize,
    /// Gradient components in Eh/Bohr: `[gx, gy, gz]`.
    pub gradient: [f64; 3],
}

/// Parsed ORCA output file.
#[derive(Debug, Clone, Default)]
pub struct OrcaOutput {
    /// Final electronic energy (Hartree).
    pub electronic_energy: Option<f64>,
    /// Dispersion correction (Hartree), if present.
    pub dispersion_correction: Option<f64>,
    /// Cartesian gradient for each atom.
    pub gradients: Vec<AtomGradient>,
    /// Excited states (TDDFT or EOM-CCSD).
    pub excited_states: Vec<OrcaExcitedState>,
    /// Whether the geometry optimization converged.
    pub opt_converged: bool,
    /// Number of SCF iterations performed.
    pub scf_iterations: usize,
    /// ORCA version string.
    pub version: Option<String>,
    /// Final geometry after optimization.
    pub final_geometry: Vec<Atom>,
}

impl OrcaOutput {
    /// Create an empty output record.
    pub fn new() -> Self {
        Self::default()
    }

    /// Parse an ORCA output file string and extract key results.
    pub fn parse(content: &str) -> Self {
        let mut out = OrcaOutput::new();

        for line in content.lines() {
            let trimmed = line.trim();

            // Version
            if trimmed.starts_with("Program Version") {
                out.version = Some(trimmed.to_string());
            }

            // Electronic energy
            if trimmed.starts_with("FINAL SINGLE POINT ENERGY")
                && let Some(val) = trimmed
                    .split_whitespace()
                    .last()
                    .and_then(|s| s.parse().ok())
            {
                out.electronic_energy = Some(val);
            }

            // Dispersion
            if trimmed.starts_with("Dispersion correction")
                && let Some(val) = trimmed
                    .split_whitespace()
                    .last()
                    .and_then(|s| s.parse().ok())
            {
                out.dispersion_correction = Some(val);
            }

            // Convergence
            if trimmed.contains("OPTIMIZATION HAS CONVERGED") {
                out.opt_converged = true;
            }

            // SCF iterations
            if trimmed.starts_with("SCF CONVERGED AFTER")
                && let Some(n) = trimmed
                    .split_whitespace()
                    .find(|t| t.parse::<usize>().is_ok())
                    .and_then(|t| t.parse().ok())
            {
                out.scf_iterations = n;
            }

            // Excited states
            if trimmed.starts_with("STATE") && trimmed.contains("EXCITATION ENERGY") {
                // e.g. "STATE  1:  E=   0.1234 au    3.36 eV  f= 0.0012"
                if let Some(ev) = parse_ev_from_state_line(trimmed) {
                    let osc = parse_f_from_state_line(trimmed).unwrap_or(0.0);
                    let state_num = out.excited_states.len() + 1;
                    out.excited_states
                        .push(OrcaExcitedState::new(state_num, ev, osc, "SINGLET"));
                }
            }
        }
        out
    }

    /// Return total energy including dispersion correction.
    pub fn total_energy(&self) -> Option<f64> {
        match (self.electronic_energy, self.dispersion_correction) {
            (Some(e), Some(d)) => Some(e + d),
            (Some(e), None) => Some(e),
            _ => None,
        }
    }

    /// Return the oscillator-strength-weighted centroid absorption wavelength (nm).
    pub fn absorption_centroid_nm(&self) -> Option<f64> {
        let total_osc: f64 = self
            .excited_states
            .iter()
            .map(|s| s.oscillator_strength)
            .sum();
        if total_osc < 1e-15 {
            return None;
        }
        let weighted: f64 = self
            .excited_states
            .iter()
            .map(|s| (1239.84 / s.energy_ev) * s.oscillator_strength)
            .sum();
        Some(weighted / total_osc)
    }
}

/// Parse eV value from ORCA excited state line.
fn parse_ev_from_state_line(line: &str) -> Option<f64> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    for (i, &tok) in tokens.iter().enumerate() {
        if tok == "eV" && i > 0 {
            return tokens[i - 1].parse().ok();
        }
    }
    None
}

/// Parse oscillator strength `f=` from ORCA state line.
fn parse_f_from_state_line(line: &str) -> Option<f64> {
    for part in line.split_whitespace() {
        if let Some(val) = part.strip_prefix("f=") {
            return val.parse().ok();
        }
    }
    None
}

// ---------------------------------------------------------------------------
// MolecularOrbital
// ---------------------------------------------------------------------------

/// Spin channel for a molecular orbital.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Spin {
    /// Alpha (or closed-shell) spin channel.
    Alpha,
    /// Beta spin channel.
    Beta,
}

/// A single molecular orbital with expansion coefficients over basis functions.
#[derive(Debug, Clone)]
pub struct MolecularOrbital {
    /// MO index (0-based).
    pub index: usize,
    /// Orbital energy (Hartree).
    pub energy: f64,
    /// Occupation number (0, 1, or 2).
    pub occupation: f64,
    /// Spin channel.
    pub spin: Spin,
    /// Symmetry label (e.g., "A1", "B2u").
    pub symmetry: String,
    /// Expansion coefficients over the AO basis functions.
    pub coefficients: Vec<f64>,
}

impl MolecularOrbital {
    /// Create a new MO record.
    pub fn new(index: usize, energy: f64, occupation: f64, spin: Spin) -> Self {
        Self {
            index,
            energy,
            occupation,
            spin,
            symmetry: String::new(),
            coefficients: Vec::new(),
        }
    }

    /// Return `true` if this is the HOMO (highest occupied MO in alpha channel).
    pub fn is_homo(&self, homo_idx: usize) -> bool {
        self.index == homo_idx
    }

    /// Return `true` if this is the LUMO.
    pub fn is_lumo(&self, lumo_idx: usize) -> bool {
        self.index == lumo_idx
    }

    /// Compute the norm of the MO coefficient vector (should be ~1 for an orthonormal MO).
    pub fn coeff_norm(&self) -> f64 {
        self.coefficients.iter().map(|c| c * c).sum::<f64>().sqrt()
    }
}

/// A set of molecular orbitals with helper functions for HOMO/LUMO gap.
#[derive(Debug, Clone, Default)]
pub struct MolecularOrbitalSet {
    /// Alpha MOs.
    pub alpha_mos: Vec<MolecularOrbital>,
    /// Beta MOs (empty for restricted calculations).
    pub beta_mos: Vec<MolecularOrbital>,
    /// Number of alpha electrons.
    pub n_alpha: usize,
    /// Number of beta electrons.
    pub n_beta: usize,
}

impl MolecularOrbitalSet {
    /// Create an empty MO set.
    pub fn new(n_alpha: usize, n_beta: usize) -> Self {
        Self {
            alpha_mos: Vec::new(),
            beta_mos: Vec::new(),
            n_alpha,
            n_beta,
        }
    }

    /// Index of the alpha HOMO (0-based).
    pub fn homo_idx(&self) -> Option<usize> {
        if self.n_alpha == 0 {
            None
        } else {
            Some(self.n_alpha - 1)
        }
    }

    /// Index of the alpha LUMO.
    pub fn lumo_idx(&self) -> Option<usize> {
        if self.n_alpha < self.alpha_mos.len() {
            Some(self.n_alpha)
        } else {
            None
        }
    }

    /// HOMO–LUMO gap in eV (None if LUMO is not present).
    pub fn homo_lumo_gap_ev(&self) -> Option<f64> {
        let homo = self.homo_idx()?;
        let lumo = self.lumo_idx()?;
        let e_homo = self.alpha_mos.get(homo)?.energy;
        let e_lumo = self.alpha_mos.get(lumo)?.energy;
        // 1 Hartree = 27.2114 eV
        Some((e_lumo - e_homo) * 27.2114)
    }

    /// Add an alpha MO.
    pub fn add_alpha(&mut self, mo: MolecularOrbital) {
        self.alpha_mos.push(mo);
    }

    /// Add a beta MO.
    pub fn add_beta(&mut self, mo: MolecularOrbital) {
        self.beta_mos.push(mo);
    }
}

// ---------------------------------------------------------------------------
// BasisSet — Gaussian-type orbitals
// ---------------------------------------------------------------------------

/// Angular momentum shell type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShellType {
    /// s-type (l=0).
    S,
    /// p-type (l=1).
    P,
    /// d-type (l=2).
    D,
    /// f-type (l=3).
    F,
    /// g-type (l=4).
    G,
    /// SP-type (shared s+p exponents, Pople style).
    SP,
}

impl fmt::Display for ShellType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ShellType::S => write!(f, "S"),
            ShellType::P => write!(f, "P"),
            ShellType::D => write!(f, "D"),
            ShellType::F => write!(f, "F"),
            ShellType::G => write!(f, "G"),
            ShellType::SP => write!(f, "SP"),
        }
    }
}

/// A single contracted Gaussian shell.
#[derive(Debug, Clone)]
pub struct GaussianShell {
    /// Shell type.
    pub shell_type: ShellType,
    /// Exponents of the primitive Gaussians.
    pub exponents: Vec<f64>,
    /// Contraction coefficients for s (or s-part of SP).
    pub coefficients_s: Vec<f64>,
    /// Contraction coefficients for p-part of SP (empty if not SP).
    pub coefficients_p: Vec<f64>,
}

impl GaussianShell {
    /// Create a pure-type shell (S, P, D, …).
    pub fn new(shell_type: ShellType, exponents: Vec<f64>, coefficients: Vec<f64>) -> Self {
        Self {
            shell_type,
            exponents,
            coefficients_s: coefficients,
            coefficients_p: Vec::new(),
        }
    }

    /// Create an SP-type shell (Pople 6-31G style).
    pub fn new_sp(exponents: Vec<f64>, coeff_s: Vec<f64>, coeff_p: Vec<f64>) -> Self {
        Self {
            shell_type: ShellType::SP,
            exponents,
            coefficients_s: coeff_s,
            coefficients_p: coeff_p,
        }
    }

    /// Number of primitive Gaussians in this shell.
    pub fn n_primitives(&self) -> usize {
        self.exponents.len()
    }

    /// Number of contracted basis functions for this shell (angular momentum degeneracy).
    pub fn n_basis_functions(&self) -> usize {
        match self.shell_type {
            ShellType::S => 1,
            ShellType::P => 3,
            ShellType::D => 6, // Cartesian; use 5 for spherical
            ShellType::F => 10,
            ShellType::G => 15,
            ShellType::SP => 4, // 1 s + 3 p
        }
    }
}

/// Basis set definition for one element.
#[derive(Debug, Clone)]
pub struct ElementBasis {
    /// Element symbol.
    pub element: String,
    /// Basis shells.
    pub shells: Vec<GaussianShell>,
}

impl ElementBasis {
    /// Create an empty element basis.
    pub fn new(element: impl Into<String>) -> Self {
        Self {
            element: element.into(),
            shells: Vec::new(),
        }
    }

    /// Add a shell.
    pub fn add_shell(&mut self, shell: GaussianShell) {
        self.shells.push(shell);
    }

    /// Total number of contracted basis functions.
    pub fn n_basis_functions(&self) -> usize {
        self.shells.iter().map(|s| s.n_basis_functions()).sum()
    }
}

/// A named basis set (e.g., "STO-3G", "6-31G*", "cc-pVDZ").
#[derive(Debug, Clone)]
pub struct BasisSet {
    /// Name of the basis set.
    pub name: String,
    /// Per-element basis data.
    pub elements: HashMap<String, ElementBasis>,
}

impl BasisSet {
    /// Create an empty basis set.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            elements: HashMap::new(),
        }
    }

    /// Add an element basis.
    pub fn add_element(&mut self, basis: ElementBasis) {
        self.elements.insert(basis.element.clone(), basis);
    }

    /// Total number of basis functions for a molecule described by `atom_symbols`.
    pub fn total_basis_functions(&self, atom_symbols: &[&str]) -> usize {
        atom_symbols
            .iter()
            .filter_map(|&sym| self.elements.get(sym))
            .map(|eb| eb.n_basis_functions())
            .sum()
    }
}

// ---------------------------------------------------------------------------
// NormalModes
// ---------------------------------------------------------------------------

/// Normal mode data for a molecule.
#[derive(Debug, Clone)]
pub struct NormalModes {
    /// Number of atoms.
    pub n_atoms: usize,
    /// Vibrational frequencies (cm⁻¹), length = 3N - 6 (or 3N - 5 for linear).
    pub frequencies: Vec<f64>,
    /// IR intensities (km/mol).
    pub ir_intensities: Vec<f64>,
    /// Raman activities (Å⁴/amu), optional.
    pub raman_activities: Vec<f64>,
    /// Mass-weighted Cartesian displacement vectors, shape \[n_modes\]\[3*n_atoms\].
    pub displacements: Vec<Vec<f64>>,
    /// Reduced masses (amu) per mode.
    pub reduced_masses: Vec<f64>,
    /// Force constants (mdyn/Å) per mode.
    pub force_constants: Vec<f64>,
}

impl NormalModes {
    /// Create a new normal modes record.
    pub fn new(n_atoms: usize) -> Self {
        Self {
            n_atoms,
            frequencies: Vec::new(),
            ir_intensities: Vec::new(),
            raman_activities: Vec::new(),
            displacements: Vec::new(),
            reduced_masses: Vec::new(),
            force_constants: Vec::new(),
        }
    }

    /// Expected number of vibrational modes for a non-linear molecule.
    pub fn expected_modes_nonlinear(&self) -> usize {
        3 * self.n_atoms - 6
    }

    /// Expected number of vibrational modes for a linear molecule.
    pub fn expected_modes_linear(&self) -> usize {
        3 * self.n_atoms - 5
    }

    /// Number of imaginary (negative) frequencies.
    pub fn n_imaginary(&self) -> usize {
        self.frequencies.iter().filter(|&&f| f < 0.0).count()
    }

    /// Largest IR intensity and the corresponding frequency.
    pub fn strongest_ir_band(&self) -> Option<(f64, f64)> {
        self.ir_intensities
            .iter()
            .copied()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, intensity)| (self.frequencies[i], intensity))
    }

    /// Zero-point vibrational energy (cm⁻¹) = Σ(ν_i / 2) over positive frequencies.
    pub fn zpve_cm(&self) -> f64 {
        self.frequencies
            .iter()
            .filter(|&&f| f > 0.0)
            .map(|&f| f * 0.5)
            .sum()
    }
}

// ---------------------------------------------------------------------------
// ElectronicStructure
// ---------------------------------------------------------------------------

/// Electronic structure thermodynamic data.
#[derive(Debug, Clone)]
pub struct ElectronicStructure {
    /// Total electronic energy (Hartree).
    pub total_energy: f64,
    /// Zero-point energy correction (Hartree).
    pub zero_point_energy: f64,
    /// Thermal correction to internal energy at temperature T (Hartree).
    pub thermal_energy: f64,
    /// Thermal correction to enthalpy (Hartree).
    pub enthalpy_correction: f64,
    /// Entropy (cal/mol/K).
    pub entropy: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Gibbs free energy = H - T*S (Hartree).
    pub gibbs_free_energy: f64,
    /// Translational entropy contribution (cal/mol/K).
    pub translational_entropy: f64,
    /// Rotational entropy contribution (cal/mol/K).
    pub rotational_entropy: f64,
    /// Vibrational entropy contribution (cal/mol/K).
    pub vibrational_entropy: f64,
}

impl ElectronicStructure {
    /// Create an electronic structure record.
    pub fn new(total_energy: f64, temperature: f64) -> Self {
        Self {
            total_energy,
            zero_point_energy: 0.0,
            thermal_energy: 0.0,
            enthalpy_correction: 0.0,
            entropy: 0.0,
            temperature,
            gibbs_free_energy: 0.0,
            translational_entropy: 0.0,
            rotational_entropy: 0.0,
            vibrational_entropy: 0.0,
        }
    }

    /// E0 + ZPE (electronic energy plus zero-point correction).
    pub fn e0_plus_zpe(&self) -> f64 {
        self.total_energy + self.zero_point_energy
    }

    /// Enthalpy H = E + ZPE + thermal + kT (Hartree).
    ///
    /// kT at temperature T in Hartree: k_B * T / E_h
    pub fn enthalpy(&self) -> f64 {
        // kT/E_h ≈ T / 315775 (Boltzmann in Hartree/K)
        let kt = self.temperature / 315_775.13;
        self.total_energy + self.zero_point_energy + self.thermal_energy + kt
    }

    /// Free energy G = H - T*S, with S in cal/mol/K converted to Hartree.
    pub fn free_energy(&self) -> f64 {
        // 1 cal = 4.184 J; 1 Hartree = 4.3597e-18 J; 1 mol = 6.022e23
        // S_hartree = S_cal_mol_K * 4.184 / (4.3597e-18 * 6.022e23) * T
        let s_hartree_per_k = self.entropy * 4.184 / (4.3597e-18 * 6.022e23);
        self.enthalpy() - self.temperature * s_hartree_per_k
    }
}

// ---------------------------------------------------------------------------
// ExcitedStates
// ---------------------------------------------------------------------------

/// A single excited state transition.
#[derive(Debug, Clone)]
pub struct ExcitedStateTransition {
    /// Excited state index (1-based).
    pub state: usize,
    /// Excitation energy (eV).
    pub energy_ev: f64,
    /// Oscillator strength (dimensionless).
    pub oscillator_strength: f64,
    /// Dominant transition MOs (from MO i to MO a, with coefficient).
    pub dominant_transitions: Vec<(usize, usize, f64)>,
    /// Transition dipole moment vector (a.u.).
    pub transition_dipole: [f64; 3],
    /// Singlet or triplet.
    pub multiplicity: String,
}

impl ExcitedStateTransition {
    /// Create a new excited state transition.
    pub fn new(state: usize, energy_ev: f64, oscillator_strength: f64) -> Self {
        Self {
            state,
            energy_ev,
            oscillator_strength,
            dominant_transitions: Vec::new(),
            transition_dipole: [0.0; 3],
            multiplicity: "Singlet".to_string(),
        }
    }

    /// Absorption wavelength in nm.
    pub fn wavelength_nm(&self) -> f64 {
        1239.84 / self.energy_ev
    }

    /// Add a dominant transition.
    pub fn add_transition(&mut self, from: usize, to: usize, coeff: f64) {
        self.dominant_transitions.push((from, to, coeff));
    }
}

/// Collection of excited states and their absorption/emission spectra.
#[derive(Debug, Clone)]
pub struct ExcitedStates {
    /// All excited state transitions.
    pub states: Vec<ExcitedStateTransition>,
    /// Method used (e.g., "TDDFT", "EOM-CCSD", "ADC(2)").
    pub method: String,
    /// Ground state energy (Hartree).
    pub ground_energy: f64,
}

impl ExcitedStates {
    /// Create a new excited states collection.
    pub fn new(method: impl Into<String>, ground_energy: f64) -> Self {
        Self {
            states: Vec::new(),
            method: method.into(),
            ground_energy,
        }
    }

    /// Add an excited state.
    pub fn add_state(&mut self, state: ExcitedStateTransition) {
        self.states.push(state);
    }

    /// Return the state with the largest oscillator strength (brightest absorption).
    pub fn brightest_state(&self) -> Option<&ExcitedStateTransition> {
        self.states.iter().max_by(|a, b| {
            a.oscillator_strength
                .partial_cmp(&b.oscillator_strength)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Simulate a Lorentzian-broadened UV/vis spectrum.
    ///
    /// Returns `(wavelengths_nm, intensities)` sampled at `n_points` between
    /// `lambda_min` and `lambda_max` nm with broadening `fwhm_nm`.
    pub fn uv_vis_spectrum(
        &self,
        lambda_min: f64,
        lambda_max: f64,
        n_points: usize,
        fwhm_nm: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        let gamma = fwhm_nm / 2.0;
        let step = (lambda_max - lambda_min) / n_points.max(1) as f64;
        let wavelengths: Vec<f64> = (0..n_points)
            .map(|i| lambda_min + i as f64 * step)
            .collect();
        let intensities: Vec<f64> = wavelengths
            .iter()
            .map(|&lam| {
                self.states
                    .iter()
                    .map(|s| {
                        let dl = lam - s.wavelength_nm();
                        s.oscillator_strength * gamma * gamma / (dl * dl + gamma * gamma)
                    })
                    .sum()
            })
            .collect();
        (wavelengths, intensities)
    }
}

// ---------------------------------------------------------------------------
// NboOutput — Natural Bond Orbital analysis
// ---------------------------------------------------------------------------

/// A single natural bond orbital with occupancy and hybrid composition.
#[derive(Debug, Clone)]
pub struct NaturalBondOrbital {
    /// NBO index (1-based).
    pub index: usize,
    /// NBO type label (e.g., "BD", "LP", "CR", "BD*").
    pub nbo_type: String,
    /// Atom indices involved (1 for lone pair/core, 2 for bond).
    pub atoms: Vec<usize>,
    /// Occupancy (electrons).
    pub occupancy: f64,
    /// NBO energy (Hartree).
    pub energy: f64,
    /// Hybrid composition: `(atom_idx, hybrid_type, coefficient)`.
    pub hybrids: Vec<(usize, String, f64)>,
}

impl NaturalBondOrbital {
    /// Create a new NBO record.
    pub fn new(index: usize, nbo_type: impl Into<String>, occupancy: f64, energy: f64) -> Self {
        Self {
            index,
            nbo_type: nbo_type.into(),
            atoms: Vec::new(),
            occupancy,
            energy,
            hybrids: Vec::new(),
        }
    }

    /// Return `true` if this is a lone pair.
    pub fn is_lone_pair(&self) -> bool {
        self.nbo_type.starts_with("LP")
    }

    /// Return `true` if this is a bonding NBO.
    pub fn is_bonding(&self) -> bool {
        self.nbo_type == "BD"
    }

    /// Return `true` if this is an antibonding NBO.
    pub fn is_antibonding(&self) -> bool {
        self.nbo_type == "BD*"
    }
}

/// Wiberg bond order matrix (upper triangular, indexed by 0-based atom pairs).
#[derive(Debug, Clone)]
pub struct WibergBondOrderMatrix {
    /// Number of atoms.
    pub n_atoms: usize,
    /// Bond orders stored as a flat upper-triangular list (i < j).
    data: Vec<f64>,
}

impl WibergBondOrderMatrix {
    /// Create an n×n zero Wiberg bond order matrix.
    pub fn new(n_atoms: usize) -> Self {
        let size = n_atoms * (n_atoms + 1) / 2;
        Self {
            n_atoms,
            data: vec![0.0; size],
        }
    }

    /// Flat index for atom pair (i, j) with i <= j.
    fn idx(&self, i: usize, j: usize) -> usize {
        let (a, b) = if i <= j { (i, j) } else { (j, i) };
        a * self.n_atoms - a * (a + 1) / 2 + b
    }

    /// Set bond order for atom pair (i, j).
    pub fn set(&mut self, i: usize, j: usize, bo: f64) {
        let idx = self.idx(i, j);
        if idx < self.data.len() {
            self.data[idx] = bo;
        }
    }

    /// Get bond order for atom pair (i, j).
    pub fn get(&self, i: usize, j: usize) -> f64 {
        let idx = self.idx(i, j);
        if idx < self.data.len() {
            self.data[idx]
        } else {
            0.0
        }
    }

    /// Valence of atom `i` = sum of all bond orders involving atom `i`.
    pub fn valence(&self, i: usize) -> f64 {
        (0..self.n_atoms)
            .filter(|&j| j != i)
            .map(|j| self.get(i, j))
            .sum()
    }
}

/// Full NBO analysis output.
#[derive(Debug, Clone)]
pub struct NboOutput {
    /// List of all NBOs.
    pub nbos: Vec<NaturalBondOrbital>,
    /// Natural atomic charges per atom.
    pub natural_charges: Vec<f64>,
    /// Wiberg bond order matrix.
    pub bond_orders: WibergBondOrderMatrix,
    /// Second-order perturbation energies: `(donor_idx, acceptor_idx, E2_kcal_mol)`.
    pub second_order_interactions: Vec<(usize, usize, f64)>,
}

impl NboOutput {
    /// Create a new NBO output for `n_atoms`.
    pub fn new(n_atoms: usize) -> Self {
        Self {
            nbos: Vec::new(),
            natural_charges: vec![0.0; n_atoms],
            bond_orders: WibergBondOrderMatrix::new(n_atoms),
            second_order_interactions: Vec::new(),
        }
    }

    /// Total charge from natural atomic charges.
    pub fn total_charge(&self) -> f64 {
        self.natural_charges.iter().sum()
    }

    /// Count lone pairs.
    pub fn n_lone_pairs(&self) -> usize {
        self.nbos.iter().filter(|n| n.is_lone_pair()).count()
    }

    /// Count bonding NBOs.
    pub fn n_bonding(&self) -> usize {
        self.nbos.iter().filter(|n| n.is_bonding()).count()
    }

    /// Strongest donor–acceptor interaction (highest E2).
    pub fn strongest_interaction(&self) -> Option<&(usize, usize, f64)> {
        self.second_order_interactions
            .iter()
            .max_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- Atom ---
    #[test]
    fn test_atom_display() {
        let a = Atom::new("C", 0.0, 0.0, 0.0);
        let s = a.to_string();
        assert!(s.contains("C"), "display should contain element: {s}");
    }

    // --- GaussianInput ---
    #[test]
    fn test_gaussian_input_render_route() {
        let route = GaussianRouteSection::new("B3LYP", "6-31G*", "Opt");
        let inp = GaussianInput::new(route, "Test molecule", 0, 1);
        let rendered = inp.render();
        assert!(rendered.contains("B3LYP/6-31G*"), "route: {rendered}");
        assert!(rendered.contains("0 1"), "charge/mult: {rendered}");
    }

    #[test]
    fn test_gaussian_input_num_atoms() {
        let route = GaussianRouteSection::new("HF", "STO-3G", "SP");
        let mut inp = GaussianInput::new(route, "H2", 0, 1);
        inp.add_atom(Atom::new("H", 0.0, 0.0, 0.0));
        inp.add_atom(Atom::new("H", 0.74, 0.0, 0.0));
        assert_eq!(inp.num_atoms(), 2);
    }

    #[test]
    fn test_gaussian_input_link0() {
        let route = GaussianRouteSection::new("B3LYP", "6-31G*", "SP");
        let mut inp = GaussianInput::new(route, "test", 0, 1);
        inp.add_link0("chk=job.chk");
        let rendered = inp.render();
        assert!(rendered.contains("%chk=job.chk"), "link0: {rendered}");
    }

    #[test]
    fn test_gaussian_route_keywords() {
        let mut route = GaussianRouteSection::new("MP2", "cc-pVDZ", "SP");
        route.add_keyword("SCF", "Tight");
        let line = route.render();
        assert!(line.contains("SCF=Tight"), "keywords: {line}");
    }

    // --- GaussianOutput ---
    #[test]
    fn test_gaussian_output_parse_scf_energy() {
        let log = "SCF Done:  E(RB3LYP) =  -78.5873 A.U. after    8 cycles\n";
        let out = GaussianOutput::parse(log);
        let e = out.scf_energy.expect("SCF energy should be parsed");
        assert!((e + 78.5873).abs() < 1e-3, "e={e}");
    }

    #[test]
    fn test_gaussian_output_parse_frequencies() {
        let log = " Frequencies --   1625.54  3050.12  3100.45\n IR Inten    --      10.20      0.50      5.30\n";
        let out = GaussianOutput::parse(log);
        assert_eq!(out.frequencies.len(), 3);
        assert!((out.frequencies[0].freq_cm - 1625.54).abs() < 0.01);
    }

    #[test]
    fn test_gaussian_output_imaginary_frequency() {
        let log = " Frequencies --   -150.3   500.0\n";
        let out = GaussianOutput::parse(log);
        assert_eq!(out.n_imaginary, 1);
    }

    #[test]
    fn test_gaussian_output_zpe() {
        let log = "Zero-point correction=  0.054321 (Hartree/Particle)\n";
        let out = GaussianOutput::parse(log);
        let zpe = out.zero_point_energy.expect("ZPE");
        assert!((zpe - 0.054321).abs() < 1e-6, "zpe={zpe}");
    }

    #[test]
    fn test_gaussian_output_opt_converged() {
        let log = "Optimization completed.\n";
        let out = GaussianOutput::parse(log);
        assert!(out.opt_converged);
    }

    #[test]
    fn test_gaussian_output_lowest_frequency() {
        let log = " Frequencies --   300.0   500.0   100.0\n";
        let out = GaussianOutput::parse(log);
        let lo = out.lowest_frequency().expect("lowest freq");
        assert!((lo - 100.0).abs() < 0.01);
    }

    // --- OrcaInput ---
    #[test]
    fn test_orca_input_render_keywords() {
        let mut inp = OrcaInput::new(0, 1);
        inp.add_keyword("B3LYP");
        inp.add_keyword("def2-TZVP");
        let rendered = inp.render();
        assert!(rendered.contains("! B3LYP def2-TZVP"), "render: {rendered}");
    }

    #[test]
    fn test_orca_input_coord_block() {
        let mut inp = OrcaInput::new(0, 1);
        inp.add_atom(Atom::new("O", 0.0, 0.0, 0.0));
        inp.add_atom(Atom::new("H", 0.96, 0.0, 0.0));
        let rendered = inp.render();
        assert!(rendered.contains("* xyz 0 1"), "coord block: {rendered}");
        assert!(rendered.contains("O"), "oxygen: {rendered}");
    }

    // --- OrcaOutput ---
    #[test]
    fn test_orca_output_parse_energy() {
        let log = "FINAL SINGLE POINT ENERGY       -76.3456789\n";
        let out = OrcaOutput::parse(log);
        let e = out.electronic_energy.expect("energy");
        assert!((e + 76.3456789).abs() < 1e-6, "e={e}");
    }

    #[test]
    fn test_orca_output_total_energy_with_dispersion() {
        let log =
            "FINAL SINGLE POINT ENERGY       -76.3456789\nDispersion correction    -0.0012345\n";
        let out = OrcaOutput::parse(log);
        let total = out.total_energy().expect("total");
        assert!((total + 76.3469134).abs() < 1e-4, "total={total}");
    }

    #[test]
    fn test_orca_output_opt_converged() {
        let log = "OPTIMIZATION HAS CONVERGED\n";
        let out = OrcaOutput::parse(log);
        assert!(out.opt_converged);
    }

    // --- MolecularOrbitalSet ---
    #[test]
    fn test_mo_homo_lumo_gap() {
        let mut mo_set = MolecularOrbitalSet::new(5, 5);
        // HOMO at index 4, energy = -0.3 Eh; LUMO at index 5, energy = 0.1 Eh
        for i in 0..5 {
            let e = -0.5 + i as f64 * 0.05;
            mo_set.add_alpha(MolecularOrbital::new(i, e, 2.0, Spin::Alpha));
        }
        mo_set.add_alpha(MolecularOrbital::new(5, 0.1, 0.0, Spin::Alpha));
        let gap = mo_set.homo_lumo_gap_ev().expect("gap");
        // HOMO energy = -0.5 + 4*0.05 = -0.3; LUMO = 0.1; gap = 0.4 Eh * 27.2114
        let expected = 0.4 * 27.2114;
        assert!((gap - expected).abs() < 0.01, "gap={gap}");
    }

    #[test]
    fn test_mo_homo_idx() {
        let mo_set = MolecularOrbitalSet::new(3, 3);
        assert_eq!(mo_set.homo_idx(), Some(2));
    }

    #[test]
    fn test_mo_coeff_norm() {
        let mut mo = MolecularOrbital::new(0, -0.5, 2.0, Spin::Alpha);
        mo.coefficients = vec![1.0, 0.0, 0.0, 0.0];
        assert!((mo.coeff_norm() - 1.0).abs() < 1e-12);
    }

    // --- BasisSet ---
    #[test]
    fn test_basis_set_total_functions() {
        let mut bs = BasisSet::new("STO-3G");
        // Carbon: 1s + 2s + 2p(×3) = 5 basis functions
        let mut c_basis = ElementBasis::new("C");
        c_basis.add_shell(GaussianShell::new(
            ShellType::S,
            vec![71.62, 13.04, 3.531],
            vec![0.1543, 0.5353, 0.4446],
        ));
        c_basis.add_shell(GaussianShell::new(
            ShellType::SP,
            vec![2.941, 0.6835, 0.2222],
            vec![-0.09996, 0.3995, 0.7001],
        ));
        bs.add_element(c_basis);
        let total = bs.total_basis_functions(&["C", "C"]);
        // Each C: 1 + 4 = 5 basis functions → 10 total
        assert_eq!(total, 10, "total={total}");
    }

    #[test]
    fn test_shell_n_primitives() {
        let shell = GaussianShell::new(ShellType::D, vec![1.0, 2.0, 3.0], vec![0.1, 0.5, 0.4]);
        assert_eq!(shell.n_primitives(), 3);
    }

    // --- NormalModes ---
    #[test]
    fn test_normal_modes_n_imaginary() {
        let mut nm = NormalModes::new(3);
        nm.frequencies = vec![
            -150.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0,
        ];
        assert_eq!(nm.n_imaginary(), 1);
    }

    #[test]
    fn test_normal_modes_zpve() {
        let mut nm = NormalModes::new(2);
        nm.frequencies = vec![1000.0, 2000.0];
        let zpve = nm.zpve_cm();
        assert!((zpve - 1500.0).abs() < 1e-6, "zpve={zpve}");
    }

    #[test]
    fn test_normal_modes_strongest_ir() {
        let mut nm = NormalModes::new(3);
        nm.frequencies = vec![500.0, 1000.0, 1500.0];
        nm.ir_intensities = vec![10.0, 250.0, 5.0];
        let (freq, intensity) = nm.strongest_ir_band().expect("IR band");
        assert!((freq - 1000.0).abs() < 1e-6);
        assert!((intensity - 250.0).abs() < 1e-6);
    }

    // --- ElectronicStructure ---
    #[test]
    fn test_electronic_structure_e0_zpe() {
        let mut es = ElectronicStructure::new(-76.4, 298.15);
        es.zero_point_energy = 0.02;
        assert!((es.e0_plus_zpe() + 76.38).abs() < 1e-6);
    }

    #[test]
    fn test_electronic_structure_enthalpy_positive() {
        let es = ElectronicStructure::new(-76.4, 298.15);
        // Enthalpy should be close to total_energy (small correction)
        let h = es.enthalpy();
        assert!(
            (h - es.total_energy).abs() < 0.01,
            "h-e={}",
            h - es.total_energy
        );
    }

    // --- ExcitedStates ---
    #[test]
    fn test_excited_states_wavelength() {
        let s = ExcitedStateTransition::new(1, 4.0, 0.3);
        // 1239.84 / 4.0 = 309.96 nm
        let wl = s.wavelength_nm();
        assert!((wl - 309.96).abs() < 0.01, "wl={wl}");
    }

    #[test]
    fn test_excited_states_brightest() {
        let mut es = ExcitedStates::new("TDDFT", -76.4);
        es.add_state(ExcitedStateTransition::new(1, 4.0, 0.05));
        es.add_state(ExcitedStateTransition::new(2, 5.0, 0.80));
        es.add_state(ExcitedStateTransition::new(3, 6.0, 0.10));
        let brightest = es.brightest_state().expect("brightest");
        assert_eq!(brightest.state, 2);
    }

    #[test]
    fn test_excited_states_uv_vis_spectrum() {
        let mut es = ExcitedStates::new("TDDFT", -76.4);
        es.add_state(ExcitedStateTransition::new(1, 4.0, 1.0));
        let (wls, ints) = es.uv_vis_spectrum(200.0, 600.0, 100, 10.0);
        assert_eq!(wls.len(), 100);
        assert_eq!(ints.len(), 100);
        // Intensity should peak near 310 nm
        let max_idx = ints
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        assert!(
            wls[max_idx] > 280.0 && wls[max_idx] < 340.0,
            "peak at {}",
            wls[max_idx]
        );
    }

    // --- NboOutput ---
    #[test]
    fn test_nbo_total_charge() {
        let mut nbo = NboOutput::new(3);
        nbo.natural_charges = vec![-0.5, 0.25, 0.25];
        assert!((nbo.total_charge()).abs() < 1e-12);
    }

    #[test]
    fn test_nbo_count_lone_pairs() {
        let mut nbo = NboOutput::new(2);
        nbo.nbos.push(NaturalBondOrbital::new(1, "LP", 1.98, -0.5));
        nbo.nbos.push(NaturalBondOrbital::new(2, "BD", 1.95, -0.4));
        nbo.nbos.push(NaturalBondOrbital::new(3, "LP", 1.97, -0.45));
        assert_eq!(nbo.n_lone_pairs(), 2);
        assert_eq!(nbo.n_bonding(), 1);
    }

    #[test]
    fn test_wiberg_bond_order() {
        let mut mat = WibergBondOrderMatrix::new(3);
        mat.set(0, 1, 1.95);
        mat.set(1, 2, 0.98);
        assert!((mat.get(0, 1) - 1.95).abs() < 1e-12);
        assert!((mat.get(1, 0) - 1.95).abs() < 1e-12); // symmetric
        let val = mat.valence(1);
        assert!((val - 2.93).abs() < 1e-6, "valence={val}");
    }

    #[test]
    fn test_nbo_strongest_interaction() {
        let mut nbo = NboOutput::new(4);
        nbo.second_order_interactions = vec![(0, 1, 5.0), (1, 2, 30.0), (2, 3, 15.0)];
        let strongest = nbo.strongest_interaction().expect("strongest");
        assert!((strongest.2 - 30.0).abs() < 1e-12);
    }
}
