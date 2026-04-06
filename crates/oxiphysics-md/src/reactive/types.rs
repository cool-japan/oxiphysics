//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Manages the current bonded topology in a reactive MD simulation.
///
/// Detects bond breaking and forming by comparing current distances against
/// cutoffs, and records breaking events.
#[allow(dead_code)]
pub struct TopologyManager {
    /// Current list of bonded pairs `(i, j)` with `i < j`.
    pub current_bonds: Vec<(usize, usize)>,
    /// History of bond-breaking events.
    pub break_events: Vec<BondBreakingEvent>,
    /// Bond-order parameters per pair.
    pub bond_params: Vec<ReaxFFBond>,
    /// Distance cutoff beyond which a bond is considered broken.
    pub break_cutoff: f64,
    /// Distance cutoff within which a new bond can form.
    pub form_cutoff: f64,
}
impl TopologyManager {
    /// Create a new topology manager.
    #[allow(dead_code)]
    pub fn new(break_cutoff: f64, form_cutoff: f64) -> Self {
        Self {
            current_bonds: Vec::new(),
            break_events: Vec::new(),
            bond_params: Vec::new(),
            break_cutoff,
            form_cutoff,
        }
    }
    /// Rebuild bond topology from current positions.
    ///
    /// Bonds whose distance exceeds `break_cutoff` are removed and logged.
    /// New bonds within `form_cutoff` are added.
    #[allow(dead_code)]
    pub fn update(&mut self, positions: &[[f64; 3]], time: f64) {
        let n = positions.len();
        let mut broken_indices = Vec::new();
        for (idx, &(i, j)) in self.current_bonds.iter().enumerate() {
            let r = dist(positions[i], positions[j]);
            if r > self.break_cutoff {
                broken_indices.push(idx);
                self.break_events.push(BondBreakingEvent {
                    atom_i: i,
                    atom_j: j,
                    time,
                    reason: BreakReason::Stretched,
                });
            }
        }
        for &idx in broken_indices.iter().rev() {
            self.current_bonds.swap_remove(idx);
            if !self.bond_params.is_empty() && idx < self.bond_params.len() {
                self.bond_params.swap_remove(idx);
            }
        }
        for i in 0..n {
            for j in (i + 1)..n {
                if self.current_bonds.contains(&(i, j)) {
                    continue;
                }
                let r = dist(positions[i], positions[j]);
                if r < self.form_cutoff {
                    self.current_bonds.push((i, j));
                }
            }
        }
    }
    /// Update topology using bond-order criterion instead of distance.
    ///
    /// Bonds with BO < `bo_threshold` are broken.
    #[allow(dead_code)]
    pub fn update_by_bond_order(
        &mut self,
        positions: &[[f64; 3]],
        bond: &ReaxFFBond,
        bo_threshold: f64,
        time: f64,
    ) {
        let mut broken = Vec::new();
        for (idx, &(i, j)) in self.current_bonds.iter().enumerate() {
            let r = dist(positions[i], positions[j]);
            let bo = compute_bond_order(r, bond);
            if bo < bo_threshold {
                broken.push(idx);
                self.break_events.push(BondBreakingEvent {
                    atom_i: i,
                    atom_j: j,
                    time,
                    reason: BreakReason::BondOrderLow,
                });
            }
        }
        for &idx in broken.iter().rev() {
            self.current_bonds.swap_remove(idx);
            if !self.bond_params.is_empty() && idx < self.bond_params.len() {
                self.bond_params.swap_remove(idx);
            }
        }
    }
    /// Number of currently active bonds.
    #[allow(dead_code)]
    pub fn num_bonds(&self) -> usize {
        self.current_bonds.len()
    }
    /// Number of breaking events recorded.
    #[allow(dead_code)]
    pub fn num_break_events(&self) -> usize {
        self.break_events.len()
    }
}
/// A log entry describing a single reaction event.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub struct ReactionEvent {
    /// Simulation step at which the event was detected.
    pub step: usize,
    /// Atom index of the donor (bond-breaking: losing partner; forming: gaining partner).
    pub atom_i: usize,
    /// Atom index of the partner.
    pub atom_j: usize,
    /// Type of event.
    pub event_type: ReactionEventType,
}
/// Collects running statistics during a reactive MD simulation.
///
/// Records step-by-step energies, temperatures, and bond-event counts.
#[allow(dead_code)]
#[derive(Debug, Default, Clone)]
pub struct ReactiveMDStats {
    /// Total potential energies at each recorded step.
    pub energies: Vec<f64>,
    /// Instantaneous temperature estimates at each recorded step.
    pub temperatures: Vec<f64>,
    /// Number of bond-breaking events up to each recorded step.
    pub break_counts: Vec<usize>,
    /// Number of bond-forming events up to each recorded step.
    pub form_counts: Vec<usize>,
}
impl ReactiveMDStats {
    /// Create a new empty statistics collector.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self::default()
    }
    /// Record a new step.
    #[allow(dead_code)]
    pub fn record(&mut self, energy: f64, temperature: f64, n_breaks: usize, n_forms: usize) {
        self.energies.push(energy);
        self.temperatures.push(temperature);
        self.break_counts.push(n_breaks);
        self.form_counts.push(n_forms);
    }
    /// Mean potential energy over all recorded steps.
    #[allow(dead_code)]
    pub fn mean_energy(&self) -> f64 {
        if self.energies.is_empty() {
            return 0.0;
        }
        self.energies.iter().sum::<f64>() / self.energies.len() as f64
    }
    /// Mean temperature over all recorded steps.
    #[allow(dead_code)]
    pub fn mean_temperature(&self) -> f64 {
        if self.temperatures.is_empty() {
            return 0.0;
        }
        self.temperatures.iter().sum::<f64>() / self.temperatures.len() as f64
    }
    /// Total number of bond-breaking events at the last recorded step.
    #[allow(dead_code)]
    pub fn total_breaks(&self) -> usize {
        self.break_counts.last().copied().unwrap_or(0)
    }
    /// Number of recorded steps.
    #[allow(dead_code)]
    pub fn n_steps(&self) -> usize {
        self.energies.len()
    }
}
/// Over-coordination parameters for the ReaxFF penalty term.
///
/// When the sum of bond orders on an atom exceeds its valence, ReaxFF applies
/// a penalty energy to restore the correct coordination.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct OverCoordParams {
    /// Expected valence of the atom (e.g. 4 for C, 2 for O).
    pub val_i: f64,
    /// Over-coordination penalty exponent p_ovun1.
    pub p_ovun1: f64,
    /// Over-coordination penalty coefficient p_ovun2 (kJ/mol per bond-order unit).
    pub p_ovun2: f64,
}
impl OverCoordParams {
    /// Standard carbon over-coordination parameters.
    #[allow(dead_code)]
    pub fn carbon() -> Self {
        Self {
            val_i: 4.0,
            p_ovun1: 0.1736,
            p_ovun2: 0.9951,
        }
    }
    /// Standard oxygen over-coordination parameters.
    #[allow(dead_code)]
    pub fn oxygen() -> Self {
        Self {
            val_i: 2.0,
            p_ovun1: 0.0918,
            p_ovun2: 1.3321,
        }
    }
    /// Standard nitrogen over-coordination parameters.
    #[allow(dead_code)]
    pub fn nitrogen() -> Self {
        Self {
            val_i: 3.0,
            p_ovun1: 0.1450,
            p_ovun2: 1.0000,
        }
    }
}
/// Tabulated bond energy as a function of bond order.
///
/// Allows look-up interpolation: given BO, returns the corresponding bond energy
/// from a user-supplied table of `(bond_order, energy)` pairs sorted by BO.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BondEnergyTable {
    /// Sorted (bond_order, energy) pairs.
    pub entries: Vec<(f64, f64)>,
}
#[allow(dead_code)]
impl BondEnergyTable {
    /// Create a new lookup table from unsorted entries.
    pub fn new(mut entries: Vec<(f64, f64)>) -> Self {
        entries.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        Self { entries }
    }
    /// Linear interpolation: return energy for the given bond order.
    ///
    /// Clamps to the table range outside the domain.
    pub fn lookup(&self, bo: f64) -> f64 {
        let n = self.entries.len();
        if n == 0 {
            return 0.0;
        }
        if bo <= self.entries[0].0 {
            return self.entries[0].1;
        }
        if bo >= self.entries[n - 1].0 {
            return self.entries[n - 1].1;
        }
        let pos = self.entries.partition_point(|&(b, _)| b <= bo);
        let lo = &self.entries[pos - 1];
        let hi = &self.entries[pos];
        let t = (bo - lo.0) / (hi.0 - lo.0);
        lo.1 + t * (hi.1 - lo.1)
    }
    /// Number of table entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// True if the table is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}
/// ReaxFF-style bond order parameters for a pair of atom types.
///
/// The bond order is computed from the ReaxFF sigma-bond formula:
/// ```text
/// BO' = exp(p_bo1 * (r / r0)^p_bo2) + exp(p_bo3 * (r / r0)^p_bo4)
/// ```
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ReaxFFBond {
    /// Equilibrium sigma-bond length.
    pub r0: f64,
    /// Dissociation energy.
    pub d_e: f64,
    /// Morse width parameter.
    pub beta: f64,
    /// Sigma-bond order coefficient.
    pub p_bo1: f64,
    /// Sigma-bond order exponent.
    pub p_bo2: f64,
    /// Pi-bond order coefficient.
    pub p_bo3: f64,
    /// Pi-bond order exponent.
    pub p_bo4: f64,
}
impl ReaxFFBond {
    /// Construct typical C-C single bond parameters.
    #[allow(dead_code)]
    pub fn carbon_carbon() -> Self {
        Self {
            r0: 1.54,
            d_e: 6.0,
            beta: 1.8,
            p_bo1: -0.0777,
            p_bo2: 6.0,
            p_bo3: -0.1000,
            p_bo4: 15.0,
        }
    }
    /// Construct typical C-H bond parameters.
    #[allow(dead_code)]
    pub fn carbon_hydrogen() -> Self {
        Self {
            r0: 1.09,
            d_e: 4.5,
            beta: 2.0,
            p_bo1: -0.05,
            p_bo2: 6.0,
            p_bo3: 0.0,
            p_bo4: 1.0,
        }
    }
    /// Construct typical O-H bond parameters.
    #[allow(dead_code)]
    pub fn oxygen_hydrogen() -> Self {
        Self {
            r0: 0.96,
            d_e: 5.0,
            beta: 2.2,
            p_bo1: -0.06,
            p_bo2: 6.0,
            p_bo3: 0.0,
            p_bo4: 1.0,
        }
    }
}
/// Oxidation state of an atom based on its formal charge and bonding environment.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct OxidationState {
    /// Atom index.
    pub atom_idx: usize,
    /// Current oxidation state (integer-like value).
    pub state: i32,
    /// Electronegativity (Pauling scale).
    pub electronegativity: f64,
}
impl OxidationState {
    /// Create a new oxidation state.
    #[allow(dead_code)]
    pub fn new(atom_idx: usize, state: i32, electronegativity: f64) -> Self {
        Self {
            atom_idx,
            state,
            electronegativity,
        }
    }
}
/// A reactive simulation container that tracks bond topology over time.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ReactiveSimulation {
    /// Atom positions.
    pub positions: Vec<[f64; 3]>,
    /// Bond-order parameters for all pairs.
    pub bond_params: ReaxFFBond,
    /// Bond-order threshold below which a bond is considered broken.
    pub bo_threshold: f64,
    /// Current list of bonded pairs.
    pub current_bonds: Vec<(usize, usize)>,
    /// Accumulated reaction event log.
    pub event_log: Vec<ReactionEvent>,
}
impl ReactiveSimulation {
    /// Create a new reactive simulation.
    #[allow(dead_code)]
    pub fn new(positions: Vec<[f64; 3]>, bond_params: ReaxFFBond, bo_threshold: f64) -> Self {
        let n = positions.len();
        let current_bonds = bond_list_from_bond_orders(&positions, &bond_params, bo_threshold);
        let _ = n;
        Self {
            positions,
            bond_params,
            bo_threshold,
            current_bonds,
            event_log: Vec::new(),
        }
    }
    /// Update atom positions and detect reaction events.
    ///
    /// Compares the new bond topology derived from `new_positions` with the
    /// stored `current_bonds`.  Broken and formed bonds are appended to
    /// `event_log` and `current_bonds` is updated.
    ///
    /// # Arguments
    /// * `new_positions` — updated atom positions
    /// * `step`          — current simulation step (for logging)
    ///
    /// Returns the number of new reaction events detected.
    #[allow(dead_code)]
    pub fn detect_reaction_events(&mut self, new_positions: &[[f64; 3]], step: usize) -> usize {
        let new_bonds =
            bond_list_from_bond_orders(new_positions, &self.bond_params, self.bo_threshold);
        let (broken, formed) = detect_topology_changes(&self.current_bonds, &new_bonds);
        let mut n_events = 0usize;
        for (i, j) in &broken {
            self.event_log.push(ReactionEvent {
                step,
                atom_i: *i,
                atom_j: *j,
                event_type: ReactionEventType::BondBroken,
            });
            n_events += 1;
        }
        for (i, j) in &formed {
            self.event_log.push(ReactionEvent {
                step,
                atom_i: *i,
                atom_j: *j,
                event_type: ReactionEventType::BondFormed,
            });
            n_events += 1;
        }
        self.current_bonds = new_bonds;
        self.positions = new_positions.to_vec();
        n_events
    }
    /// Total number of bond-breaking events recorded.
    #[allow(dead_code)]
    pub fn count_breaks(&self) -> usize {
        self.event_log
            .iter()
            .filter(|e| e.event_type == ReactionEventType::BondBroken)
            .count()
    }
    /// Total number of bond-forming events recorded.
    #[allow(dead_code)]
    pub fn count_forms(&self) -> usize {
        self.event_log
            .iter()
            .filter(|e| e.event_type == ReactionEventType::BondFormed)
            .count()
    }
}
/// A bond-order potential that includes σ and π bond contributions.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BondOrderPotential {
    /// σ-bond Morse-like parameters.
    pub sigma: ReaxFFBond,
    /// π-bond parameters.
    pub pi: PiBondParams,
}
impl BondOrderPotential {
    /// Create a C=C double bond potential.
    #[allow(dead_code)]
    pub fn carbon_double() -> Self {
        Self {
            sigma: ReaxFFBond::carbon_carbon(),
            pi: PiBondParams::carbon_carbon_double(),
        }
    }
    /// Compute the π bond contribution to the bond energy.
    ///
    /// ```text
    /// E_π = D_e^π * BO_pi * exp(p_be1 * (1 − BO_pi^p_be2))
    /// ```
    ///
    /// # Arguments
    /// * `bo_pi` — the π component of the bond order (0 ≤ bo_pi ≤ 1)
    ///
    /// Returns energy in kJ/mol.
    #[allow(dead_code)]
    pub fn compute_pi_bond_contribution(&self, bo_pi: f64) -> f64 {
        if bo_pi <= 0.0 {
            return 0.0;
        }
        let p = &self.pi;
        p.d_e_pi * bo_pi * (p.p_be1 * (1.0 - bo_pi.powf(p.p_be2))).exp()
    }
    /// Compute the total bond energy: σ-bond Morse + π-bond contribution.
    ///
    /// # Arguments
    /// * `r`     — current bond distance (Å)
    /// * `bo_pi` — π component of the bond order
    #[allow(dead_code)]
    pub fn compute_total_bond_energy(&self, r: f64, bo_pi: f64) -> f64 {
        let e_sigma = compute_reaxff_bond_energy(r, &self.sigma);
        let e_pi = self.compute_pi_bond_contribution(bo_pi);
        e_sigma + e_pi
    }
}
/// Atom-type parameters used in a reactive force field.
///
/// Stores standard element-level properties used for multiple reactive
/// force field flavours (ReaxFF, REBO, etc.).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ReactiveAtomParams {
    /// Element symbol.
    pub element: &'static str,
    /// Atomic mass (g/mol).
    pub mass: f64,
    /// Pauling electronegativity.
    pub electronegativity: f64,
    /// Chemical hardness (eV).
    pub hardness: f64,
    /// Atomic radius (Å).
    pub radius: f64,
    /// Valence (expected coordination number).
    pub valence: f64,
}
impl ReactiveAtomParams {
    /// Standard C parameters.
    #[allow(dead_code)]
    pub fn carbon() -> Self {
        Self {
            element: "C",
            mass: 12.011,
            electronegativity: 2.55,
            hardness: 5.0,
            radius: 0.77,
            valence: 4.0,
        }
    }
    /// Standard H parameters.
    #[allow(dead_code)]
    pub fn hydrogen() -> Self {
        Self {
            element: "H",
            mass: 1.008,
            electronegativity: 2.20,
            hardness: 6.9,
            radius: 0.31,
            valence: 1.0,
        }
    }
    /// Standard O parameters.
    #[allow(dead_code)]
    pub fn oxygen() -> Self {
        Self {
            element: "O",
            mass: 15.999,
            electronegativity: 3.44,
            hardness: 8.0,
            radius: 0.73,
            valence: 2.0,
        }
    }
    /// Standard N parameters.
    #[allow(dead_code)]
    pub fn nitrogen() -> Self {
        Self {
            element: "N",
            mass: 14.007,
            electronegativity: 3.04,
            hardness: 7.3,
            radius: 0.75,
            valence: 3.0,
        }
    }
}
/// π-bond (pi) contribution parameters for a bond-order potential.
///
/// The pi bond energy modifies the total bond energy based on the pi bond order:
///
/// ```text
/// E_π = D_e^π * BO_pi * exp(p_be1 * (1 − BO_pi^p_be2))
/// ```
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PiBondParams {
    /// π-bond dissociation energy (kJ/mol).
    pub d_e_pi: f64,
    /// π-bond order scaling exponent p_be1.
    pub p_be1: f64,
    /// π-bond order exponent p_be2.
    pub p_be2: f64,
}
impl PiBondParams {
    /// Typical C=C π bond parameters.
    #[allow(dead_code)]
    pub fn carbon_carbon_double() -> Self {
        Self {
            d_e_pi: 134.0,
            p_be1: -0.7738,
            p_be2: 1.0,
        }
    }
    /// Typical C≡C π bond parameters (cumulative).
    #[allow(dead_code)]
    pub fn carbon_carbon_triple() -> Self {
        Self {
            d_e_pi: 180.0,
            p_be1: -0.90,
            p_be2: 1.5,
        }
    }
}
/// Type of a reaction event.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub enum ReactionEventType {
    /// A bond was broken (BO dropped below threshold).
    BondBroken,
    /// A bond was formed (BO rose above threshold).
    BondFormed,
}
/// Records a bond-breaking event during a reactive MD trajectory.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BondBreakingEvent {
    /// Index of the first atom in the bond.
    pub atom_i: usize,
    /// Index of the second atom in the bond.
    pub atom_j: usize,
    /// Simulation time at which the break was detected.
    pub time: f64,
    /// Reason for the break.
    pub reason: BreakReason,
}
impl BondBreakingEvent {
    /// Create a new bond-breaking event.
    #[allow(dead_code)]
    pub fn new(atom_i: usize, atom_j: usize, time: f64, reason: BreakReason) -> Self {
        Self {
            atom_i,
            atom_j,
            time,
            reason,
        }
    }
}
/// Morse-like bond-order parameters.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BondOrder {
    /// Equilibrium bond length.
    pub r0: f64,
    /// Dissociation energy (well depth).
    pub d_e: f64,
    /// Range parameter (controls well width).
    pub beta: f64,
    /// Smoothing / scaling parameter (reserved for bond-order weighting).
    pub s: f64,
}
/// Reason a bond was considered broken during reactive MD.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub enum BreakReason {
    /// Bond length exceeded a threshold distance.
    Stretched,
    /// Total energy exceeded the dissociation energy.
    EnergyExceeded,
    /// Force along the bond exceeded a threshold.
    ForceExceeded,
    /// Bond order dropped below threshold.
    BondOrderLow,
}
/// A reactive site: a central atom with its bonded neighbors and bond orders.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ReactiveSite {
    /// Index of the central atom.
    pub central_atom: usize,
    /// Indices of the neighbor atoms.
    pub neighbor_atoms: Vec<usize>,
    /// Bond orders to each neighbor (same length as `neighbor_atoms`).
    pub bond_orders: Vec<f64>,
}
impl ReactiveSite {
    /// Create a new reactive site.
    #[allow(dead_code)]
    pub fn new(central_atom: usize) -> Self {
        Self {
            central_atom,
            neighbor_atoms: Vec::new(),
            bond_orders: Vec::new(),
        }
    }
    /// Add a neighbor with the given bond order.
    #[allow(dead_code)]
    pub fn add_neighbor(&mut self, neighbor: usize, bond_order: f64) {
        self.neighbor_atoms.push(neighbor);
        self.bond_orders.push(bond_order);
    }
    /// Total bond order (sum over all neighbors).
    #[allow(dead_code)]
    pub fn total_bond_order(&self) -> f64 {
        self.bond_orders.iter().sum()
    }
    /// Number of neighbors.
    #[allow(dead_code)]
    pub fn coordination_number(&self) -> usize {
        self.neighbor_atoms.len()
    }
    /// Maximum bond order among all neighbors.
    #[allow(dead_code)]
    pub fn max_bond_order(&self) -> f64 {
        self.bond_orders.iter().cloned().fold(0.0_f64, f64::max)
    }
    /// Remove neighbors with bond order below threshold.
    #[allow(dead_code)]
    pub fn prune_weak_bonds(&mut self, threshold: f64) {
        let mut i = 0;
        while i < self.bond_orders.len() {
            if self.bond_orders[i] < threshold {
                self.bond_orders.swap_remove(i);
                self.neighbor_atoms.swap_remove(i);
            } else {
                i += 1;
            }
        }
    }
}
/// A minimal reactive force-field system.
///
/// Holds atom positions, simulation box, Morse bond parameters, and the list
/// of bonded pairs.  Forces are accumulated per atom using Newton's third law.
#[allow(dead_code)]
pub struct ReaxFFSimple {
    /// Atom positions in Cartesian coordinates.
    pub atoms: Vec<[f64; 3]>,
    /// Simulation box lengths along x, y, z (PBC not applied here).
    pub box_lengths: [f64; 3],
    /// Bond-order parameter sets (one entry per bonded pair in `pairs`).
    pub bond_params: Vec<BondOrder>,
    /// Bonded pairs as `(atom_i, atom_j)` indices into `atoms`.
    pub pairs: Vec<(usize, usize)>,
}
impl ReaxFFSimple {
    /// Total Morse bond energy summed over all listed pairs.
    #[allow(dead_code)]
    pub fn compute_bond_energy(&self) -> f64 {
        self.pairs
            .iter()
            .zip(self.bond_params.iter())
            .map(|(&(i, j), params)| {
                let r = dist(self.atoms[i], self.atoms[j]);
                morse_potential(r, params)
            })
            .sum()
    }
    /// Per-atom force vectors from Morse bond interactions.
    ///
    /// Returns a `Vec<[f64;3]>` of length `atoms.len()`.
    #[allow(dead_code)]
    pub fn compute_bond_forces(&self) -> Vec<[f64; 3]> {
        let n = self.atoms.len();
        let mut forces = vec![[0.0f64; 3]; n];
        for (&(i, j), params) in self.pairs.iter().zip(self.bond_params.iter()) {
            let r_vec = [
                self.atoms[j][0] - self.atoms[i][0],
                self.atoms[j][1] - self.atoms[i][1],
                self.atoms[j][2] - self.atoms[i][2],
            ];
            let r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();
            if r < 1e-14 {
                continue;
            }
            let dv_dr = morse_force(r, params);
            let r_hat = [r_vec[0] / r, r_vec[1] / r, r_vec[2] / r];
            for k in 0..3 {
                forces[i][k] -= dv_dr * r_hat[k];
                forces[j][k] += dv_dr * r_hat[k];
            }
        }
        forces
    }
    /// Compute per-bond bond orders using ReaxFF formula.
    #[allow(dead_code)]
    pub fn compute_bond_orders(&self, reaxff_bond: &ReaxFFBond) -> Vec<f64> {
        self.pairs
            .iter()
            .map(|&(i, j)| {
                let r = dist(self.atoms[i], self.atoms[j]);
                compute_bond_order(r, reaxff_bond)
            })
            .collect()
    }
}
