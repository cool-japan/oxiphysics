//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Isotropic Parrinello-Rahman barostat.
///
/// Extended-Lagrangian pressure control: the cubic-equivalent box edge
/// `L = V^(1/3)` is a dynamical variable with inertia, obeying
///
/// ```text
/// L̈ = 3·V·(P − P₀) / W
/// ```
///
/// where `P` is the instantaneous pressure, `P₀` the target pressure and `W`
/// the barostat mass.  The box velocity `L̇` ([`ParrinelloRahmanBarostat::box_velocity`])
/// is persistent state: it is accumulated over successive steps, so the box
/// oscillates around the equilibrium volume instead of creeping towards it.
/// Positions are rescaled together with the box and, because the canonical
/// momentum `m·L²·ṡ` of the scaled coordinates is conserved, the peculiar
/// velocities scale inversely with the box edge (an expanding box cools the
/// system).
#[derive(Debug, Clone)]
pub struct ParrinelloRahmanBarostat {
    /// Target pressure.
    pub target_pressure: f64,
    /// Coupling time constant τ_P.
    ///
    /// Used to derive the barostat mass when `mass` is not positive: the mass
    /// is then chosen so that the box coupled to an ideal gas with the
    /// instantaneous kinetic energy of the system oscillates with period τ_P.
    pub tau_p: f64,
    /// Barostat mass `W` in the convention `L̈ = 3·V·(P − P₀) / W`.
    ///
    /// A value that is not finite and positive selects the mass derived from
    /// `tau_p` (see [`ParrinelloRahmanBarostat::effective_mass`]).
    pub mass: f64,
    /// Rate of change `dL/dt` of the box edge `L = V^(1/3)` (initially zero).
    pub box_velocity: f64,
}
impl ParrinelloRahmanBarostat {
    /// Create a new Parrinello-Rahman barostat with the box at rest.
    pub fn new(target_pressure: f64, tau_p: f64, mass: f64) -> Self {
        Self {
            target_pressure,
            tau_p,
            mass,
            box_velocity: 0.0,
        }
    }
    /// Barostat mass `W` actually used for the box equation of motion.
    ///
    /// Returns `mass` if it is finite and positive.  Otherwise the mass is
    /// derived from `tau_p`: linearising `L̈ = 3·V·(P − P₀)/W` around equilibrium
    /// for an ideal gas (`P = N·k_B·T/V`, `N·k_B·T = 2·KE/3`) gives the angular
    /// frequency `ω² = 9·N·k_B·T / (W·L)`, and requiring a period `τ_P = 2π/ω`
    ///
    /// ```text
    /// W = 3·KE·τ_P² / (2·π²·L)
    /// ```
    ///
    /// Returns `None` when neither source yields a positive finite mass (no
    /// `tau_p`, no kinetic energy, or a degenerate box): the box is then left
    /// untouched.
    pub fn effective_mass(&self, kinetic_energy: f64, box_edge: f64) -> Option<f64> {
        if self.mass.is_finite() && self.mass > 0.0 {
            return Some(self.mass);
        }
        let derived = 3.0 * kinetic_energy * self.tau_p * self.tau_p
            / (2.0 * std::f64::consts::PI * std::f64::consts::PI * box_edge);
        (derived.is_finite() && derived > 0.0).then_some(derived)
    }
    /// Compute the derivative of the enthalpy with respect to the cell matrix h.
    ///
    /// In the Parrinello-Rahman barostat the cell matrix `h` obeys an equation
    /// of motion driven by the pressure mismatch.  The enthalpy is:
    ///
    /// ```text
    /// H = U + P₀ · V(h)
    /// ```
    ///
    /// so the force on the cell matrix is:
    ///
    /// ```text
    /// dH/dh_ab = dU/dh_ab + P₀ · dV/dh_ab
    /// ```
    ///
    /// For an isotropic (cubic) box with edge length `L` and volume `V = L³`
    /// we have `dV/dL = 3·L²`.  This simplified implementation returns the
    /// scalar force on the (cubic) edge length:
    ///
    /// ```text
    /// dH/dL = -(P - P₀) · 3·V
    /// ```
    ///
    /// (positive = cell wants to expand; negative = cell wants to contract).
    ///
    /// # Arguments
    /// * `current_pressure` – instantaneous pressure P.
    /// * `volume`           – current box volume V.
    ///
    /// # Returns
    /// Scalar force on the cell edge length.
    pub fn compute_cell_derivative(&self, current_pressure: f64, volume: f64) -> f64 {
        -(current_pressure - self.target_pressure) * 3.0 * volume
    }
    /// Advance the box velocity by one time step.
    ///
    /// v_box += 3 · V · (P − P₀) / W · dt
    ///
    /// # Arguments
    /// * `box_vel`          – current box velocity (rate of change of box edge).
    /// * `current_pressure` – instantaneous pressure.
    /// * `volume`           – current box volume.
    /// * `dt`               – time step.
    pub fn box_velocity_update(
        &self,
        box_vel: f64,
        current_pressure: f64,
        volume: f64,
        dt: f64,
    ) -> f64 {
        box_vel + 3.0 * volume * (current_pressure - self.target_pressure) / self.mass * dt
    }
    /// Perform one full Parrinello-Rahman integration step with caller-owned
    /// box velocity.
    ///
    /// Updates the box velocity, rescales the box and all positions.  Unlike
    /// [`Barostat::apply`](super::functions::Barostat::apply) this variant
    /// keeps the box velocity in the caller's `box_vel`, uses the configured
    /// `mass` and does not touch particle velocities.
    ///
    /// # Arguments
    /// * `positions`        – mutable atomic positions.
    /// * `box_lengths`      – mutable box dimensions `[Lx, Ly, Lz]`.
    /// * `box_vel`          – mutable box velocity.
    /// * `current_pressure` – instantaneous pressure.
    /// * `volume`           – current box volume.
    /// * `dt`               – time step.
    pub fn apply_step(
        &self,
        positions: &mut [[f64; 3]],
        box_lengths: &mut [f64; 3],
        box_vel: &mut f64,
        current_pressure: f64,
        volume: f64,
        dt: f64,
    ) {
        *box_vel = self.box_velocity_update(*box_vel, current_pressure, volume, dt);
        let l = volume.cbrt();
        let dl = *box_vel * dt;
        let new_l = l + dl;
        if new_l > 0.0 {
            let scale = new_l / l;
            for pos in positions.iter_mut() {
                pos[0] *= scale;
                pos[1] *= scale;
                pos[2] *= scale;
            }
            box_lengths[0] *= scale;
            box_lengths[1] *= scale;
            box_lengths[2] *= scale;
        }
    }
}
/// Pressure coupling mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PressureCouplingMode {
    /// No pressure coupling (NVT/NVE).
    None,
    /// Isotropic: all box dimensions coupled uniformly.
    Isotropic,
    /// Semi-isotropic: x–y coupled together, z independent.
    SemiIsotropic,
    /// Anisotropic: all three axes coupled independently.
    Anisotropic,
}
/// Anisotropic Berendsen barostat with independent coupling per axis.
///
/// Allows different scale factors for x, y, z directions, useful for
/// non-cubic simulation boxes or anisotropic systems.
#[derive(Debug, Clone)]
pub struct AnisotropicBerendsenBarostat {
    /// Target pressure per axis \[P_x, P_y, P_z\].
    pub target_pressure: [f64; 3],
    /// Coupling time constant τ_P.
    pub tau_p: f64,
    /// Compressibility per axis \[β_x, β_y, β_z\].
    pub compressibility: [f64; 3],
}
impl AnisotropicBerendsenBarostat {
    /// Create a new anisotropic Berendsen barostat with uniform target pressure.
    pub fn new(target_pressure: f64, tau_p: f64, compressibility: f64) -> Self {
        Self {
            target_pressure: [target_pressure; 3],
            tau_p,
            compressibility: [compressibility; 3],
        }
    }
    /// Create with independent axis parameters.
    pub fn anisotropic(target_pressure: [f64; 3], tau_p: f64, compressibility: [f64; 3]) -> Self {
        Self {
            target_pressure,
            tau_p,
            compressibility,
        }
    }
    /// Per-axis scale factors given current diagonal pressures.
    pub fn scale_factors(&self, current_pressure: [f64; 3], dt: f64) -> [f64; 3] {
        let mut mu = [0.0f64; 3];
        for k in 0..3 {
            let mu_k = 1.0
                - self.compressibility[k] * dt / self.tau_p
                    * (self.target_pressure[k] - current_pressure[k]);
            mu[k] = mu_k.cbrt();
        }
        mu
    }
    /// Apply anisotropic scaling to positions and box lengths.
    pub fn apply_anisotropic(
        &self,
        positions: &mut [[f64; 3]],
        box_lengths: &mut [f64; 3],
        current_pressure: [f64; 3],
        dt: f64,
    ) {
        let mu = self.scale_factors(current_pressure, dt);
        for pos in positions.iter_mut() {
            pos[0] *= mu[0];
            pos[1] *= mu[1];
            pos[2] *= mu[2];
        }
        box_lengths[0] *= mu[0];
        box_lengths[1] *= mu[1];
        box_lengths[2] *= mu[2];
    }
}
/// Full 3×3 pressure tensor for anisotropic analysis.
///
/// The pressure tensor P_ab = (1/V) * \[ Σ_i m_i v_i_a v_i_b + W_ab \]
/// where W_ab = Σ_i r_i_a f_i_b is the virial tensor.
#[derive(Debug, Clone, Copy)]
pub struct PressureTensor {
    /// The 3×3 tensor stored as `tensor[row][col]`.
    pub tensor: [[f64; 3]; 3],
}
impl PressureTensor {
    /// Compute the full pressure tensor from velocities, masses, positions, and forces.
    ///
    /// # Arguments
    /// * `volume`    – simulation box volume.
    /// * `velocities` – `[vx, vy, vz]` per atom.
    /// * `masses`    – mass per atom.
    /// * `positions` – `[x, y, z]` per atom.
    /// * `forces`    – `[fx, fy, fz]` per atom.
    pub fn compute(
        volume: f64,
        velocities: &[[f64; 3]],
        masses: &[f64],
        positions: &[[f64; 3]],
        forces: &[[f64; 3]],
    ) -> Self {
        let n = velocities.len();
        let mut tensor = [[0.0f64; 3]; 3];
        for i in 0..n {
            let m = masses[i];
            for a in 0..3 {
                for b in 0..3 {
                    tensor[a][b] += m * velocities[i][a] * velocities[i][b];
                    tensor[a][b] += positions[i][a] * forces[i][b];
                }
            }
        }
        if volume > 1e-30 {
            let inv_v = 1.0 / volume;
            for row in &mut tensor {
                for v in row.iter_mut() {
                    *v *= inv_v;
                }
            }
        }
        Self { tensor }
    }
    /// Scalar (isotropic) pressure = (P_xx + P_yy + P_zz) / 3.
    pub fn scalar_pressure(&self) -> f64 {
        (self.tensor[0][0] + self.tensor[1][1] + self.tensor[2][2]) / 3.0
    }
    /// Diagonal components \[P_xx, P_yy, P_zz\].
    pub fn diagonal(&self) -> [f64; 3] {
        [self.tensor[0][0], self.tensor[1][1], self.tensor[2][2]]
    }
    /// Off-diagonal shear stresses \[P_xy, P_xz, P_yz\].
    pub fn off_diagonal(&self) -> [f64; 3] {
        [self.tensor[0][1], self.tensor[0][2], self.tensor[1][2]]
    }
    /// Pressure anisotropy = max(P_ii) - min(P_ii).
    pub fn anisotropy(&self) -> f64 {
        let d = self.diagonal();
        let max = d[0].max(d[1]).max(d[2]);
        let min = d[0].min(d[1]).min(d[2]);
        max - min
    }
    /// Compute the virial tensor W_ab = Σ_i r_i_a f_i_b.
    pub fn virial_tensor(positions: &[[f64; 3]], forces: &[[f64; 3]]) -> [[f64; 3]; 3] {
        let n = positions.len();
        let mut w = [[0.0f64; 3]; 3];
        for i in 0..n {
            for a in 0..3 {
                for b in 0..3 {
                    w[a][b] += positions[i][a] * forces[i][b];
                }
            }
        }
        w
    }
}
/// Full 3×3 box tensor for Parrinello-Rahman dynamics.
///
/// Tracks the simulation-cell matrix H and its velocity Hdot.  Useful for
/// anisotropic pressure coupling (e.g., for crystals or stretched systems).
///
/// The box is represented as a 3×3 matrix `h` where each column is a
/// cell vector.
#[derive(Debug, Clone)]
pub struct BoxTensor {
    /// Cell matrix H (rows = cell vectors or columns, depending on convention).
    /// Here: h\[i\]\[j\] is element (i, j).
    pub h: [[f64; 3]; 3],
    /// Velocity of box matrix Hdot.
    pub hdot: [[f64; 3]; 3],
    /// Target pressure tensor.
    pub p_ref: [[f64; 3]; 3],
    /// Barostat mass W.
    pub mass: f64,
}
impl BoxTensor {
    /// Create a cubic BoxTensor of edge length `l`.
    pub fn cubic(l: f64) -> Self {
        let h = [[l, 0.0, 0.0], [0.0, l, 0.0], [0.0, 0.0, l]];
        let p_ref = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        Self {
            h,
            hdot: [[0.0; 3]; 3],
            p_ref,
            mass: 100.0,
        }
    }
    /// Create with explicit parameters.
    pub fn new(h: [[f64; 3]; 3], p_ref: [[f64; 3]; 3], mass: f64) -> Self {
        Self {
            h,
            hdot: [[0.0; 3]; 3],
            p_ref,
            mass,
        }
    }
    /// Volume = det(H).
    pub fn volume(&self) -> f64 {
        let h = &self.h;
        h[0][0] * (h[1][1] * h[2][2] - h[1][2] * h[2][1])
            - h[0][1] * (h[1][0] * h[2][2] - h[1][2] * h[2][0])
            + h[0][2] * (h[1][0] * h[2][1] - h[1][1] * h[2][0])
    }
    /// Diagonal (lattice parameter lengths) \[a, b, c\].
    pub fn cell_lengths(&self) -> [f64; 3] {
        [
            (self.h[0][0] * self.h[0][0]
                + self.h[0][1] * self.h[0][1]
                + self.h[0][2] * self.h[0][2])
                .sqrt(),
            (self.h[1][0] * self.h[1][0]
                + self.h[1][1] * self.h[1][1]
                + self.h[1][2] * self.h[1][2])
                .sqrt(),
            (self.h[2][0] * self.h[2][0]
                + self.h[2][1] * self.h[2][1]
                + self.h[2][2] * self.h[2][2])
                .sqrt(),
        ]
    }
    /// Advance box dynamics by time step dt (simple Euler update).
    ///
    /// The box force on H is proportional to V * (P_current - P_ref).
    /// Hdot is updated and H is advanced.
    pub fn advance(&mut self, p_current: &[[f64; 3]; 3], dt: f64) {
        let vol = self.volume();
        for (i, p_row) in p_current.iter().enumerate() {
            for (j, &p_val) in p_row.iter().enumerate() {
                let g = vol * (p_val - self.p_ref[i][j]);
                self.hdot[i][j] += g / self.mass * dt;
                self.h[i][j] += self.hdot[i][j] * dt;
            }
        }
    }
    /// Isotropic pressure from trace: P = (P_xx + P_yy + P_zz) / 3.
    pub fn isotropic_pressure(p_tensor: &[[f64; 3]; 3]) -> f64 {
        (p_tensor[0][0] + p_tensor[1][1] + p_tensor[2][2]) / 3.0
    }
}
/// Natural logarithm of the Metropolis acceptance ratio of an isotropic volume
/// move `V_old → V_new` that was proposed uniformly in `ln V`:
///
/// `ln(acc) = −β·(ΔU + P₀·(V_new − V_old)) + (N + 1)·ln(V_new/V_old)`
///
/// with `β = 1/kt`.  Returns `−∞` (always reject) if either volume is not
/// positive and finite, or if `kt` is not positive and finite.
fn volume_move_log_acceptance(
    target_pressure: f64,
    delta_energy: f64,
    v_old: f64,
    v_new: f64,
    n_atoms: usize,
    kt: f64,
) -> f64 {
    let valid_volumes = v_old.is_finite() && v_new.is_finite() && v_old > 0.0 && v_new > 0.0;
    let valid_thermal_energy = kt.is_finite() && kt > 0.0;
    if !(valid_volumes && valid_thermal_energy) {
        return f64::NEG_INFINITY;
    }
    let beta = 1.0 / kt;
    -beta * (delta_energy + target_pressure * (v_new - v_old))
        + (n_atoms as f64 + 1.0) * (v_new / v_old).ln()
}
/// Extended Monte Carlo barostat with an internal LCG random number generator.
///
/// Can be used directly in a simulation loop without an external RNG.  A move
/// consists of three steps that the caller can drive separately:
///
/// 1. [`McBarostatRng::propose_volume`] draws the trial volume `V'`;
/// 2. the caller rescales its configuration to `V'` and evaluates the
///    potential-energy change `ΔU = U(V') − U(V)` **at that volume**;
/// 3. [`McBarostatRng::accept`] applies the Metropolis rule.
///
/// [`McBarostatRng::try_volume_move`] runs the three steps in one call, taking
/// the energy evaluation as a closure of the proposed volume.
#[derive(Debug, Clone)]
pub struct McBarostatRng {
    /// Target pressure P₀.
    pub target_pressure: f64,
    /// Width of the uniform proposal for `ln V`: the trial volume is
    /// `V·exp(δ)` with `δ` uniform in `[−move_size/2, +move_size/2)`.
    pub move_size: f64,
    /// Accepted move counter.
    pub n_accepted: u64,
    /// Total attempt counter.
    pub n_attempted: u64,
    /// LCG RNG state.
    pub(super) rng_state: u64,
}
impl McBarostatRng {
    /// Create a new MC barostat with default move size 0.01.
    pub fn new(target_pressure: f64) -> Self {
        Self {
            target_pressure,
            move_size: 0.01,
            n_accepted: 0,
            n_attempted: 0,
            rng_state: 55555,
        }
    }
    /// Create with explicit seed.
    pub fn with_seed(target_pressure: f64, move_size: f64, seed: u64) -> Self {
        Self {
            target_pressure,
            move_size,
            n_accepted: 0,
            n_attempted: 0,
            rng_state: seed,
        }
    }
    fn next_uniform(&mut self) -> f64 {
        self.rng_state = self
            .rng_state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.rng_state >> 11) as f64 / (1u64 << 53) as f64
    }
    /// Acceptance fraction so far.
    pub fn acceptance_rate(&self) -> f64 {
        if self.n_attempted == 0 {
            return 0.0;
        }
        self.n_accepted as f64 / self.n_attempted as f64
    }
    /// Draw a trial volume `V' = V·exp((ξ − ½)·move_size)` with `ξ` uniform in
    /// `[0, 1)`.
    ///
    /// The proposal is symmetric in `ln V` (not in `V`), which is why the
    /// acceptance rule of [`McBarostatRng::accept`] carries `N + 1` in its
    /// volume Jacobian.  Consumes one random number.  The caller must evaluate
    /// the energy change of the returned volume before calling `accept`.
    pub fn propose_volume(&mut self, current_volume: f64) -> f64 {
        let xi = self.next_uniform();
        current_volume * ((xi - 0.5) * self.move_size).exp()
    }
    /// Metropolis decision for the volume move `v_old → v_new`.
    ///
    /// Accepts with probability `min(1, exp(ln_acc))` where
    ///
    /// `ln_acc = −(ΔU + P₀·(V' − V))/kt + (N + 1)·ln(V'/V)`.
    ///
    /// # Arguments
    /// * `v_old`        – volume `V` before the move.
    /// * `v_new`        – proposed volume `V'` (from [`McBarostatRng::propose_volume`]).
    /// * `delta_energy` – potential-energy change `ΔU = U(V') − U(V)` of the
    ///   configuration rescaled to `V'`.
    /// * `n_atoms`      – number of atoms `N`.
    /// * `kt`           – thermal energy `k_B·T` in the units of `delta_energy`.
    ///
    /// A move with a non-positive or non-finite volume, a non-positive `kt` or
    /// a NaN energy change is rejected.  Updates the attempt and acceptance
    /// counters and consumes one random number.
    pub fn accept(
        &mut self,
        v_old: f64,
        v_new: f64,
        delta_energy: f64,
        n_atoms: usize,
        kt: f64,
    ) -> bool {
        self.n_attempted += 1;
        let log_acc = volume_move_log_acceptance(
            self.target_pressure,
            delta_energy,
            v_old,
            v_new,
            n_atoms,
            kt,
        );
        let threshold = self.next_uniform();
        // `log_acc` is NaN for a NaN energy change; both comparisons are then
        // false and the move is rejected.
        let accepted = log_acc >= 0.0 || threshold < log_acc.exp();
        if accepted {
            self.n_accepted += 1;
        }
        accepted
    }
    /// Propose a volume move, evaluate its energy cost and decide acceptance.
    ///
    /// Returns `Some(new_volume)` if accepted, `None` if rejected.
    ///
    /// `energy_change` is called once with the **proposed** volume `V'` and
    /// must return the potential-energy change `ΔU = U(V') − U(V)` of the
    /// configuration rescaled from the current volume to `V'`.  It is not
    /// called when the current or the proposed volume is not positive and
    /// finite; such a move is rejected and counted as an attempt.
    ///
    /// # Arguments
    /// * `current_volume` – current box volume.
    /// * `energy_change`  – closure returning `ΔU` for a trial volume.
    /// * `n_atoms`        – number of atoms in the system.
    /// * `temp`           – temperature T.
    /// * `kb`             – Boltzmann constant.
    pub fn try_volume_move(
        &mut self,
        current_volume: f64,
        mut energy_change: impl FnMut(f64) -> f64,
        n_atoms: usize,
        temp: f64,
        kb: f64,
    ) -> Option<f64> {
        let new_volume = self.propose_volume(current_volume);
        let usable = |volume: f64| volume.is_finite() && volume > 0.0;
        if !usable(current_volume) || !usable(new_volume) {
            self.n_attempted += 1;
            return None;
        }
        let delta_energy = energy_change(new_volume);
        self.accept(current_volume, new_volume, delta_energy, n_atoms, kb * temp)
            .then_some(new_volume)
    }
    /// Reset acceptance statistics.
    pub fn reset_statistics(&mut self) {
        self.n_accepted = 0;
        self.n_attempted = 0;
    }
}
/// Running statistics for box volume fluctuations.
///
/// Tracks mean, variance, min, max of the volume over time.  Useful for
/// equilibration monitoring and isothermal compressibility estimation.
#[derive(Debug, Clone)]
pub struct VolumeTracker {
    /// Number of recorded samples.
    pub count: u64,
    /// Running sum of volumes.
    pub(super) sum: f64,
    /// Running sum of volume^2.
    pub(super) sum_sq: f64,
    /// Minimum volume recorded.
    pub min_volume: f64,
    /// Maximum volume recorded.
    pub max_volume: f64,
    /// History of all recorded volumes (optional, bounded).
    pub(super) history: Vec<f64>,
    /// Maximum number of history entries to store.
    pub(super) max_history: usize,
}
impl VolumeTracker {
    /// Create a new volume tracker with bounded history.
    pub fn new(max_history: usize) -> Self {
        Self {
            count: 0,
            sum: 0.0,
            sum_sq: 0.0,
            min_volume: f64::MAX,
            max_volume: f64::MIN,
            history: Vec::new(),
            max_history,
        }
    }
    /// Record a volume sample.
    pub fn record(&mut self, volume: f64) {
        self.count += 1;
        self.sum += volume;
        self.sum_sq += volume * volume;
        if volume < self.min_volume {
            self.min_volume = volume;
        }
        if volume > self.max_volume {
            self.max_volume = volume;
        }
        if self.history.len() < self.max_history {
            self.history.push(volume);
        }
    }
    /// Mean volume.
    pub fn mean(&self) -> f64 {
        if self.count == 0 {
            return 0.0;
        }
        self.sum / self.count as f64
    }
    /// Variance of volume.
    pub fn variance(&self) -> f64 {
        if self.count < 2 {
            return 0.0;
        }
        let n = self.count as f64;
        (self.sum_sq / n) - (self.sum / n).powi(2)
    }
    /// Standard deviation of volume.
    pub fn std_dev(&self) -> f64 {
        self.variance().max(0.0).sqrt()
    }
    /// Relative volume fluctuation (σ_V / `V`).
    pub fn relative_fluctuation(&self) -> f64 {
        let m = self.mean();
        if m.abs() < 1e-30 {
            return 0.0;
        }
        self.std_dev() / m
    }
    /// Estimate isothermal compressibility from volume fluctuations.
    ///
    /// κ_T = `δV²` / (`V` · k_B · T)
    pub fn compressibility_estimate(&self, kb_t: f64) -> f64 {
        let m = self.mean();
        if m.abs() < 1e-30 || kb_t.abs() < 1e-30 {
            return 0.0;
        }
        self.variance() / (m * kb_t)
    }
    /// Return the stored volume history.
    pub fn history(&self) -> &[f64] {
        &self.history
    }
    /// Reset all accumulated statistics.
    pub fn reset(&mut self) {
        self.count = 0;
        self.sum = 0.0;
        self.sum_sq = 0.0;
        self.min_volume = f64::MAX;
        self.max_volume = f64::MIN;
        self.history.clear();
    }
}
/// Berendsen barostat for isotropic pressure coupling.
///
/// Rescales box dimensions and positions by
/// `mu = (1 − β·dt/τ·(P₀ − P))^(1/3)`.
#[derive(Debug, Clone)]
pub struct BerendsenBarostat {
    /// Target (reference) pressure.
    pub target_pressure: f64,
    /// Pressure coupling time constant τ_P.
    pub tau_p: f64,
    /// Isothermal compressibility β.
    pub compressibility: f64,
}
impl BerendsenBarostat {
    /// Create a new Berendsen barostat.
    ///
    /// # Arguments
    /// * `target_pressure`  – desired pressure P₀.
    /// * `tau_p`            – coupling time constant.
    /// * `compressibility`  – isothermal compressibility β.
    pub fn new(target_pressure: f64, tau_p: f64, compressibility: f64) -> Self {
        Self {
            target_pressure,
            tau_p,
            compressibility,
        }
    }
    /// Compute the isotropic scale factor μ for one time step.
    ///
    /// μ = (1 − β·dt/τ·(P₀ − P))^(1/3)
    pub fn scale_factor(&self, current_pressure: f64, dt: f64) -> f64 {
        self.scale_factor_to(current_pressure, self.target_pressure, dt)
    }
    /// Compute the isotropic scale factor μ for an explicit target pressure.
    ///
    /// μ = (1 − β·dt/τ·(P_target − P))^(1/3)
    ///
    /// [`BerendsenBarostat::scale_factor`] is this with the barostat's own
    /// `target_pressure`; the [`Barostat`](super::functions::Barostat) impl uses
    /// the target handed over by the simulation driver.
    pub fn scale_factor_to(&self, current_pressure: f64, target_pressure: f64, dt: f64) -> f64 {
        let mu_cube =
            1.0 - self.compressibility * dt / self.tau_p * (target_pressure - current_pressure);
        mu_cube.cbrt()
    }
    /// Compute the isotropic pressure-scaling matrix μ·I for one time step.
    ///
    /// In the Berendsen barostat the box is rescaled by the scalar factor
    /// `μ = (1 − β·dt/τ·(P₀ − P))^(1/3)`.  For anisotropic extensions each
    /// diagonal element of the cell matrix `h` is multiplied by the
    /// corresponding axis scale factor.  This helper returns the diagonal
    /// `[μ_x, μ_y, μ_z]` scaling vector for an isotropic (cubic) box, where
    /// all three elements are equal to `μ`.
    ///
    /// # Arguments
    /// * `current_pressure` – instantaneous scalar pressure P.
    /// * `dt`               – integration time step.
    ///
    /// # Returns
    /// `[μ, μ, μ]` – diagonal elements of the scaling matrix.
    pub fn compute_scaling_matrix(&self, current_pressure: f64, dt: f64) -> [f64; 3] {
        let mu = self.scale_factor(current_pressure, dt);
        [mu, mu, mu]
    }
    /// Scale all positions and box lengths by the Berendsen factor.
    ///
    /// # Arguments
    /// * `positions`        – mutable slice of `[x,y,z]` atomic positions.
    /// * `box_lengths`      – mutable `[Lx, Ly, Lz]` box dimensions.
    /// * `current_pressure` – instantaneous pressure.
    /// * `dt`               – integration time step.
    pub fn apply(
        &self,
        positions: &mut [[f64; 3]],
        box_lengths: &mut [f64; 3],
        current_pressure: f64,
        dt: f64,
    ) {
        let mu = self.scale_factor(current_pressure, dt);
        for pos in positions.iter_mut() {
            pos[0] *= mu;
            pos[1] *= mu;
            pos[2] *= mu;
        }
        box_lengths[0] *= mu;
        box_lengths[1] *= mu;
        box_lengths[2] *= mu;
    }
}
/// A no-op barostat (NVT/NVE ensemble).
#[derive(Debug, Clone, Default)]
pub struct NoBarostat;
impl NoBarostat {
    /// Create a no-op barostat.
    pub fn new() -> Self {
        Self
    }
}
/// Martyna-Tobias-Klein (MTK) barostat.
///
/// Couples box dynamics to a Nosé-Hoover chain thermostat for correct
/// NPT ensemble sampling.  This simplified implementation stores the
/// box velocity and a single thermostat variable for the barostat DOF.
#[derive(Debug, Clone)]
pub struct MtkBarostat {
    /// Target pressure P₀.
    pub target_pressure: f64,
    /// Barostat mass W (controls oscillation frequency).
    pub mass: f64,
    /// Box velocity (rate of change of ln(V)/3).
    pub epsilon_dot: f64,
    /// Thermostat variable for the barostat DOF.
    pub xi_baro: f64,
    /// Thermostat mass for the barostat DOF.
    pub q_baro: f64,
    /// Reference temperature (for the barostat thermostat).
    pub temperature: f64,
    /// Number of atoms (needed for pressure offset).
    pub n_atoms: usize,
}
impl MtkBarostat {
    /// Create a new MTK barostat.
    pub fn new(
        target_pressure: f64,
        mass: f64,
        q_baro: f64,
        temperature: f64,
        n_atoms: usize,
    ) -> Self {
        Self {
            target_pressure,
            mass,
            epsilon_dot: 0.0,
            xi_baro: 0.0,
            q_baro,
            temperature,
            n_atoms,
        }
    }
    /// Compute the barostat force:
    ///
    /// G_epsilon = 3V(P - P_ref) + (1/N_f) * 2*KE
    ///
    /// where N_f is the number of degrees of freedom.
    pub fn barostat_force(&self, volume: f64, current_pressure: f64, kinetic_energy: f64) -> f64 {
        let n_f = (3 * self.n_atoms) as f64;
        3.0 * volume * (current_pressure - self.target_pressure) + (2.0 * kinetic_energy) / n_f
    }
    /// Update epsilon_dot (box velocity) by dt.
    pub fn update_epsilon_dot(
        &mut self,
        volume: f64,
        current_pressure: f64,
        kinetic_energy: f64,
        dt: f64,
    ) {
        let g_eps = self.barostat_force(volume, current_pressure, kinetic_energy);
        let accel = g_eps / self.mass - self.xi_baro * self.epsilon_dot;
        self.epsilon_dot += accel * dt;
    }
    /// Update the barostat thermostat variable xi_baro.
    pub fn update_thermostat(&mut self, dt: f64) {
        let kb_t = self.temperature;
        let g_xi = self.mass * self.epsilon_dot * self.epsilon_dot - kb_t;
        self.xi_baro += (g_xi / self.q_baro) * dt;
    }
    /// Compute the isotropic scale factor for one step.
    pub fn scale_factor(&self, dt: f64) -> f64 {
        (self.epsilon_dot * dt).exp()
    }
    /// Apply MTK barostat step to positions and box.
    pub fn apply_step(
        &mut self,
        positions: &mut [[f64; 3]],
        box_lengths: &mut [f64; 3],
        volume: f64,
        current_pressure: f64,
        kinetic_energy: f64,
        dt: f64,
    ) {
        self.update_thermostat(dt * 0.5);
        self.update_epsilon_dot(volume, current_pressure, kinetic_energy, dt);
        self.update_thermostat(dt * 0.5);
        let scale = self.scale_factor(dt);
        for pos in positions.iter_mut() {
            pos[0] *= scale;
            pos[1] *= scale;
            pos[2] *= scale;
        }
        box_lengths[0] *= scale;
        box_lengths[1] *= scale;
        box_lengths[2] *= scale;
    }
    /// Barostat kinetic energy: 0.5 * W * epsilon_dot^2.
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * self.epsilon_dot * self.epsilon_dot
    }
    /// Thermostat energy contribution.
    pub fn thermostat_energy(&self) -> f64 {
        self.temperature * self.xi_baro
    }
}
/// Semi-isotropic Berendsen barostat.
///
/// Couples the lateral (x–y) dimensions isotropically and the z-direction
/// independently.  This is the standard approach for membrane simulations
/// where the bilayer surface tension is controlled separately from the
/// normal pressure.
#[derive(Debug, Clone)]
pub struct SemiIsotropicBerendsenBarostat {
    /// Target lateral pressure (P_xy).
    pub target_p_xy: f64,
    /// Target normal pressure (P_z).
    pub target_p_z: f64,
    /// Pressure coupling time constant τ_P.
    pub tau_p: f64,
    /// Lateral compressibility β_xy.
    pub compressibility_xy: f64,
    /// Normal compressibility β_z.
    pub compressibility_z: f64,
}
impl SemiIsotropicBerendsenBarostat {
    /// Create a new semi-isotropic barostat.
    pub fn new(
        target_p_xy: f64,
        target_p_z: f64,
        tau_p: f64,
        compressibility_xy: f64,
        compressibility_z: f64,
    ) -> Self {
        Self {
            target_p_xy,
            target_p_z,
            tau_p,
            compressibility_xy,
            compressibility_z,
        }
    }
    /// Compute the lateral scale factor μ_xy.
    ///
    /// μ_xy = (1 − β_xy · dt/τ · (P₀_xy − P_xy))^(1/3)
    pub fn scale_xy(&self, current_p_xy: f64, dt: f64) -> f64 {
        let mu =
            1.0 - self.compressibility_xy * dt / self.tau_p * (self.target_p_xy - current_p_xy);
        mu.cbrt()
    }
    /// Compute the normal scale factor μ_z.
    ///
    /// μ_z = (1 − β_z · dt/τ · (P₀_z − P_z))^(1/3)
    pub fn scale_z(&self, current_p_z: f64, dt: f64) -> f64 {
        let mu = 1.0 - self.compressibility_z * dt / self.tau_p * (self.target_p_z - current_p_z);
        mu.cbrt()
    }
    /// Apply semi-isotropic scaling to positions and box lengths.
    ///
    /// x and y coordinates scaled by `mu_xy`, z by `mu_z`.
    pub fn apply_scaling(
        &self,
        positions: &mut [[f64; 3]],
        box_lengths: &mut [f64; 3],
        current_p_xy: f64,
        current_p_z: f64,
        dt: f64,
    ) {
        let mu_xy = self.scale_xy(current_p_xy, dt);
        let mu_z = self.scale_z(current_p_z, dt);
        for pos in positions.iter_mut() {
            pos[0] *= mu_xy;
            pos[1] *= mu_xy;
            pos[2] *= mu_z;
        }
        box_lengths[0] *= mu_xy;
        box_lengths[1] *= mu_xy;
        box_lengths[2] *= mu_z;
    }
}
/// Monte Carlo barostat for the NPT ensemble.
///
/// Proposes isotropic volume moves that are uniform in `ln V` and
/// accepts/rejects them via a Metropolis criterion that includes the P·ΔV work
/// and the Jacobian `(N + 1)·ln(V_new/V_old)` (the `+ 1` comes from sampling
/// `ln V` instead of `V`; the `N` from the scaled particle coordinates).
#[derive(Debug, Clone)]
pub struct MonteCarloBarostat {
    /// Target pressure P₀.
    pub target_pressure: f64,
    /// Relative move size (fraction of ln-volume step).
    pub move_size: f64,
}
impl MonteCarloBarostat {
    /// Create a Monte Carlo barostat with default `move_size = 0.01`.
    pub fn new(target_pressure: f64) -> Self {
        Self {
            target_pressure,
            move_size: 0.01,
        }
    }
    /// Propose a new volume by a random logarithmic displacement.
    ///
    /// V_new = V_old * exp((ξ − 0.5) * move_size)
    ///
    /// # Arguments
    /// * `current_volume` – current box volume V_old.
    /// * `rng_uniform`    – uniform random number in \[0, 1).
    pub fn propose_volume_change(&self, current_volume: f64, rng_uniform: f64) -> f64 {
        current_volume * ((rng_uniform - 0.5) * self.move_size).exp()
    }
    /// Log-probability of accepting a volume move (Metropolis criterion).
    ///
    /// ln(acc) = −β·(ΔE + P₀·(V_new − V_old)) + (N + 1)·ln(V_new/V_old)
    ///
    /// where β = 1/(k_B·T).  The Jacobian uses `N + 1` because
    /// [`MonteCarloBarostat::propose_volume_change`] draws the new volume
    /// uniformly in `ln V`; for an ideal gas the sampled volume distribution is
    /// then `p(V) ∝ V^N·exp(−β·P₀·V)`.
    ///
    /// Returns `-∞` (always reject) if either volume is not positive and
    /// finite, because `ln(V_new/V_old)` does not exist, or if `k_B·T` is not
    /// positive and finite.
    ///
    /// # Arguments
    /// * `delta_energy`  – potential-energy change ΔE = E_new − E_old.
    /// * `v_old`         – volume V_old before the move.
    /// * `v_new`         – proposed volume V_new.
    /// * `n_atoms`       – number of atoms N.
    /// * `temp`          – temperature T.
    /// * `kb`            – Boltzmann constant k_B.
    pub fn acceptance_log_prob(
        &self,
        delta_energy: f64,
        v_old: f64,
        v_new: f64,
        n_atoms: usize,
        temp: f64,
        kb: f64,
    ) -> f64 {
        volume_move_log_acceptance(
            self.target_pressure,
            delta_energy,
            v_old,
            v_new,
            n_atoms,
            kb * temp,
        )
    }
    /// Compute the change in potential energy due to a volume move.
    ///
    /// In the NPT Monte Carlo acceptance criterion the relevant energy change
    /// is the sum of the *internal* potential-energy change and the *P·ΔV*
    /// pressure-volume work:
    ///
    /// ```text
    /// ΔU_total = ΔU_pot + P₀ · ΔV − (N + 1) · k_B · T · ln(V_new / V_old)
    /// ```
    ///
    /// The last term is the Jacobian for the Cartesian-coordinate change of
    /// variables when the box is rescaled (`N`) together with the measure
    /// factor of the log-volume proposal (`+ 1`).  This method returns
    /// `ΔU_total` so that the caller can perform a simple Metropolis test
    /// `min(1, exp(−β · ΔU_total))`.
    ///
    /// # Arguments
    /// * `delta_pot`     – direct change in potential energy ΔU_pot.
    /// * `v_old`         – old box volume V_old.
    /// * `v_new`         – new box volume V_new.
    /// * `n_atoms`       – number of atoms N.
    /// * `temp`          – temperature T.
    /// * `boltzmann_k`   – Boltzmann constant k_B.
    pub fn compute_volume_change_energy(
        &self,
        delta_pot: f64,
        v_old: f64,
        v_new: f64,
        n_atoms: usize,
        temp: f64,
        boltzmann_k: f64,
    ) -> f64 {
        let delta_v = v_new - v_old;
        let jacobian = if v_old > 1e-300 && v_new > 1e-300 {
            -(n_atoms as f64 + 1.0) * boltzmann_k * temp * (v_new / v_old).ln()
        } else {
            0.0
        };
        delta_pot + self.target_pressure * delta_v + jacobian
    }
}
