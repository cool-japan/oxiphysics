//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use std::collections::HashMap;

/// Result of an auction (winner index and price paid).
#[derive(Debug, Clone, PartialEq)]
pub struct AuctionResult {
    /// Index of the winning bidder.
    pub winner: usize,
    /// Price the winner pays.
    pub price: f64,
}
/// First-price sealed-bid auction.
///
/// The highest bidder wins and pays their own bid.
#[derive(Debug, Clone, Default)]
pub struct FirstPriceAuction;
impl FirstPriceAuction {
    /// Construct a new first-price auction solver.
    pub fn new() -> Self {
        Self
    }
    /// Run the auction with the given bids.
    ///
    /// Returns `None` if `bids` is empty.
    pub fn run(&self, bids: &[f64]) -> Option<AuctionResult> {
        if bids.is_empty() {
            return None;
        }
        let (winner, &price) = bids
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))?;
        Some(AuctionResult { winner, price })
    }
    /// Compute the symmetric Bayes-Nash equilibrium bid for a player with
    /// value `v` in a first-price auction with `n` bidders whose values are
    /// uniform on `[0, 1]`.
    ///
    /// Equilibrium bid: `b(v) = v * (n-1)/n`.
    pub fn equilibrium_bid(v: f64, n: usize) -> f64 {
        if n <= 1 {
            return v;
        }
        v * (n - 1) as f64 / n as f64
    }
}
/// English (ascending-bid) auction simulation.
#[derive(Debug, Clone)]
pub struct EnglishAuction {
    /// Bidders' true values.
    pub values: Vec<f64>,
}
impl EnglishAuction {
    /// Create a new English auction with the given bidder values.
    pub fn new(values: Vec<f64>) -> Self {
        Self { values }
    }
    /// Run the English auction (each bidder bids truthfully up to their value).
    ///
    /// Returns the winner index and clearing price.
    pub fn run(&self) -> Option<AuctionResult> {
        if self.values.is_empty() {
            return None;
        }
        let mut sorted: Vec<(usize, f64)> = self.values.iter().cloned().enumerate().collect();
        sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let winner = sorted[0].0;
        let price = if sorted.len() > 1 { sorted[1].1 } else { 0.0 };
        Some(AuctionResult { winner, price })
    }
}
/// A two-player action in a repeated Prisoner's Dilemma.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PdAction {
    /// Cooperate.
    Cooperate,
    /// Defect.
    Defect,
}
/// Tit-for-Tat strategy: cooperate on round 0, then copy opponent's last move.
#[derive(Debug, Clone, Default)]
pub struct TitForTat;
impl TitForTat {
    /// Create a new Tit-for-Tat agent.
    pub fn new() -> Self {
        Self
    }
    /// Choose the next action given the history of the opponent's moves.
    pub fn choose(&self, opponent_history: &[PdAction]) -> PdAction {
        opponent_history
            .last()
            .copied()
            .unwrap_or(PdAction::Cooperate)
    }
}
/// A congestion game where resources (e.g. roads) have load-dependent costs.
#[derive(Debug, Clone)]
pub struct CongestionGame {
    /// Number of resources.
    pub n_resources: usize,
    /// Cost function for each resource: `cost[r](load)` = cost for each user
    /// when `load` users share resource `r`.
    ///
    /// Stored as a vector of per-resource cost polynomials (coefficients in
    /// ascending order: `c0 + c1*x + c2*x^2 + …`).
    pub cost_poly: Vec<Vec<f64>>,
}
impl CongestionGame {
    /// Create a new congestion game with the given per-resource cost polynomials.
    pub fn new(cost_poly: Vec<Vec<f64>>) -> Self {
        let n = cost_poly.len();
        Self {
            n_resources: n,
            cost_poly,
        }
    }
    /// Evaluate the cost for resource `r` at load `x`.
    pub fn cost(&self, r: usize, x: f64) -> f64 {
        self.cost_poly[r]
            .iter()
            .enumerate()
            .map(|(k, &c)| c * x.powi(k as i32))
            .sum()
    }
    /// Compute the private cost for a player who uses a set of resources.
    ///
    /// `loads[r]` is the total number of users on resource `r`.
    pub fn private_cost(&self, resources: &[usize], loads: &[f64]) -> f64 {
        resources.iter().map(|&r| self.cost(r, loads[r])).sum()
    }
    /// Compute the social cost (sum over all users of their private costs).
    ///
    /// Equivalent to `Σ_r load[r] * cost(r, load[r])`.
    pub fn social_cost(&self, loads: &[f64]) -> f64 {
        loads
            .iter()
            .enumerate()
            .map(|(r, &l)| l * self.cost(r, l))
            .sum()
    }
    /// Compute the potential function value.
    ///
    /// For a congestion game with cost polynomials, the Rosenthal potential is:
    /// `Φ(σ) = Σ_r Σ_{k=1}^{load[r]} cost(r, k)`.
    pub fn rosenthal_potential(&self, loads: &[f64]) -> f64 {
        loads
            .iter()
            .enumerate()
            .map(|(r, &l)| {
                let n = l.round() as usize;
                (1..=n).map(|k| self.cost(r, k as f64)).sum::<f64>()
            })
            .sum()
    }
    /// Compute the Price of Anarchy (PoA) given a Nash equilibrium load vector
    /// and the social optimum load vector.
    ///
    /// `PoA = social_cost(nash) / social_cost(opt)`.
    pub fn price_of_anarchy(&self, nash_loads: &[f64], opt_loads: &[f64]) -> f64 {
        let sc_nash = self.social_cost(nash_loads);
        let sc_opt = self.social_cost(opt_loads);
        if sc_opt.abs() < 1e-30 {
            return 1.0;
        }
        sc_nash / sc_opt
    }
}
/// Replicator dynamics for evolutionary game theory.
///
/// Given a symmetric payoff matrix `A` where `A[i][j]` is the payoff to
/// strategy `i` when playing against strategy `j`, the replicator equation is:
///
/// `dx_i/dt = x_i ( f_i(x) - f̄(x) )`
///
/// where `f_i(x) = Σ_j A[i][j] x_j` and `f̄(x) = Σ_i x_i f_i(x)`.
#[derive(Debug, Clone)]
pub struct ReplicatorDynamics {
    /// Symmetric payoff matrix.
    pub payoff: Vec<Vec<f64>>,
}
impl ReplicatorDynamics {
    /// Create a new replicator dynamics model from a payoff matrix.
    pub fn new(payoff: Vec<Vec<f64>>) -> Self {
        Self { payoff }
    }
    /// Number of strategies.
    pub fn n_strategies(&self) -> usize {
        self.payoff.len()
    }
    /// Compute the fitness of each strategy given population state `x`.
    pub fn fitness(&self, x: &[f64]) -> Vec<f64> {
        let n = self.n_strategies();
        (0..n)
            .map(|i| (0..n).map(|j| self.payoff[i][j] * x[j]).sum())
            .collect()
    }
    /// Compute the mean fitness of the population.
    pub fn mean_fitness(&self, x: &[f64]) -> f64 {
        let f = self.fitness(x);
        x.iter().zip(f.iter()).map(|(xi, fi)| xi * fi).sum()
    }
    /// Compute `dx/dt` for the current state `x`.
    pub fn derivative(&self, x: &[f64]) -> Vec<f64> {
        let f = self.fitness(x);
        let mean = self.mean_fitness(x);
        x.iter()
            .zip(f.iter())
            .map(|(xi, fi)| xi * (fi - mean))
            .collect()
    }
    /// Advance the state by one Euler step of size `dt`.
    pub fn step_euler(&self, x: &[f64], dt: f64) -> Vec<f64> {
        let dx = self.derivative(x);
        let raw: Vec<f64> = x
            .iter()
            .zip(dx.iter())
            .map(|(xi, dxi)| xi + dt * dxi)
            .collect();
        normalise_simplex(&raw)
    }
    /// Run the replicator dynamics for `steps` Euler steps of size `dt`.
    ///
    /// Returns the final population state.
    pub fn run(&self, x0: &[f64], dt: f64, steps: usize) -> Vec<f64> {
        let mut x = x0.to_vec();
        for _ in 0..steps {
            x = self.step_euler(&x, dt);
        }
        x
    }
    /// Check whether `x` is an evolutionarily stable strategy (ESS) numerically.
    ///
    /// Uses the standard definition: strategy `i` (a pure strategy at index
    /// `i`) is an ESS if `A[i][i] > A[j][i]` for all `j ≠ i`.
    pub fn is_pure_ess(&self, strategy: usize) -> bool {
        let n = self.n_strategies();
        let aii = self.payoff[strategy][strategy];
        (0..n)
            .filter(|&j| j != strategy)
            .all(|j| aii > self.payoff[j][strategy])
    }
}
/// A direct revelation mechanism specifying an allocation rule and payment
/// rule for `n` agents.
#[derive(Debug, Clone)]
pub struct DirectMechanism {
    /// Number of agents.
    pub n_agents: usize,
}
impl DirectMechanism {
    /// Construct a new direct mechanism for `n_agents` agents.
    pub fn new(n_agents: usize) -> Self {
        Self { n_agents }
    }
    /// Check individual rationality (IR) for all agents.
    ///
    /// For each agent `i`, their expected utility `u_i = v_i * q_i - p_i`
    /// must be non-negative.
    ///
    /// `values[i]` is agent `i`'s reported value, `allocations[i]` the
    /// probability of winning, `payments[i]` the expected payment.
    pub fn is_individually_rational(
        &self,
        values: &[f64],
        allocations: &[f64],
        payments: &[f64],
    ) -> bool {
        (0..self.n_agents).all(|i| values[i] * allocations[i] - payments[i] >= -1e-9)
    }
    /// Check incentive compatibility (IC) for a single agent `i`.
    ///
    /// Given that all others report truthfully, agent `i` cannot gain by
    /// deviating.  This simplified check tests truthful vs. reporting `v_hat`.
    pub fn is_incentive_compatible_agent(
        &self,
        agent: usize,
        v_true: f64,
        q_truth: f64,
        p_truth: f64,
        q_dev: f64,
        p_dev: f64,
    ) -> bool {
        let _ = agent;
        let u_truth = v_true * q_truth - p_truth;
        let u_dev = v_true * q_dev - p_dev;
        u_truth >= u_dev - 1e-9
    }
    /// Myerson's revenue-optimal payment for a single-item auction with
    /// bidder `i` having value uniformly distributed on `[0, 1]`.
    ///
    /// Virtual value: `ψ(v) = v - (1 - F(v)) / f(v) = 2v - 1` for uniform.
    pub fn myerson_virtual_value_uniform(v: f64) -> f64 {
        2.0 * v - 1.0
    }
    /// Compute the revenue-maximising reserve price for a single-item auction
    /// with uniform values on `[0, 1]`.
    ///
    /// Set `ψ(r) = 0 → r = 0.5`.
    pub fn optimal_reserve_price_uniform() -> f64 {
        0.5
    }
}
/// Social welfare criteria.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WelfareCriterion {
    /// Utilitarian: sum of utilities.
    Utilitarian,
    /// Rawlsian: maximin (maximise the minimum utility).
    Rawlsian,
    /// Nash bargaining: product of utility gains from a disagreement point.
    NashBargaining,
}
/// Second-price sealed-bid (Vickrey) auction.
///
/// The highest bidder wins but pays the second-highest bid.
#[derive(Debug, Clone, Default)]
pub struct VickreyAuction;
impl VickreyAuction {
    /// Construct a new Vickrey auction solver.
    pub fn new() -> Self {
        Self
    }
    /// Run the Vickrey auction with the given bids.
    ///
    /// Returns `None` if `bids` is empty.
    pub fn run(&self, bids: &[f64]) -> Option<AuctionResult> {
        if bids.is_empty() {
            return None;
        }
        if bids.len() == 1 {
            return Some(AuctionResult {
                winner: 0,
                price: 0.0,
            });
        }
        let mut sorted_bids: Vec<(usize, f64)> = bids.iter().cloned().enumerate().collect();
        sorted_bids.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let winner = sorted_bids[0].0;
        let price = sorted_bids[1].1;
        Some(AuctionResult { winner, price })
    }
    /// In a Vickrey auction, truthful bidding is a dominant strategy.
    ///
    /// This function verifies that bidding `v` (true value) is optimal
    /// regardless of others' bids for the given set of competing bids.
    pub fn truthful_dominates(&self, v: f64, others: &[f64]) -> bool {
        let max_other = others.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let utility = |bid: f64| -> f64 { if bid > max_other { v - max_other } else { 0.0 } };
        let truth_utility = utility(v);
        let deviations = [
            0.0,
            0.5 * v,
            1.2 * v,
            max_other - 0.1,
            max_other + 0.1,
            2.0 * max_other,
        ];
        deviations
            .iter()
            .all(|&bid| truth_utility >= utility(bid) - 1e-9)
    }
}
/// Social welfare optimiser for a set of outcome alternatives.
#[derive(Debug, Clone)]
pub struct SocialWelfareOptimiser {
    /// Matrix of utilities: `utilities[outcome][agent]`.
    pub utilities: Vec<Vec<f64>>,
    /// Disagreement point for Nash bargaining.
    pub disagreement: Vec<f64>,
}
impl SocialWelfareOptimiser {
    /// Create a new social welfare optimiser.
    pub fn new(utilities: Vec<Vec<f64>>, disagreement: Vec<f64>) -> Self {
        Self {
            utilities,
            disagreement,
        }
    }
    /// Compute the welfare of a given outcome under the specified criterion.
    pub fn welfare(&self, outcome: usize, criterion: WelfareCriterion) -> f64 {
        let u = &self.utilities[outcome];
        match criterion {
            WelfareCriterion::Utilitarian => u.iter().sum(),
            WelfareCriterion::Rawlsian => u.iter().cloned().fold(f64::INFINITY, f64::min),
            WelfareCriterion::NashBargaining => u
                .iter()
                .zip(self.disagreement.iter())
                .map(|(&ui, &di)| (ui - di).max(0.0))
                .product(),
        }
    }
    /// Find the outcome maximising the given welfare criterion.
    pub fn optimal_outcome(&self, criterion: WelfareCriterion) -> usize {
        (0..self.utilities.len())
            .max_by(|&a, &b| {
                self.welfare(a, criterion)
                    .partial_cmp(&self.welfare(b, criterion))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(0)
    }
    /// Compute the Gini coefficient for the welfare distribution of the
    /// optimal utilitarian outcome.
    pub fn gini_coefficient(&self, outcome: usize) -> f64 {
        let mut u: Vec<f64> = self.utilities[outcome].to_vec();
        u.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let n = u.len() as f64;
        let total: f64 = u.iter().sum();
        if total.abs() < 1e-30 {
            return 0.0;
        }
        let numerator: f64 = u
            .iter()
            .enumerate()
            .map(|(i, &ui)| (2.0 * (i as f64 + 1.0) - n - 1.0) * ui)
            .sum();
        numerator / (n * total)
    }
}
/// Regret-matching algorithm for finding Nash equilibria.
///
/// Each player maintains cumulative regret for each action and mixes
/// proportional to positive cumulative regret.
#[derive(Debug, Clone)]
pub struct RegretMatcher {
    /// Number of actions.
    pub n_actions: usize,
    /// Cumulative regret for each action.
    pub regret: Vec<f64>,
    /// Cumulative strategy (for computing average strategy).
    pub strategy_sum: Vec<f64>,
}
impl RegretMatcher {
    /// Create a new regret matcher with the given number of actions.
    pub fn new(n_actions: usize) -> Self {
        Self {
            n_actions,
            regret: vec![0.0; n_actions],
            strategy_sum: vec![0.0; n_actions],
        }
    }
    /// Compute the current mixed strategy proportional to positive regret.
    pub fn current_strategy(&self) -> Vec<f64> {
        let pos: Vec<f64> = self.regret.iter().map(|&r| r.max(0.0)).collect();
        let total: f64 = pos.iter().sum();
        if total < 1e-12 {
            return vec![1.0 / self.n_actions as f64; self.n_actions];
        }
        pos.iter().map(|&p| p / total).collect()
    }
    /// Get the average strategy (for Nash equilibrium approximation).
    pub fn average_strategy(&self) -> Vec<f64> {
        let total: f64 = self.strategy_sum.iter().sum();
        if total < 1e-12 {
            return vec![1.0 / self.n_actions as f64; self.n_actions];
        }
        self.strategy_sum.iter().map(|&s| s / total).collect()
    }
    /// Update regret after observing action utilities `utilities[a]` when the
    /// current strategy was `strategy`.
    pub fn update_regret(&mut self, strategy: &[f64], utilities: &[f64]) {
        let realised: f64 = strategy
            .iter()
            .zip(utilities.iter())
            .map(|(s, u)| s * u)
            .sum();
        for a in 0..self.n_actions {
            self.regret[a] += utilities[a] - realised;
            self.strategy_sum[a] += strategy[a];
        }
    }
    /// Run `iterations` steps of regret matching for a 2-player zero-sum game.
    ///
    /// Returns the approximate Nash equilibrium strategies for both players.
    pub fn run_cfr_2p_zero_sum(payoff: &[Vec<f64>], iterations: usize) -> (Vec<f64>, Vec<f64>) {
        let m = payoff.len();
        if m == 0 {
            return (vec![], vec![]);
        }
        let n = payoff[0].len();
        let mut rm_a = RegretMatcher::new(m);
        let mut rm_b = RegretMatcher::new(n);
        for _ in 0..iterations {
            let sa = rm_a.current_strategy();
            let sb = rm_b.current_strategy();
            let ua: Vec<f64> = (0..m)
                .map(|i| (0..n).map(|j| sb[j] * payoff[i][j]).sum())
                .collect();
            let ub: Vec<f64> = (0..n)
                .map(|j| (0..m).map(|i| sa[i] * (-payoff[i][j])).sum())
                .collect();
            rm_a.update_regret(&sa, &ua);
            rm_b.update_regret(&sb, &ub);
        }
        (rm_a.average_strategy(), rm_b.average_strategy())
    }
}
/// Payoff parameters for the Prisoner's Dilemma.
///
/// Standard ordering: `T > R > P > S` and `2R > T + S`.
#[derive(Debug, Clone)]
pub struct PrisonersDilemma {
    /// Temptation payoff (defect vs cooperate).
    pub t: f64,
    /// Reward payoff (cooperate vs cooperate).
    pub r: f64,
    /// Punishment payoff (defect vs defect).
    pub p: f64,
    /// Sucker payoff (cooperate vs defect).
    pub s: f64,
}
impl PrisonersDilemma {
    /// Create a new Prisoner's Dilemma with the given payoffs.
    pub fn new(t: f64, r: f64, p: f64, s: f64) -> Self {
        Self { t, r, p, s }
    }
    /// Return the payoff for the row player given (row_action, col_action).
    pub fn payoff(&self, row: PdAction, col: PdAction) -> f64 {
        match (row, col) {
            (PdAction::Cooperate, PdAction::Cooperate) => self.r,
            (PdAction::Cooperate, PdAction::Defect) => self.s,
            (PdAction::Defect, PdAction::Cooperate) => self.t,
            (PdAction::Defect, PdAction::Defect) => self.p,
        }
    }
}
/// A cooperative game in characteristic-function form.
///
/// `v[S]` is the worth of coalition `S` (represented as a bitmask over `n`
/// players).
#[derive(Debug, Clone)]
pub struct CooperativeGame {
    /// Number of players.
    pub n: usize,
    /// Characteristic function mapping coalition bitmask → worth.
    pub v: HashMap<u32, f64>,
}
impl CooperativeGame {
    /// Create a new cooperative game.
    pub fn new(n: usize) -> Self {
        Self {
            n,
            v: HashMap::new(),
        }
    }
    /// Set the worth of a coalition given by its bitmask.
    pub fn set_worth(&mut self, coalition: u32, worth: f64) {
        self.v.insert(coalition, worth);
    }
    /// Get the worth of a coalition.
    pub fn worth(&self, coalition: u32) -> f64 {
        *self.v.get(&coalition).unwrap_or(&0.0)
    }
    /// Compute the Shapley value for all `n` players.
    ///
    /// The Shapley value is the unique value allocation satisfying efficiency,
    /// symmetry, and the dummy axiom.  It is computed by averaging a player's
    /// marginal contribution over all permutations of the players.
    ///
    /// Complexity: `O(n! * n)` — practical for n ≤ 10.
    pub fn shapley_values(&self) -> Vec<f64> {
        let n = self.n;
        let mut phi = vec![0.0; n];
        let total_perms = factorial(n) as f64;
        for perm in permutations(n) {
            let mut coalition: u32 = 0;
            for &player in &perm {
                let marginal = self.worth(coalition | (1 << player)) - self.worth(coalition);
                phi[player] += marginal / total_perms;
                coalition |= 1 << player;
            }
        }
        phi
    }
    /// Check whether the given `payoffs` vector is in the core.
    ///
    /// The core requires: (i) efficiency: `Σ payoffs = v(grand coalition)`, and
    /// (ii) coalitional rationality: for every coalition `S`, `Σ_{i∈S} payoffs[i] >= v(S)`.
    pub fn is_in_core(&self, payoffs: &[f64]) -> bool {
        let n = self.n;
        let total: f64 = payoffs.iter().sum();
        let grand = self.worth((1 << n) - 1);
        if (total - grand).abs() > 1e-9 {
            return false;
        }
        for mask in 0..((1u32 << n) - 1) {
            let coalition_val: f64 = (0..n)
                .filter(|&i| mask & (1 << i) != 0)
                .map(|i| payoffs[i])
                .sum();
            if coalition_val < self.worth(mask) - 1e-9 {
                return false;
            }
        }
        true
    }
    /// Compute the Banzhaf power index for all players.
    ///
    /// The raw Banzhaf index for player `i` is the number of coalitions where
    /// `i` is a swing voter.  This returns the normalised version.
    pub fn banzhaf_index(&self) -> Vec<f64> {
        let n = self.n;
        let mut swings = vec![0.0; n];
        for mask in 0u32..(1 << n) {
            for (i, swing) in swings.iter_mut().enumerate() {
                if mask & (1 << i) != 0 {
                    let without = mask & !(1 << i);
                    if self.worth(mask) > self.worth(without) + 1e-9 {
                        *swing += 1.0;
                    }
                }
            }
        }
        let total: f64 = swings.iter().sum();
        if total < 1e-30 {
            return vec![1.0 / n as f64; n];
        }
        swings.iter().map(|&s| s / total).collect()
    }
}
/// A two-player normal-form game stored as payoff matrices.
///
/// `payoff_a[i][j]` is player A's payoff when A plays row `i` and B plays
/// column `j`. `payoff_b` carries the same shape for player B.
#[derive(Debug, Clone)]
pub struct NormalFormGame {
    /// Payoff matrix for player A (rows × columns).
    pub payoff_a: Vec<Vec<f64>>,
    /// Payoff matrix for player B (rows × columns).
    pub payoff_b: Vec<Vec<f64>>,
}
impl NormalFormGame {
    /// Construct a new game from two payoff matrices.
    ///
    /// # Panics
    /// Panics if the matrices have different shapes.
    pub fn new(payoff_a: Vec<Vec<f64>>, payoff_b: Vec<Vec<f64>>) -> Self {
        assert_eq!(
            payoff_a.len(),
            payoff_b.len(),
            "payoff matrices must have the same number of rows"
        );
        if !payoff_a.is_empty() {
            assert_eq!(
                payoff_a[0].len(),
                payoff_b[0].len(),
                "payoff matrices must have the same number of columns"
            );
        }
        Self { payoff_a, payoff_b }
    }
    /// Return the number of rows (player A's strategies).
    pub fn rows(&self) -> usize {
        self.payoff_a.len()
    }
    /// Return the number of columns (player B's strategies).
    pub fn cols(&self) -> usize {
        if self.payoff_a.is_empty() {
            0
        } else {
            self.payoff_a[0].len()
        }
    }
    /// Find all pure-strategy Nash equilibria.
    ///
    /// A profile `(i*, j*)` is a Nash equilibrium when:
    /// * `payoff_a[i*][j*] >= payoff_a[i][j*]` for all i (A has no profitable deviation), and
    /// * `payoff_b[i*][j*] >= payoff_b[i*][j]` for all j (B has no profitable deviation).
    pub fn pure_nash_equilibria(&self) -> Vec<(usize, usize)> {
        let rows = self.rows();
        let cols = self.cols();
        let mut result = Vec::new();
        for (i, pa_row) in self.payoff_a.iter().enumerate() {
            for (j, &pa_ij) in pa_row.iter().enumerate() {
                let a_best = (0..rows).all(|k| pa_ij >= self.payoff_a[k][j]);
                let b_best = (0..cols).all(|l| self.payoff_b[i][j] >= self.payoff_b[i][l]);
                if a_best && b_best {
                    result.push((i, j));
                }
            }
        }
        result
    }
    /// Compute the expected payoff for player A given mixed strategies.
    ///
    /// `mixed_a[i]` is the probability A plays row `i`;
    /// `mixed_b[j]` is the probability B plays column `j`.
    pub fn expected_payoff_a(&self, mixed_a: &[f64], mixed_b: &[f64]) -> f64 {
        self.payoff_a
            .iter()
            .zip(mixed_a.iter())
            .map(|(row, &mai)| {
                row.iter()
                    .zip(mixed_b.iter())
                    .map(|(&pa_ij, &mbj)| mai * mbj * pa_ij)
                    .sum::<f64>()
            })
            .sum()
    }
    /// Compute the expected payoff for player B given mixed strategies.
    pub fn expected_payoff_b(&self, mixed_a: &[f64], mixed_b: &[f64]) -> f64 {
        self.payoff_b
            .iter()
            .zip(mixed_a.iter())
            .map(|(row, &mai)| {
                row.iter()
                    .zip(mixed_b.iter())
                    .map(|(&pb_ij, &mbj)| mai * mbj * pb_ij)
                    .sum::<f64>()
            })
            .sum()
    }
    /// Compute the best-response set for player A against a fixed `mixed_b`.
    ///
    /// Returns the set of row indices that maximise player A's expected payoff.
    pub fn best_response_a(&self, mixed_b: &[f64]) -> Vec<usize> {
        let expected: Vec<f64> = (0..self.rows())
            .map(|i| {
                (0..self.cols())
                    .map(|j| mixed_b[j] * self.payoff_a[i][j])
                    .sum()
            })
            .collect();
        let best = expected.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        expected
            .iter()
            .enumerate()
            .filter(|&(_, v)| (v - best).abs() < 1e-9)
            .map(|(i, _)| i)
            .collect()
    }
    /// Compute the best-response set for player B against a fixed `mixed_a`.
    pub fn best_response_b(&self, mixed_a: &[f64]) -> Vec<usize> {
        let expected: Vec<f64> = (0..self.cols())
            .map(|j| {
                (0..self.rows())
                    .map(|i| mixed_a[i] * self.payoff_b[i][j])
                    .sum()
            })
            .collect();
        let best = expected.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        expected
            .iter()
            .enumerate()
            .filter(|&(_, v)| (v - best).abs() < 1e-9)
            .map(|(j, _)| j)
            .collect()
    }
    /// Check whether a strategy profile is a Nash equilibrium (pure or mixed).
    ///
    /// Returns `true` if both players are playing best responses.
    pub fn is_nash_equilibrium(&self, mixed_a: &[f64], mixed_b: &[f64]) -> bool {
        let bra = self.best_response_a(mixed_b);
        let brb = self.best_response_b(mixed_a);
        let support_a_ok = mixed_a
            .iter()
            .enumerate()
            .filter(|&(_, p)| *p > 1e-9)
            .all(|(i, _)| bra.contains(&i));
        let support_b_ok = mixed_b
            .iter()
            .enumerate()
            .filter(|&(_, p)| *p > 1e-9)
            .all(|(j, _)| brb.contains(&j));
        support_a_ok && support_b_ok
    }
    /// Solve for the mixed-strategy Nash equilibrium of a 2×2 game.
    ///
    /// Returns `(p, q)` where `p` is the probability player A plays row 0 and
    /// `q` is the probability player B plays column 0. Returns `None` when no
    /// fully-mixed equilibrium exists.
    pub fn mixed_nash_2x2(&self) -> Option<(f64, f64)> {
        if self.rows() != 2 || self.cols() != 2 {
            return None;
        }
        let a00 = self.payoff_a[0][0];
        let a01 = self.payoff_a[0][1];
        let a10 = self.payoff_a[1][0];
        let a11 = self.payoff_a[1][1];
        let denom_q = (a00 - a01 - a10 + a11).abs();
        if denom_q < 1e-12 {
            return None;
        }
        let q = (a11 - a10) / (a00 - a01 - a10 + a11);
        let b00 = self.payoff_b[0][0];
        let b01 = self.payoff_b[0][1];
        let b10 = self.payoff_b[1][0];
        let b11 = self.payoff_b[1][1];
        let denom_p = (b00 - b01 - b10 + b11).abs();
        if denom_p < 1e-12 {
            return None;
        }
        let p = (b11 - b01) / (b00 - b01 - b10 + b11);
        if !(0.0..=1.0).contains(&p) || !(0.0..=1.0).contains(&q) {
            return None;
        }
        Some((p, q))
    }
    /// Compute the social welfare (sum of both players' payoffs) at a profile.
    pub fn social_welfare(&self, i: usize, j: usize) -> f64 {
        self.payoff_a[i][j] + self.payoff_b[i][j]
    }
    /// Find the social optimum: strategy profile maximising total welfare.
    pub fn social_optimum(&self) -> (usize, usize) {
        let mut best = (0, 0);
        let mut best_val = f64::NEG_INFINITY;
        for (i, _) in self.payoff_a.iter().enumerate() {
            for j in 0..self.cols() {
                let w = self.social_welfare(i, j);
                if w > best_val {
                    best_val = w;
                    best = (i, j);
                }
            }
        }
        best
    }
}
/// Minimax / maximin solver for two-player zero-sum games.
///
/// The payoff matrix stores player A's payoffs; player B's payoffs are the
/// negatives.
#[derive(Debug, Clone)]
pub struct ZeroSumGame {
    /// Player A's payoff matrix (player B gets -payoff_a).
    pub payoff: Vec<Vec<f64>>,
}
impl ZeroSumGame {
    /// Construct a zero-sum game from player A's payoff matrix.
    pub fn new(payoff: Vec<Vec<f64>>) -> Self {
        Self { payoff }
    }
    /// Return the number of rows.
    pub fn rows(&self) -> usize {
        self.payoff.len()
    }
    /// Return the number of columns.
    pub fn cols(&self) -> usize {
        if self.payoff.is_empty() {
            0
        } else {
            self.payoff[0].len()
        }
    }
    /// Compute the maximin value for player A (pure strategies).
    ///
    /// Player A maximises the minimum gain: `max_i min_j payoff[i][j]`.
    pub fn maximin(&self) -> f64 {
        (0..self.rows())
            .map(|i| {
                (0..self.cols())
                    .map(|j| self.payoff[i][j])
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Compute the minimax value for player B (pure strategies).
    ///
    /// Player B minimises the maximum loss: `min_j max_i payoff[i][j]`.
    pub fn minimax(&self) -> f64 {
        (0..self.cols())
            .map(|j| {
                (0..self.rows())
                    .map(|i| self.payoff[i][j])
                    .fold(f64::NEG_INFINITY, f64::max)
            })
            .fold(f64::INFINITY, f64::min)
    }
    /// Find a saddle point (pure-strategy Nash equilibrium) if one exists.
    ///
    /// A saddle point is an entry that is simultaneously a row-minimum and a
    /// column-maximum.  Returns `Some((row, col))` or `None`.
    pub fn saddle_point(&self) -> Option<(usize, usize)> {
        for (i, row) in self.payoff.iter().enumerate() {
            for (j, &pij) in row.iter().enumerate() {
                let is_row_min = row.iter().all(|&pl| pij <= pl);
                let is_col_max = self.payoff.iter().all(|r| pij >= r[j]);
                if is_row_min && is_col_max {
                    return Some((i, j));
                }
            }
        }
        None
    }
    /// Check whether the game has a saddle point.
    pub fn has_saddle_point(&self) -> bool {
        self.saddle_point().is_some()
    }
    /// Compute the value of the game under mixed strategies using linear
    /// programming (approximate iterative method).
    ///
    /// Uses the fictitious-play / Brown-Robinson iteration for `iterations`
    /// steps.  Returns the approximate value.
    pub fn mixed_value_fictitious_play(&self, iterations: usize) -> f64 {
        let m = self.rows();
        let n = self.cols();
        if m == 0 || n == 0 {
            return 0.0;
        }
        let mut count_a = vec![0u64; m];
        let mut count_b = vec![0u64; n];
        count_a[0] = 1;
        count_b[0] = 1;
        for t in 1..iterations {
            let br_a = (0..m)
                .max_by(|&i1, &i2| {
                    let v1: f64 = (0..n).map(|j| count_b[j] as f64 * self.payoff[i1][j]).sum();
                    let v2: f64 = (0..n).map(|j| count_b[j] as f64 * self.payoff[i2][j]).sum();
                    v1.partial_cmp(&v2).unwrap_or(std::cmp::Ordering::Equal)
                })
                .unwrap_or(0);
            let br_b = (0..n)
                .min_by(|&j1, &j2| {
                    let v1: f64 = (0..m).map(|i| count_a[i] as f64 * self.payoff[i][j1]).sum();
                    let v2: f64 = (0..m).map(|i| count_a[i] as f64 * self.payoff[i][j2]).sum();
                    v1.partial_cmp(&v2).unwrap_or(std::cmp::Ordering::Equal)
                })
                .unwrap_or(0);
            count_a[br_a] += 1;
            count_b[br_b] += 1;
            let _ = t;
        }
        let total = iterations as f64;
        let mixed_a: Vec<f64> = count_a.iter().map(|&c| c as f64 / total).collect();
        let mixed_b: Vec<f64> = count_b.iter().map(|&c| c as f64 / total).collect();
        let value: f64 = self
            .payoff
            .iter()
            .zip(mixed_a.iter())
            .map(|(row, &mai)| {
                row.iter()
                    .zip(mixed_b.iter())
                    .map(|(&pij, &mbj)| mai * mbj * pij)
                    .sum::<f64>()
            })
            .sum();
        value
    }
    /// Compute the dominant strategy for player A if one exists.
    ///
    /// Returns `Some(row)` if row strictly dominates all others.
    pub fn dominant_strategy_a(&self) -> Option<usize> {
        let m = self.rows();
        let n = self.cols();
        for i in 0..m {
            let dominates = (0..m)
                .filter(|&k| k != i)
                .all(|k| (0..n).all(|j| self.payoff[i][j] > self.payoff[k][j]));
            if dominates {
                return Some(i);
            }
        }
        None
    }
}
/// Grim Trigger strategy: cooperate until the opponent defects once, then
/// defect forever.
#[derive(Debug, Clone)]
pub struct GrimTrigger {
    pub(super) triggered: bool,
}
impl GrimTrigger {
    /// Create a new Grim Trigger agent.
    pub fn new() -> Self {
        Self { triggered: false }
    }
    /// Update the state based on the opponent's latest action.
    pub fn update(&mut self, opponent_action: PdAction) {
        if opponent_action == PdAction::Defect {
            self.triggered = true;
        }
    }
    /// Choose the next action.
    pub fn choose(&self) -> PdAction {
        if self.triggered {
            PdAction::Defect
        } else {
            PdAction::Cooperate
        }
    }
}
/// Two-player Nash bargaining problem.
///
/// Finds the utility allocation that maximises the Nash product
/// `(u1 - d1) * (u2 - d2)` subject to the feasibility set.
#[derive(Debug, Clone)]
pub struct NashBargaining {
    /// Feasible utility pairs as `(u1, u2)`.
    pub feasible: Vec<(f64, f64)>,
    /// Disagreement point `(d1, d2)`.
    pub disagreement: (f64, f64),
}
impl NashBargaining {
    /// Create a new Nash bargaining problem.
    pub fn new(feasible: Vec<(f64, f64)>, disagreement: (f64, f64)) -> Self {
        Self {
            feasible,
            disagreement,
        }
    }
    /// Find the Nash bargaining solution from the finite feasible set.
    ///
    /// Returns the feasible pair that maximises the Nash product, or `None`
    /// if no individually-rational feasible point exists.
    pub fn solution(&self) -> Option<(f64, f64)> {
        let (d1, d2) = self.disagreement;
        self.feasible
            .iter()
            .filter(|&&(u1, u2)| u1 >= d1 && u2 >= d2)
            .max_by(|&&(u1a, u2a), &&(u1b, u2b)| {
                let pa = (u1a - d1) * (u2a - d2);
                let pb = (u1b - d1) * (u2b - d2);
                pa.partial_cmp(&pb).unwrap_or(std::cmp::Ordering::Equal)
            })
            .copied()
    }
    /// Kalai-Smorodinsky bargaining solution (from a two-point feasible set).
    ///
    /// Returns the point on the Pareto frontier that equalises relative gains.
    pub fn kalai_smorodinsky(&self) -> Option<(f64, f64)> {
        let (d1, d2) = self.disagreement;
        let max1 = self
            .feasible
            .iter()
            .map(|&(u1, _)| u1)
            .fold(f64::NEG_INFINITY, f64::max);
        let max2 = self
            .feasible
            .iter()
            .map(|&(_, u2)| u2)
            .fold(f64::NEG_INFINITY, f64::max);
        if max1 <= d1 || max2 <= d2 {
            return None;
        }
        let target_ratio = (max1 - d1) / (max2 - d2);
        self.feasible
            .iter()
            .filter(|&&(u1, u2)| u1 >= d1 && u2 >= d2)
            .min_by(|&&(u1a, u2a), &&(u1b, u2b)| {
                let ra = (u1a - d1) / ((u2a - d2).max(1e-30));
                let rb = (u1b - d1) / ((u2b - d2).max(1e-30));
                (ra - target_ratio)
                    .abs()
                    .partial_cmp(&(rb - target_ratio).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .copied()
    }
}
/// A node in an extensive-form game tree.
#[derive(Debug, Clone)]
pub struct GameNode {
    /// Index of the player to move at this node (`0` = player 1, `1` = player 2, etc.).
    /// `usize::MAX` marks a terminal node.
    pub player: usize,
    /// Payoff vector for all players at this terminal node (only used when
    /// `player == usize::MAX`).
    pub payoffs: Vec<f64>,
    /// Child nodes (one per available action).
    pub children: Vec<GameNode>,
    /// Label for this node (optional, for display/debugging).
    pub label: String,
}
impl GameNode {
    /// Create a new decision node for `player` with given children.
    pub fn decision(player: usize, label: impl Into<String>, children: Vec<GameNode>) -> Self {
        Self {
            player,
            payoffs: vec![],
            children,
            label: label.into(),
        }
    }
    /// Create a terminal node with the given payoff vector.
    pub fn terminal(payoffs: Vec<f64>, label: impl Into<String>) -> Self {
        Self {
            player: usize::MAX,
            payoffs,
            children: vec![],
            label: label.into(),
        }
    }
    /// Return `true` if this is a terminal node.
    pub fn is_terminal(&self) -> bool {
        self.player == usize::MAX
    }
    /// Perform backward induction and return the optimal payoff vector
    /// at this subtree root.
    ///
    /// Each player maximises their own payoff component.
    pub fn backward_induction(&self) -> Vec<f64> {
        if self.is_terminal() {
            return self.payoffs.clone();
        }
        let child_payoffs: Vec<Vec<f64>> = self
            .children
            .iter()
            .map(|c| c.backward_induction())
            .collect();
        let mover = self.player;
        child_payoffs
            .into_iter()
            .max_by(|a, b| {
                let va = if mover < a.len() {
                    a[mover]
                } else {
                    f64::NEG_INFINITY
                };
                let vb = if mover < b.len() {
                    b[mover]
                } else {
                    f64::NEG_INFINITY
                };
                va.partial_cmp(&vb).unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or_default()
    }
    /// Count the total number of terminal nodes in this subtree.
    pub fn terminal_count(&self) -> usize {
        if self.is_terminal() {
            return 1;
        }
        self.children.iter().map(|c| c.terminal_count()).sum()
    }
    /// Compute the depth of this subtree.
    pub fn depth(&self) -> usize {
        if self.is_terminal() {
            return 0;
        }
        1 + self.children.iter().map(|c| c.depth()).max().unwrap_or(0)
    }
}
/// A strategy in a game with imperfect information, represented as a mapping
/// from information-set identifiers to action distributions.
#[derive(Debug, Clone)]
pub struct BehaviouralStrategy {
    /// Maps each information set id to a probability distribution over actions.
    pub policy: HashMap<u32, Vec<f64>>,
}
impl BehaviouralStrategy {
    /// Create a new behavioural strategy.
    pub fn new() -> Self {
        Self {
            policy: HashMap::new(),
        }
    }
    /// Set the action distribution for an information set.
    pub fn set(&mut self, info_set: u32, distribution: Vec<f64>) {
        self.policy.insert(info_set, distribution);
    }
    /// Get the probability of action `a` at information set `info_set`.
    pub fn prob(&self, info_set: u32, a: usize) -> f64 {
        self.policy
            .get(&info_set)
            .and_then(|d| d.get(a))
            .copied()
            .unwrap_or(0.0)
    }
    /// Check whether the strategy is valid (all distributions sum to ~1).
    pub fn is_valid(&self) -> bool {
        self.policy
            .values()
            .all(|d| (d.iter().sum::<f64>() - 1.0).abs() < 1e-9)
    }
}
/// A correlated equilibrium is a probability distribution over strategy
/// profiles such that no player benefits from deviating conditional on their
/// recommended action.
#[derive(Debug, Clone)]
pub struct CorrelatedEquilibrium {
    /// Distribution over strategy profiles (indexed as flat row-major for a
    /// 2-player game with `rows` × `cols` strategies).
    pub distribution: Vec<f64>,
    /// Number of rows (player A strategies).
    pub rows: usize,
    /// Number of columns (player B strategies).
    pub cols: usize,
}
impl CorrelatedEquilibrium {
    /// Create a new correlated equilibrium distribution.
    pub fn new(distribution: Vec<f64>, rows: usize, cols: usize) -> Self {
        Self {
            distribution,
            rows,
            cols,
        }
    }
    /// Return the probability of profile `(i, j)`.
    pub fn prob(&self, i: usize, j: usize) -> f64 {
        self.distribution[i * self.cols + j]
    }
    /// Verify the correlated equilibrium conditions for a 2-player game.
    ///
    /// For each player and each recommended action, verify that following the
    /// recommendation is at least as good as any deviation.
    pub fn verify(&self, payoff_a: &[Vec<f64>], payoff_b: &[Vec<f64>]) -> bool {
        let rows = self.rows;
        let cols = self.cols;
        for i in 0..rows {
            let prob_i: f64 = (0..cols).map(|j| self.prob(i, j)).sum();
            if prob_i < 1e-12 {
                continue;
            }
            for i_dev in 0..rows {
                if i_dev == i {
                    continue;
                }
                let gain: f64 = (0..cols)
                    .map(|j| self.prob(i, j) * (payoff_a[i_dev][j] - payoff_a[i][j]))
                    .sum();
                if gain > 1e-9 {
                    return false;
                }
            }
        }
        for j in 0..cols {
            let prob_j: f64 = (0..rows).map(|i| self.prob(i, j)).sum();
            if prob_j < 1e-12 {
                continue;
            }
            for j_dev in 0..cols {
                if j_dev == j {
                    continue;
                }
                let gain: f64 = (0..rows)
                    .map(|i| self.prob(i, j) * (payoff_b[i][j_dev] - payoff_b[i][j]))
                    .sum();
                if gain > 1e-9 {
                    return false;
                }
            }
        }
        true
    }
}
