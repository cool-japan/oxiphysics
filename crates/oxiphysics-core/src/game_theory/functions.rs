//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{GrimTrigger, PdAction, PrisonersDilemma, TitForTat};

/// Normalise a vector to sum to 1 (project onto the probability simplex).
pub fn normalise_simplex(x: &[f64]) -> Vec<f64> {
    let total: f64 = x.iter().map(|&v| v.max(0.0)).sum();
    if total < 1e-30 {
        let n = x.len();
        return vec![1.0 / n as f64; n];
    }
    x.iter().map(|&v| v.max(0.0) / total).collect()
}
/// Compute `n!` for small `n`.
pub fn factorial(n: usize) -> usize {
    (1..=n).product()
}
/// Generate all permutations of `{0, 1, ..., n-1}`.
pub fn permutations(n: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut perm: Vec<usize> = (0..n).collect();
    permutations_helper(&mut perm, 0, &mut result);
    result
}
/// Helper for generating permutations recursively.
pub fn permutations_helper(perm: &mut Vec<usize>, start: usize, result: &mut Vec<Vec<usize>>) {
    if start == perm.len() {
        result.push(perm.clone());
        return;
    }
    for i in start..perm.len() {
        perm.swap(start, i);
        permutations_helper(perm, start + 1, result);
        perm.swap(start, i);
    }
}
/// Simulate a repeated Prisoner's Dilemma between Tit-for-Tat and Grim Trigger.
///
/// Returns `(payoff_tft, payoff_grim)` after `rounds` rounds.
pub fn simulate_tft_vs_grim(pd: &PrisonersDilemma, rounds: usize, discount: f64) -> (f64, f64) {
    let tft = TitForTat::new();
    let mut grim = GrimTrigger::new();
    let mut history_tft: Vec<PdAction> = vec![];
    let mut history_grim: Vec<PdAction> = vec![];
    let mut total_tft = 0.0;
    let mut total_grim = 0.0;
    let mut weight = 1.0;
    for _ in 0..rounds {
        let a_tft = tft.choose(&history_grim);
        let a_grim = grim.choose();
        total_tft += weight * pd.payoff(a_tft, a_grim);
        total_grim += weight * pd.payoff(a_grim, a_tft);
        history_tft.push(a_tft);
        history_grim.push(a_grim);
        grim.update(a_tft);
        weight *= discount;
    }
    (total_tft, total_grim)
}
/// Folk theorem check: whether mutual cooperation is sustainable in an
/// infinitely repeated game.
///
/// Returns `true` if the discount factor `delta` is high enough that
/// cooperation is a Nash equilibrium of the repeated game.
///
/// Condition: `delta >= (T - R) / (T - P)`.
pub fn folk_theorem_cooperation_sustainable(pd: &PrisonersDilemma, delta: f64) -> bool {
    let threshold = (pd.t - pd.r) / (pd.t - pd.p);
    delta >= threshold
}
/// Compute the Price of Stability given a list of Nash equilibrium social
/// costs and the social optimum cost.
///
/// `PoS = min_{NE} social_cost(NE) / social_cost(OPT)`.
pub fn price_of_stability(nash_costs: &[f64], opt_cost: f64) -> f64 {
    if opt_cost.abs() < 1e-30 {
        return 1.0;
    }
    let min_ne = nash_costs.iter().cloned().fold(f64::INFINITY, f64::min);
    min_ne / opt_cost
}
/// Compute the Price of Anarchy given a list of Nash equilibrium social costs
/// and the social optimum cost.
///
/// `PoA = max_{NE} social_cost(NE) / social_cost(OPT)`.
pub fn price_of_anarchy_global(nash_costs: &[f64], opt_cost: f64) -> f64 {
    if opt_cost.abs() < 1e-30 {
        return 1.0;
    }
    let max_ne = nash_costs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    max_ne / opt_cost
}
/// Compute the efficiency ratio: `social_welfare(NE) / social_welfare(OPT)`.
///
/// Ranges from 0 to 1 (higher is better).
pub fn efficiency_ratio(ne_welfare: f64, opt_welfare: f64) -> f64 {
    if opt_welfare.abs() < 1e-30 {
        return 1.0;
    }
    (ne_welfare / opt_welfare).min(1.0)
}
#[cfg(test)]
mod tests_normal_form {

    use crate::game_theory::NormalFormGame;

    fn prisoners_dilemma_game() -> NormalFormGame {
        NormalFormGame::new(
            vec![vec![3.0, 0.0], vec![5.0, 1.0]],
            vec![vec![3.0, 5.0], vec![0.0, 1.0]],
        )
    }
    #[test]
    fn test_pd_pure_nash() {
        let g = prisoners_dilemma_game();
        let ne = g.pure_nash_equilibria();
        assert_eq!(ne, vec![(1, 1)], "PD's unique NE is (Defect, Defect)");
    }
    #[test]
    fn test_no_pure_nash_matching_pennies() {
        let g = NormalFormGame::new(
            vec![vec![1.0, -1.0], vec![-1.0, 1.0]],
            vec![vec![-1.0, 1.0], vec![1.0, -1.0]],
        );
        assert!(g.pure_nash_equilibria().is_empty());
    }
    #[test]
    fn test_expected_payoff_uniform_mix() {
        let g = prisoners_dilemma_game();
        let mix = vec![0.5, 0.5];
        let ea = g.expected_payoff_a(&mix, &mix);
        assert!((ea - 2.25).abs() < 1e-9);
    }
    #[test]
    fn test_best_response_a_pure() {
        let g = prisoners_dilemma_game();
        let br = g.best_response_a(&[0.0, 1.0]);
        assert_eq!(br, vec![1]);
    }
    #[test]
    fn test_mixed_nash_2x2_matching_pennies() {
        let g = NormalFormGame::new(
            vec![vec![1.0, -1.0], vec![-1.0, 1.0]],
            vec![vec![-1.0, 1.0], vec![1.0, -1.0]],
        );
        let ne = g.mixed_nash_2x2();
        assert!(ne.is_some());
        let (p, q) = ne.unwrap();
        assert!((p - 0.5).abs() < 1e-9);
        assert!((q - 0.5).abs() < 1e-9);
    }
    #[test]
    fn test_social_optimum_pd() {
        let g = prisoners_dilemma_game();
        let opt = g.social_optimum();
        assert_eq!(opt, (0, 0), "PD social optimum is (Cooperate, Cooperate)");
    }
    #[test]
    fn test_game_rows_cols() {
        let g = prisoners_dilemma_game();
        assert_eq!(g.rows(), 2);
        assert_eq!(g.cols(), 2);
    }
    #[test]
    fn test_is_nash_equilibrium_valid() {
        let g = prisoners_dilemma_game();
        assert!(g.is_nash_equilibrium(&[0.0, 1.0], &[0.0, 1.0]));
    }
    #[test]
    fn test_social_welfare_computation() {
        let g = prisoners_dilemma_game();
        assert!((g.social_welfare(0, 0) - 6.0).abs() < 1e-9);
        assert!((g.social_welfare(1, 1) - 2.0).abs() < 1e-9);
    }
}
#[cfg(test)]
mod tests_zero_sum {

    use crate::game_theory::ZeroSumGame;

    fn rock_paper_scissors() -> ZeroSumGame {
        ZeroSumGame::new(vec![
            vec![0.0, -1.0, 1.0],
            vec![1.0, 0.0, -1.0],
            vec![-1.0, 1.0, 0.0],
        ])
    }
    #[test]
    fn test_saddle_point_simple() {
        let g = ZeroSumGame::new(vec![vec![3.0, 2.0], vec![4.0, 1.0]]);
        let sp = g.saddle_point();
        assert!(sp.is_some());
    }
    #[test]
    fn test_no_saddle_point_rps() {
        let g = rock_paper_scissors();
        assert!(!g.has_saddle_point());
    }
    #[test]
    fn test_minimax_ge_maximin() {
        let g = rock_paper_scissors();
        assert!(g.minimax() >= g.maximin() - 1e-9);
    }
    #[test]
    fn test_rps_maximin() {
        let g = rock_paper_scissors();
        assert!((g.maximin() - (-1.0)).abs() < 1e-9);
    }
    #[test]
    fn test_rps_minimax() {
        let g = rock_paper_scissors();
        assert!((g.minimax() - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_fictitious_play_converges() {
        let g = rock_paper_scissors();
        let v = g.mixed_value_fictitious_play(5000);
        assert!(v.abs() < 0.15, "RPS value should be near 0, got {v}");
    }
    #[test]
    fn test_dominant_strategy_none_rps() {
        let g = rock_paper_scissors();
        assert!(g.dominant_strategy_a().is_none());
    }
}
#[cfg(test)]
mod tests_extensive_form {

    use crate::game_theory::GameNode;

    fn simple_sequential_game() -> GameNode {
        GameNode::decision(
            0,
            "root",
            vec![
                GameNode::decision(
                    1,
                    "L",
                    vec![
                        GameNode::terminal(vec![3.0, 2.0], "LL"),
                        GameNode::terminal(vec![1.0, 4.0], "LR"),
                    ],
                ),
                GameNode::decision(
                    1,
                    "R",
                    vec![
                        GameNode::terminal(vec![2.0, 1.0], "RL"),
                        GameNode::terminal(vec![4.0, 3.0], "RR"),
                    ],
                ),
            ],
        )
    }
    #[test]
    fn test_backward_induction_depth2() {
        let game = simple_sequential_game();
        let val = game.backward_induction();
        assert_eq!(val, vec![4.0, 3.0]);
    }
    #[test]
    fn test_terminal_count() {
        let game = simple_sequential_game();
        assert_eq!(game.terminal_count(), 4);
    }
    #[test]
    fn test_game_depth() {
        let game = simple_sequential_game();
        assert_eq!(game.depth(), 2);
    }
    #[test]
    fn test_terminal_node() {
        let t = GameNode::terminal(vec![1.0, 2.0], "t");
        assert!(t.is_terminal());
        assert_eq!(t.terminal_count(), 1);
        assert_eq!(t.depth(), 0);
    }
}
#[cfg(test)]
mod tests_evolutionary {

    use crate::game_theory::ReplicatorDynamics;

    use crate::game_theory::normalise_simplex;

    fn hawk_dove_payoff() -> Vec<Vec<f64>> {
        vec![vec![-0.5, 2.0], vec![0.0, 1.0]]
    }
    #[test]
    fn test_replicator_step_sums_to_one() {
        let rd = ReplicatorDynamics::new(hawk_dove_payoff());
        let x = vec![0.6, 0.4];
        let x1 = rd.step_euler(&x, 0.01);
        let sum: f64 = x1.iter().sum();
        assert!((sum - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_replicator_converges() {
        let rd = ReplicatorDynamics::new(hawk_dove_payoff());
        let x = rd.run(&[0.5, 0.5], 0.01, 5000);
        assert!(
            (x[0] - 2.0 / 3.0).abs() < 0.05,
            "Hawks ESS proportion should be ~2/3, got {}",
            x[0]
        );
    }
    #[test]
    fn test_mean_fitness_dominance() {
        let rd = ReplicatorDynamics::new(vec![vec![2.0, 0.0], vec![0.0, 1.0]]);
        let x = vec![0.8, 0.2];
        let f = rd.fitness(&x);
        let mean = rd.mean_fitness(&x);
        assert!(f[0] > f[1]);
        assert!(mean > 0.0);
    }
    #[test]
    fn test_pure_ess_stag_hunt() {
        let rd = ReplicatorDynamics::new(vec![vec![3.0, 0.0], vec![2.0, 2.0]]);
        assert!(rd.is_pure_ess(0), "Stag is ESS");
        assert!(rd.is_pure_ess(1), "Hare is ESS");
    }
    #[test]
    fn test_normalise_simplex() {
        let x = vec![-0.1, 0.5, 0.3];
        let n = normalise_simplex(&x);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < 1e-9);
        assert!(n[0] < 1e-12, "negative component should be 0");
    }
}
#[cfg(test)]
mod tests_auction {

    use crate::game_theory::DirectMechanism;
    use crate::game_theory::EnglishAuction;
    use crate::game_theory::FirstPriceAuction;

    use crate::game_theory::VickreyAuction;

    #[test]
    fn test_first_price_winner() {
        let a = FirstPriceAuction::new();
        let r = a.run(&[0.3, 0.9, 0.5]).unwrap();
        assert_eq!(r.winner, 1);
        assert!((r.price - 0.9).abs() < 1e-9);
    }
    #[test]
    fn test_first_price_empty() {
        let a = FirstPriceAuction::new();
        assert!(a.run(&[]).is_none());
    }
    #[test]
    fn test_first_price_equilibrium_bid() {
        let b = FirstPriceAuction::equilibrium_bid(0.8, 2);
        assert!((b - 0.4).abs() < 1e-9);
    }
    #[test]
    fn test_vickrey_pays_second_price() {
        let v = VickreyAuction::new();
        let r = v.run(&[0.2, 0.9, 0.7]).unwrap();
        assert_eq!(r.winner, 1);
        assert!((r.price - 0.7).abs() < 1e-9);
    }
    #[test]
    fn test_vickrey_single_bidder() {
        let v = VickreyAuction::new();
        let r = v.run(&[0.5]).unwrap();
        assert_eq!(r.winner, 0);
        assert!((r.price - 0.0).abs() < 1e-9);
    }
    #[test]
    fn test_vickrey_truthful_dominant() {
        let v = VickreyAuction::new();
        assert!(v.truthful_dominates(0.8, &[0.3, 0.6]));
    }
    #[test]
    fn test_english_auction_outcome() {
        let e = EnglishAuction::new(vec![0.2, 0.9, 0.7]);
        let r = e.run().unwrap();
        assert_eq!(r.winner, 1);
        assert!((r.price - 0.7).abs() < 1e-9);
    }
    #[test]
    fn test_myerson_virtual_value() {
        assert!((DirectMechanism::myerson_virtual_value_uniform(0.5)).abs() < 1e-9);
        assert!((DirectMechanism::myerson_virtual_value_uniform(0.75) - 0.5).abs() < 1e-9);
    }
    #[test]
    fn test_optimal_reserve_price() {
        assert!((DirectMechanism::optimal_reserve_price_uniform() - 0.5).abs() < 1e-9);
    }
}
#[cfg(test)]
mod tests_cooperative {

    use crate::game_theory::CooperativeGame;

    fn three_player_majority_game() -> CooperativeGame {
        let mut g = CooperativeGame::new(3);
        g.set_worth(0b000, 0.0);
        g.set_worth(0b001, 0.0);
        g.set_worth(0b010, 0.0);
        g.set_worth(0b100, 0.0);
        g.set_worth(0b011, 1.0);
        g.set_worth(0b101, 1.0);
        g.set_worth(0b110, 1.0);
        g.set_worth(0b111, 1.0);
        g
    }
    #[test]
    fn test_shapley_value_symmetric() {
        let g = three_player_majority_game();
        let sv = g.shapley_values();
        assert!((sv[0] - 1.0 / 3.0).abs() < 1e-9);
        assert!((sv[1] - 1.0 / 3.0).abs() < 1e-9);
        assert!((sv[2] - 1.0 / 3.0).abs() < 1e-9);
    }
    #[test]
    fn test_shapley_efficiency() {
        let g = three_player_majority_game();
        let sv = g.shapley_values();
        let total: f64 = sv.iter().sum();
        assert!((total - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_core_majority_game_empty() {
        let g = three_player_majority_game();
        let split = vec![1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
        assert!(!g.is_in_core(&split));
    }
    #[test]
    fn test_banzhaf_index_symmetric() {
        let g = three_player_majority_game();
        let bz = g.banzhaf_index();
        let total: f64 = bz.iter().sum();
        assert!((total - 1.0).abs() < 1e-9);
        assert!((bz[0] - bz[1]).abs() < 1e-9);
        assert!((bz[1] - bz[2]).abs() < 1e-9);
    }
    #[test]
    fn test_worth_coalitions() {
        let g = three_player_majority_game();
        assert!((g.worth(0b001) - 0.0).abs() < 1e-9);
        assert!((g.worth(0b011) - 1.0).abs() < 1e-9);
        assert!((g.worth(0b111) - 1.0).abs() < 1e-9);
    }
}
#[cfg(test)]
mod tests_repeated {

    use crate::game_theory::GrimTrigger;

    use crate::game_theory::PdAction;
    use crate::game_theory::PrisonersDilemma;
    use crate::game_theory::TitForTat;

    use crate::game_theory::folk_theorem_cooperation_sustainable;

    use crate::game_theory::simulate_tft_vs_grim;
    fn pd() -> PrisonersDilemma {
        PrisonersDilemma::new(5.0, 3.0, 1.0, 0.0)
    }
    #[test]
    fn test_tft_cooperates_first() {
        let tft = TitForTat::new();
        assert_eq!(tft.choose(&[]), PdAction::Cooperate);
    }
    #[test]
    fn test_tft_copies_defect() {
        let tft = TitForTat::new();
        assert_eq!(tft.choose(&[PdAction::Defect]), PdAction::Defect);
    }
    #[test]
    fn test_tft_copies_cooperate() {
        let tft = TitForTat::new();
        assert_eq!(
            tft.choose(&[PdAction::Cooperate, PdAction::Cooperate]),
            PdAction::Cooperate
        );
    }
    #[test]
    fn test_grim_trigger_initial() {
        let gt = GrimTrigger::new();
        assert_eq!(gt.choose(), PdAction::Cooperate);
    }
    #[test]
    fn test_grim_trigger_after_defect() {
        let mut gt = GrimTrigger::new();
        gt.update(PdAction::Defect);
        assert_eq!(gt.choose(), PdAction::Defect);
    }
    #[test]
    fn test_grim_trigger_stays_triggered() {
        let mut gt = GrimTrigger::new();
        gt.update(PdAction::Defect);
        gt.update(PdAction::Cooperate);
        assert_eq!(gt.choose(), PdAction::Defect);
    }
    #[test]
    fn test_pd_payoffs() {
        let p = pd();
        assert!((p.payoff(PdAction::Cooperate, PdAction::Cooperate) - 3.0).abs() < 1e-9);
        assert!((p.payoff(PdAction::Defect, PdAction::Cooperate) - 5.0).abs() < 1e-9);
        assert!((p.payoff(PdAction::Cooperate, PdAction::Defect) - 0.0).abs() < 1e-9);
        assert!((p.payoff(PdAction::Defect, PdAction::Defect) - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_simulate_tft_vs_grim_mutual_cooperate() {
        let p = pd();
        let (pt, pg) = simulate_tft_vs_grim(&p, 10, 1.0);
        assert!((pt - 30.0).abs() < 1e-9);
        assert!((pg - 30.0).abs() < 1e-9);
    }
    #[test]
    fn test_folk_theorem_threshold() {
        let p = pd();
        assert!(!folk_theorem_cooperation_sustainable(&p, 0.4));
        assert!(folk_theorem_cooperation_sustainable(&p, 0.5));
        assert!(folk_theorem_cooperation_sustainable(&p, 0.9));
    }
}
#[cfg(test)]
mod tests_congestion {

    use crate::game_theory::CongestionGame;

    fn linear_congestion_game() -> CongestionGame {
        CongestionGame::new(vec![vec![0.0, 1.0], vec![2.0]])
    }
    #[test]
    fn test_cost_evaluation() {
        let g = linear_congestion_game();
        assert!((g.cost(0, 3.0) - 3.0).abs() < 1e-9);
        assert!((g.cost(1, 100.0) - 2.0).abs() < 1e-9);
    }
    #[test]
    fn test_social_cost() {
        let g = linear_congestion_game();
        let sc = g.social_cost(&[3.0, 1.0]);
        assert!((sc - 11.0).abs() < 1e-9);
    }
    #[test]
    fn test_poa_one_if_equal() {
        let g = linear_congestion_game();
        let loads = vec![2.0, 2.0];
        assert!((g.price_of_anarchy(&loads, &loads) - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_private_cost() {
        let g = linear_congestion_game();
        let loads = vec![3.0, 1.0];
        let pc = g.private_cost(&[0], &loads);
        assert!((pc - 3.0).abs() < 1e-9);
    }
    #[test]
    fn test_rosenthal_potential() {
        let g = linear_congestion_game();
        let pot = g.rosenthal_potential(&[2.0, 1.0]);
        assert!((pot - 5.0).abs() < 1e-9);
    }
}
#[cfg(test)]
mod tests_welfare {

    use crate::game_theory::SocialWelfareOptimiser;

    use crate::game_theory::WelfareCriterion;

    fn test_optimiser() -> SocialWelfareOptimiser {
        SocialWelfareOptimiser::new(
            vec![vec![4.0, 2.0], vec![1.0, 5.0], vec![3.0, 3.0]],
            vec![0.0, 0.0],
        )
    }
    #[test]
    fn test_utilitarian_welfare() {
        let opt = test_optimiser();
        assert!((opt.welfare(0, WelfareCriterion::Utilitarian) - 6.0).abs() < 1e-9);
        assert!((opt.welfare(2, WelfareCriterion::Utilitarian) - 6.0).abs() < 1e-9);
    }
    #[test]
    fn test_rawlsian_welfare() {
        let opt = test_optimiser();
        let idx = opt.optimal_outcome(WelfareCriterion::Rawlsian);
        assert_eq!(idx, 2);
    }
    #[test]
    fn test_nash_bargaining_welfare() {
        let opt = test_optimiser();
        let idx = opt.optimal_outcome(WelfareCriterion::NashBargaining);
        assert_eq!(idx, 2);
    }
    #[test]
    fn test_gini_coefficient_equal() {
        let opt = SocialWelfareOptimiser::new(vec![vec![2.0, 2.0]], vec![0.0, 0.0]);
        let g = opt.gini_coefficient(0);
        assert!(g.abs() < 1e-9, "equal utilities → Gini=0");
    }
}
#[cfg(test)]
mod tests_correlated_equilibrium {

    use crate::game_theory::CorrelatedEquilibrium;

    #[test]
    fn test_correlated_eq_pure_ne() {
        let payoff_a = vec![vec![3.0, 0.0], vec![5.0, 1.0]];
        let payoff_b = vec![vec![3.0, 5.0], vec![0.0, 1.0]];
        let dist = vec![0.0, 0.0, 0.0, 1.0];
        let ce = CorrelatedEquilibrium::new(dist, 2, 2);
        assert!(ce.verify(&payoff_a, &payoff_b));
    }
    #[test]
    fn test_correlated_eq_prob() {
        let dist = vec![0.25, 0.25, 0.25, 0.25];
        let ce = CorrelatedEquilibrium::new(dist, 2, 2);
        assert!((ce.prob(0, 0) - 0.25).abs() < 1e-9);
        assert!((ce.prob(1, 1) - 0.25).abs() < 1e-9);
    }
}
#[cfg(test)]
mod tests_bargaining {

    use crate::game_theory::NashBargaining;

    #[test]
    fn test_nash_bargaining_solution() {
        let nb = NashBargaining::new(
            vec![(0.0, 4.0), (2.0, 2.0), (4.0, 0.0), (1.0, 3.0), (3.0, 1.0)],
            (0.0, 0.0),
        );
        let sol = nb.solution().unwrap();
        assert_eq!(sol, (2.0, 2.0));
    }
    #[test]
    fn test_nash_bargaining_no_ir_solution() {
        let nb = NashBargaining::new(vec![(0.5, 0.5)], (1.0, 1.0));
        assert!(nb.solution().is_none());
    }
    #[test]
    fn test_kalai_smorodinsky() {
        let nb = NashBargaining::new(
            vec![(4.0, 0.0), (2.0, 2.0), (0.0, 4.0), (3.0, 1.0), (1.0, 3.0)],
            (0.0, 0.0),
        );
        let ks = nb.kalai_smorodinsky();
        assert!(ks.is_some());
    }
}
#[cfg(test)]
mod tests_regret {

    use crate::game_theory::RegretMatcher;

    #[test]
    fn test_regret_matcher_initial_uniform() {
        let rm = RegretMatcher::new(3);
        let s = rm.current_strategy();
        for &p in &s {
            assert!((p - 1.0 / 3.0).abs() < 1e-9);
        }
    }
    #[test]
    fn test_cfr_rps_converges_to_uniform() {
        let payoff = vec![
            vec![0.0, -1.0, 1.0],
            vec![1.0, 0.0, -1.0],
            vec![-1.0, 1.0, 0.0],
        ];
        let (sa, sb) = RegretMatcher::run_cfr_2p_zero_sum(&payoff, 10000);
        for &p in &sa {
            assert!((p - 1.0 / 3.0).abs() < 0.05, "sa={p}");
        }
        for &p in &sb {
            assert!((p - 1.0 / 3.0).abs() < 0.05, "sb={p}");
        }
    }
    #[test]
    fn test_regret_update_increases_winning_action() {
        let mut rm = RegretMatcher::new(2);
        let strategy = vec![0.5, 0.5];
        let utilities = vec![2.0, 0.0];
        rm.update_regret(&strategy, &utilities);
        assert!(rm.regret[0] > 0.0);
    }
}
#[cfg(test)]
mod tests_mechanism_design {

    use crate::game_theory::DirectMechanism;

    #[test]
    fn test_ir_check() {
        let m = DirectMechanism::new(2);
        let values = vec![0.8, 0.6];
        let alloc = vec![1.0, 0.0];
        let pay = vec![0.6, 0.0];
        assert!(m.is_individually_rational(&values, &alloc, &pay));
    }
    #[test]
    fn test_ir_violation() {
        let m = DirectMechanism::new(2);
        let values = vec![0.3, 0.6];
        let alloc = vec![1.0, 0.0];
        let pay = vec![0.5, 0.0];
        assert!(!m.is_individually_rational(&values, &alloc, &pay));
    }
    #[test]
    fn test_ic_check() {
        let m = DirectMechanism::new(2);
        assert!(m.is_incentive_compatible_agent(0, 0.8, 1.0, 0.6, 0.0, 0.0));
    }
}
#[cfg(test)]
mod tests_behavioural_strategy {

    use crate::game_theory::BehaviouralStrategy;

    #[test]
    fn test_behavioural_strategy_valid() {
        let mut s = BehaviouralStrategy::new();
        s.set(0, vec![0.3, 0.7]);
        s.set(1, vec![0.5, 0.5]);
        assert!(s.is_valid());
    }
    #[test]
    fn test_behavioural_strategy_invalid() {
        let mut s = BehaviouralStrategy::new();
        s.set(0, vec![0.4, 0.4]);
        assert!(!s.is_valid());
    }
    #[test]
    fn test_behavioural_strategy_prob() {
        let mut s = BehaviouralStrategy::new();
        s.set(0, vec![0.3, 0.7]);
        assert!((s.prob(0, 1) - 0.7).abs() < 1e-9);
        assert!((s.prob(99, 0) - 0.0).abs() < 1e-9);
    }
}
#[cfg(test)]
mod tests_efficiency {

    use crate::game_theory::efficiency_ratio;

    use crate::game_theory::price_of_anarchy_global;
    use crate::game_theory::price_of_stability;

    #[test]
    fn test_price_of_stability() {
        let pos = price_of_stability(&[3.0, 4.0, 5.0], 2.0);
        assert!((pos - 1.5).abs() < 1e-9);
    }
    #[test]
    fn test_price_of_anarchy_global() {
        let poa = price_of_anarchy_global(&[3.0, 4.0, 5.0], 2.0);
        assert!((poa - 2.5).abs() < 1e-9);
    }
    #[test]
    fn test_efficiency_ratio() {
        let r = efficiency_ratio(8.0, 10.0);
        assert!((r - 0.8).abs() < 1e-9);
    }
    #[test]
    fn test_efficiency_ratio_capped() {
        let r = efficiency_ratio(12.0, 10.0);
        assert!((r - 1.0).abs() < 1e-9);
    }
}
