//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use super::super::*;
    fn quad(x: &[f64]) -> f64 {
        (x[0] - 2.0).powi(2)
    }
    fn quad_grad(x: &[f64]) -> Vec<f64> {
        vec![2.0 * (x[0] - 2.0)]
    }
    fn bowl(x: &[f64]) -> f64 {
        x.iter()
            .enumerate()
            .map(|(i, xi)| (xi - i as f64).powi(2))
            .sum()
    }
    fn bowl_grad(x: &[f64]) -> Vec<f64> {
        x.iter()
            .enumerate()
            .map(|(i, xi)| 2.0 * (xi - i as f64))
            .collect()
    }
    fn rosenbrock(x: &[f64]) -> f64 {
        (1.0 - x[0]).powi(2) + 100.0 * (x[1] - x[0].powi(2)).powi(2)
    }
    fn rosenbrock_grad(x: &[f64]) -> Vec<f64> {
        let dx0 = -2.0 * (1.0 - x[0]) - 400.0 * x[0] * (x[1] - x[0].powi(2));
        let dx1 = 200.0 * (x[1] - x[0].powi(2));
        vec![dx0, dx1]
    }
    #[test]
    fn test_gd_converges_quad() {
        let res = gradient_descent(quad, quad_grad, vec![0.0], 0.1, 1e-6, 1000);
        assert!(res.converged);
        assert!((res.x[0] - 2.0).abs() < 1e-4);
    }
    #[test]
    fn test_gd_fval_at_min() {
        let res = gradient_descent(quad, quad_grad, vec![0.0], 0.1, 1e-8, 5000);
        assert!(res.f_val < 1e-8);
    }
    #[test]
    fn test_gd_bowl_2d() {
        let res = gradient_descent(bowl, bowl_grad, vec![5.0, 5.0], 0.1, 1e-6, 5000);
        assert!(res.converged);
        assert!((res.x[0] - 0.0).abs() < 1e-3);
        assert!((res.x[1] - 1.0).abs() < 1e-3);
    }
    #[test]
    fn test_gd_n_iter_positive() {
        let res = gradient_descent(quad, quad_grad, vec![0.0], 0.1, 1e-6, 1000);
        assert!(res.n_iter > 0);
    }
    #[test]
    fn test_gd_no_convergence_low_iter() {
        let res = gradient_descent(quad, quad_grad, vec![100.0], 0.01, 1e-12, 2);
        assert!(!res.converged);
        assert_eq!(res.n_iter, 2);
    }
    #[test]
    fn test_adam_converges_quad() {
        let res = adam(
            quad,
            quad_grad,
            vec![0.0],
            AdamConfig {
                lr: 0.1,
                beta1: 0.9,
                beta2: 0.999,
                eps: 1e-8,
                max_iter: 5000,
            },
        );
        assert!((res.x[0] - 2.0).abs() < 1e-3);
    }
    #[test]
    fn test_adam_bowl_2d() {
        let res = adam(
            bowl,
            bowl_grad,
            vec![5.0, 5.0],
            AdamConfig {
                lr: 0.1,
                beta1: 0.9,
                beta2: 0.999,
                eps: 1e-8,
                max_iter: 5000,
            },
        );
        assert!((res.x[0] - 0.0).abs() < 1e-2);
        assert!((res.x[1] - 1.0).abs() < 1e-2);
    }
    #[test]
    fn test_adam_fval_reasonable() {
        let res = adam(
            quad,
            quad_grad,
            vec![0.0],
            AdamConfig {
                lr: 0.1,
                beta1: 0.9,
                beta2: 0.999,
                eps: 1e-8,
                max_iter: 5000,
            },
        );
        assert!(res.f_val < 1e-4);
    }
    #[test]
    fn test_adam_n_iter_positive() {
        let res = adam(
            quad,
            quad_grad,
            vec![0.0],
            AdamConfig {
                lr: 0.1,
                beta1: 0.9,
                beta2: 0.999,
                eps: 1e-8,
                max_iter: 100,
            },
        );
        assert!(res.n_iter > 0 && res.n_iter <= 100);
    }
    #[test]
    fn test_lbfgs_quad() {
        let res = lbfgs(quad, quad_grad, vec![0.0], 5, 1e-8, 200);
        assert!(res.converged);
        assert!((res.x[0] - 2.0).abs() < 1e-6);
    }
    #[test]
    fn test_lbfgs_bowl() {
        let res = lbfgs(bowl, bowl_grad, vec![5.0, 5.0], 5, 1e-8, 500);
        assert!(res.converged);
        assert!((res.x[0] - 0.0).abs() < 1e-5);
        assert!((res.x[1] - 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_lbfgs_rosenbrock() {
        let res = lbfgs(rosenbrock, rosenbrock_grad, vec![-1.0, 1.0], 10, 1e-6, 2000);
        assert!((res.x[0] - 1.0).abs() < 1e-3);
        assert!((res.x[1] - 1.0).abs() < 1e-3);
    }
    #[test]
    fn test_lbfgs_fval() {
        let res = lbfgs(quad, quad_grad, vec![0.0], 5, 1e-8, 200);
        assert!(res.f_val < 1e-12);
    }
    #[test]
    fn test_nelder_mead_1d() {
        let res = nelder_mead(|x| (x[0] - 3.0).powi(2), vec![0.0], 1.0, 1e-8, 1000);
        assert!((res.x[0] - 3.0).abs() < 1e-4);
    }
    #[test]
    fn test_nelder_mead_2d() {
        let res = nelder_mead(bowl, vec![5.0, 5.0], 1.0, 1e-8, 10000);
        assert!((res.x[0] - 0.0).abs() < 1e-2, "x[0]={}", res.x[0]);
        assert!((res.x[1] - 1.0).abs() < 1e-2, "x[1]={}", res.x[1]);
    }
    #[test]
    fn test_nelder_mead_rosenbrock() {
        let res = nelder_mead(rosenbrock, vec![0.0, 0.0], 0.5, 1e-8, 10000);
        assert!(res.f_val < 1e-4);
    }
    #[test]
    fn test_nelder_mead_converged_flag() {
        let res = nelder_mead(|x| (x[0] - 1.0).powi(2), vec![0.0], 1.0, 1e-8, 2000);
        assert!(res.converged);
    }
    #[test]
    fn test_golden_section_parabola() {
        let (xm, fm) = golden_section(|x| (x - 3.0).powi(2), 0.0, 6.0, 1e-9);
        assert!((xm - 3.0).abs() < 1e-6);
        assert!(fm < 1e-10);
    }
    #[test]
    fn test_golden_section_negative_side() {
        let (xm, _) = golden_section(|x| (x + 1.0).powi(2), -5.0, 5.0, 1e-9);
        assert!((xm + 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_golden_section_returns_fmin() {
        let (_, fm) = golden_section(|x| x * x + 1.0, -2.0, 2.0, 1e-9);
        assert!((fm - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_brent_parabola() {
        let (xm, fm) = brent(|x| (x - 2.0).powi(2), 0.0, 2.0, 4.0, 1e-9);
        assert!((xm - 2.0).abs() < 1e-6);
        assert!(fm < 1e-10);
    }
    #[test]
    fn test_brent_cosine() {
        let (xm, fm) = brent(
            |x| x.cos(),
            0.0,
            std::f64::consts::PI,
            2.0 * std::f64::consts::PI,
            1e-9,
        );
        assert!((xm - std::f64::consts::PI).abs() < 1e-6);
        assert!((fm + 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_bisection_simple() {
        let root = bisection(|x| x - 1.5, 0.0, 3.0, 1e-9);
        assert!(root.is_some());
        assert!((root.unwrap() - 1.5).abs() < 1e-7);
    }
    #[test]
    fn test_bisection_sqrt2() {
        let root = bisection(|x| x * x - 2.0, 1.0, 2.0, 1e-10);
        assert!(root.is_some());
        assert!((root.unwrap() - 2.0_f64.sqrt()).abs() < 1e-8);
    }
    #[test]
    fn test_bisection_no_bracket() {
        let root = bisection(|x| x * x + 1.0, -2.0, 2.0, 1e-9);
        assert!(root.is_none());
    }
    #[test]
    fn test_numerical_gradient_quad() {
        let g = numerical_gradient(quad, &[3.0], 1e-5);
        assert!((g[0] - 2.0).abs() < 1e-6);
    }
    #[test]
    fn test_numerical_gradient_bowl() {
        let x = vec![2.0, 3.0];
        let g = numerical_gradient(bowl, &x, 1e-5);
        assert!((g[0] - 4.0).abs() < 1e-6);
        assert!((g[1] - 4.0).abs() < 1e-6);
    }
    #[test]
    fn test_numerical_hessian_quad() {
        let h = numerical_hessian(quad, &[2.0], 1e-4);
        assert!((h[0][0] - 2.0).abs() < 1e-4);
    }
    #[test]
    fn test_numerical_hessian_bowl_2d() {
        let h = numerical_hessian(bowl, &[0.0, 1.0], 1e-4);
        assert!((h[0][0] - 2.0).abs() < 1e-4);
        assert!((h[1][1] - 2.0).abs() < 1e-4);
        assert!(h[0][1].abs() < 1e-4);
    }
    #[test]
    fn test_backtracking_decreases_f() {
        let x = vec![0.0];
        let dir = vec![1.0];
        let f0 = quad(&x);
        let gd = -2.0 * 2.0;
        let alpha = backtracking_line_search(quad, &x, &dir, f0, gd, BacktrackingConfig::default());
        let x_new: Vec<f64> = x
            .iter()
            .zip(dir.iter())
            .map(|(xi, di)| xi + alpha * di)
            .collect();
        assert!(quad(&x_new) <= f0 + 1e-4 * alpha * gd);
    }
    #[test]
    fn test_backtracking_returns_positive_alpha() {
        let x = vec![5.0];
        let dir = vec![-1.0];
        let f0 = quad(&x);
        let gd = -6.0;
        let alpha = backtracking_line_search(quad, &x, &dir, f0, gd, BacktrackingConfig::default());
        assert!(alpha > 0.0);
    }
    #[test]
    fn test_constrained_unconstrained_minimum_inside() {
        let res = constrained_min_box(quad, quad_grad, vec![0.0], &[0.0], &[5.0], 1e-7, 1000);
        assert!((res.x[0] - 2.0).abs() < 1e-4);
    }
    #[test]
    fn test_constrained_minimum_at_boundary() {
        let f = |x: &[f64]| (x[0] - 5.0).powi(2);
        let g = |x: &[f64]| vec![2.0 * (x[0] - 5.0)];
        let res = constrained_min_box(f, g, vec![0.0], &[0.0], &[3.0], 1e-7, 1000);
        assert!((res.x[0] - 3.0).abs() < 1e-4);
    }
    #[test]
    fn test_momentum_converges_quad() {
        let res = gradient_descent_momentum(quad, quad_grad, vec![0.0], 0.05, 0.9, 1e-6, 5000);
        assert!(res.converged || res.f_val < 1e-6, "f_val={}", res.f_val);
        assert!((res.x[0] - 2.0).abs() < 1e-3, "x={}", res.x[0]);
    }
    #[test]
    fn test_momentum_bowl_2d() {
        let res = gradient_descent_momentum(bowl, bowl_grad, vec![5.0, 5.0], 0.05, 0.9, 1e-6, 5000);
        assert!((res.x[0] - 0.0).abs() < 1e-2, "x[0]={}", res.x[0]);
        assert!((res.x[1] - 1.0).abs() < 1e-2, "x[1]={}", res.x[1]);
    }
    #[test]
    fn test_momentum_n_iter_positive() {
        let res = gradient_descent_momentum(quad, quad_grad, vec![0.0], 0.05, 0.9, 1e-6, 100);
        assert!(res.n_iter > 0 && res.n_iter <= 100);
    }
    #[test]
    fn test_momentum_zero_at_minimum() {
        let res = gradient_descent_momentum(quad, quad_grad, vec![2.0], 0.05, 0.9, 1e-8, 100);
        assert!(res.f_val < 1e-8);
    }
    #[test]
    fn test_sa_quad_finds_minimum() {
        let res = simulated_annealing(quad, vec![0.0], 10.0, 0.95, 0.5, 10000);
        assert!((res.x[0] - 2.0).abs() < 0.5, "sa x[0]={}", res.x[0]);
    }
    #[test]
    fn test_sa_returns_finite() {
        let res = simulated_annealing(bowl, vec![5.0, 5.0], 5.0, 0.99, 0.5, 5000);
        assert!(res.f_val.is_finite());
        assert!(res.x.iter().all(|x| x.is_finite()));
    }
    #[test]
    fn test_sa_converges_flag_at_low_temp() {
        let res = simulated_annealing(quad, vec![2.0], 1.0, 0.5, 0.1, 1000);
        assert!(res.converged);
    }
    #[test]
    fn test_sa_f_val_matches_returned_x() {
        let res = simulated_annealing(quad, vec![5.0], 2.0, 0.9, 0.3, 2000);
        let expected = quad(&res.x);
        assert!((res.f_val - expected).abs() < 1e-12);
    }
    #[test]
    fn test_pso_quad() {
        let res = particle_swarm(
            quad,
            &[-10.0],
            &[10.0],
            PsoConfig {
                n_particles: 20,
                w: 0.7,
                c1: 2.0,
                c2: 2.0,
                max_iter: 200,
            },
        );
        assert!((res.x[0] - 2.0).abs() < 0.5, "pso x[0]={}", res.x[0]);
    }
    #[test]
    fn test_pso_bowl_2d() {
        let res = particle_swarm(
            bowl,
            &[-10.0, -10.0],
            &[10.0, 10.0],
            PsoConfig {
                n_particles: 30,
                w: 0.7,
                c1: 2.0,
                c2: 2.0,
                max_iter: 300,
            },
        );
        assert!(res.f_val < 2.0, "pso f_val={}", res.f_val);
    }
    #[test]
    fn test_pso_respects_bounds() {
        let res = particle_swarm(
            quad,
            &[0.0],
            &[1.0],
            PsoConfig {
                n_particles: 20,
                w: 0.7,
                c1: 2.0,
                c2: 2.0,
                max_iter: 200,
            },
        );
        assert!(
            res.x[0] >= 0.0 && res.x[0] <= 1.0,
            "out of bounds: x={}",
            res.x[0]
        );
    }
    #[test]
    fn test_pso_n_iter_equals_max() {
        let res = particle_swarm(
            quad,
            &[-5.0],
            &[5.0],
            PsoConfig {
                n_particles: 10,
                w: 0.7,
                c1: 2.0,
                c2: 2.0,
                max_iter: 50,
            },
        );
        assert_eq!(res.n_iter, 50);
    }
    #[test]
    fn test_ga_quad_finds_minimum() {
        let res = genetic_algorithm(
            quad,
            &[-10.0],
            &[10.0],
            GeneticAlgorithmConfig {
                pop_size: 50,
                crossover_rate: 0.8,
                mutation_rate: 0.1,
                mutation_scale: 0.5,
                max_generations: 200,
            },
        );
        assert!((res.x[0] - 2.0).abs() < 1.0, "ga x[0]={}", res.x[0]);
    }
    #[test]
    fn test_ga_bowl_2d() {
        let res = genetic_algorithm(
            bowl,
            &[-5.0, -5.0],
            &[5.0, 5.0],
            GeneticAlgorithmConfig {
                pop_size: 80,
                crossover_rate: 0.8,
                mutation_rate: 0.2,
                mutation_scale: 0.5,
                max_generations: 500,
            },
        );
        assert!(res.f_val < 3.0, "ga f_val={}", res.f_val);
    }
    #[test]
    fn test_ga_respects_bounds() {
        let res = genetic_algorithm(
            quad,
            &[0.0],
            &[1.5],
            GeneticAlgorithmConfig {
                pop_size: 30,
                crossover_rate: 0.8,
                mutation_rate: 0.1,
                mutation_scale: 0.2,
                max_generations: 100,
            },
        );
        assert!(
            res.x[0] >= 0.0 && res.x[0] <= 1.5,
            "out of bounds: x={}",
            res.x[0]
        );
    }
    #[test]
    fn test_ga_n_iter_equals_max_generations() {
        let res = genetic_algorithm(
            quad,
            &[-5.0],
            &[5.0],
            GeneticAlgorithmConfig {
                pop_size: 10,
                crossover_rate: 0.8,
                mutation_rate: 0.1,
                mutation_scale: 0.5,
                max_generations: 20,
            },
        );
        assert_eq!(res.n_iter, 20);
    }
    #[test]
    fn test_lbfgsb_quad_inside_bounds() {
        let res = lbfgsb(
            quad,
            quad_grad,
            vec![0.0],
            &[0.0],
            &[5.0],
            LbfgsbConfig {
                m: 5,
                tol: 1e-8,
                max_iter: 200,
            },
        );
        assert!(res.converged);
        assert!((res.x[0] - 2.0).abs() < 1e-5, "x={}", res.x[0]);
    }
    #[test]
    fn test_lbfgsb_quad_min_outside_bounds() {
        let f = |x: &[f64]| (x[0] - 5.0).powi(2);
        let g = |x: &[f64]| vec![2.0 * (x[0] - 5.0)];
        let res = lbfgsb(
            f,
            g,
            vec![0.0],
            &[-3.0],
            &[3.0],
            LbfgsbConfig {
                m: 5,
                tol: 1e-8,
                max_iter: 200,
            },
        );
        assert!((res.x[0] - 3.0).abs() < 1e-4, "x={}", res.x[0]);
    }
    #[test]
    fn test_lbfgsb_bowl_2d() {
        let res = lbfgsb(
            bowl,
            bowl_grad,
            vec![4.0, 4.0],
            &[-10.0, -10.0],
            &[10.0, 10.0],
            LbfgsbConfig {
                m: 5,
                tol: 1e-8,
                max_iter: 500,
            },
        );
        assert!(res.converged, "lbfgsb should converge on bowl");
        assert!((res.x[0] - 0.0).abs() < 1e-4, "x[0]={}", res.x[0]);
        assert!((res.x[1] - 1.0).abs() < 1e-4, "x[1]={}", res.x[1]);
    }
    #[test]
    fn test_lbfgsb_rosenbrock() {
        let res = lbfgsb(
            rosenbrock,
            rosenbrock_grad,
            vec![-1.0, 1.0],
            &[-5.0, -5.0],
            &[5.0, 5.0],
            LbfgsbConfig {
                m: 10,
                tol: 1e-6,
                max_iter: 2000,
            },
        );
        assert!((res.x[0] - 1.0).abs() < 0.01, "x[0]={}", res.x[0]);
        assert!((res.x[1] - 1.0).abs() < 0.01, "x[1]={}", res.x[1]);
    }
    #[test]
    fn test_lbfgsb_respects_bounds() {
        let res = lbfgsb(
            quad,
            quad_grad,
            vec![0.0],
            &[0.0],
            &[1.0],
            LbfgsbConfig {
                m: 5,
                tol: 1e-8,
                max_iter: 200,
            },
        );
        assert!(
            res.x[0] >= 0.0 && res.x[0] <= 1.0 + 1e-10,
            "out of bounds: x={}",
            res.x[0]
        );
    }
    #[test]
    fn test_nelder_mead_struct_minimizes_quadratic() {
        let mut nm = NelderMead::new(vec![0.0], 1.0);
        let best = nm.minimize(|x| (x[0] - 3.0).powi(2), 1000, 1e-8);
        assert!(
            (best[0] - 3.0).abs() < 1e-4,
            "NelderMead struct: x={}",
            best[0]
        );
    }
    #[test]
    fn test_nelder_mead_struct_2d() {
        let mut nm = NelderMead::new(vec![5.0, 5.0], 1.0);
        let best = nm.minimize(bowl, 5000, 1e-8);
        assert!((best[0] - 0.0).abs() < 1e-2, "x[0]={}", best[0]);
        assert!((best[1] - 1.0).abs() < 1e-2, "x[1]={}", best[1]);
    }
    #[test]
    fn test_nelder_mead_struct_simplex_size() {
        let nm = NelderMead::new(vec![1.0, 2.0], 0.5);
        assert_eq!(nm.simplex.len(), 3);
        assert_eq!(nm.n, 2);
    }
    #[test]
    fn test_de_evolves_toward_minimum() {
        let mut rng = rand::rng();
        let bounds = vec![(-5.0_f64, 5.0_f64)];
        let mut de = DifferentialEvolution::new(20, 1, &bounds, &mut rng);
        for _ in 0..100 {
            de.step(quad, &mut rng);
        }
        let best = de.best(quad);
        assert!((best[0] - 2.0).abs() < 1.0, "DE best: x={}", best[0]);
    }
    #[test]
    fn test_de_population_size() {
        let mut rng = rand::rng();
        let bounds = vec![(-5.0, 5.0), (-5.0, 5.0)];
        let de = DifferentialEvolution::new(15, 2, &bounds, &mut rng);
        assert_eq!(de.pop.len(), 15);
        assert_eq!(de.pop[0].len(), 2);
    }
    #[test]
    fn test_de_fitness_improves() {
        let mut rng = rand::rng();
        let bounds = vec![(-5.0_f64, 5.0_f64)];
        let mut de = DifferentialEvolution::new(20, 1, &bounds, &mut rng);
        let init_best = de.best(quad);
        let f_init = quad(&init_best);
        for _ in 0..50 {
            de.step(quad, &mut rng);
        }
        let final_best = de.best(quad);
        let f_final = quad(&final_best);
        assert!(
            f_final <= f_init + 1e-10,
            "DE should not worsen: f_init={f_init}, f_final={f_final}"
        );
    }
    #[test]
    fn test_trust_region_reduces_residual() {
        let mut tr = TrustRegion::new(1.0, 0.1);
        let mut x = vec![5.0_f64];
        let f_init = quad(&x);
        for _ in 0..50 {
            x = tr.step(quad, quad_grad, &x);
        }
        let f_final = quad(&x);
        assert!(
            f_final < f_init,
            "TrustRegion: f did not decrease; f_init={f_init}, f_final={f_final}"
        );
    }
    #[test]
    fn test_trust_region_converges_quad() {
        let mut tr = TrustRegion::new(1.0, 0.1);
        let mut x = vec![0.0_f64];
        for _ in 0..200 {
            x = tr.step(quad, quad_grad, &x);
        }
        assert!((x[0] - 2.0).abs() < 0.5, "TrustRegion: x={}", x[0]);
    }
    #[test]
    fn test_cg_minimizes_quadratic() {
        let x = conjugate_gradient_minimize(quad, quad_grad, vec![0.0], 500, 1e-8);
        assert!((x[0] - 2.0).abs() < 1e-4, "CG: x={}", x[0]);
    }
    #[test]
    fn test_cg_minimizes_bowl_2d() {
        let x = conjugate_gradient_minimize(bowl, bowl_grad, vec![5.0, 5.0], 2000, 1e-8);
        assert!((x[0] - 0.0).abs() < 1e-3, "CG bowl: x[0]={}", x[0]);
        assert!((x[1] - 1.0).abs() < 1e-3, "CG bowl: x[1]={}", x[1]);
    }
    #[test]
    fn test_cg_already_at_minimum() {
        let x = conjugate_gradient_minimize(quad, quad_grad, vec![2.0], 100, 1e-12);
        assert!((x[0] - 2.0).abs() < 1e-10, "CG at minimum: x={}", x[0]);
    }
    #[test]
    fn test_wolfe_returns_positive_alpha() {
        let x = vec![5.0_f64];
        let d = vec![-1.0_f64];
        let f0 = quad(&x);
        let g0d: f64 = quad_grad(&x)
            .iter()
            .zip(d.iter())
            .map(|(g, di)| g * di)
            .sum();
        let alpha = wolfe_line_search(quad, quad_grad, &x, &d, f0, g0d, 1e-4, 0.9, 50);
        assert!(alpha > 0.0, "alpha={alpha}");
    }
    #[test]
    fn test_wolfe_armijo_satisfied() {
        let x = vec![5.0_f64];
        let d = vec![-1.0_f64];
        let f0 = quad(&x);
        let g0d: f64 = quad_grad(&x)
            .iter()
            .zip(d.iter())
            .map(|(g, di)| g * di)
            .sum();
        let c1 = 1e-4;
        let alpha = wolfe_line_search(quad, quad_grad, &x, &d, f0, g0d, c1, 0.9, 50);
        let x_new = vec![x[0] + alpha * d[0]];
        assert!(quad(&x_new) <= f0 + c1 * alpha * g0d, "Armijo violated");
    }
    #[test]
    fn test_wolfe_bowl_2d() {
        let x = vec![5.0, 5.0];
        let g = bowl_grad(&x);
        let d: Vec<f64> = g.iter().map(|gi| -*gi).collect();
        let f0 = bowl(&x);
        let g0d: f64 = g.iter().zip(d.iter()).map(|(gi, di)| gi * di).sum();
        let alpha = wolfe_line_search(bowl, bowl_grad, &x, &d, f0, g0d, 1e-4, 0.9, 50);
        let x_new: Vec<f64> = x
            .iter()
            .zip(d.iter())
            .map(|(xi, di)| xi + alpha * di)
            .collect();
        assert!(bowl(&x_new) < f0, "Wolfe step should decrease f");
    }
    #[test]
    fn test_aug_lagrangian_equality_constraint() {
        let f = |x: &[f64]| x[0] * x[0];
        let gf = |x: &[f64]| vec![2.0 * x[0]];
        let c = |x: &[f64]| x[0] - 3.0;
        let dc = |_x: &[f64]| vec![1.0_f64];
        let res = augmented_lagrangian(f, gf, &[c], &[dc], vec![0.0], 1.0, 20, 50, 1e-6);
        assert!((res.x[0] - 3.0).abs() < 0.1, "aug_lag x={}", res.x[0]);
    }
    #[test]
    fn test_aug_lagrangian_unconstrained_reduces_to_min() {
        let f = |x: &[f64]| x[0] * x[0];
        let gf = |x: &[f64]| vec![2.0 * x[0]];
        let constraints: Vec<fn(&[f64]) -> f64> = vec![];
        let grad_constraints = Vec::<fn(&[f64]) -> Vec<f64>>::new();
        let res = augmented_lagrangian(
            f,
            gf,
            &constraints,
            &grad_constraints,
            vec![5.0],
            1.0,
            20,
            100,
            1e-6,
        );
        assert!(res.x[0].abs() < 1.0, "x={}", res.x[0]);
    }
    #[test]
    fn test_coordinate_descent_quadratic() {
        let res = coordinate_descent(quad, vec![5.0], 8.0, 1e-6, 500);
        assert!((res.x[0] - 2.0).abs() < 0.1, "coord_desc x={}", res.x[0]);
    }
    #[test]
    fn test_coordinate_descent_bowl_2d() {
        let res = coordinate_descent(bowl, vec![5.0, 5.0], 10.0, 1e-5, 1000);
        assert!(res.f_val < 1.0, "coord_desc f_val={}", res.f_val);
    }
    #[test]
    fn test_coordinate_descent_returns_finite() {
        let res = coordinate_descent(
            |x| (x[0] - 1.0).powi(2) + (x[1] + 2.0).powi(2),
            vec![0.0, 0.0],
            5.0,
            1e-4,
            200,
        );
        assert!(res.x.iter().all(|v| v.is_finite()));
    }
    #[test]
    fn test_powell_quadratic() {
        let res = powell(quad, vec![5.0], 1e-6, 100);
        assert!((res.x[0] - 2.0).abs() < 0.5, "powell x={}", res.x[0]);
    }
    #[test]
    fn test_powell_bowl_2d() {
        let res = powell(bowl, vec![5.0, 5.0], 1e-5, 500);
        assert!(res.f_val < 1.0, "powell f_val={}", res.f_val);
    }
    #[test]
    fn test_powell_already_at_min() {
        let res = powell(quad, vec![2.0], 1e-8, 50);
        assert!(res.f_val < 1e-8, "powell at min: f_val={}", res.f_val);
    }
    #[test]
    fn test_sgd_cosine_reduces_f() {
        let f_init = quad(&[5.0]);
        let res = sgd_cosine_annealing(quad, quad_grad, vec![5.0], 0.1, 1e-4, 200);
        assert!(
            res.f_val < f_init,
            "SGD cosine should reduce f: f_init={f_init} f_final={}",
            res.f_val
        );
    }
    #[test]
    fn test_sgd_cosine_n_iter() {
        let res = sgd_cosine_annealing(quad, quad_grad, vec![0.0], 0.1, 1e-4, 100);
        assert_eq!(res.n_iter, 100);
    }
    #[test]
    fn test_sgd_cosine_finite() {
        let res = sgd_cosine_annealing(bowl, bowl_grad, vec![5.0, 5.0], 0.05, 1e-5, 500);
        assert!(res.x.iter().all(|v| v.is_finite()));
    }
    #[test]
    fn test_bfgs_quadratic() {
        let res = bfgs(quad, quad_grad, vec![0.0], 1e-8, 200);
        assert!(res.converged, "BFGS should converge on quadratic");
        assert!((res.x[0] - 2.0).abs() < 1e-5, "x={}", res.x[0]);
    }
    #[test]
    fn test_bfgs_bowl_2d() {
        let res = bfgs(bowl, bowl_grad, vec![5.0, 5.0], 1e-8, 500);
        assert!(res.converged, "BFGS should converge on bowl");
        assert!((res.x[0] - 0.0).abs() < 1e-4, "x[0]={}", res.x[0]);
        assert!((res.x[1] - 1.0).abs() < 1e-4, "x[1]={}", res.x[1]);
    }
    #[test]
    fn test_bfgs_rosenbrock() {
        let res = bfgs(rosenbrock, rosenbrock_grad, vec![-1.0, 1.0], 1e-6, 2000);
        assert!((res.x[0] - 1.0).abs() < 0.01, "x[0]={}", res.x[0]);
        assert!((res.x[1] - 1.0).abs() < 0.01, "x[1]={}", res.x[1]);
    }
    #[test]
    fn test_bfgs_f_val_at_min() {
        let res = bfgs(quad, quad_grad, vec![0.0], 1e-8, 200);
        assert!(res.f_val < 1e-12, "f_val={}", res.f_val);
    }
    #[test]
    fn test_clip_gradient_clips_large() {
        let mut g = vec![3.0_f64, 4.0_f64];
        clip_gradient(&mut g, 1.0);
        let norm: f64 = g.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!((norm - 1.0).abs() < 1e-10, "norm={norm}");
    }
    #[test]
    fn test_clip_gradient_no_op_when_small() {
        let mut g = vec![0.1_f64, 0.1_f64];
        let g_orig = g.clone();
        clip_gradient(&mut g, 1.0);
        for (a, b) in g.iter().zip(g_orig.iter()) {
            assert!((a - b).abs() < 1e-15);
        }
    }
    #[test]
    fn test_gradient_descent_clipped_converges() {
        let res = gradient_descent_clipped(quad, quad_grad, vec![0.0], 0.1, 5.0, 1e-6, 1000);
        assert!(res.converged || res.f_val < 1e-6);
        assert!((res.x[0] - 2.0).abs() < 1e-3, "clipped GD: x={}", res.x[0]);
    }
    #[test]
    fn test_gradient_descent_clipped_respects_clip() {
        let res = gradient_descent_clipped(quad, quad_grad, vec![1000.0], 0.01, 1.0, 1e-4, 5000);
        assert!(res.x[0].is_finite(), "clipped GD: x should be finite");
    }
    #[test]
    fn test_nelder_mead_already_at_minimum() {
        let res = nelder_mead(quad, vec![2.0], 0.1, 1e-10, 100);
        assert!(res.converged, "NM at minimum should converge");
    }
    #[test]
    fn test_de_bowl_2d() {
        let mut rng = rand::rng();
        let bounds = vec![(-5.0_f64, 5.0_f64), (-5.0_f64, 5.0_f64)];
        let mut de = DifferentialEvolution::new(20, 2, &bounds, &mut rng);
        for _ in 0..200 {
            de.step(bowl, &mut rng);
        }
        let best = de.best(bowl);
        assert!(bowl(&best) < 1.0, "DE bowl: f_val={}", bowl(&best));
    }
    #[test]
    fn test_de_f_cr_defaults() {
        let mut rng = rand::rng();
        let bounds = vec![(-1.0_f64, 1.0_f64)];
        let de = DifferentialEvolution::new(10, 1, &bounds, &mut rng);
        assert!((de.f_weight - 0.8).abs() < 1e-12);
        assert!((de.cr - 0.9).abs() < 1e-12);
    }
    #[test]
    fn test_de_population_in_bounds() {
        let mut rng = rand::rng();
        let bounds = vec![(-3.0_f64, 3.0_f64), (0.0_f64, 5.0_f64)];
        let de = DifferentialEvolution::new(15, 2, &bounds, &mut rng);
        for individual in &de.pop {
            assert!(individual[0] >= -3.0 && individual[0] <= 3.0);
            assert!(individual[1] >= 0.0 && individual[1] <= 5.0);
        }
    }
    #[test]
    fn test_pso_rosenbrock() {
        let res = particle_swarm(
            rosenbrock,
            &[-2.0, -2.0],
            &[2.0, 2.0],
            PsoConfig {
                n_particles: 30,
                w: 0.7,
                c1: 2.0,
                c2: 2.0,
                max_iter: 500,
            },
        );
        assert!(res.f_val < 2.0, "PSO rosenbrock f_val={}", res.f_val);
    }
    #[test]
    fn test_pso_1d_minimum() {
        let res = particle_swarm(
            |x: &[f64]| (x[0] - 3.0).powi(2),
            &[0.0],
            &[6.0],
            PsoConfig {
                n_particles: 20,
                w: 0.7,
                c1: 2.0,
                c2: 2.0,
                max_iter: 300,
            },
        );
        assert!((res.x[0] - 3.0).abs() < 0.5, "PSO 1D: x={}", res.x[0]);
    }
    #[test]
    fn test_cg_minimize_rosenbrock() {
        let x =
            conjugate_gradient_minimize(rosenbrock, rosenbrock_grad, vec![-1.0, 1.0], 5000, 1e-6);
        assert!((x[0] - 1.0).abs() < 0.05, "CG rosenbrock x[0]={}", x[0]);
        assert!((x[1] - 1.0).abs() < 0.05, "CG rosenbrock x[1]={}", x[1]);
    }
    #[test]
    fn test_lbfgs_converges_3d() {
        let f = |x: &[f64]| x.iter().map(|xi| xi * xi).sum::<f64>();
        let g = |x: &[f64]| x.iter().map(|xi| 2.0 * xi).collect::<Vec<_>>();
        let res = lbfgs(f, g, vec![3.0, 3.0, 3.0], 5, 1e-10, 200);
        assert!(res.converged);
        for xi in &res.x {
            assert!(xi.abs() < 1e-5, "x should be near 0, got {xi}");
        }
    }
    #[test]
    fn test_trust_region_at_minimum_stable() {
        let mut tr = TrustRegion::new(1.0, 0.1);
        let x0 = vec![2.0];
        let x1 = tr.step(quad, quad_grad, &x0);
        assert!(
            (x1[0] - 2.0).abs() < 1e-6,
            "TR at minimum should stay: x={}",
            x1[0]
        );
    }
    #[test]
    fn test_golden_section_narrow_bracket() {
        let (x, _) = golden_section(|x| x * x, -0.1, 0.1, 1e-12);
        assert!(x.abs() < 1e-8, "golden section near 0: x={x}");
    }
    #[test]
    fn test_golden_section_sine_minimum() {
        use std::f64::consts::PI;
        let (x, fx) = golden_section(|x| x.sin(), -PI, 0.0, 1e-8);
        assert!((x - (-PI / 2.0)).abs() < 1e-4, "golden section sin: x={x}");
        assert!((fx - (-1.0)).abs() < 1e-6, "golden section sin f(x)={fx}");
    }
    #[test]
    fn test_bisection_root_at_boundary() {
        let root = bisection(|x| x - 1.0, 0.0, 2.0, 1e-10).expect("root at x=1");
        assert!((root - 1.0).abs() < 1e-8, "root={root}");
    }
    #[test]
    fn test_numerical_gradient_rosenbrock() {
        let x = vec![1.0, 1.0];
        let g = numerical_gradient(rosenbrock, &x, 1e-5);
        let g_exact = rosenbrock_grad(&x);
        for i in 0..2 {
            assert!(
                (g[i] - g_exact[i]).abs() < 1e-4,
                "numerical gradient at {i}: {} vs {}",
                g[i],
                g_exact[i]
            );
        }
    }
    #[test]
    fn test_optimizer_de_finds_quadratic_minimum() {
        let opt = Optimizer::new();
        let bounds = vec![(-5.0_f64, 5.0_f64)];
        let result =
            opt.differential_evolution(|x| (x[0] - 3.0).powi(2), &bounds, 20, 0.8, 0.9, 200);
        assert!(
            (result.x[0] - 3.0).abs() < 0.2,
            "DE should find min near x=3, got {}",
            result.x[0]
        );
    }
    #[test]
    fn test_optimizer_de_2d_bowl() {
        let opt = Optimizer::new();
        let bounds = vec![(-5.0_f64, 5.0_f64), (-5.0_f64, 5.0_f64)];
        let result =
            opt.differential_evolution(|x| x[0] * x[0] + x[1] * x[1], &bounds, 30, 0.8, 0.9, 300);
        assert!(result.f_val < 0.5, "DE 2D bowl: f_val={}", result.f_val);
    }
    #[test]
    fn test_optimizer_de_returns_opt_result() {
        let opt = Optimizer::new();
        let bounds = vec![(0.0_f64, 1.0_f64)];
        let result = opt.differential_evolution(|x| x[0], &bounds, 10, 0.5, 0.7, 50);
        assert!(
            result.f_val.is_finite(),
            "DE result f_val should be finite: {}",
            result.f_val
        );
        assert_eq!(result.n_iter, 50);
    }
    #[test]
    fn test_optimizer_sa_finds_minimum() {
        let opt = Optimizer::new();
        let result = opt.simulated_annealing(
            |x| (x[0] - 2.0).powi(2),
            vec![0.0_f64],
            10.0,
            0.99,
            0.5,
            5000,
        );
        assert!(
            (result.x[0] - 2.0).abs() < 0.5,
            "SA should find min near x=2, got {}",
            result.x[0]
        );
    }
    #[test]
    fn test_optimizer_sa_2d_sphere() {
        let opt = Optimizer::new();
        let result = opt.simulated_annealing(
            |x| x[0] * x[0] + x[1] * x[1],
            vec![5.0_f64, 5.0_f64],
            20.0,
            0.995,
            0.3,
            5000,
        );
        assert!(result.f_val < 5.0, "SA 2D sphere: f_val={}", result.f_val);
    }
    #[test]
    fn test_cmaes_step_reduces_mean_distance_to_optimum() {
        let mut cmaes = CmaEs::new(vec![5.0_f64], 1.0);
        let f = |x: &[f64]| (x[0] - 0.0).powi(2);
        for _ in 0..30 {
            cmaes.step(&f);
        }
        assert!(
            cmaes.mean[0].abs() < 3.0,
            "CMA-ES mean should converge toward 0, got {}",
            cmaes.mean[0]
        );
    }
    #[test]
    fn test_cmaes_generation_counter_increments() {
        let mut cmaes = CmaEs::new(vec![1.0_f64, 1.0_f64], 0.5);
        let f = |x: &[f64]| x[0] * x[0] + x[1] * x[1];
        assert_eq!(cmaes.generation, 0);
        cmaes.step(&f);
        assert_eq!(cmaes.generation, 1);
        cmaes.step(&f);
        assert_eq!(cmaes.generation, 2);
    }
    #[test]
    fn test_cmaes_sigma_stays_positive() {
        let mut cmaes = CmaEs::new(vec![3.0_f64, -3.0_f64], 1.0);
        let f = |x: &[f64]| x.iter().map(|v| v * v).sum::<f64>();
        for _ in 0..10 {
            cmaes.step(&f);
            assert!(
                cmaes.sigma > 0.0,
                "sigma must remain positive: {}",
                cmaes.sigma
            );
        }
    }
    #[test]
    fn test_optimizer_cmaes_step_via_facade() {
        let opt = Optimizer::new();
        let mut cmaes = CmaEs::new(vec![4.0_f64], 2.0);
        let f = |x: &[f64]| (x[0] - 1.0).powi(2);
        let best_f = opt.cmaes_step(&mut cmaes, &f);
        assert!(best_f.is_finite(), "cmaes_step should return finite f_val");
    }
    #[test]
    fn test_cmaes_covariance_remains_symmetric() {
        let mut cmaes = CmaEs::new(vec![2.0_f64, -1.0_f64], 0.5);
        let f = |x: &[f64]| x[0] * x[0] + 2.0 * x[1] * x[1];
        let n = 2;
        for _ in 0..5 {
            cmaes.step(&f);
        }
        for i in 0..n {
            for j in 0..n {
                let diff = (cmaes.cov[i * n + j] - cmaes.cov[j * n + i]).abs();
                assert!(diff < 1e-10, "cov[{i},{j}] != cov[{j},{i}]: diff={diff}");
            }
        }
    }
}
