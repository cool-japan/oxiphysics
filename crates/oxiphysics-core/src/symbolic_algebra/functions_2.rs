//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[cfg(test)]
mod tests {
    use super::*;
    use crate::symbolic_algebra::Expr;
    use crate::symbolic_algebra::PartialFractionTerm;
    use crate::symbolic_algebra::Polynomial;
    use crate::symbolic_algebra::RationalFunction;
    use std::collections::HashMap;
    #[test]
    fn test_const_eval() {
        let e = cst(42.0);
        let env = HashMap::new();
        assert!((eval_expr(&e, &env).unwrap() - 42.0).abs() < 1e-15);
    }
    #[test]
    fn test_var_eval() {
        let e = var("x");
        let mut env = HashMap::new();
        env.insert("x".into(), 7.0);
        assert!((eval_expr(&e, &env).unwrap() - 7.0).abs() < 1e-15);
    }
    #[test]
    fn test_add_eval() {
        let e = cst(3.0).add_expr(cst(4.0));
        let env = HashMap::new();
        assert!((eval_expr(&e, &env).unwrap() - 7.0).abs() < 1e-15);
    }
    #[test]
    fn test_mul_eval() {
        let e = cst(3.0).mul_expr(cst(5.0));
        let env = HashMap::new();
        assert!((eval_expr(&e, &env).unwrap() - 15.0).abs() < 1e-15);
    }
    #[test]
    fn test_div_eval() {
        let e = cst(10.0).div_expr(cst(4.0));
        let env = HashMap::new();
        assert!((eval_expr(&e, &env).unwrap() - 2.5).abs() < 1e-15);
    }
    #[test]
    fn test_div_by_zero() {
        let e = cst(1.0).div_expr(cst(0.0));
        let env = HashMap::new();
        assert!(eval_expr(&e, &env).is_none());
    }
    #[test]
    fn test_pow_eval() {
        let e = cst(2.0).pow_expr(cst(10.0));
        let env = HashMap::new();
        assert!((eval_expr(&e, &env).unwrap() - 1024.0).abs() < 1e-10);
    }
    #[test]
    fn test_diff_const() {
        let d = diff(&cst(5.0), "x");
        let env = HashMap::new();
        assert!((eval_expr(&simplify(&d), &env).unwrap()).abs() < 1e-15);
    }
    #[test]
    fn test_diff_var() {
        let d = diff(&var("x"), "x");
        let env = HashMap::new();
        assert!((eval_expr(&simplify(&d), &env).unwrap() - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_diff_x_squared() {
        let e = var("x").pow_expr(cst(2.0));
        let d = simplify(&diff(&e, "x"));
        let mut env = HashMap::new();
        env.insert("x".into(), 3.0);
        assert!((eval_expr(&d, &env).unwrap() - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_diff_product_rule() {
        let e = var("x").mul_expr(var("x").sin_expr());
        let d = simplify(&diff(&e, "x"));
        let mut env = HashMap::new();
        env.insert("x".into(), 1.0);
        let expected = 1.0_f64.sin() + 1.0_f64.cos();
        assert!((eval_expr(&d, &env).unwrap() - expected).abs() < 1e-8);
    }
    #[test]
    fn test_simplify_zero_add() {
        let e = cst(0.0).add_expr(var("x"));
        let s = simplify(&e);
        assert_eq!(s, var("x"));
    }
    #[test]
    fn test_simplify_one_mul() {
        let e = cst(1.0).mul_expr(var("x"));
        let s = simplify(&e);
        assert_eq!(s, var("x"));
    }
    #[test]
    fn test_simplify_zero_mul() {
        let e = cst(0.0).mul_expr(var("x"));
        let s = simplify(&e);
        assert!(s.is_zero());
    }
    #[test]
    fn test_simplify_pow_one() {
        let e = var("x").pow_expr(cst(1.0));
        let s = simplify(&e);
        assert_eq!(s, var("x"));
    }
    #[test]
    fn test_simplify_exp_ln() {
        let e = var("x").ln_expr().exp_expr();
        let s = simplify(&e);
        assert_eq!(s, var("x"));
    }
    #[test]
    fn test_polynomial_eval() {
        let p = Polynomial::new(vec![2.0, 3.0, 1.0]);
        assert!((p.eval(2.0) - 12.0).abs() < 1e-12);
    }
    #[test]
    fn test_polynomial_add() {
        let a = Polynomial::new(vec![1.0, 2.0]);
        let b = Polynomial::new(vec![3.0, 4.0, 5.0]);
        let c = a.add(&b);
        assert!((c.eval(1.0) - 15.0).abs() < 1e-12);
    }
    #[test]
    fn test_polynomial_mul() {
        let a = Polynomial::new(vec![1.0, 1.0]);
        let b = Polynomial::new(vec![1.0, 1.0]);
        let c = a.mul(&b);
        assert!((c.eval(2.0) - 9.0).abs() < 1e-12);
    }
    #[test]
    fn test_polynomial_derivative() {
        let p = Polynomial::new(vec![1.0, 0.0, 3.0]);
        let dp = p.derivative();
        assert!((dp.eval(2.0) - 12.0).abs() < 1e-12);
    }
    #[test]
    fn test_polynomial_integral() {
        let p = Polynomial::new(vec![6.0]);
        let ip = p.integral();
        assert!((ip.eval(3.0) - 18.0).abs() < 1e-12);
    }
    #[test]
    fn test_polynomial_div_rem() {
        let a = Polynomial::new(vec![-1.0, 0.0, 1.0]);
        let b = Polynomial::new(vec![-1.0, 1.0]);
        let (q, r) = a.div_rem(&b);
        assert!((q.eval(3.0) - 4.0).abs() < 1e-10);
        assert!(r.is_zero());
    }
    #[test]
    fn test_polynomial_gcd() {
        let a = Polynomial::new(vec![-1.0, 0.0, 1.0]);
        let b = Polynomial::new(vec![-1.0, 1.0]);
        let g = a.gcd(&b);
        assert!((g.eval(1.0)).abs() < 1e-10);
        assert!(g.degree() == 1);
    }
    #[test]
    fn test_taylor_exp() {
        let e = var("x").exp_expr();
        let t = taylor_expand(&e, "x", 0.0, 5);
        assert!((t.eval(1.0) - std::f64::consts::E).abs() < 0.01);
    }
    #[test]
    fn test_taylor_sin() {
        let e = var("x").sin_expr();
        let t = taylor_expand(&e, "x", 0.0, 7);
        assert!((t.eval(0.5) - 0.5_f64.sin()).abs() < 1e-6);
    }
    #[test]
    fn test_partial_fractions() {
        let num = Polynomial::constant(1.0);
        let roots = vec![1.0, 2.0];
        let terms = partial_fraction_decompose(&num, &roots);
        assert_eq!(terms.len(), 2);
        let val: f64 = terms.iter().map(|t| t.coefficient / (3.0 - t.root)).sum();
        assert!((val - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_cse_basic() {
        let sub = var("x").add_expr(cst(1.0));
        let e = sub.clone().mul_expr(sub);
        let result = common_subexpression_elimination(&e);
        assert!(!result.bindings.is_empty());
    }
    #[test]
    fn test_substitute() {
        let e = var("x").add_expr(var("y"));
        let s = substitute(&e, "x", &cst(5.0));
        let mut env = HashMap::new();
        env.insert("y".into(), 3.0);
        assert!((eval_expr(&s, &env).unwrap() - 8.0).abs() < 1e-15);
    }
    #[test]
    fn test_substitute_many() {
        let e = var("x").mul_expr(var("y"));
        let mut subs = HashMap::new();
        subs.insert("x".into(), cst(2.0));
        subs.insert("y".into(), cst(3.0));
        let s = substitute_many(&e, &subs);
        let env = HashMap::new();
        assert!((eval_expr(&simplify(&s), &env).unwrap() - 6.0).abs() < 1e-15);
    }
    #[test]
    fn test_display_const() {
        assert_eq!(format!("{}", cst(3.0)), "3");
    }
    #[test]
    fn test_display_add() {
        let e = var("x").add_expr(cst(1.0));
        assert_eq!(format!("{}", e), "(x + 1)");
    }
    #[test]
    fn test_variables() {
        let e = var("x").add_expr(var("y").mul_expr(var("x")));
        let v = e.variables();
        assert_eq!(v, vec!["x", "y"]);
    }
    #[test]
    fn test_node_count() {
        let e = var("x").add_expr(cst(1.0));
        assert_eq!(e.node_count(), 3);
    }
    #[test]
    fn test_depth() {
        let e = var("x").add_expr(cst(1.0));
        assert_eq!(e.depth(), 2);
    }
    #[test]
    fn test_diff_n() {
        let e = var("x").pow_expr(cst(3.0));
        let d2 = diff_n(&e, "x", 2);
        let mut env = HashMap::new();
        env.insert("x".into(), 2.0);
        assert!((eval_expr(&d2, &env).unwrap() - 12.0).abs() < 1e-8);
    }
    #[test]
    fn test_gradient() {
        let f = var("x")
            .pow_expr(cst(2.0))
            .add_expr(var("y").pow_expr(cst(2.0)));
        let g = gradient(&f, &["x", "y"]);
        let mut env = HashMap::new();
        env.insert("x".into(), 3.0);
        env.insert("y".into(), 4.0);
        assert!((eval_expr(&g[0], &env).unwrap() - 6.0).abs() < 1e-8);
        assert!((eval_expr(&g[1], &env).unwrap() - 8.0).abs() < 1e-8);
    }
    #[test]
    fn test_hessian() {
        let f = var("x").pow_expr(cst(2.0)).mul_expr(var("y"));
        let h = hessian(&f, &["x", "y"]);
        let mut env = HashMap::new();
        env.insert("x".into(), 1.0);
        env.insert("y".into(), 2.0);
        assert!((eval_expr(&h[0][0], &env).unwrap() - 4.0).abs() < 1e-8);
    }
    #[test]
    fn test_jacobian() {
        let f1 = var("x").mul_expr(var("y"));
        let f2 = var("x").add_expr(var("y"));
        let j = jacobian(&[f1, f2], &["x", "y"]);
        let mut env = HashMap::new();
        env.insert("x".into(), 2.0);
        env.insert("y".into(), 3.0);
        assert!((eval_expr(&j[0][0], &env).unwrap() - 3.0).abs() < 1e-8);
        assert!((eval_expr(&j[1][1], &env).unwrap() - 1.0).abs() < 1e-8);
    }
    #[test]
    fn test_gradient_check() {
        let e = var("x").pow_expr(cst(3.0));
        let err = gradient_check(&e, "x", &[1.0, 2.0, 3.0], 1e-6);
        assert!(err < 1e-4);
    }
    #[test]
    fn test_rational_function_eval() {
        let rf = RationalFunction::new(
            Polynomial::new(vec![1.0, 1.0]),
            Polynomial::new(vec![-1.0, 1.0]),
        );
        assert!((rf.eval(2.0).unwrap() - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_rational_function_add() {
        let r1 = RationalFunction::new(Polynomial::constant(1.0), Polynomial::x());
        let r2 = RationalFunction::new(Polynomial::constant(1.0), Polynomial::x());
        let r3 = r1.add(&r2);
        assert!((r3.eval(3.0).unwrap() - 2.0 / 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_lagrange_interpolation() {
        let p = lagrange_interpolation(&[(0.0, 1.0), (1.0, 3.0), (2.0, 7.0)]);
        assert!((p.eval(3.0) - 13.0).abs() < 1e-10);
    }
    #[test]
    fn test_chebyshev_t() {
        let t2 = chebyshev_t(2);
        assert!((t2.eval(0.5) - (-0.5)).abs() < 1e-12);
    }
    #[test]
    fn test_hermite_he() {
        let h2 = hermite_he(2);
        assert!((h2.eval(2.0) - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_legendre_p() {
        let p2 = legendre_p(2);
        assert!((p2.eval(1.0) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_find_root_newton() {
        let p = Polynomial::new(vec![-4.0, 0.0, 1.0]);
        let r = find_root_newton(&p, 3.0, 1e-12, 100).unwrap();
        assert!((r - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_polynomial_composition() {
        let p = Polynomial::new(vec![1.0, 1.0]);
        let q = Polynomial::new(vec![0.0, 2.0]);
        let c = p.compose(&q);
        assert!((c.eval(3.0) - 7.0).abs() < 1e-12);
    }
    #[test]
    fn test_sturm_root_count() {
        let p = Polynomial::new(vec![-1.0, 0.0, 1.0]);
        assert_eq!(count_roots_in_interval(&p, -2.0, 2.0), 2);
    }
    #[test]
    fn test_expand_product() {
        let e = var("x")
            .add_expr(cst(1.0))
            .mul_expr(var("x").add_expr(cst(2.0)));
        let ex = expand(&e);
        let mut env = HashMap::new();
        env.insert("x".into(), 5.0);
        assert!((eval_expr(&simplify(&ex), &env).unwrap() - 42.0).abs() < 1e-10);
    }
    #[test]
    fn test_is_polynomial_expr() {
        let e = var("x").pow_expr(cst(2.0)).add_expr(cst(1.0));
        assert!(is_polynomial_expr(&e));
        let e2 = var("x").sin_expr();
        assert!(!is_polynomial_expr(&e2));
    }
    #[test]
    fn test_polynomial_display() {
        let p = Polynomial::new(vec![1.0, 0.0, 3.0]);
        let s = format!("{}", p);
        assert!(s.contains("x^2"));
    }
    #[test]
    fn test_polynomial_sign_changes() {
        let p = Polynomial::new(vec![1.0, -1.0, 1.0]);
        assert_eq!(p.sign_changes(), 2);
    }
    #[test]
    fn test_expr_to_polynomial() {
        let e = var("x")
            .pow_expr(cst(2.0))
            .add_expr(cst(3.0).mul_expr(var("x")))
            .add_expr(cst(2.0));
        let p = expr_to_polynomial(&e, "x").unwrap();
        assert!((p.eval(1.0) - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_collect_terms() {
        let e = var("x").add_expr(var("x")).add_expr(cst(1.0));
        let c = collect_terms(&e, "x");
        let mut env = HashMap::new();
        env.insert("x".into(), 5.0);
        assert!((eval_expr(&c, &env).unwrap() - 11.0).abs() < 1e-10);
    }
    #[test]
    fn test_map_expr() {
        let e = var("x").add_expr(cst(5.0));
        let mapped = map_expr(&e, &|ex| match ex {
            Expr::Const(_) => cst(0.0),
            other => other.clone(),
        });
        let mut env = HashMap::new();
        env.insert("x".into(), 3.0);
        assert!((eval_expr(&mapped, &env).unwrap() - 3.0).abs() < 1e-15);
    }
    #[test]
    fn test_simplify_double_neg() {
        let e = var("x").neg_expr().neg_expr();
        let s = simplify(&e);
        assert_eq!(s, var("x"));
    }
    #[test]
    fn test_polynomial_leading_coeff() {
        let p = Polynomial::new(vec![1.0, 2.0, 5.0]);
        assert!((p.leading_coeff() - 5.0).abs() < 1e-15);
    }
    #[test]
    fn test_partial_fractions_to_expr() {
        let terms = vec![
            PartialFractionTerm {
                coefficient: -1.0,
                root: 1.0,
            },
            PartialFractionTerm {
                coefficient: 1.0,
                root: 2.0,
            },
        ];
        let e = partial_fractions_to_expr(&terms, "x");
        let mut env = HashMap::new();
        env.insert("x".into(), 3.0);
        assert!((eval_expr(&e, &env).unwrap() - 0.5).abs() < 1e-12);
    }
}
