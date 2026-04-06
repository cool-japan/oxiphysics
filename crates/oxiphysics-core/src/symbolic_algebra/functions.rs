//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;

use super::types::{CseResult, Expr, PartialFractionTerm, Polynomial};

/// Create a variable expression.
pub fn var(name: &str) -> Expr {
    Expr::Var(name.to_string())
}
/// Create a constant expression.
pub fn cst(v: f64) -> Expr {
    Expr::Const(v)
}
/// Evaluate an expression given a variable-binding environment.
///
/// Returns `None` when an unbound variable is encountered or a numerical
/// error occurs (e.g., `ln(0)`).
pub fn eval_expr(expr: &Expr, env: &HashMap<String, f64>) -> Option<f64> {
    match expr {
        Expr::Const(v) => Some(*v),
        Expr::Var(s) => env.get(s).copied(),
        Expr::Add(a, b) => Some(eval_expr(a, env)? + eval_expr(b, env)?),
        Expr::Mul(a, b) => Some(eval_expr(a, env)? * eval_expr(b, env)?),
        Expr::Div(a, b) => {
            let denom = eval_expr(b, env)?;
            if denom.abs() < 1e-300 {
                None
            } else {
                Some(eval_expr(a, env)? / denom)
            }
        }
        Expr::Pow(a, b) => Some(eval_expr(a, env)?.powf(eval_expr(b, env)?)),
        Expr::Neg(a) => Some(-eval_expr(a, env)?),
        Expr::Sin(a) => Some(eval_expr(a, env)?.sin()),
        Expr::Cos(a) => Some(eval_expr(a, env)?.cos()),
        Expr::Exp(a) => Some(eval_expr(a, env)?.exp()),
        Expr::Ln(a) => {
            let v = eval_expr(a, env)?;
            if v <= 0.0 { None } else { Some(v.ln()) }
        }
    }
}
/// Substitute every occurrence of `name` with `replacement` in `expr`.
pub fn substitute(expr: &Expr, name: &str, replacement: &Expr) -> Expr {
    match expr {
        Expr::Const(v) => Expr::Const(*v),
        Expr::Var(s) => {
            if s == name {
                replacement.clone()
            } else {
                Expr::Var(s.clone())
            }
        }
        Expr::Add(a, b) => Expr::Add(
            Box::new(substitute(a, name, replacement)),
            Box::new(substitute(b, name, replacement)),
        ),
        Expr::Mul(a, b) => Expr::Mul(
            Box::new(substitute(a, name, replacement)),
            Box::new(substitute(b, name, replacement)),
        ),
        Expr::Div(a, b) => Expr::Div(
            Box::new(substitute(a, name, replacement)),
            Box::new(substitute(b, name, replacement)),
        ),
        Expr::Pow(a, b) => Expr::Pow(
            Box::new(substitute(a, name, replacement)),
            Box::new(substitute(b, name, replacement)),
        ),
        Expr::Neg(a) => Expr::Neg(Box::new(substitute(a, name, replacement))),
        Expr::Sin(a) => Expr::Sin(Box::new(substitute(a, name, replacement))),
        Expr::Cos(a) => Expr::Cos(Box::new(substitute(a, name, replacement))),
        Expr::Exp(a) => Expr::Exp(Box::new(substitute(a, name, replacement))),
        Expr::Ln(a) => Expr::Ln(Box::new(substitute(a, name, replacement))),
    }
}
/// Symbolically differentiate `expr` with respect to variable `var`.
pub fn diff(expr: &Expr, v: &str) -> Expr {
    match expr {
        Expr::Const(_) => cst(0.0),
        Expr::Var(s) => {
            if s == v {
                cst(1.0)
            } else {
                cst(0.0)
            }
        }
        Expr::Add(a, b) => diff(a, v).add_expr(diff(b, v)),
        Expr::Mul(a, b) => diff(a, v)
            .mul_expr((**b).clone())
            .add_expr((**a).clone().mul_expr(diff(b, v))),
        Expr::Div(a, b) => {
            let num = diff(a, v)
                .mul_expr((**b).clone())
                .sub_expr((**a).clone().mul_expr(diff(b, v)));
            let den = (**b).clone().pow_expr(cst(2.0));
            num.div_expr(den)
        }
        Expr::Pow(base, exp) => {
            let f = base;
            let g = exp;
            let fp = diff(f, v);
            let gp = diff(g, v);
            let term1 = gp.mul_expr((**f).clone().ln_expr());
            let term2 = (**g).clone().mul_expr(fp.div_expr((**f).clone()));
            expr.clone().mul_expr(term1.add_expr(term2))
        }
        Expr::Neg(a) => Expr::Neg(Box::new(diff(a, v))),
        Expr::Sin(a) => diff(a, v).mul_expr((**a).clone().cos_expr()),
        Expr::Cos(a) => diff(a, v).mul_expr((**a).clone().sin_expr()).neg_expr(),
        Expr::Exp(a) => diff(a, v).mul_expr(expr.clone()),
        Expr::Ln(a) => diff(a, v).div_expr((**a).clone()),
    }
}
/// Apply algebraic simplification rules repeatedly until a fixed point.
pub fn simplify(expr: &Expr) -> Expr {
    let mut cur = expr.clone();
    for _ in 0..20 {
        let next = simplify_once(&cur);
        if next == cur {
            break;
        }
        cur = next;
    }
    cur
}
/// One pass of simplification rules.
#[allow(dead_code)]
pub fn simplify_once(expr: &Expr) -> Expr {
    match expr {
        Expr::Const(v) => Expr::Const(*v),
        Expr::Var(s) => Expr::Var(s.clone()),
        Expr::Add(a, b) => {
            let a = simplify_once(a);
            let b = simplify_once(b);
            if a.is_zero() {
                return b;
            }
            if b.is_zero() {
                return a;
            }
            if let (Expr::Const(va), Expr::Const(vb)) = (&a, &b) {
                return Expr::Const(va + vb);
            }
            if a == b {
                return cst(2.0).mul_expr(a);
            }
            Expr::Add(Box::new(a), Box::new(b))
        }
        Expr::Mul(a, b) => {
            let a = simplify_once(a);
            let b = simplify_once(b);
            if a.is_zero() || b.is_zero() {
                return cst(0.0);
            }
            if a.is_one() {
                return b;
            }
            if b.is_one() {
                return a;
            }
            if let (Expr::Const(va), Expr::Const(vb)) = (&a, &b) {
                return Expr::Const(va * vb);
            }
            if a == b {
                return a.pow_expr(cst(2.0));
            }
            Expr::Mul(Box::new(a), Box::new(b))
        }
        Expr::Div(a, b) => {
            let a = simplify_once(a);
            let b = simplify_once(b);
            if a.is_zero() {
                return cst(0.0);
            }
            if b.is_one() {
                return a;
            }
            if a == b {
                return cst(1.0);
            }
            if let (Expr::Const(va), Expr::Const(vb)) = (&a, &b)
                && vb.abs() > 1e-300
            {
                return Expr::Const(va / vb);
            }
            Expr::Div(Box::new(a), Box::new(b))
        }
        Expr::Pow(a, b) => {
            let a = simplify_once(a);
            let b = simplify_once(b);
            if b.is_zero() {
                return cst(1.0);
            }
            if b.is_one() {
                return a;
            }
            if a.is_zero() {
                return cst(0.0);
            }
            if a.is_one() {
                return cst(1.0);
            }
            if let (Expr::Const(va), Expr::Const(vb)) = (&a, &b) {
                return Expr::Const(va.powf(*vb));
            }
            Expr::Pow(Box::new(a), Box::new(b))
        }
        Expr::Neg(a) => {
            let a = simplify_once(a);
            if let Expr::Const(v) = &a {
                return Expr::Const(-v);
            }
            if let Expr::Neg(inner) = a {
                return *inner;
            }
            Expr::Neg(Box::new(a))
        }
        Expr::Sin(a) => {
            let a = simplify_once(a);
            if let Expr::Const(v) = &a {
                return Expr::Const(v.sin());
            }
            Expr::Sin(Box::new(a))
        }
        Expr::Cos(a) => {
            let a = simplify_once(a);
            if let Expr::Const(v) = &a {
                return Expr::Const(v.cos());
            }
            Expr::Cos(Box::new(a))
        }
        Expr::Exp(a) => {
            let a = simplify_once(a);
            if let Expr::Const(v) = &a {
                return Expr::Const(v.exp());
            }
            if let Expr::Ln(inner) = &a {
                return *inner.clone();
            }
            Expr::Exp(Box::new(a))
        }
        Expr::Ln(a) => {
            let a = simplify_once(a);
            if let Expr::Const(v) = &a
                && *v > 0.0
            {
                return Expr::Const(v.ln());
            }
            if let Expr::Exp(inner) = &a {
                return *inner.clone();
            }
            Expr::Ln(Box::new(a))
        }
    }
}
/// Compute the Taylor expansion of `expr` about `x0` up to order `n`.
///
/// Returns a [`Polynomial`] whose variable is `(x - x0)`, i.e. you
/// should evaluate the result at `x - x0`.
pub fn taylor_expand(expr: &Expr, var_name: &str, x0: f64, n: usize) -> Polynomial {
    let mut env = HashMap::new();
    env.insert(var_name.to_string(), x0);
    let mut coeffs = Vec::with_capacity(n + 1);
    let mut current = expr.clone();
    let mut factorial = 1.0_f64;
    for k in 0..=n {
        if k > 0 {
            factorial *= k as f64;
        }
        let val = eval_expr(&simplify(&current), &env).unwrap_or(0.0);
        coeffs.push(val / factorial);
        current = diff(&current, var_name);
    }
    Polynomial::new(coeffs)
}
/// Perform simple common subexpression elimination.
///
/// Finds subexpressions that appear more than once and extracts them
/// into named temporaries `_cse0`, `_cse1`, ...
pub fn common_subexpression_elimination(expr: &Expr) -> CseResult {
    let mut counts: HashMap<String, usize> = HashMap::new();
    count_subexpressions(expr, &mut counts);
    let mut bindings: Vec<(String, Expr)> = Vec::new();
    let mut mapping: HashMap<String, String> = HashMap::new();
    let reduced = cse_rewrite(expr, &counts, &mut bindings, &mut mapping);
    CseResult { bindings, reduced }
}
/// Count occurrences of each subexpression by its Display string.
pub fn count_subexpressions(expr: &Expr, counts: &mut HashMap<String, usize>) {
    let key = format!("{}", expr);
    *counts.entry(key).or_insert(0) += 1;
    match expr {
        Expr::Const(_) | Expr::Var(_) => {}
        Expr::Add(a, b) | Expr::Mul(a, b) | Expr::Pow(a, b) | Expr::Div(a, b) => {
            count_subexpressions(a, counts);
            count_subexpressions(b, counts);
        }
        Expr::Neg(a) | Expr::Sin(a) | Expr::Cos(a) | Expr::Exp(a) | Expr::Ln(a) => {
            count_subexpressions(a, counts);
        }
    }
}
/// Rewrite expression, replacing repeated subexpressions with CSE vars.
pub fn cse_rewrite(
    expr: &Expr,
    counts: &HashMap<String, usize>,
    bindings: &mut Vec<(String, Expr)>,
    mapping: &mut HashMap<String, String>,
) -> Expr {
    match expr {
        Expr::Const(_) | Expr::Var(_) => return expr.clone(),
        _ => {}
    }
    let key = format!("{}", expr);
    if let Some(&c) = counts.get(&key)
        && c > 1
    {
        if let Some(name) = mapping.get(&key) {
            return Expr::Var(name.clone());
        }
        let name = format!("_cse{}", bindings.len());
        mapping.insert(key, name.clone());
        let rewritten = cse_rewrite_children(expr, counts, bindings, mapping);
        bindings.push((name.clone(), rewritten));
        return Expr::Var(name);
    }
    cse_rewrite_children(expr, counts, bindings, mapping)
}
/// Rewrite only the children of an expression.
pub fn cse_rewrite_children(
    expr: &Expr,
    counts: &HashMap<String, usize>,
    bindings: &mut Vec<(String, Expr)>,
    mapping: &mut HashMap<String, String>,
) -> Expr {
    match expr {
        Expr::Const(_) | Expr::Var(_) => expr.clone(),
        Expr::Add(a, b) => Expr::Add(
            Box::new(cse_rewrite(a, counts, bindings, mapping)),
            Box::new(cse_rewrite(b, counts, bindings, mapping)),
        ),
        Expr::Mul(a, b) => Expr::Mul(
            Box::new(cse_rewrite(a, counts, bindings, mapping)),
            Box::new(cse_rewrite(b, counts, bindings, mapping)),
        ),
        Expr::Div(a, b) => Expr::Div(
            Box::new(cse_rewrite(a, counts, bindings, mapping)),
            Box::new(cse_rewrite(b, counts, bindings, mapping)),
        ),
        Expr::Pow(a, b) => Expr::Pow(
            Box::new(cse_rewrite(a, counts, bindings, mapping)),
            Box::new(cse_rewrite(b, counts, bindings, mapping)),
        ),
        Expr::Neg(a) => Expr::Neg(Box::new(cse_rewrite(a, counts, bindings, mapping))),
        Expr::Sin(a) => Expr::Sin(Box::new(cse_rewrite(a, counts, bindings, mapping))),
        Expr::Cos(a) => Expr::Cos(Box::new(cse_rewrite(a, counts, bindings, mapping))),
        Expr::Exp(a) => Expr::Exp(Box::new(cse_rewrite(a, counts, bindings, mapping))),
        Expr::Ln(a) => Expr::Ln(Box::new(cse_rewrite(a, counts, bindings, mapping))),
    }
}
/// Decompose `numerator / denominator` into partial fractions given the
/// known roots of the denominator.
///
/// Assumes the denominator has simple real roots and the degree of the
/// numerator is less than the degree of the denominator.
///
/// Uses the cover-up method (Heaviside's method): for root r_i,
/// A_i = numerator(r_i) / product_{j != i} (r_i - r_j).
pub fn partial_fraction_decompose(
    numerator: &Polynomial,
    roots: &[f64],
) -> Vec<PartialFractionTerm> {
    let n = roots.len();
    let mut terms = Vec::with_capacity(n);
    for i in 0..n {
        let ri = roots[i];
        let num_val = numerator.eval(ri);
        let mut den_product = 1.0;
        for (j, &rj) in roots.iter().enumerate() {
            if j != i {
                den_product *= ri - rj;
            }
        }
        terms.push(PartialFractionTerm {
            coefficient: num_val / den_product,
            root: ri,
        });
    }
    terms
}
/// Reconstruct an expression from partial fraction terms.
///
/// Returns sum of `A_i / (x - r_i)` terms as an `Expr`.
pub fn partial_fractions_to_expr(terms: &[PartialFractionTerm], var_name: &str) -> Expr {
    let v = var(var_name);
    let mut acc: Option<Expr> = None;
    for t in terms {
        let term = cst(t.coefficient).div_expr(v.clone().sub_expr(cst(t.root)));
        acc = Some(match acc {
            None => term,
            Some(a) => a.add_expr(term),
        });
    }
    acc.unwrap_or(cst(0.0))
}
/// Try to convert an expression tree to a polynomial in a single variable.
///
/// Returns `None` if the expression contains non-polynomial operations
/// (sin, cos, exp, ln) or involves more than one variable.
pub fn expr_to_polynomial(expr: &Expr, var_name: &str) -> Option<Polynomial> {
    match expr {
        Expr::Const(v) => Some(Polynomial::constant(*v)),
        Expr::Var(s) => {
            if s == var_name {
                Some(Polynomial::x())
            } else {
                None
            }
        }
        Expr::Add(a, b) => {
            let pa = expr_to_polynomial(a, var_name)?;
            let pb = expr_to_polynomial(b, var_name)?;
            Some(pa.add(&pb))
        }
        Expr::Mul(a, b) => {
            let pa = expr_to_polynomial(a, var_name)?;
            let pb = expr_to_polynomial(b, var_name)?;
            Some(pa.mul(&pb))
        }
        Expr::Neg(a) => {
            let pa = expr_to_polynomial(a, var_name)?;
            Some(pa.scale(-1.0))
        }
        Expr::Pow(base, exp) => {
            if let Expr::Const(n) = exp.as_ref() {
                let ni = *n as usize;
                if (*n - ni as f64).abs() < 1e-12 && ni <= 20 {
                    let pb = expr_to_polynomial(base, var_name)?;
                    let mut result = Polynomial::constant(1.0);
                    for _ in 0..ni {
                        result = result.mul(&pb);
                    }
                    return Some(result);
                }
            }
            None
        }
        Expr::Div(_, _) => None,
        _ => None,
    }
}
/// Check structural equality of two expressions.
#[allow(dead_code)]
pub fn expr_equal(a: &Expr, b: &Expr) -> bool {
    a == b
}
/// Compute the n-th derivative of `expr` with respect to `var`.
pub fn diff_n(expr: &Expr, v: &str, n: usize) -> Expr {
    let mut result = expr.clone();
    for _ in 0..n {
        result = diff(&result, v);
    }
    simplify(&result)
}
/// Compute the gradient: partial derivatives with respect to each variable.
pub fn gradient(expr: &Expr, vars: &[&str]) -> Vec<Expr> {
    vars.iter().map(|v| simplify(&diff(expr, v))).collect()
}
/// Compute the Hessian matrix (second partial derivatives).
///
/// Returns a row-major flattened matrix of size `vars.len() x vars.len()`.
pub fn hessian(expr: &Expr, vars: &[&str]) -> Vec<Vec<Expr>> {
    let grad = gradient(expr, vars);
    grad.iter()
        .map(|gi| vars.iter().map(|v| simplify(&diff(gi, v))).collect())
        .collect()
}
/// Compute the Jacobian matrix for a system of expressions.
///
/// `exprs[i]` is differentiated with respect to `vars[j]`.
pub fn jacobian(exprs: &[Expr], vars: &[&str]) -> Vec<Vec<Expr>> {
    exprs
        .iter()
        .map(|e| vars.iter().map(|v| simplify(&diff(e, v))).collect())
        .collect()
}
/// Find a root of polynomial `p` near `x0` using Newton's method.
///
/// Returns `None` if the method does not converge in `max_iter` steps.
pub fn find_root_newton(p: &Polynomial, x0: f64, tol: f64, max_iter: usize) -> Option<f64> {
    let dp = p.derivative();
    let mut x = x0;
    for _ in 0..max_iter {
        let fx = p.eval(x);
        if fx.abs() < tol {
            return Some(x);
        }
        let dfx = dp.eval(x);
        if dfx.abs() < 1e-300 {
            return None;
        }
        x -= fx / dfx;
    }
    if p.eval(x).abs() < tol * 100.0 {
        Some(x)
    } else {
        None
    }
}
/// Construct the Lagrange interpolating polynomial through the given points.
pub fn lagrange_interpolation(points: &[(f64, f64)]) -> Polynomial {
    let n = points.len();
    let mut result = Polynomial::zero();
    for i in 0..n {
        let mut basis = Polynomial::constant(1.0);
        for j in 0..n {
            if j == i {
                continue;
            }
            let scale = 1.0 / (points[i].0 - points[j].0);
            let factor = Polynomial::new(vec![-points[j].0 * scale, scale]);
            basis = basis.mul(&factor);
        }
        result = result.add(&basis.scale(points[i].1));
    }
    result
}
/// Generate the n-th Chebyshev polynomial of the first kind T_n(x).
pub fn chebyshev_t(n: usize) -> Polynomial {
    if n == 0 {
        return Polynomial::constant(1.0);
    }
    if n == 1 {
        return Polynomial::x();
    }
    let two_x = Polynomial::new(vec![0.0, 2.0]);
    let mut t_prev2 = Polynomial::constant(1.0);
    let mut t_prev1 = Polynomial::x();
    for _ in 2..=n {
        let t_next = two_x.mul(&t_prev1).sub(&t_prev2);
        t_prev2 = t_prev1;
        t_prev1 = t_next;
    }
    t_prev1
}
/// Generate the n-th (probabilist's) Hermite polynomial He_n(x).
pub fn hermite_he(n: usize) -> Polynomial {
    if n == 0 {
        return Polynomial::constant(1.0);
    }
    if n == 1 {
        return Polynomial::x();
    }
    let x_poly = Polynomial::x();
    let mut h_prev2 = Polynomial::constant(1.0);
    let mut h_prev1 = Polynomial::x();
    for k in 2..=n {
        let h_next = x_poly.mul(&h_prev1).sub(&h_prev2.scale((k - 1) as f64));
        h_prev2 = h_prev1;
        h_prev1 = h_next;
    }
    h_prev1
}
/// Generate the n-th Legendre polynomial P_n(x) via Bonnet's recursion.
pub fn legendre_p(n: usize) -> Polynomial {
    if n == 0 {
        return Polynomial::constant(1.0);
    }
    if n == 1 {
        return Polynomial::x();
    }
    let x_poly = Polynomial::x();
    let mut p_prev2 = Polynomial::constant(1.0);
    let mut p_prev1 = Polynomial::x();
    for k in 2..=n {
        let kf = k as f64;
        let term1 = x_poly.mul(&p_prev1).scale(2.0 * kf - 1.0);
        let term2 = p_prev2.scale(kf - 1.0);
        let p_next = term1.sub(&term2).scale(1.0 / kf);
        p_prev2 = p_prev1;
        p_prev1 = p_next;
    }
    p_prev1
}
/// Compute an approximate complexity score for an expression.
///
/// This can be used to compare different algebraically equivalent forms
/// and pick the simplest.
#[allow(dead_code)]
pub fn complexity(expr: &Expr) -> usize {
    match expr {
        Expr::Const(_) | Expr::Var(_) => 1,
        Expr::Add(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => 1 + complexity(a) + complexity(b),
        Expr::Pow(a, b) => 2 + complexity(a) + complexity(b),
        Expr::Neg(a) => 1 + complexity(a),
        Expr::Sin(a) | Expr::Cos(a) | Expr::Exp(a) | Expr::Ln(a) => 2 + complexity(a),
    }
}
/// Apply a transformation function to every node in an expression
/// (bottom-up).
pub fn map_expr<F>(expr: &Expr, f: &F) -> Expr
where
    F: Fn(&Expr) -> Expr,
{
    let mapped = match expr {
        Expr::Const(_) | Expr::Var(_) => expr.clone(),
        Expr::Add(a, b) => Expr::Add(Box::new(map_expr(a, f)), Box::new(map_expr(b, f))),
        Expr::Mul(a, b) => Expr::Mul(Box::new(map_expr(a, f)), Box::new(map_expr(b, f))),
        Expr::Div(a, b) => Expr::Div(Box::new(map_expr(a, f)), Box::new(map_expr(b, f))),
        Expr::Pow(a, b) => Expr::Pow(Box::new(map_expr(a, f)), Box::new(map_expr(b, f))),
        Expr::Neg(a) => Expr::Neg(Box::new(map_expr(a, f))),
        Expr::Sin(a) => Expr::Sin(Box::new(map_expr(a, f))),
        Expr::Cos(a) => Expr::Cos(Box::new(map_expr(a, f))),
        Expr::Exp(a) => Expr::Exp(Box::new(map_expr(a, f))),
        Expr::Ln(a) => Expr::Ln(Box::new(map_expr(a, f))),
    };
    f(&mapped)
}
/// Algebraically expand an expression (distribute `*` over `+`).
pub fn expand(expr: &Expr) -> Expr {
    match expr {
        Expr::Const(_) | Expr::Var(_) => expr.clone(),
        Expr::Add(a, b) => expand(a).add_expr(expand(b)),
        Expr::Neg(a) => expand(a).neg_expr(),
        Expr::Mul(a, b) => {
            let ea = expand(a);
            let eb = expand(b);
            distribute(&ea, &eb)
        }
        Expr::Div(a, b) => expand(a).div_expr(expand(b)),
        Expr::Pow(base, exp) => {
            if let Expr::Const(n) = exp.as_ref() {
                let ni = *n as usize;
                if (*n - ni as f64).abs() < 1e-12 && (2..=8).contains(&ni) {
                    let eb = expand(base);
                    let mut result = eb.clone();
                    for _ in 1..ni {
                        result = distribute(&result, &eb);
                    }
                    return result;
                }
            }
            expand(base).pow_expr(expand(exp))
        }
        Expr::Sin(a) => Expr::Sin(Box::new(expand(a))),
        Expr::Cos(a) => Expr::Cos(Box::new(expand(a))),
        Expr::Exp(a) => Expr::Exp(Box::new(expand(a))),
        Expr::Ln(a) => Expr::Ln(Box::new(expand(a))),
    }
}
/// Distribute multiplication: (a+b)*c = a*c + b*c, etc.
pub fn distribute(a: &Expr, b: &Expr) -> Expr {
    match a {
        Expr::Add(a1, a2) => distribute(a1, b).add_expr(distribute(a2, b)),
        _ => match b {
            Expr::Add(b1, b2) => distribute(a, b1).add_expr(distribute(a, b2)),
            _ => a.clone().mul_expr(b.clone()),
        },
    }
}
/// Collect like terms in a polynomial expression.
///
/// Attempts to convert to a polynomial, then converts back to a
/// canonical expression form.  Returns the original if conversion fails.
pub fn collect_terms(expr: &Expr, var_name: &str) -> Expr {
    match expr_to_polynomial(expr, var_name) {
        Some(p) => p.to_expr(var_name),
        None => expr.clone(),
    }
}
/// Numerically check symbolic differentiation using finite differences.
///
/// Returns the maximum absolute difference between the symbolic derivative
/// evaluated at `test_points` and the finite-difference approximation.
pub fn gradient_check(expr: &Expr, var_name: &str, test_points: &[f64], h: f64) -> f64 {
    let d = simplify(&diff(expr, var_name));
    let mut max_err = 0.0_f64;
    for &x in test_points {
        let mut env = HashMap::new();
        env.insert(var_name.to_string(), x);
        let symbolic_val = eval_expr(&d, &env).unwrap_or(f64::NAN);
        env.insert(var_name.to_string(), x + h);
        let fp = eval_expr(expr, &env).unwrap_or(f64::NAN);
        env.insert(var_name.to_string(), x - h);
        let fm = eval_expr(expr, &env).unwrap_or(f64::NAN);
        let numerical = (fp - fm) / (2.0 * h);
        let err = (symbolic_val - numerical).abs();
        if err.is_finite() {
            max_err = max_err.max(err);
        }
    }
    max_err
}
/// Substitute multiple variables at once.
pub fn substitute_many(expr: &Expr, subs: &HashMap<String, Expr>) -> Expr {
    match expr {
        Expr::Const(v) => Expr::Const(*v),
        Expr::Var(s) => subs.get(s).cloned().unwrap_or_else(|| Expr::Var(s.clone())),
        Expr::Add(a, b) => Expr::Add(
            Box::new(substitute_many(a, subs)),
            Box::new(substitute_many(b, subs)),
        ),
        Expr::Mul(a, b) => Expr::Mul(
            Box::new(substitute_many(a, subs)),
            Box::new(substitute_many(b, subs)),
        ),
        Expr::Div(a, b) => Expr::Div(
            Box::new(substitute_many(a, subs)),
            Box::new(substitute_many(b, subs)),
        ),
        Expr::Pow(a, b) => Expr::Pow(
            Box::new(substitute_many(a, subs)),
            Box::new(substitute_many(b, subs)),
        ),
        Expr::Neg(a) => Expr::Neg(Box::new(substitute_many(a, subs))),
        Expr::Sin(a) => Expr::Sin(Box::new(substitute_many(a, subs))),
        Expr::Cos(a) => Expr::Cos(Box::new(substitute_many(a, subs))),
        Expr::Exp(a) => Expr::Exp(Box::new(substitute_many(a, subs))),
        Expr::Ln(a) => Expr::Ln(Box::new(substitute_many(a, subs))),
    }
}
/// Compute the Sturm sequence for polynomial `p`.
pub fn sturm_sequence(p: &Polynomial) -> Vec<Polynomial> {
    let mut seq = vec![p.clone(), p.derivative()];
    loop {
        let n = seq.len();
        let (_q, r) = seq[n - 2].div_rem(&seq[n - 1]);
        let neg_r = r.scale(-1.0);
        if neg_r.is_zero() {
            break;
        }
        seq.push(neg_r);
    }
    seq
}
/// Count the number of real roots of `p` in the interval `(a, b]`.
///
/// Uses Sturm's theorem.
pub fn count_roots_in_interval(p: &Polynomial, a: f64, b: f64) -> usize {
    let seq = sturm_sequence(p);
    let sign_changes_at = |x: f64| -> usize {
        let vals: Vec<f64> = seq.iter().map(|s| s.eval(x)).collect();
        let nonzero: Vec<f64> = vals.into_iter().filter(|v| v.abs() > 1e-15).collect();
        let mut c = 0;
        for w in nonzero.windows(2) {
            if w[0] * w[1] < 0.0 {
                c += 1;
            }
        }
        c
    };
    let sa = sign_changes_at(a);
    let sb = sign_changes_at(b);
    sa.saturating_sub(sb)
}
/// Check if an expression is a polynomial (no transcendentals).
pub fn is_polynomial_expr(expr: &Expr) -> bool {
    match expr {
        Expr::Const(_) | Expr::Var(_) => true,
        Expr::Add(a, b) | Expr::Mul(a, b) => is_polynomial_expr(a) && is_polynomial_expr(b),
        Expr::Pow(base, exp) => {
            if let Expr::Const(n) = exp.as_ref() {
                let ni = *n as i64;
                (*n - ni as f64).abs() < 1e-12 && ni >= 0 && is_polynomial_expr(base)
            } else {
                false
            }
        }
        Expr::Neg(a) => is_polynomial_expr(a),
        Expr::Div(_, _) | Expr::Sin(_) | Expr::Cos(_) | Expr::Exp(_) | Expr::Ln(_) => false,
    }
}
/// Check if an expression is a rational function (ratio of polynomials).
#[allow(dead_code)]
pub fn is_rational_expr(expr: &Expr) -> bool {
    match expr {
        Expr::Div(a, b) => is_polynomial_expr(a) && is_polynomial_expr(b),
        _ => is_polynomial_expr(expr),
    }
}
