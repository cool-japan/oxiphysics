// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Symbolic mathematics: expression trees, evaluation, differentiation, and simplification.
//!
//! # Quick start
//!
//! ```no_run
//! use oxiphysics_core::symbolic_math::{cst, var, diff, simplify, eval};
//! use std::collections::HashMap;
//!
//! // f(x) = x^2 + 1
//! let x = var("x");
//! let expr = x.pow(cst(2.0)).add_expr(cst(1.0));
//! let df = simplify(&diff(&expr, "x")); // 2*x
//! let mut vars = HashMap::new();
//! vars.insert("x".to_string(), 3.0);
//! assert!((eval(&df, &vars).unwrap() - 6.0).abs() < 1e-12);
//! ```

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Expression tree
// ---------------------------------------------------------------------------

/// A symbolic mathematical expression tree.
///
/// Each variant represents one kind of mathematical operation or atom.
#[derive(Debug, Clone, PartialEq)]
pub enum Expr {
    /// A numeric constant.
    Const(f64),
    /// A named variable.
    Var(String),
    /// Addition of two sub-expressions.
    Add(Box<Expr>, Box<Expr>),
    /// Multiplication of two sub-expressions.
    Mul(Box<Expr>, Box<Expr>),
    /// Exponentiation: base raised to exponent.
    Pow(Box<Expr>, Box<Expr>),
    /// Unary negation.
    Neg(Box<Expr>),
    /// Sine function.
    Sin(Box<Expr>),
    /// Cosine function.
    Cos(Box<Expr>),
    /// Natural exponential function.
    Exp(Box<Expr>),
    /// Natural logarithm.
    Ln(Box<Expr>),
}

// ---------------------------------------------------------------------------
// Convenience constructors (free functions)
// ---------------------------------------------------------------------------

/// Create a named-variable expression.
pub fn var(name: &str) -> Expr {
    Expr::Var(name.to_string())
}

/// Create a numeric constant expression.
pub fn cst(v: f64) -> Expr {
    Expr::Const(v)
}

// ---------------------------------------------------------------------------
// Expr: helper methods
// ---------------------------------------------------------------------------

impl Expr {
    /// Add another expression to `self`: `self + rhs`.
    pub fn add_expr(self, rhs: Expr) -> Expr {
        Expr::Add(Box::new(self), Box::new(rhs))
    }

    /// Subtract another expression from `self`: `self + (-rhs)`.
    pub fn sub_expr(self, rhs: Expr) -> Expr {
        Expr::Add(Box::new(self), Box::new(Expr::Neg(Box::new(rhs))))
    }

    /// Multiply `self` by another expression: `self * rhs`.
    pub fn mul_expr(self, rhs: Expr) -> Expr {
        Expr::Mul(Box::new(self), Box::new(rhs))
    }

    /// Raise `self` to the power of `exp`.
    pub fn pow(self, exp: Expr) -> Expr {
        Expr::Pow(Box::new(self), Box::new(exp))
    }

    /// Negate `self`: `-self`.
    ///
    /// Prefer using the unary `-` operator via `std::ops::Neg` instead of
    /// calling this method directly.
    fn neg_expr(self) -> Expr {
        Expr::Neg(Box::new(self))
    }

    /// Apply sine to `self`.
    pub fn sin(self) -> Expr {
        Expr::Sin(Box::new(self))
    }

    /// Apply cosine to `self`.
    pub fn cos(self) -> Expr {
        Expr::Cos(Box::new(self))
    }

    /// Apply the natural exponential to `self`.
    pub fn exp(self) -> Expr {
        Expr::Exp(Box::new(self))
    }

    /// Apply the natural logarithm to `self`.
    pub fn ln(self) -> Expr {
        Expr::Ln(Box::new(self))
    }
}

impl std::ops::Neg for Expr {
    type Output = Expr;
    fn neg(self) -> Expr {
        self.neg_expr()
    }
}

// ---------------------------------------------------------------------------
// Evaluation
// ---------------------------------------------------------------------------

/// Evaluate `expr` by substituting variable values from `vars`.
///
/// Returns `Err` if a variable is encountered that is not present in `vars`,
/// or if the expression is undefined (e.g., `ln(0)`).
pub fn eval(expr: &Expr, vars: &HashMap<String, f64>) -> Result<f64, String> {
    match expr {
        Expr::Const(c) => Ok(*c),
        Expr::Var(name) => vars
            .get(name)
            .copied()
            .ok_or_else(|| format!("undefined variable: {name}")),
        Expr::Add(a, b) => Ok(eval(a, vars)? + eval(b, vars)?),
        Expr::Mul(a, b) => Ok(eval(a, vars)? * eval(b, vars)?),
        Expr::Pow(base, exp) => Ok(eval(base, vars)?.powf(eval(exp, vars)?)),
        Expr::Neg(inner) => Ok(-eval(inner, vars)?),
        Expr::Sin(inner) => Ok(eval(inner, vars)?.sin()),
        Expr::Cos(inner) => Ok(eval(inner, vars)?.cos()),
        Expr::Exp(inner) => Ok(eval(inner, vars)?.exp()),
        Expr::Ln(inner) => {
            let v = eval(inner, vars)?;
            if v <= 0.0 {
                Err(format!("ln of non-positive value: {v}"))
            } else {
                Ok(v.ln())
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Symbolic differentiation
// ---------------------------------------------------------------------------

/// Symbolically differentiate `expr` with respect to the variable named `var`.
///
/// Returns a new expression for the derivative (not yet simplified).
pub fn diff(expr: &Expr, var: &str) -> Expr {
    match expr {
        Expr::Const(_) => cst(0.0),
        Expr::Var(name) => {
            if name == var {
                cst(1.0)
            } else {
                cst(0.0)
            }
        }
        // (f + g)' = f' + g'
        Expr::Add(f, g) => Expr::Add(Box::new(diff(f, var)), Box::new(diff(g, var))),
        // (f * g)' = f' * g + f * g'
        Expr::Mul(f, g) => Expr::Add(
            Box::new(Expr::Mul(Box::new(diff(f, var)), g.clone())),
            Box::new(Expr::Mul(f.clone(), Box::new(diff(g, var)))),
        ),
        // (f^n)' = n * f^(n-1) * f'   — general power rule via chain rule
        Expr::Pow(base, exp) => {
            // d/dx [f^g] = f^g * (g' * ln(f) + g * f'/f)
            // For the common case where exp is a Const, simplify:
            if let Expr::Const(n) = exp.as_ref() {
                let n = *n;
                // n * base^(n-1) * base'
                Expr::Mul(
                    Box::new(Expr::Mul(
                        Box::new(cst(n)),
                        Box::new(Expr::Pow(base.clone(), Box::new(cst(n - 1.0)))),
                    )),
                    Box::new(diff(base, var)),
                )
            } else {
                // General: (f^g)' = f^g * (g' * ln(f) + g * f'/f)
                let f = base.as_ref();
                let g = exp.as_ref();
                let fg = Expr::Pow(base.clone(), exp.clone());
                let term1 = Expr::Mul(Box::new(diff(g, var)), Box::new(Expr::Ln(base.clone())));
                let term2 = Expr::Mul(
                    g.clone().into(),
                    Box::new(Expr::Mul(
                        Box::new(diff(f, var)),
                        Box::new(Expr::Pow(base.clone(), Box::new(cst(-1.0)))),
                    )),
                );
                Expr::Mul(
                    Box::new(fg),
                    Box::new(Expr::Add(Box::new(term1), Box::new(term2))),
                )
            }
        }
        // (-f)' = -(f')
        Expr::Neg(f) => Expr::Neg(Box::new(diff(f, var))),
        // sin(f)' = cos(f) * f'
        Expr::Sin(f) => Expr::Mul(Box::new(Expr::Cos(f.clone())), Box::new(diff(f, var))),
        // cos(f)' = -sin(f) * f'
        Expr::Cos(f) => Expr::Neg(Box::new(Expr::Mul(
            Box::new(Expr::Sin(f.clone())),
            Box::new(diff(f, var)),
        ))),
        // exp(f)' = exp(f) * f'
        Expr::Exp(f) => Expr::Mul(Box::new(Expr::Exp(f.clone())), Box::new(diff(f, var))),
        // ln(f)' = f' / f  = f' * f^(-1)
        Expr::Ln(f) => Expr::Mul(
            Box::new(diff(f, var)),
            Box::new(Expr::Pow(f.clone(), Box::new(cst(-1.0)))),
        ),
    }
}

// ---------------------------------------------------------------------------
// Simplification
// ---------------------------------------------------------------------------

/// Apply basic algebraic simplifications to `expr`.
///
/// Rules applied (recursively, bottom-up):
/// - `0 + x = x`, `x + 0 = x`
/// - `0 * x = 0`, `x * 0 = 0`, `1 * x = x`, `x * 1 = x`
/// - `x ^ 0 = 1`, `x ^ 1 = x`
/// - `--x = x` (double negation)
/// - `-0 = 0`
/// - Fold constant sub-expressions to a single `Const`.
pub fn simplify(expr: &Expr) -> Expr {
    match expr {
        // Atoms are already in simplest form.
        Expr::Const(_) | Expr::Var(_) => expr.clone(),

        Expr::Add(a, b) => {
            let a = simplify(a);
            let b = simplify(b);
            // Constant folding
            if let (Expr::Const(x), Expr::Const(y)) = (&a, &b) {
                return cst(x + y);
            }
            // 0 + b = b
            if matches!(a, Expr::Const(x) if x == 0.0) {
                return b;
            }
            // a + 0 = a
            if matches!(b, Expr::Const(x) if x == 0.0) {
                return a;
            }
            Expr::Add(Box::new(a), Box::new(b))
        }

        Expr::Mul(a, b) => {
            let a = simplify(a);
            let b = simplify(b);
            // Constant folding
            if let (Expr::Const(x), Expr::Const(y)) = (&a, &b) {
                return cst(x * y);
            }
            // 0 * b = 0  or  a * 0 = 0
            if matches!(a, Expr::Const(x) if x == 0.0) {
                return cst(0.0);
            }
            if matches!(b, Expr::Const(x) if x == 0.0) {
                return cst(0.0);
            }
            // 1 * b = b
            if matches!(a, Expr::Const(x) if x == 1.0) {
                return b;
            }
            // a * 1 = a
            if matches!(b, Expr::Const(x) if x == 1.0) {
                return a;
            }
            Expr::Mul(Box::new(a), Box::new(b))
        }

        Expr::Pow(base, exp) => {
            let base = simplify(base);
            let exp = simplify(exp);
            // Constant folding
            if let (Expr::Const(b), Expr::Const(e)) = (&base, &exp) {
                return cst(b.powf(*e));
            }
            // x^0 = 1
            if matches!(exp, Expr::Const(e) if e == 0.0) {
                return cst(1.0);
            }
            // x^1 = x
            if matches!(exp, Expr::Const(e) if e == 1.0) {
                return base;
            }
            Expr::Pow(Box::new(base), Box::new(exp))
        }

        Expr::Neg(inner) => {
            let inner = simplify(inner);
            // Constant folding
            if let Expr::Const(c) = &inner {
                return cst(-c);
            }
            // Double negation: -(-x) = x
            if let Expr::Neg(x) = inner {
                return *x;
            }
            Expr::Neg(Box::new(inner))
        }

        Expr::Sin(inner) => Expr::Sin(Box::new(simplify(inner))),
        Expr::Cos(inner) => Expr::Cos(Box::new(simplify(inner))),
        Expr::Exp(inner) => {
            let inner = simplify(inner);
            if let Expr::Const(c) = &inner {
                return cst(c.exp());
            }
            Expr::Exp(Box::new(inner))
        }
        Expr::Ln(inner) => {
            let inner = simplify(inner);
            if let Expr::Const(c) = &inner
                && *c > 0.0
            {
                return cst(c.ln());
            }
            Expr::Ln(Box::new(inner))
        }
    }
}

// ---------------------------------------------------------------------------
// Pretty-print
// ---------------------------------------------------------------------------

/// Convert `expr` to a human-readable infix string.
pub fn to_string(expr: &Expr) -> String {
    expr_to_str(expr)
}

fn expr_to_str(expr: &Expr) -> String {
    match expr {
        Expr::Const(c) => {
            // Print integers without a decimal point for readability
            if c.fract() == 0.0 && c.abs() < 1e15 {
                format!("{}", *c as i64)
            } else {
                format!("{c}")
            }
        }
        Expr::Var(name) => name.clone(),
        Expr::Add(a, b) => {
            let bs = expr_to_str(b);
            // If b is a negation, display as subtraction
            if let Expr::Neg(inner) = b.as_ref() {
                format!("({} - {})", expr_to_str(a), expr_to_str(inner))
            } else if let Some(bs_stripped) = bs.strip_prefix('-') {
                format!("({} - {})", expr_to_str(a), bs_stripped)
            } else {
                format!("({} + {})", expr_to_str(a), bs)
            }
        }
        Expr::Mul(a, b) => format!("({} * {})", expr_to_str(a), expr_to_str(b)),
        Expr::Pow(base, exp) => format!("({}^{})", expr_to_str(base), expr_to_str(exp)),
        Expr::Neg(inner) => format!("(-{})", expr_to_str(inner)),
        Expr::Sin(inner) => format!("sin({})", expr_to_str(inner)),
        Expr::Cos(inner) => format!("cos({})", expr_to_str(inner)),
        Expr::Exp(inner) => format!("exp({})", expr_to_str(inner)),
        Expr::Ln(inner) => format!("ln({})", expr_to_str(inner)),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn vars(bindings: &[(&str, f64)]) -> HashMap<String, f64> {
        bindings.iter().map(|(k, v)| (k.to_string(), *v)).collect()
    }

    // --- eval ---

    #[test]
    fn eval_const() {
        assert_eq!(eval(&cst(3.125), &HashMap::new()).unwrap(), 3.125);
    }

    #[test]
    fn eval_var_found() {
        let e = var("x");
        assert_eq!(eval(&e, &vars(&[("x", 5.0)])).unwrap(), 5.0);
    }

    #[test]
    fn eval_var_missing_returns_err() {
        let e = var("y");
        assert!(eval(&e, &HashMap::new()).is_err());
    }

    #[test]
    fn eval_add() {
        let e = var("x").add_expr(cst(1.0));
        assert_eq!(eval(&e, &vars(&[("x", 4.0)])).unwrap(), 5.0);
    }

    #[test]
    fn eval_mul() {
        let e = var("x").mul_expr(cst(3.0));
        assert_eq!(eval(&e, &vars(&[("x", 2.0)])).unwrap(), 6.0);
    }

    #[test]
    fn eval_pow() {
        let e = var("x").pow(cst(3.0));
        assert!((eval(&e, &vars(&[("x", 2.0)])).unwrap() - 8.0).abs() < 1e-12);
    }

    #[test]
    fn eval_neg() {
        let e = -var("x");
        assert_eq!(eval(&e, &vars(&[("x", 7.0)])).unwrap(), -7.0);
    }

    #[test]
    fn eval_sin() {
        let e = var("x").sin();
        let got = eval(&e, &vars(&[("x", 0.0)])).unwrap();
        assert!(got.abs() < 1e-12);
    }

    #[test]
    fn eval_cos() {
        let e = var("x").cos();
        let got = eval(&e, &vars(&[("x", 0.0)])).unwrap();
        assert!((got - 1.0).abs() < 1e-12);
    }

    #[test]
    fn eval_exp() {
        let e = var("x").exp();
        let got = eval(&e, &vars(&[("x", 0.0)])).unwrap();
        assert!((got - 1.0).abs() < 1e-12);
    }

    #[test]
    fn eval_ln() {
        let e = var("x").ln();
        let got = eval(&e, &vars(&[("x", 1.0)])).unwrap();
        assert!(got.abs() < 1e-12);
    }

    #[test]
    fn eval_ln_nonpositive_returns_err() {
        let e = var("x").ln();
        assert!(eval(&e, &vars(&[("x", 0.0)])).is_err());
        assert!(eval(&e, &vars(&[("x", -1.0)])).is_err());
    }

    #[test]
    fn eval_complex_poly() {
        // 3x^2 + 2x + 1 at x=2 → 3*4 + 4 + 1 = 17
        let x = var("x");
        let e = cst(3.0)
            .mul_expr(x.clone().pow(cst(2.0)))
            .add_expr(cst(2.0).mul_expr(x.clone()))
            .add_expr(cst(1.0));
        let got = eval(&e, &vars(&[("x", 2.0)])).unwrap();
        assert!((got - 17.0).abs() < 1e-12);
    }

    // --- diff ---

    #[test]
    fn diff_const_is_zero() {
        let e = diff(&cst(42.0), "x");
        assert_eq!(simplify(&e), cst(0.0));
    }

    #[test]
    fn diff_var_self_is_one() {
        let e = diff(&var("x"), "x");
        assert_eq!(simplify(&e), cst(1.0));
    }

    #[test]
    fn diff_var_other_is_zero() {
        let e = diff(&var("y"), "x");
        assert_eq!(simplify(&e), cst(0.0));
    }

    #[test]
    fn diff_linear() {
        // d/dx (3x) = 3
        let e = cst(3.0).mul_expr(var("x"));
        let d = simplify(&diff(&e, "x"));
        let got = eval(&d, &HashMap::new()).unwrap();
        assert!((got - 3.0).abs() < 1e-12);
    }

    #[test]
    fn diff_quadratic() {
        // d/dx (x^2) = 2x  → at x=3: 6
        let e = var("x").pow(cst(2.0));
        let d = simplify(&diff(&e, "x"));
        let got = eval(&d, &vars(&[("x", 3.0)])).unwrap();
        assert!((got - 6.0).abs() < 1e-12);
    }

    #[test]
    fn diff_cubic() {
        // d/dx (x^3) = 3x^2 → at x=2: 12
        let e = var("x").pow(cst(3.0));
        let d = simplify(&diff(&e, "x"));
        let got = eval(&d, &vars(&[("x", 2.0)])).unwrap();
        assert!((got - 12.0).abs() < 1e-12);
    }

    #[test]
    fn diff_sin() {
        // d/dx sin(x) = cos(x) → at x=0: 1
        let e = var("x").sin();
        let d = simplify(&diff(&e, "x"));
        let got = eval(&d, &vars(&[("x", 0.0)])).unwrap();
        assert!((got - 1.0).abs() < 1e-12);
    }

    #[test]
    fn diff_cos() {
        // d/dx cos(x) = -sin(x) → at x=0: 0
        let e = var("x").cos();
        let d = simplify(&diff(&e, "x"));
        let got = eval(&d, &vars(&[("x", 0.0)])).unwrap();
        assert!(got.abs() < 1e-12);
    }

    #[test]
    fn diff_exp() {
        // d/dx exp(x) = exp(x) → at x=0: 1
        let e = var("x").exp();
        let d = simplify(&diff(&e, "x"));
        let got = eval(&d, &vars(&[("x", 0.0)])).unwrap();
        assert!((got - 1.0).abs() < 1e-12);
    }

    #[test]
    fn diff_ln() {
        // d/dx ln(x) = 1/x → at x=2: 0.5
        let e = var("x").ln();
        let d = simplify(&diff(&e, "x"));
        let got = eval(&d, &vars(&[("x", 2.0)])).unwrap();
        assert!((got - 0.5).abs() < 1e-12);
    }

    #[test]
    fn diff_product_rule() {
        // d/dx (x * sin(x)) = sin(x) + x*cos(x) → at x=0: 0
        let e = var("x").mul_expr(var("x").sin());
        let d = simplify(&diff(&e, "x"));
        let got = eval(&d, &vars(&[("x", 0.0)])).unwrap();
        assert!(got.abs() < 1e-12);
    }

    #[test]
    fn diff_chain_sin_of_poly() {
        // d/dx sin(x^2) = cos(x^2) * 2x → at x=0: 0
        let e = var("x").pow(cst(2.0)).sin();
        let d = simplify(&diff(&e, "x"));
        let got = eval(&d, &vars(&[("x", 0.0)])).unwrap();
        assert!(got.abs() < 1e-12);
    }

    #[test]
    fn diff_neg() {
        // d/dx (-x) = -1
        let e = -var("x");
        let d = simplify(&diff(&e, "x"));
        let got = eval(&d, &HashMap::new()).unwrap();
        assert!((got + 1.0).abs() < 1e-12);
    }

    // --- simplify ---

    #[test]
    fn simplify_zero_plus_x() {
        let e = cst(0.0).add_expr(var("x"));
        assert_eq!(simplify(&e), var("x"));
    }

    #[test]
    fn simplify_x_plus_zero() {
        let e = var("x").add_expr(cst(0.0));
        assert_eq!(simplify(&e), var("x"));
    }

    #[test]
    fn simplify_zero_times_x() {
        let e = cst(0.0).mul_expr(var("x"));
        assert_eq!(simplify(&e), cst(0.0));
    }

    #[test]
    fn simplify_x_times_zero() {
        let e = var("x").mul_expr(cst(0.0));
        assert_eq!(simplify(&e), cst(0.0));
    }

    #[test]
    fn simplify_one_times_x() {
        let e = cst(1.0).mul_expr(var("x"));
        assert_eq!(simplify(&e), var("x"));
    }

    #[test]
    fn simplify_x_times_one() {
        let e = var("x").mul_expr(cst(1.0));
        assert_eq!(simplify(&e), var("x"));
    }

    #[test]
    fn simplify_x_pow_zero() {
        let e = var("x").pow(cst(0.0));
        assert_eq!(simplify(&e), cst(1.0));
    }

    #[test]
    fn simplify_x_pow_one() {
        let e = var("x").pow(cst(1.0));
        assert_eq!(simplify(&e), var("x"));
    }

    #[test]
    fn simplify_double_neg() {
        let e = -(-var("x"));
        assert_eq!(simplify(&e), var("x"));
    }

    #[test]
    fn simplify_const_fold_add() {
        let e = cst(3.0).add_expr(cst(4.0));
        assert_eq!(simplify(&e), cst(7.0));
    }

    #[test]
    fn simplify_const_fold_mul() {
        let e = cst(3.0).mul_expr(cst(4.0));
        assert_eq!(simplify(&e), cst(12.0));
    }

    #[test]
    fn simplify_const_fold_pow() {
        let e = cst(2.0).pow(cst(10.0));
        assert_eq!(simplify(&e), cst(1024.0));
    }

    // --- to_string ---

    #[test]
    fn to_string_const() {
        assert_eq!(to_string(&cst(3.0)), "3");
    }

    #[test]
    fn to_string_var() {
        assert_eq!(to_string(&var("theta")), "theta");
    }

    #[test]
    fn to_string_add() {
        let e = var("x").add_expr(cst(1.0));
        let s = to_string(&e);
        assert!(s.contains("x") && s.contains("1") && s.contains("+"));
    }

    #[test]
    fn to_string_mul() {
        let e = var("a").mul_expr(var("b"));
        let s = to_string(&e);
        assert!(s.contains("a") && s.contains("b") && s.contains("*"));
    }

    #[test]
    fn to_string_pow() {
        let e = var("x").pow(cst(2.0));
        let s = to_string(&e);
        assert!(s.contains("x") && s.contains("2") && s.contains("^"));
    }

    #[test]
    fn to_string_sin() {
        let s = to_string(&var("x").sin());
        assert!(s.starts_with("sin("));
    }

    #[test]
    fn to_string_cos() {
        let s = to_string(&var("x").cos());
        assert!(s.starts_with("cos("));
    }

    #[test]
    fn to_string_exp() {
        let s = to_string(&var("x").exp());
        assert!(s.starts_with("exp("));
    }

    #[test]
    fn to_string_ln() {
        let s = to_string(&var("x").ln());
        assert!(s.starts_with("ln("));
    }

    #[test]
    fn to_string_neg() {
        let s = to_string(&(-var("x")));
        assert!(s.contains("x") && s.contains('-'));
    }

    // --- combined ---

    #[test]
    fn diff_poly_numeric_check() {
        // d/dx (x^4 - 3x^2 + 2) at x=1: 4*1 - 6*1 = -2
        let x = var("x");
        let poly = x
            .clone()
            .pow(cst(4.0))
            .sub_expr(cst(3.0).mul_expr(x.clone().pow(cst(2.0))))
            .add_expr(cst(2.0));
        let d = simplify(&diff(&poly, "x"));
        let got = eval(&d, &vars(&[("x", 1.0)])).unwrap();
        assert!((got - (-2.0)).abs() < 1e-10);
    }

    #[test]
    fn diff_exp_of_linear() {
        // d/dx exp(3x) = 3*exp(3x) → at x=0: 3
        let e = cst(3.0).mul_expr(var("x")).exp();
        let d = simplify(&diff(&e, "x"));
        let got = eval(&d, &vars(&[("x", 0.0)])).unwrap();
        assert!((got - 3.0).abs() < 1e-12);
    }
}
