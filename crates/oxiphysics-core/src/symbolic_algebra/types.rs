//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
/// A single term in a partial fraction decomposition: `coeff / (x - root)`.
#[derive(Debug, Clone)]
pub struct PartialFractionTerm {
    /// The coefficient (residue) A_i.
    pub coefficient: f64,
    /// The pole location r_i.
    pub root: f64,
}
/// A rational function P(x) / Q(x) represented by two polynomials.
#[derive(Debug, Clone)]
pub struct RationalFunction {
    /// Numerator polynomial.
    pub numerator: Polynomial,
    /// Denominator polynomial.
    pub denominator: Polynomial,
}
impl RationalFunction {
    /// Create a new rational function.
    pub fn new(num: Polynomial, den: Polynomial) -> Self {
        RationalFunction {
            numerator: num,
            denominator: den,
        }
    }
    /// Evaluate the rational function at x.
    pub fn eval(&self, x: f64) -> Option<f64> {
        let d = self.denominator.eval(x);
        if d.abs() < 1e-300 {
            None
        } else {
            Some(self.numerator.eval(x) / d)
        }
    }
    /// Add two rational functions: a/b + c/d = (ad + bc) / (bd).
    pub fn add(&self, other: &RationalFunction) -> RationalFunction {
        let num = self
            .numerator
            .mul(&other.denominator)
            .add(&other.numerator.mul(&self.denominator));
        let den = self.denominator.mul(&other.denominator);
        RationalFunction::new(num, den)
    }
    /// Multiply two rational functions: (a/b) * (c/d) = (ac) / (bd).
    pub fn mul(&self, other: &RationalFunction) -> RationalFunction {
        let num = self.numerator.mul(&other.numerator);
        let den = self.denominator.mul(&other.denominator);
        RationalFunction::new(num, den)
    }
    /// Convert to an `Expr`.
    pub fn to_expr(&self, var_name: &str) -> Expr {
        self.numerator
            .to_expr(var_name)
            .div_expr(self.denominator.to_expr(var_name))
    }
}
/// A univariate polynomial stored as a coefficient vector (index = degree).
///
/// `coeffs[i]` is the coefficient of x^i.  The vector may have trailing
/// zeros after arithmetic; call [`Polynomial::trim`] to normalise.
#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial {
    /// Coefficient vector: `coeffs[i]` is the coefficient of x^i.
    pub coeffs: Vec<f64>,
}
impl Polynomial {
    /// Create a polynomial from coefficients (lowest degree first).
    pub fn new(coeffs: Vec<f64>) -> Self {
        let mut p = Polynomial { coeffs };
        p.trim();
        p
    }
    /// Remove trailing zero coefficients.
    pub fn trim(&mut self) {
        while self.coeffs.len() > 1 && self.coeffs.last().is_some_and(|c| c.abs() < 1e-15) {
            self.coeffs.pop();
        }
    }
    /// Degree of the polynomial (0 for the zero polynomial).
    pub fn degree(&self) -> usize {
        if self.coeffs.is_empty() {
            0
        } else {
            self.coeffs.len().saturating_sub(1)
        }
    }
    /// Evaluate p(x) using Horner's method.
    pub fn eval(&self, x: f64) -> f64 {
        let mut result = 0.0;
        for c in self.coeffs.iter().rev() {
            result = result * x + c;
        }
        result
    }
    /// Return the zero polynomial.
    pub fn zero() -> Self {
        Polynomial { coeffs: vec![0.0] }
    }
    /// Return the constant polynomial `c`.
    pub fn constant(c: f64) -> Self {
        Polynomial { coeffs: vec![c] }
    }
    /// Return the monomial `x` (degree 1, leading coefficient 1).
    pub fn x() -> Self {
        Polynomial {
            coeffs: vec![0.0, 1.0],
        }
    }
    /// Formal derivative p'(x).
    pub fn derivative(&self) -> Self {
        if self.coeffs.len() <= 1 {
            return Polynomial::zero();
        }
        let c: Vec<f64> = self
            .coeffs
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, &c)| c * i as f64)
            .collect();
        Polynomial::new(c)
    }
    /// Formal integral with constant of integration = 0.
    pub fn integral(&self) -> Self {
        let mut c = vec![0.0];
        for (i, &coeff) in self.coeffs.iter().enumerate() {
            c.push(coeff / (i as f64 + 1.0));
        }
        Polynomial::new(c)
    }
    /// Add two polynomials.
    pub fn add(&self, other: &Polynomial) -> Polynomial {
        let n = self.coeffs.len().max(other.coeffs.len());
        let mut c = vec![0.0; n];
        for (i, &v) in self.coeffs.iter().enumerate() {
            c[i] += v;
        }
        for (i, &v) in other.coeffs.iter().enumerate() {
            c[i] += v;
        }
        Polynomial::new(c)
    }
    /// Subtract: self - other.
    pub fn sub(&self, other: &Polynomial) -> Polynomial {
        let n = self.coeffs.len().max(other.coeffs.len());
        let mut c = vec![0.0; n];
        for (i, &v) in self.coeffs.iter().enumerate() {
            c[i] += v;
        }
        for (i, &v) in other.coeffs.iter().enumerate() {
            c[i] -= v;
        }
        Polynomial::new(c)
    }
    /// Multiply two polynomials.
    pub fn mul(&self, other: &Polynomial) -> Polynomial {
        if self.coeffs.is_empty() || other.coeffs.is_empty() {
            return Polynomial::zero();
        }
        let n = self.coeffs.len() + other.coeffs.len() - 1;
        let mut c = vec![0.0; n];
        for (i, &a) in self.coeffs.iter().enumerate() {
            for (j, &b) in other.coeffs.iter().enumerate() {
                c[i + j] += a * b;
            }
        }
        Polynomial::new(c)
    }
    /// Scalar multiplication.
    pub fn scale(&self, s: f64) -> Polynomial {
        Polynomial::new(self.coeffs.iter().map(|c| c * s).collect())
    }
    /// Polynomial long division: returns `(quotient, remainder)`.
    pub fn div_rem(&self, divisor: &Polynomial) -> (Polynomial, Polynomial) {
        assert!(
            !(divisor.coeffs.is_empty()
                || divisor.coeffs.len() == 1 && divisor.coeffs[0].abs() < 1e-15),
            "Division by zero polynomial"
        );
        if self.degree() < divisor.degree() {
            return (Polynomial::zero(), self.clone());
        }
        let mut rem = self.coeffs.clone();
        let dlen = divisor.coeffs.len();
        let lead = *divisor.coeffs.last().expect("divisor has non-zero degree");
        let mut quot = vec![0.0; rem.len() - dlen + 1];
        for i in (0..quot.len()).rev() {
            let coeff = rem[i + dlen - 1] / lead;
            quot[i] = coeff;
            for (j, &d) in divisor.coeffs.iter().enumerate() {
                rem[i + j] -= coeff * d;
            }
        }
        (Polynomial::new(quot), Polynomial::new(rem))
    }
    /// GCD of two polynomials via the Euclidean algorithm.
    pub fn gcd(&self, other: &Polynomial) -> Polynomial {
        let mut a = self.clone();
        let mut b = other.clone();
        for _ in 0..200 {
            b.trim();
            if b.coeffs.len() == 1 && b.coeffs[0].abs() < 1e-12 {
                break;
            }
            let (_q, r) = a.div_rem(&b);
            a = b;
            b = r;
        }
        a.trim();
        let lead = *a.coeffs.last().unwrap_or(&1.0);
        if lead.abs() > 1e-15 {
            a = a.scale(1.0 / lead);
        }
        a
    }
    /// Convert to an `Expr` tree using the given variable name.
    pub fn to_expr(&self, var_name: &str) -> Expr {
        let v = Expr::Var(var_name.to_string());
        let mut terms: Vec<Expr> = Vec::new();
        for (i, &c) in self.coeffs.iter().enumerate() {
            if c.abs() < 1e-15 {
                continue;
            }
            let term = if i == 0 {
                cst(c)
            } else if i == 1 {
                cst(c).mul_expr(v.clone())
            } else {
                cst(c).mul_expr(v.clone().pow_expr(cst(i as f64)))
            };
            terms.push(term);
        }
        if terms.is_empty() {
            cst(0.0)
        } else {
            let mut acc = terms.remove(0);
            for t in terms {
                acc = acc.add_expr(t);
            }
            acc
        }
    }
}
impl Polynomial {
    /// Compose: compute p(q(x)).
    pub fn compose(&self, q: &Polynomial) -> Polynomial {
        let mut result = Polynomial::zero();
        let mut q_power = Polynomial::constant(1.0);
        for &c in &self.coeffs {
            result = result.add(&q_power.scale(c));
            q_power = q_power.mul(q);
        }
        result
    }
    /// Return the leading coefficient.
    pub fn leading_coeff(&self) -> f64 {
        *self.coeffs.last().unwrap_or(&0.0)
    }
    /// Check if this is the zero polynomial.
    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|c| c.abs() < 1e-15)
    }
    /// Return the number of sign changes (for Descartes' rule).
    pub fn sign_changes(&self) -> usize {
        let nonzero: Vec<f64> = self
            .coeffs
            .iter()
            .copied()
            .filter(|c| c.abs() > 1e-15)
            .collect();
        let mut count = 0;
        for w in nonzero.windows(2) {
            if w[0] * w[1] < 0.0 {
                count += 1;
            }
        }
        count
    }
}
/// Result of common subexpression elimination.
#[derive(Debug, Clone)]
pub struct CseResult {
    /// List of `(name, expression)` for each extracted subexpression.
    pub bindings: Vec<(String, Expr)>,
    /// The rewritten top-level expression (uses the introduced names).
    pub reduced: Expr,
}
/// A symbolic expression node.
///
/// This is an algebraic expression tree supporting arithmetic, elementary
/// transcendental functions, and power expressions.
#[derive(Debug, Clone, PartialEq)]
pub enum Expr {
    /// Numeric constant.
    Const(f64),
    /// Named variable.
    Var(String),
    /// `lhs + rhs`.
    Add(Box<Expr>, Box<Expr>),
    /// `lhs * rhs`.
    Mul(Box<Expr>, Box<Expr>),
    /// `base ^ exponent`.
    Pow(Box<Expr>, Box<Expr>),
    /// Unary negation `-e`.
    Neg(Box<Expr>),
    /// `sin(e)`.
    Sin(Box<Expr>),
    /// `cos(e)`.
    Cos(Box<Expr>),
    /// `exp(e)`.
    Exp(Box<Expr>),
    /// `ln(e)` (natural logarithm).
    Ln(Box<Expr>),
    /// `lhs / rhs` (division kept explicit for partial-fraction work).
    Div(Box<Expr>, Box<Expr>),
}
impl Expr {
    /// `self + rhs`
    pub fn add_expr(self, rhs: Expr) -> Expr {
        Expr::Add(Box::new(self), Box::new(rhs))
    }
    /// `self - rhs`
    pub fn sub_expr(self, rhs: Expr) -> Expr {
        Expr::Add(Box::new(self), Box::new(Expr::Neg(Box::new(rhs))))
    }
    /// `self * rhs`
    pub fn mul_expr(self, rhs: Expr) -> Expr {
        Expr::Mul(Box::new(self), Box::new(rhs))
    }
    /// `self / rhs`
    pub fn div_expr(self, rhs: Expr) -> Expr {
        Expr::Div(Box::new(self), Box::new(rhs))
    }
    /// `self ^ rhs`
    pub fn pow_expr(self, rhs: Expr) -> Expr {
        Expr::Pow(Box::new(self), Box::new(rhs))
    }
    /// `sin(self)`
    pub fn sin_expr(self) -> Expr {
        Expr::Sin(Box::new(self))
    }
    /// `cos(self)`
    pub fn cos_expr(self) -> Expr {
        Expr::Cos(Box::new(self))
    }
    /// `exp(self)`
    pub fn exp_expr(self) -> Expr {
        Expr::Exp(Box::new(self))
    }
    /// `ln(self)`
    pub fn ln_expr(self) -> Expr {
        Expr::Ln(Box::new(self))
    }
    /// Negate `self`.
    pub fn neg_expr(self) -> Expr {
        Expr::Neg(Box::new(self))
    }
    /// Return true if the expression is the constant zero.
    pub fn is_zero(&self) -> bool {
        matches!(self, Expr::Const(v) if * v == 0.0)
    }
    /// Return true if the expression is the constant one.
    pub fn is_one(&self) -> bool {
        matches!(self, Expr::Const(v) if (* v - 1.0).abs() < 1e-15)
    }
    /// Collect all variable names appearing in the expression.
    pub fn variables(&self) -> Vec<String> {
        let mut out = Vec::new();
        self.collect_vars(&mut out);
        out.sort();
        out.dedup();
        out
    }
    /// Helper to recursively collect variable names.
    fn collect_vars(&self, out: &mut Vec<String>) {
        match self {
            Expr::Const(_) => {}
            Expr::Var(s) => out.push(s.clone()),
            Expr::Add(a, b) | Expr::Mul(a, b) | Expr::Pow(a, b) | Expr::Div(a, b) => {
                a.collect_vars(out);
                b.collect_vars(out);
            }
            Expr::Neg(a) | Expr::Sin(a) | Expr::Cos(a) | Expr::Exp(a) | Expr::Ln(a) => {
                a.collect_vars(out);
            }
        }
    }
    /// Count total number of nodes in the expression tree.
    pub fn node_count(&self) -> usize {
        match self {
            Expr::Const(_) | Expr::Var(_) => 1,
            Expr::Add(a, b) | Expr::Mul(a, b) | Expr::Pow(a, b) | Expr::Div(a, b) => {
                1 + a.node_count() + b.node_count()
            }
            Expr::Neg(a) | Expr::Sin(a) | Expr::Cos(a) | Expr::Exp(a) | Expr::Ln(a) => {
                1 + a.node_count()
            }
        }
    }
    /// Maximum depth of the expression tree.
    pub fn depth(&self) -> usize {
        match self {
            Expr::Const(_) | Expr::Var(_) => 1,
            Expr::Add(a, b) | Expr::Mul(a, b) | Expr::Pow(a, b) | Expr::Div(a, b) => {
                1 + a.depth().max(b.depth())
            }
            Expr::Neg(a) | Expr::Sin(a) | Expr::Cos(a) | Expr::Exp(a) | Expr::Ln(a) => {
                1 + a.depth()
            }
        }
    }
}
