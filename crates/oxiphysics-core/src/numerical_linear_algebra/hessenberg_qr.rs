//! Eigenvalues of a real upper-Hessenberg matrix by the Francis double-shift QR
//! algorithm.
//!
//! This is the classical EISPACK routine `hqr` (Martin, Peters & Wilkinson, 1970;
//! see also Numerical Recipes, section 11.6) for the eigenvalues only.  The
//! implicit double shift keeps the arithmetic real while converging to complex
//! conjugate pairs, a 2x2 block that splits off the active window is solved in
//! closed form, and exceptional shifts after 10 and 20 sweeps break the cycles
//! (for example cyclic permutation matrices) that ordinary shifts leave unchanged.
//!
//! The working matrix uses 1-based indices, exactly as in the reference
//! algorithm, so that the index arithmetic of the original can be followed
//! line by line.

/// Sweeps allowed per eigenvalue before the iteration gives up.
const MAX_SWEEPS_PER_EIGENVALUE: usize = 30;

/// `|a|` carrying the sign of `b` (positive for `b >= 0`).
fn sign(a: f64, b: f64) -> f64 {
    if b >= 0.0 { a.abs() } else { -a.abs() }
}

/// Eigenvalues `(re, im)` of the row-major upper-Hessenberg matrix `h`
/// (`n x n`), sorted by real part and then imaginary part.
///
/// Real eigenvalues have `im == 0`; a complex-conjugate pair is returned as
/// `(re, -im)` followed by `(re, +im)`.  Entries of `h` below the sub-diagonal
/// are ignored.  If an eigenvalue has not converged after 30 sweeps the remaining
/// diagonal of the unresolved window is returned as a real estimate instead of
/// iterating forever.
pub(super) fn francis_eigenvalues(h: &[f64], n: usize) -> Vec<(f64, f64)> {
    if n == 0 {
        return Vec::new();
    }
    let mut a = vec![vec![0.0_f64; n + 1]; n + 1];
    for i in 1..=n {
        for j in 1..=n {
            a[i][j] = h[(i - 1) * n + (j - 1)];
        }
    }
    let mut wr = vec![0.0_f64; n + 1];
    let mut wi = vec![0.0_f64; n + 1];

    // 1-norm-like size of the Hessenberg part, used as a deflation scale.
    let anorm: f64 = a
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, row)| {
            row.iter()
                .skip(i.saturating_sub(1).max(1))
                .map(|v| v.abs())
                .sum::<f64>()
        })
        .sum();
    let mut nn = n;
    let mut t = 0.0_f64;
    'eigenvalues: while nn >= 1 {
        let mut its = 0usize;
        loop {
            // Look for a single small sub-diagonal element.
            let mut l = nn;
            while l >= 2 {
                let mut s = a[l - 1][l - 1].abs() + a[l][l].abs();
                if s == 0.0 {
                    s = anorm;
                }
                if a[l][l - 1].abs() + s == s {
                    a[l][l - 1] = 0.0;
                    break;
                }
                l -= 1;
            }
            let mut x = a[nn][nn];
            if l == nn {
                // One root found.
                wr[nn] = x + t;
                wi[nn] = 0.0;
                nn -= 1;
                break;
            }
            let mut y = a[nn - 1][nn - 1];
            let mut w = a[nn][nn - 1] * a[nn - 1][nn];
            if l == nn - 1 {
                // Two roots found: solve the trailing 2x2 block in closed form.
                let p = 0.5 * (y - x);
                let q = p * p + w;
                let z = q.abs().sqrt();
                x += t;
                if q >= 0.0 {
                    // A real pair.
                    let z = p + sign(z, p);
                    wr[nn - 1] = x + z;
                    wr[nn] = if z != 0.0 { x - w / z } else { x + z };
                    wi[nn - 1] = 0.0;
                    wi[nn] = 0.0;
                } else {
                    // A complex-conjugate pair.
                    wr[nn - 1] = x + p;
                    wr[nn] = x + p;
                    wi[nn - 1] = -z;
                    wi[nn] = z;
                }
                nn -= 2;
                break;
            }
            if its == MAX_SWEEPS_PER_EIGENVALUE {
                // No convergence: report the diagonal of the unresolved window.
                for i in 1..=nn {
                    wr[i] = a[i][i] + t;
                    wi[i] = 0.0;
                }
                break 'eigenvalues;
            }
            if its == 10 || its == 20 {
                // Exceptional shift.
                t += x;
                for (i, row) in a.iter_mut().enumerate().take(nn + 1).skip(1) {
                    row[i] -= x;
                }
                let s = a[nn][nn - 1].abs() + a[nn - 1][nn - 2].abs();
                x = 0.75 * s;
                y = x;
                w = -0.4375 * s * s;
            }
            its += 1;

            // Form the first column of (H - s1 I)(H - s2 I) from the rows m, m+1,
            // m+2 and stop as soon as two consecutive sub-diagonals are small.
            let mut m = nn - 2;
            let (mut p, mut q, mut r);
            loop {
                let z = a[m][m];
                let shift_r = x - z;
                let shift_s = y - z;
                p = (shift_r * shift_s - w) / a[m + 1][m] + a[m][m + 1];
                q = a[m + 1][m + 1] - z - shift_r - shift_s;
                r = a[m + 2][m + 1];
                let scale = p.abs() + q.abs() + r.abs();
                if scale != 0.0 {
                    p /= scale;
                    q /= scale;
                    r /= scale;
                }
                if m == l {
                    break;
                }
                let u = a[m][m - 1].abs() * (q.abs() + r.abs());
                let v = p.abs() * (a[m - 1][m - 1].abs() + z.abs() + a[m + 1][m + 1].abs());
                if u + v == v {
                    break;
                }
                m -= 1;
            }
            for i in (m + 2)..=nn {
                a[i][i - 2] = 0.0;
                if i != m + 2 {
                    a[i][i - 3] = 0.0;
                }
            }
            // Double QR step on rows l..nn and columns m..nn.
            for k in m..nn {
                if k != m {
                    p = a[k][k - 1];
                    q = a[k + 1][k - 1];
                    r = if k != nn - 1 { a[k + 2][k - 1] } else { 0.0 };
                    x = p.abs() + q.abs() + r.abs();
                    if x != 0.0 {
                        p /= x;
                        q /= x;
                        r /= x;
                    }
                }
                let s = sign((p * p + q * q + r * r).sqrt(), p);
                if s == 0.0 {
                    continue;
                }
                if k == m {
                    if l != m {
                        a[k][k - 1] = -a[k][k - 1];
                    }
                } else {
                    a[k][k - 1] = -s * x;
                }
                p += s;
                x = p / s;
                y = q / s;
                let z = r / s;
                q /= p;
                r /= p;
                // Row modification of rows k, k+1 (and k+2) over columns k..=nn.
                let (upper, lower) = a.split_at_mut(k + 1);
                let row_k = &mut upper[k][k..=nn];
                let (first_lower, second_lower) = lower.split_at_mut(1);
                let row_k1 = &mut first_lower[0][k..=nn];
                if k != nn - 1 {
                    let row_k2 = &mut second_lower[0][k..=nn];
                    for ((v0, v1), v2) in row_k.iter_mut().zip(row_k1.iter_mut()).zip(row_k2) {
                        let sum = *v0 + q * *v1 + r * *v2;
                        *v2 -= sum * z;
                        *v1 -= sum * y;
                        *v0 -= sum * x;
                    }
                } else {
                    for (v0, v1) in row_k.iter_mut().zip(row_k1.iter_mut()) {
                        let sum = *v0 + q * *v1;
                        *v1 -= sum * y;
                        *v0 -= sum * x;
                    }
                }
                // Column modification of columns k, k+1 (and k+2) over rows l..=last_row.
                let last_row = nn.min(k + 3);
                for row in a.iter_mut().take(last_row + 1).skip(l) {
                    let mut sum = x * row[k] + y * row[k + 1];
                    if k != nn - 1 {
                        sum += z * row[k + 2];
                        row[k + 2] -= sum * r;
                    }
                    row[k + 1] -= sum * q;
                    row[k] -= sum;
                }
            }
        }
    }
    let mut eigenvalues: Vec<(f64, f64)> = (1..=n).map(|i| (wr[i], wi[i])).collect();
    eigenvalues.sort_by(|u, v| u.0.total_cmp(&v.0).then(u.1.total_cmp(&v.1)));
    eigenvalues
}
