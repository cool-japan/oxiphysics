//! Eigenvalues of a general real matrix.
//!
//! The matrix is balanced (Parlett–Reinsch), reduced to upper-Hessenberg form
//! by Householder reflections and handed to the Francis double-shift QR
//! iteration (the classical EISPACK routine `hqr`; Martin, Peters & Wilkinson,
//! 1970).  The implicit double shift keeps the arithmetic real while the
//! iteration converges to complex-conjugate pairs, and exceptional shifts after
//! 10 and 20 sweeps break the cycles that ordinary shifts leave unchanged.

/// Sweeps allowed per eigenvalue before the iteration gives up.
const MAX_SWEEPS_PER_EIGENVALUE: usize = 40;

/// Scale rows and columns of `a` (row-major `n x n`) by powers of two so that
/// their off-diagonal norms are comparable.  This is a similarity transform, so
/// the eigenvalues are unchanged, but the QR iteration becomes far more
/// accurate for badly scaled matrices.  The entries must be finite.
fn balance(a: &mut [f64], n: usize) {
    const RADIX: f64 = 2.0;
    const RADIX_SQ: f64 = RADIX * RADIX;
    let mut balanced = false;
    while !balanced {
        balanced = true;
        for i in 0..n {
            let mut col_norm = 0.0;
            let mut row_norm = 0.0;
            for j in 0..n {
                if j != i {
                    col_norm += a[j * n + i].abs();
                    row_norm += a[i * n + j].abs();
                }
            }
            if col_norm == 0.0 || row_norm == 0.0 {
                continue;
            }
            let total = col_norm + row_norm;
            let mut factor = 1.0;
            let mut c = col_norm;
            let mut g = row_norm / RADIX;
            while c < g {
                factor *= RADIX;
                c *= RADIX_SQ;
            }
            g = row_norm * RADIX;
            while c > g {
                factor /= RADIX;
                c /= RADIX_SQ;
            }
            if (c + row_norm) / factor < 0.95 * total {
                balanced = false;
                let inv = 1.0 / factor;
                for j in 0..n {
                    a[i * n + j] *= inv;
                }
                for j in 0..n {
                    a[j * n + i] *= factor;
                }
            }
        }
    }
}

/// Reduce `a` (row-major `n x n`) to upper-Hessenberg form by an orthogonal
/// similarity transform built from Householder reflections.
fn reduce_to_hessenberg(a: &mut [f64], n: usize) {
    if n < 3 {
        return;
    }
    let mut v = vec![0.0; n];
    for k in 0..(n - 2) {
        let len = n - k - 1;
        let tail_norm = (0..len)
            .map(|i| a[(k + 1 + i) * n + k].powi(2))
            .sum::<f64>()
            .sqrt();
        if tail_norm == 0.0 {
            continue;
        }
        let head = a[(k + 1) * n + k];
        let alpha = if head > 0.0 { -tail_norm } else { tail_norm };
        for i in 0..len {
            v[i] = a[(k + 1 + i) * n + k];
        }
        v[0] -= alpha;
        let v_norm_sq: f64 = v[..len].iter().map(|x| x * x).sum();
        if v_norm_sq == 0.0 {
            continue;
        }
        // Left reflection: rows k+1..n, columns k..n.
        for j in k..n {
            let mut s = 0.0;
            for i in 0..len {
                s += v[i] * a[(k + 1 + i) * n + j];
            }
            let f = 2.0 * s / v_norm_sq;
            for i in 0..len {
                a[(k + 1 + i) * n + j] -= f * v[i];
            }
        }
        // Right reflection: all rows, columns k+1..n.
        for i in 0..n {
            let mut s = 0.0;
            for j in 0..len {
                s += a[i * n + k + 1 + j] * v[j];
            }
            let f = 2.0 * s / v_norm_sq;
            for j in 0..len {
                a[i * n + k + 1 + j] -= f * v[j];
            }
        }
        for i in (k + 2)..n {
            a[i * n + k] = 0.0;
        }
    }
}

/// `|a|` carrying the sign of `b` (positive for `b >= 0`).
fn sign(a: f64, b: f64) -> f64 {
    if b >= 0.0 { a.abs() } else { -a.abs() }
}

/// Eigenvalues `(re, im)` of the row-major upper-Hessenberg matrix `h`
/// (`n x n`) by the Francis double-shift QR iteration; `None` if an eigenvalue
/// does not converge.  A complex-conjugate pair is returned as `(re, -im)`
/// followed by `(re, +im)`.  The working copy uses 1-based indices exactly as
/// the reference algorithm does, so the index arithmetic can be followed line by
/// line.
fn francis_qr(h: &[f64], n: usize) -> Option<Vec<(f64, f64)>> {
    let mut a = vec![vec![0.0_f64; n + 1]; n + 1];
    for (i, row) in a.iter_mut().enumerate().skip(1) {
        row[1..=n].copy_from_slice(&h[(i - 1) * n..i * n]);
    }
    let mut wr = vec![0.0_f64; n + 1];
    let mut wi = vec![0.0_f64; n + 1];
    let mut anorm = 0.0;
    for (i, row) in a.iter().enumerate().skip(1) {
        let first = i.saturating_sub(1).max(1);
        anorm += row[first..=n].iter().map(|v| v.abs()).sum::<f64>();
    }
    let mut nn = n;
    let mut t = 0.0_f64;
    while nn >= 1 {
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
                    let z = p + sign(z, p);
                    wr[nn - 1] = x + z;
                    wr[nn] = if z != 0.0 { x - w / z } else { x + z };
                    wi[nn - 1] = 0.0;
                    wi[nn] = 0.0;
                } else {
                    wr[nn - 1] = x + p;
                    wr[nn] = x + p;
                    wi[nn - 1] = -z;
                    wi[nn] = z;
                }
                nn -= 2;
                break;
            }
            if its >= MAX_SWEEPS_PER_EIGENVALUE {
                return None;
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
            // Form the first column of (H - s1 I)(H - s2 I) from rows m, m+1,
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
                let (above, below) = a.split_at_mut(k + 1);
                let row_k = &mut above[k][k..=nn];
                let (next, rest) = below.split_at_mut(1);
                let row_k1 = &mut next[0][k..=nn];
                if k != nn - 1 {
                    let row_k2 = &mut rest[0][k..=nn];
                    for ((a_k, a_k1), a_k2) in row_k.iter_mut().zip(row_k1).zip(row_k2) {
                        let sum = *a_k + q * *a_k1 + r * *a_k2;
                        *a_k2 -= sum * z;
                        *a_k1 -= sum * y;
                        *a_k -= sum * x;
                    }
                } else {
                    for (a_k, a_k1) in row_k.iter_mut().zip(row_k1) {
                        let sum = *a_k + q * *a_k1;
                        *a_k1 -= sum * y;
                        *a_k -= sum * x;
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
    Some(eigenvalues)
}

/// Eigenvalues `(re, im)` of the general real row-major `n x n` matrix `a`,
/// sorted by real part and then imaginary part.
///
/// A complex-conjugate pair appears as `(re, -im)` followed by `(re, +im)`.
/// Returns `None` when `a` has the wrong length or non-finite entries, or when
/// the QR iteration fails to converge (which is very rare).
pub(crate) fn real_eigenvalues(a: &[f64], n: usize) -> Option<Vec<(f64, f64)>> {
    if n == 0 {
        return Some(Vec::new());
    }
    if a.len() != n * n || a.iter().any(|v| !v.is_finite()) {
        return None;
    }
    let mut work = a.to_vec();
    balance(&mut work, n);
    reduce_to_hessenberg(&mut work, n);
    francis_qr(&work, n)
}
