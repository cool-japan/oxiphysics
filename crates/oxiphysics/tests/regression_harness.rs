// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// Regression harness for Phase 21 validation tests.
// Shared across validation_rigid.rs, validation_sph.rs, validation_lbm.rs,
// validation_fem.rs, validation_md_lj.rs, validation_ewald_madelung.rs
// via `#[path = "regression_harness.rs"] mod harness;`.

#![allow(missing_docs)]

use std::fs;
use std::io::{self, Read, Write};
use std::path::PathBuf;

/// Directory holding per-test baseline JSON files.
pub const BASELINE_DIR: &str = "crates/oxiphysics/tests/regression_baselines";

/// Scalar regression baseline loaded from a JSON file.
#[derive(Debug, Clone)]
pub struct Baseline {
    pub name: String,
    pub expected: f64,
    pub tolerance_abs: f64,
    pub tolerance_rel: f64,
}

impl Baseline {
    /// Check whether `actual` passes the tolerance gate.
    pub fn is_close(&self, actual: f64) -> bool {
        let diff = (actual - self.expected).abs();
        let tol = self.tolerance_abs + self.tolerance_rel * self.expected.abs();
        diff <= tol
    }

    /// Rich diagnostic message on mismatch.
    pub fn diagnostic(&self, actual: f64) -> String {
        let diff = actual - self.expected;
        let abs_diff = diff.abs();
        let rel_diff = if self.expected.abs() > 1e-30 {
            abs_diff / self.expected.abs()
        } else {
            f64::INFINITY
        };
        let tol = self.tolerance_abs + self.tolerance_rel * self.expected.abs();
        format!(
            "baseline `{}` FAILED: measured = {:.9e}, expected = {:.9e}, \
             abs_diff = {:.3e}, rel_diff = {:.3e}, tolerance = {:.3e} \
             (tol_abs = {:.3e}, tol_rel = {:.3e})",
            self.name,
            actual,
            self.expected,
            abs_diff,
            rel_diff,
            tol,
            self.tolerance_abs,
            self.tolerance_rel,
        )
    }
}

/// Resolve the project root by walking up from CARGO_MANIFEST_DIR
/// of the `oxiphysics` crate. Returns the workspace root.
pub fn workspace_root() -> PathBuf {
    // CARGO_MANIFEST_DIR for integration tests points at the crate dir:
    //   /media/.../oxiphysics/crates/oxiphysics
    // The workspace root is two levels up.
    let crate_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    crate_dir
        .parent()
        .and_then(|p| p.parent())
        .map(|p| p.to_path_buf())
        .unwrap_or(crate_dir)
}

/// Load a scalar baseline from `{workspace_root}/{BASELINE_DIR}/{name}.json`.
pub fn load_baseline(name: &str) -> io::Result<Baseline> {
    let path = workspace_root()
        .join(BASELINE_DIR)
        .join(format!("{name}.json"));
    let mut f = fs::File::open(&path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!(
                "failed to open baseline `{name}` at {}: {e}",
                path.display()
            ),
        )
    })?;
    let mut src = String::new();
    f.read_to_string(&mut src)?;
    parse_scalar_baseline(name, &src).map_err(|msg| io::Error::new(io::ErrorKind::InvalidData, msg))
}

/// Parse a flat JSON object with keys `expected`, `tolerance_abs`, `tolerance_rel`.
/// Hand-rolled to avoid adding `serde`/`serde_json` deps.
pub fn parse_scalar_baseline(name: &str, src: &str) -> Result<Baseline, String> {
    let expected = extract_number(src, "expected")?;
    let tolerance_abs = extract_number(src, "tolerance_abs").unwrap_or(0.0);
    let tolerance_rel = extract_number(src, "tolerance_rel").unwrap_or(0.0);
    Ok(Baseline {
        name: name.to_string(),
        expected,
        tolerance_abs,
        tolerance_rel,
    })
}

/// Extract a numeric JSON field `"key": N.NNe+EE`. Ignores whitespace.
pub fn extract_number(src: &str, key: &str) -> Result<f64, String> {
    let needle = format!("\"{key}\"");
    let idx = src
        .find(&needle)
        .ok_or_else(|| format!("baseline missing required key `{key}`"))?;
    let rest = &src[idx + needle.len()..];
    // skip ':' and whitespace
    let after_colon = rest
        .find(':')
        .ok_or_else(|| format!("key `{key}` has no colon"))?;
    let tail = &rest[after_colon + 1..];
    let tail = tail.trim_start();
    // Read until comma, }, or whitespace
    let end = tail.find([',', '}', '\n', '\r']).unwrap_or(tail.len());
    let token = tail[..end].trim();
    token
        .parse::<f64>()
        .map_err(|e| format!("key `{key}` value `{token}` is not a number: {e}"))
}

/// Write the measured value back to the baseline JSON (used when
/// `OXI_UPDATE_BASELINES=1` is set). Keeps `tolerance_abs` / `tolerance_rel`.
pub fn update_baseline(baseline: &Baseline, new_expected: f64) -> io::Result<()> {
    let path = workspace_root()
        .join(BASELINE_DIR)
        .join(format!("{}.json", baseline.name));
    let body = format!(
        "{{\n  \"name\": \"{}\",\n  \"expected\": {:.17e},\n  \"tolerance_abs\": {:.17e},\n  \"tolerance_rel\": {:.17e}\n}}\n",
        baseline.name, new_expected, baseline.tolerance_abs, baseline.tolerance_rel,
    );
    let mut f = fs::File::create(&path)?;
    f.write_all(body.as_bytes())?;
    Ok(())
}

/// Returns true when the `OXI_UPDATE_BASELINES` env var is set to "1".
pub fn update_mode() -> bool {
    std::env::var("OXI_UPDATE_BASELINES").ok().as_deref() == Some("1")
}

/// Assert a scalar measurement is close to a `Baseline`, or (in update mode) rewrite the baseline.
#[macro_export]
macro_rules! assert_close {
    ($actual:expr, $baseline:expr) => {{
        let actual: f64 = $actual;
        let baseline: &$crate::harness::Baseline = &$baseline;
        if $crate::harness::update_mode() {
            $crate::harness::update_baseline(baseline, actual)
                .expect("OXI_UPDATE_BASELINES: rewrite failed");
            eprintln!(
                "[OXI_UPDATE_BASELINES] rewrote `{}` to {:.9e}",
                baseline.name, actual,
            );
        } else if !baseline.is_close(actual) {
            panic!("{}", baseline.diagnostic(actual));
        }
    }};
}

// ---------------------------------------------------------------
// Self-tests (run when regression_harness.rs is compiled as its own
// integration-test binary; may also run when included via #[path] from
// another validation_*.rs, which is harmless).
// ---------------------------------------------------------------

#[cfg(test)]
mod self_tests {
    use super::*;

    fn test_baseline() -> Baseline {
        Baseline {
            name: "self_test".to_string(),
            expected: 1.0,
            tolerance_abs: 0.0,
            tolerance_rel: 0.05,
        }
    }

    #[test]
    fn harness_accepts_close_values() {
        let b = test_baseline();
        assert!(b.is_close(1.04));
        assert!(b.is_close(0.96));
        assert!(b.is_close(1.0));
    }

    #[test]
    fn harness_rejects_far_values() {
        let b = test_baseline();
        assert!(!b.is_close(1.06));
        assert!(!b.is_close(0.94));
        assert!(!b.is_close(-1.0));
    }

    #[test]
    fn harness_update_env_writes_json() {
        // Exercise parse_scalar_baseline() round-trip on a temp file.
        let tmp = std::env::temp_dir().join("oxi_harness_self_test.json");
        // Seed a baseline with known content.
        let body =
            b"{\n  \"expected\": 1.0,\n  \"tolerance_abs\": 0.0,\n  \"tolerance_rel\": 0.05\n}\n";
        fs::write(&tmp, body).expect("write seed");
        // Parse round-trip via the hand-rolled parser.
        let src = fs::read_to_string(&tmp).expect("read seed");
        let parsed =
            parse_scalar_baseline("oxi_harness_self_test", &src).expect("parse seed must succeed");
        assert_eq!(parsed.name, "oxi_harness_self_test");
        assert!((parsed.expected - 1.0).abs() < 1e-12);
        assert!((parsed.tolerance_rel - 0.05).abs() < 1e-12);
        // Cleanup (best-effort).
        let _ = fs::remove_file(&tmp);
    }
}
