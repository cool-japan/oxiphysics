# Contributing to OxiPhysics

Thank you for your interest in contributing to **OxiPhysics** — a unified,
pure-Rust physics engine targeting the same problem domains as Bullet (rigid
body), OpenFOAM (CFD), LAMMPS (molecular dynamics), and CalculiX (FEM).

OxiPhysics is part of the [COOLJAPAN](https://github.com/cool-japan) ecosystem.
The repository lives at <https://github.com/cool-japan/oxiphysics>.

The project is licensed under **Apache-2.0** (© 2026 COOLJAPAN OU / Team
KitaSan). By submitting a contribution you agree that it will be distributed
under those same terms.

Please also read our [Code of Conduct](CODE_OF_CONDUCT.md) — it applies to all
project spaces and interactions.

---

## Development environment

- **Rust:** stable toolchain (install via [rustup](https://rustup.rs/)).
- **Edition:** the workspace is on **Rust edition 2024**, so a reasonably recent
  stable compiler is required.
- **Layout:** a Cargo **workspace** of 19 crates under `crates/`, with the
  top-level `oxiphysics` crate re-exporting every domain.

Build the whole workspace:

```bash
cargo build
```

For an optimized build (and to mirror what CI-equivalent local checks use):

```bash
cargo build --release
```

---

## Quality gates — run these locally before opening a PR

OxiPhysics intentionally ships **no public PR/CI gate**. By COOLJAPAN policy the
only GitHub workflows in this repository are `pypi-publish.yml` and
`npm-publish.yml` (release plumbing for the Python/WASM bindings). There is **no
automated check that runs on your pull request**, which means the burden is on
**you, the contributor**, to verify the code is clean before requesting review.

Run all four of the following and make sure each one passes:

1. **Formatting** — the tree must already be formatted:

   ```bash
   cargo fmt --all
   ```

2. **Lints** — Clippy must be completely **silent** with warnings denied:

   ```bash
   cargo clippy --workspace --all-features --all-targets -- -D warnings
   ```

3. **Tests** — the full suite runs via [`cargo nextest`](https://nexte.st/)
   (the workspace currently has **60,115 tests** and **zero stubs**):

   ```bash
   cargo nextest run --workspace --all-features
   ```

4. **Docs** — Rustdoc must build with warnings denied (no missing docs,
   no broken intra-doc links):

   ```bash
   RUSTDOCFLAGS='-D warnings' cargo doc --all-features --no-deps
   ```

If any of these fail, fix the issue before opening the pull request. A PR that
does not pass the gates above will be asked to do so before review.

---

## COOLJAPAN engineering policies

Contributions are expected to follow the same policies the rest of the codebase
upholds:

- **No `.unwrap()` in production code.** Propagate errors with `?`, return a
  `Result`, or handle the `None`/`Err` case explicitly. (`.unwrap()`/`.expect()`
  are acceptable in tests and examples where a panic is the intended failure
  mode.)
- **Zero warnings.** The compiler and Clippy must emit no warnings at all — see
  the gate above. We do **not** add blanket `#[allow(...)]` attributes to silence
  lints; fix the underlying cause.
- **Pure Rust by default.** Default features must be 100% Pure Rust with **no
  C/C++/Fortran** build-time dependencies. Any such dependency (e.g. a CUDA
  backend) must be **optional and feature-gated**, never enabled by default.
- **Latest crate versions.** Prefer the latest versions published on crates.io
  when adding or bumping dependencies.
- **Workspace-inherited dependencies.** Declare each dependency **once** at the
  workspace root and inherit it in member crates via `*.workspace = true`. Do
  **not** pin versions per crate. The only per-crate exceptions are `keywords`
  and `categories`, which legitimately differ between sub-crates.
- **snake_case naming.** Follow standard Rust naming conventions
  (`snake_case` for items/variables, `CamelCase` for types, etc.).
- **Keep files small.** Keep each source file under ~2000 lines; split large
  modules rather than letting a single file grow unbounded.

---

## Branch policy

- `master` receives **only** "Availability of X.Y.Z" release commits.
- **All development happens on `0.x.y` version branches**, never directly on
  `master`.

Create your feature/fix branch from the current development branch, do your work
there, and open the pull request against that development branch.

---

## Pull request checklist

Before requesting review, confirm:

- [ ] The PR has a clear **description** of what changed and why.
- [ ] New behavior is covered by **tests** (added or updated).
- [ ] **Docs are updated** — rustdoc comments, and the README/CHANGELOG where
      user-facing behavior changed.
- [ ] `cargo fmt --all` produces no diff.
- [ ] `cargo clippy --workspace --all-features --all-targets -- -D warnings` is
      silent — **no new warnings**.
- [ ] `cargo nextest run --workspace --all-features` passes.
- [ ] `RUSTDOCFLAGS='-D warnings' cargo doc --all-features --no-deps` passes.
- [ ] **No new `.unwrap()`** in production code.
- [ ] The branch follows the branch policy above (not targeting `master`).

---

## Code of Conduct

Participation in this project is governed by our
[Code of Conduct](CODE_OF_CONDUCT.md). Please read it before contributing.
