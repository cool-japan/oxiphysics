# oxiphysics-io TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 128,821 SLoC | 4,879 tests

File-format I/O subcrate. All completed items below are production code, tested.

> **Constraints:** Pure-Rust parsers only. Compression exclusively via oxiarc-* (no zip/flate2/zstd-c/bzip2/lz4/tar). Fuzzing runs as a cargo-fuzz local harness script — no new workflow yamls. Honest-markers philosophy: "simplified / stub / subset" doc banners stay in the code until the format is either upgraded (roadmap items below) or explicitly designated Tier-2 in the v1.0 conformance matrix — never silently relabeled.

## Completed (v0.1.0 – v0.1.3)

### Phase 1: Foundation
- [x] Core types (trait-based reader/writer abstractions, common error type)
- [x] Basic error handling
- [x] Unit tests

### Phase 2: Mesh / geometry formats
- [x] VTK legacy / XML (VTU, PVTU if present)
- [x] glTF 2.0 (scene export for visualization)
- [x] OBJ (Wavefront)
- [x] STL (binary + ASCII)
- [x] Gmsh .msh

### Phase 3: Solver-specific formats
- [x] CalculiX (.inp / .frd) — FEM solver compat
- [x] LAMMPS (data / dump / restart) — MD solver compat
- [x] OpenFOAM case directory
- [x] Abaqus .inp subset

### Phase 4: Scientific / trajectory
- [x] HDF5 containers
- [x] NetCDF
- [x] XDMF (goes with VTK HDF5 backend)
- [x] CSV / JSON time-series

### Phase 5: Scene / state
- [x] Scene JSON round-trip (via oxiphysics scene module)
- [x] Snapshot binary format (fast save / restore)

### v0.1.2 fixes
- [x] (2026-06-01) Re-exported `particle_formats` public types at crate root (`DcdWriter`, `DcdReader`, `XyzWriter`, `XyzReader`, `ParticleFrame`, `ParticleTrajectory`, `TrajectoryStats`, `BinaryFrameReader`, `BinaryFrameWriter`, `GroReader`, `GroWriter`, `DcdHeader`). Fixes failing `DcdWriter` doctest.

### Known in-code limitation markers (audited 2026-06-11)

Literal `TODO` count in `src/` is 0. There are ~52 "simplified / stub / subset" doc markers — honest labels on intentionally reduced implementations, not missing code. Audit method: ripgrep over `src/` for TODO/FIXME/simplified/stub/subset/mock banners (2026-06-11). Representative audited sample below; forward roadmap items reference these as "(closes marker N)". Anything left unclosed gets an explicit Tier-2 designation at v1.0.

1. `src/hpc_io.rs:305` — "Parallel netCDF stub for distributed array I/O"
2. `src/crystallography_io.rs:8` — "VASP POSCAR/CONTCAR (stub)"
3. `src/xtc_dcd.rs:1–6` — "Simplified... Real XTC uses lossy compression (XDR)... mimicked structurally" (not GROMACS-readable)
4. `src/hdf5_io/mod.rs:1` — "Pure-Rust in-memory HDF5 mock... entirely in memory" (not spec-compliant on disk)
5. `src/foam_io.rs:4` — "Simplified OpenFOAM polyMesh format I/O"
6. `src/cgns_format.rs:4–6` — "CGNS... text-based subset"
7. `src/restart_io.rs:28–32` — "Simplified HDF5-like layout", "Simplified MessagePack-like binary encoding"
8. `src/machine_learning_io.rs:8` — "ONNX-like simplified format"
9. `src/simulation_io.rs:312` — "simplified parser for the format produced by PhysicsSceneWriter"
10. `src/geospatial_io.rs:723,1214` — "simplified Mercator", "simplified in-memory DBF"
11. `src/scientific_formats.rs:8` — "simplified parser for Gaussian fchk"
12. `src/binary_formats.rs:206,332` — header "placeholder; patched in finalize" (benign two-pass writer)
13. `src/robotics_io/types.rs:179–463` — URDF data model (`UrdfJoint`, `UrdfVisualElement`, `UrdfGeometry`) exists; NO XML parser found
14. `src/cad_io/` — IGES entity reader real; STEP exists only as `stepparser_traits.rs` (partial)
15. `src/checkpoint_io/types.rs:591,1208` — real `oxiarc_zstd::compress_with_level`/`decompress` already wired

## v0.2.0 — Robotics interop + trustworthy parsers

Headline: the URDF XML parser (the data model already exists; the parser does not). Sequencing: `.expect()` sweep → fuzzing harness; URDF XML layer → MJCF; checkpoint_io zstd precedent → trajectory containers. Exit criteria: KUKA iiwa URDF round-trip green; zero `.expect()` reachable from public parse paths; 24h fuzz run clean; compressed + checksummed + schema-versioned containers shipped.

Definition of done for every parser/codec item from here on:
1. Typed errors only — `#[non_exhaustive]` error enums; zero `.expect()` on public parse/read paths.
2. Round-trip test (or third-party-tool verification where round-trip is asymmetric).
3. A fuzz target exists before the format can claim Tier-1 at v1.0.
4. Reference fixtures committed (small + redistributable; licensing checked); generated bulk data goes to `std::env::temp_dir()`, never the repo.
5. Honest markers: "simplified" banners are removed only by the test that proves otherwise.

Fixture policy: prefer official/public-domain references (NIST STEP files, OpenFOAM tutorial cases, RCSB structures); tool-generated fixtures (gmx, h5py) are committed together with a regeneration script so they stay auditable.

### Robotics interop (headline)
- [ ] URDF XML parser → `ArticulatedModel` import — (closes marker 13; feeds root Phase 29 robotics round-trip demos)
  - **Goal:** round-trips 3 reference robots incl. KUKA iiwa with joint limits/inertials/mesh refs; imported model's ABA gravity torques match analytic 2-link within 1e-6.
  - **Note:** data model exists (`robotics_io/types.rs` — `UrdfJoint`, `UrdfVisualElement`, `UrdfGeometry`); the XML parse + articulated mapping are new; pure-Rust XML layer.
  - **Files:** `src/robotics_io/` (XML parser NEW next to `types.rs`, articulated mapping NEW); reference URDF fixtures committed under test data.
  - **Tests:** KUKA iiwa + 2 further reference robots round-trip; analytic 2-link ABA gravity-torque check at 1e-6.
  - **Step plan:** (a) pure-Rust XML layer → (b) `<robot>/<link>/<joint>` parse into the existing types → (c) kinematic-tree validation (single root, no cycles) → (d) inertial/limit mapping into `ArticulatedModel` → (e) mesh-ref resolution relative to the URDF path → (f) writer for the round-trip half.
  - **Risk:** mesh-ref resolution (`package://` URIs, relative paths) — resolve relative to the URDF file plus a caller-supplied search-path list; never a hardcoded absolute path.
- [ ] MJCF (MuJoCo XML) import subset (body/joint/geom/actuator/default) (dep: URDF XML layer)
  - **Goal:** loads `ant.xml` + `cartpole.xml` into articulated+rigid; feeds the python gym envs.
  - **Files:** `src/robotics_io/` (MJCF reader NEW, reusing the URDF XML layer).
  - **Tests:** `ant.xml` + `cartpole.xml` load; `<default>` class resolution unit-tested separately (it is the subtle half of the subset).
  - **Risk:** MuJoCo defaults-cascade semantics — document the supported subset precisely; unsupported attributes are collected into a typed parse report, never silently dropped.

### Parser trustworthiness
- [ ] `.expect()` → `Result` sweep in all parse paths (audited 113 production sites; pairs with root 23.4) — prerequisite for the fuzzing item
  - **Goal:** zero `.expect()` reachable from any public `parse`/`read` fn; error enums `#[non_exhaustive]`.
  - **Design:** convert genuinely-fallible sites to typed errors; keep invariant-protected ones with `// INVARIANT:` comments; no-unwrap-style local script as the gate.
  - **Files:** all `src/**` parse paths (audited 113 sites); the gate script lives as a local script, not a workflow.
- [ ] Parser fuzzing harness (cargo-fuzz, local script — no workflow yaml)
  - **Goal:** fuzz targets for OBJ/STL/Gmsh/VTK/CalculiX-inp/LAMMPS/PDB/GRO/snapshot; 24h run with zero panics/OOM; corpus committed under `fuzz/corpus/`.
  - **Files:** `fuzz/` (NEW cargo-fuzz dir) with nine targets:
    `fuzz_obj`, `fuzz_stl`, `fuzz_gmsh`, `fuzz_vtk`, `fuzz_calculix_inp`, `fuzz_lammps`, `fuzz_pdb`, `fuzz_gro`, `fuzz_snapshot`.
  - **Design:** seed corpora from existing test fixtures; runner is a local script (documented invocation), never a CI yaml.
  - **Risk:** corpus bloat in-repo — minimize with `cargo fuzz cmin` before committing and cap per-target corpus size.

### Containers & durability
- [ ] oxiarc-zstd compressed trajectory containers — note: partially exists (checkpoint_io already compresses via oxiarc-zstd — marker 15)
  - **Remaining:** extend to `trajectory/` writers + binary snapshots, documented container layout.
  - **Goal:** 10x size reduction on a 1k-frame rigid trajectory; round-trip hash-identical.
  - **Files:** `src/trajectory/` writers, snapshot writer; layout doc in module rustdoc (precedent: `checkpoint_io/types.rs:591,1208`).
  - **Tests:** 1k-frame round-trip hash-equality + compression-ratio assertion on the reference trajectory.
- [ ] Checksummed integrity (xxhash, pure Rust)
  - **Goal:** per-frame checksums in snapshot/checkpoint/trajectory containers; an injected bit-flip yields a typed `ChecksumMismatch` error, never garbage data.
  - **Files:** snapshot/checkpoint/trajectory writers gain a per-frame checksum field; xxhash via a pure-Rust implementation.
  - **Tests:** bit-flip injection test per container type.
  - **Risk:** checksum placement must not break 0.1.x readers — couple it to the schema-version field so there is exactly one migration.
- [ ] Schema-versioned snapshot migration
  - **Goal:** version field + migration registry; 0.1.x snapshots load in 0.2.0 with an explicit migration log; v1→v2 test fixture committed.
  - **Files:** snapshot format module (version header) + migration registry keyed `(from, to)`.
  - **Tests:** the committed v1 fixture loads through the v1→v2 migration with the migration log asserted.
- [ ] Streaming/chunked readers for >RAM files — note: `streaming_io/` module exists — verify scope first
  - **Goal:** 50 GB synthetic LAMMPS dump + VTU appended-data streamed frame-by-frame under 512 MB RSS via iterator readers.
  - **Files:** `src/streaming_io/` (extend after scope audit), LAMMPS dump + VTU appended-data iterator readers.
  - **Design:** RSS assertion via a local script harness (synthetic file generated into `std::env::temp_dir()`, never committed).
  - **Tests:** RSS-bounded full iteration over the generated dump — env-gated (`OXIPHYSICS_BIGFILE=1` pattern) so default test runs stay fast.

### Cross-crate consumers (for sequencing)
- URDF/MJCF → oxiphysics-articulated `ArticulatedModel`, python gym envs (oxiphysics-python v0.2.0), root Phase 29 interop demos
- Compressed + checksummed containers → umbrella snapshot/rollback; root Phase 30 out-of-core goal (100M-frame trajectory under 512 MB RSS)
- Schema-versioned snapshots → root Phase 24 replay golden corpus (determinism program)
- HDF5 on-disk subset (v0.3.0) → XDMF backend, python-side h5py verification
- glTF KHR physics export (v0.3.0) → viz scene export, root Phase 29
- STEP import (v0.3.0) → CAD-to-sim collision-mesh pipeline (root Phase 29)

## v0.3.0 — Scientific-format depth

Theme: replace the honest mimics with spec-faithful codecs, verified against third-party tooling (gmx, h5py, ParaView, usdchecker). Exit criteria: gmx-readable XTC; h5py-readable HDF5 subset; ParaView-openable OpenFOAM case; STEP/OpenVDB/USD subsets land with reference-file tests.

- [ ] Real XTC/TRR GROMACS codecs (pure-Rust XDR + xdrfile lossy 3D coordinate compression port) — (closes marker 3)
  - **Goal:** byte-exact decode of a gmx-generated reference `.xtc`; OxiPhysics-encoded files readable by `gmx check`; replaces the simplified mimic in `xtc_dcd.rs`.
  - **Files:** `src/xtc_dcd.rs` (mimic replaced), pure-Rust XDR layer (NEW); gmx reference fixtures committed.
  - **Step plan:** (a) XDR primitive codec → (b) XTC header/magic → (c) port of the xdrfile 3D coordinate compressor (the hard part — bit-exact integer coding) → (d) TRR (simpler, uncompressed) → (e) writer + `gmx check` verification.
  - **Risk:** bit-exactness of the lossy coordinate coder — land the decoder against gmx-generated files first; only then attempt the encoder.
- [ ] HDF5 on-disk subset (superblock v2 + contiguous datasets, pure Rust) — (closes marker 4)
  - **Goal:** `h5py` opens an OxiPhysics-written file and reads the positions dataset (asserted from the python crate test suite); in-memory mock retained for the API layer.
  - **Files:** `src/hdf5_io/` (on-disk writer NEW beside the in-memory mock).
  - **Tests:** pure-Rust write→read round-trip here; the h5py interop assertion lives in the python crate suite.
  - **Risk:** scope creep — superblock v2 + contiguous datasets ONLY; chunking/filters/extended groups stay out of scope (binary CGNS in Deferred waits on exactly this layer).
- [ ] OpenFOAM polyMesh full-fidelity upgrade — (closes marker 5; honest-markers philosophy applied to `foam_io.rs`)
  - **Goal:** a real OpenFOAM tutorial case (polyMesh `points`/`faces`/`owner`/`neighbour`/`boundary` + volScalarField/volVectorField, ASCII first) round-trips and re-opens in ParaView; the "Simplified" banner comes off only when that test is green.
  - **Files:** `src/foam_io.rs` (subset → full-fidelity).
  - **Tests:** tutorial-case round-trip in-repo; ParaView re-open verification skip-not-fail when the tool is absent.
- [ ] glTF KHR physics extensions export (`KHR_physics_rigid_bodies` + `KHR_implicit_shapes`)
  - **Goal:** exported scene re-imports with an identical collider set; JSON validates against the extension schemas (gltf module exists with KHR_lights_punctual precedent).
  - **Files:** existing gltf module (extension serializers NEW); schema-validation fixtures committed.
  - **Tests:** export → re-import collider-set equality; JSON-schema validation against committed extension schemas.
- [ ] USD subset — usda text export first, import second
  - **Goal:** usdchecker-clean usda for a rigid scene (meshes + xforms + hierarchy).
  - **Files:** usd module (NEW — usda tokenizer/printer only; binary usdc is out of scope).
  - **Tests:** usda golden text comparison in-repo; usdchecker verification skip-not-fail.
  - **Risk:** usdchecker is a local-tool dependency — the verification test is skip-not-fail when the tool is absent.
- [ ] STEP AP203/214 import (pure Rust) building on existing `cad_io` STEP traits — (advances marker 14)
  - **Goal:** 3 NIST STEP test files imported to watertight triangle meshes usable as collision geometry.
  - **Files:** `src/cad_io/` (implementation behind `stepparser_traits.rs`; IGES reader is the in-crate precedent).
  - **Tests:** 3 NIST files → closed-2-manifold (watertight) assertion + collision-geometry usability smoke test.
  - **Risk:** entity-coverage explosion — restrict to the entity set needed for B-rep → tessellation of the 3 NIST files; unsupported entities surface as typed errors, not panics.
- [ ] OpenVDB read subset (VDB tree topology + level sets, pure Rust)
  - **Goal:** reads a reference `sphere.vdb` into `VoxelGrid`/SDF with max error <1 voxel.
  - **Files:** vdb reader module (NEW); `sphere.vdb` fixture committed.
  - **Risk:** tree-configuration generality — support the standard 5-4-3 tree first; other configurations are typed errors until demanded.
- [ ] mmCIF full dictionary + EnSight Gold completion
  - **Goal:** parses 4 reference mmCIF structures incl. multi-model (pdb.rs is PDB-only today); EnSight Gold case re-opens in ParaView.
  - **Files:** extend `pdb.rs` (mmCIF dictionary path) + EnSight Gold writer completion.
  - **Tests:** 4 RCSB reference structures (incl. a multi-model NMR entry); EnSight ParaView re-open skip-not-fail.

## v1.0 — Production & validation

Exit criteria: the README format table is generated, not hand-maintained; every Tier-1 format has a fuzz target plus a round-trip or third-party verification test; every Tier-2 format's banner comment links back to its entry here.

- [ ] Format conformance matrix + stability tiers — (resolves every remaining "simplified/subset" marker, incl. 2, 6–12, by explicit tiering)
  - **Goal:** every format documented Tier-1 (round-trip guaranteed, fuzz-covered, semver-frozen) or Tier-2 (best-effort/simplified, marker comments retained intentionally); README table generated by xtask.
  - **Files:** conformance-matrix generator as an xtask subcommand (the workspace `xtask` crate is root item 23.6).
  - **Provisional tier map at freeze time (finalized by the matrix itself):**

    | Format | Target tier | Basis |
    |---|---|---|
    | VTK/VTU, OBJ, STL, Gmsh | 1 | fuzz-covered + round-trip tests |
    | LAMMPS, CalculiX, PDB/GRO/DCD/XYZ | 1 | fuzz-covered + solver-compat tests |
    | Snapshot / checkpoint / trajectory containers | 1 | checksummed + schema-versioned (v0.2.0) |
    | URDF / MJCF | 1 | reference-robot round-trips (v0.2.0) |
    | XTC/TRR, HDF5 subset, OpenFOAM | 1 | third-party-tool verification (v0.3.0) |
    | glTF (+KHR physics), USD subset, STEP, OpenVDB | 1 or 2 | per v0.3.0 outcome |
    | Gaussian fchk (11), geospatial (10), PhysicsSceneWriter parser (9), restart encodings (7) | 2 | documented best-effort, markers retained |
    | VASP (2), ONNX-like (8), CGNS text (6) | 2 unless promoted | see Deferred |

## Deferred / research track

- [~] Parallel netCDF for distributed array I/O (marker 1) — Ready when: a pure-Rust distributed transport exists (root Phase 30 large-scale program); until then the v0.2.0 single-node streaming readers cover the scale story.
- [~] Binary CGNS (ADF/HDF5-backed) beyond the text subset (marker 6) — Ready when: the HDF5 on-disk subset (v0.3.0) ships, since binary CGNS rides on that container layer.
- [~] Real ONNX protobuf codec replacing the "ONNX-like simplified format" (marker 8) — Ready when: ML-interop demand materializes and a pure-Rust protobuf layer is approved as a dependency.
- [~] VASP POSCAR/CONTCAR full reader (marker 2) — Ready when: a crystallography/MD-interop user story lands; otherwise documented Tier-2 at v1.0.

---
Conventions: `[x]` done / `[ ]` planned / `[~]` deferred until the stated condition. Dates are YYYY-MM-DD. **Goal:**/**Design:** wording is the per-item contract; acceptance criteria are testable (third-party-tool checks skip-not-fail when the tool is absent), not aspirations.
