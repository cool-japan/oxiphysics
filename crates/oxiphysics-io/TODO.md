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
- [x] URDF XML parser → `ArticulatedModel` import — (closes marker 13; feeds root Phase 29 robotics round-trip demos)
  - **Impl files (2026-06-13):** NEW src/robotics_io/{xml.rs, urdf_parser.rs, articulated_map.rs, urdf_writer.rs, urdf_error.rs}; MODIFY Cargo.toml, src/robotics_io/mod.rs; NEW tests/urdf_roundtrip.rs, tests/fixtures/two_link_pendulum.urdf
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

<!-- AUDIT-2026-10-12 BEGIN -->
## Audit 2026-10-12 — verified findings (Phase 32 remediation)

Verified: 387 (critical 0, high 109, medium 224, low 54); 61 further raw findings pending second-round verification. IDs refer to the audit ledger; `[x]` = fixed.

### High (109)
- [x] **F1976** (high, bug, medium/medium) `abaqus_format.rs:362` — AbaqusWriter emits invalid INP (section data, undefined sets, *Plastic, *End Part)
- [x] **F1980** (high, bug, easy/small) `amber/functions/functions_impl.rs:83` — build_angles drops the third atom of every angle
- [x] **F1981** (high, bug, easy/small) `amber/functions/functions_impl.rs:89` — build_angles never reads the third atom; AmberAngle lacks atom k
- [x] **F1982** (high, fabrication, medium/small) `amber/functions/functions_impl.rs:115` — extract_lj_parameters labels the packed pair table as diagonal (i,i) types
- [x] **F1983** (high, fabrication, easy/small) `amber/functions/functions_impl.rs:132` — extract_dihedrals hardcodes zero parameters and misreads the type index
- [x] **F1984** (high, bug, easy/small) `amber/functions/functions_impl.rs:519` — parse_mdcrd re-slices quadratically and panics on a boxless final frame
- [x] **F1988** (high, fabrication, medium/medium) `ambermd_io.rs:403` — PrmtopFile write/read drops bond atoms, fabricates (2i,2i+1) pairs
- [x] **F1989** (high, bug, easy/small) `ambermd_io.rs:532` — Rst7File::read reads the box as velocities; fails on touching 12.7 fields
- [x] **F1990** (high, bug, easy/small) `ambermd_io.rs:649` — mdcrd readers break on touching 10F8.3 fields (<= -100 or >= 1000)
- [x] **F1995** (high, panic, easy/small) `amr_io/types.rs:351` — AmrCheckpoint::read panics on truncated files before mapping the error
- [x] **F1999** (high, bug, medium/medium) `animation_io.rs:509` — BVH ignores CHANNELS rotation order; export Euler convention mismatches import
- [x] **F2003** (high, bug, medium/medium) `binary_format.rs:677` — OxiFile decoder: unbounded group recursion, wrapped bounds check
- [x] **F2009** (high, panic, easy/small) `binary_formats.rs:429` — DcdReader panics on 8 zero bytes and short frames
- [x] **F2014** (high, bug, medium/medium) `cad_io/functions.rs:65` — tessellate_brep ignores edge direction; create_box mesh is broken
- [x] **F2015** (high, bug, medium/medium) `cad_io/types.rs:79` — IGES parser keeps the entity type in parameters and matches records by order
- [x] **F2020** (high, bug, medium/small) `calculix_format.rs:428` — CalculiX *BOUNDARY loses DOF ranges, pre-step BCs and node-set names
- [x] **F2027** (high, bug, easy/small) `checkpoint_io/types.rs:196` — RestartFile::load parses scalars from undecompressed bytes
- [x] **F2034** (high, bug, easy/small) `checkpoint_io/types.rs:1205` — CheckpointReader never verifies the CRC footer written by finalize
- [x] **F2035** (high, bug, easy/small) `checkpoint_io/types.rs:1275` — Checkpoint footer CRC written but never verified on read
- [x] **F2036** (high, panic, easy/small) `checkpoint_io/types.rs:1317` — Unchecked count*N cursor math and with_capacity on untrusted checkpoint bytes
- [x] **F2037** (high, bug, easy/small) `checkpoint_io/types.rs:1323` — Checkpoint block-skip arithmetic wraps; hang or capacity panic
- [x] **F2038** (high, bug, medium/medium) `crystallography_io.rs:91` — CIF parser keeps only the last loop, never leaves loop mode, stores U as B
- [x] **F2042** (high, bug, medium/small) `crystallography_io.rs:798` — POSCAR Cartesian→fractional assumes orthogonal cell, ignores scale, mislabels mode
- [x] **F2043** (high, bug, easy/small) `crystallography_io.rs:965` — simulate_xrd: Miller indices capped at ±5, and wrong LP factor
- [x] **F2044** (high, bug, medium/small) `crystallography_io.rs:1004` — simulate_xrd: wrong Lorentz-polarization factor, hkl capped at 5, width in grid points
- [x] **F2063** (high, bug, easy/trivial) `exodus/functions.rs:821` — refine_tet4_mesh emits four inverted child tetrahedra
- [ ] **F2074** (high, bug, medium/medium) `finite_element_io/types.rs:166` — MAPDL NBLOCK/EBLOCK read from wrong columns; real CDB misparsed
- [ ] **F2075** (high, fabrication, medium/medium) `finite_element_io/types.rs:382` — GmshMesh::parse v4 branch is fake; v4 $Elements panics
- [ ] **F2077** (high, panic, medium/small) `finite_element_io/types.rs:707` — parse_chexa panics on 10-field lines; continuations never joined
- [ ] **F2079** (high, bug, easy/trivial) `finite_element_io/types.rs:766` — NastranBulkData::write silently drops all CTETRA/CHEXA elements
- [ ] **F2080** (high, stub, medium/medium) `finite_element_io/types.rs:1029` — Abaqus *BOUNDARY/*CLOAD documented but ignored; GENERATE/continuations misparsed
- [ ] **F2083** (high, bug, medium/medium) `fluent_format.rs:486` — Fluent reader cannot parse standard Fluent .msh node/face layout
- [ ] **F2088** (high, bug, medium/small) `foam_io.rs:733` — OpenFOAM vector field parser swallows boundaryField vectors
- [ ] **F2098** (high, bug, medium/medium) `gltf/functions.rs:1132` — scene_to_gltf_json emits empty node graph; one mesh per primitive
- [ ] **F2102** (high, panic, medium/medium) `gmsh_format.rs:233` — Gmsh reader panics on n_tags; Gmsh 4.1 files crash it
- [ ] **F2122** (high, stub, easy/small) `hdf5_io/trajectory.rs:347` — CheckpointManager eviction forgets label but never deletes checkpoint data
- [ ] **F2130** (high, panic, easy/small) `hdf5_simple/types.rs:1005` — ShdfFile::from_bytes reserves untrusted lengths: capacity-overflow panic
- [ ] **F2145** (high, panic, easy/small) `json_io.rs:558` — JSON parser: stack overflow on nesting, char-boundary panics
- [ ] **F2146** (high, panic, easy/small) `json_io.rs:652` — JSON parser panics on UTF-8 byte slicing; unbounded recursion
- [ ] **F2148** (high, bug, medium/medium) `lammps/types.rs:533` — LammpsDataReader assumes atomic style; Pair Coeffs aborts; box fabricated
- [ ] **F2149** (high, bug, easy/small) `lammps/types.rs:972` — lammps::LammpsDumpReader drops <8-column atoms, returns wrong timestep
- [ ] **F2151** (high, panic, easy/small) `lammps/types.rs:1156` — Binary dump/restart readers allocate from untrusted atom counts
- [ ] **F2154** (high, bug, medium/small) `lammps_dump.rs:304` — LAMMPS dump columns read by position, ignoring header names
- [ ] **F2157** (high, bug, easy/small) `lammps_dump.rs:669` — RDF counts pairs as N_b and ignores PBC; g(r) off by N_a
- [ ] **F2171** (high, panic, easy/small) `material_db_io.rs:528` — Material JSON parser mixes char/byte indices: panics on non-ASCII
- [ ] **F2172** (high, panic, easy/trivial) `material_db_io.rs:591` — extract_json_str uses char count as byte index: panic or truncation
- [ ] **F2175** (high, fabrication, easy/small) `medical_imaging_io/functions.rs:28` — NIfTI reader fabricates all-zero volume for unsupported datatypes
- [ ] **F2182** (high, bug, medium/medium) `medical_imaging_io/types.rs:1065` — DICOM parser silently stops at first undefined-length element
- [ ] **F2184** (high, fabrication, medium/medium) `medical_io/types.rs:773` — NiftiWriter emits a mock header that real readers cannot use
- [ ] **F2189** (high, bug, hard/medium) `mesh_export.rs:1001` — QEM simplify uses stale targets and corrupts output vertex indices
- [ ] **F2190** (high, bug, hard/medium) `mesh_export.rs:1030` — Mesh simplify corrupts output vertices; QEM is midpoint/static
- [ ] **F2197** (high, fabrication, medium/small) `mesh_io.rs:1208` — MeshValidator counts 'inverted' triangles by normal.y < 0
- [ ] **F2199** (high, perf, medium/small) `mesh_quality/functions.rs:241` — Duplicate-vertex removal and suggest_improvements are O(V^2) with sqrt
- [ ] **F2204** (high, bug, easy/small) `molecular_visualization_io/functions.rs:142` — CIF parser zeroes coordinates with s.u. and never leaves loops
- [ ] **F2205** (high, bug, easy/trivial) `molecular_visualization_io/functions.rs:343` — write_poscar emits both 'Direct' and 'Cartesian' mode lines
- [x] **F2206** (high, bug, easy/small) `molecular_visualization_io/types.rs:275` — CHGCAR density() indexes z-fastest but grid stored x-fastest
- [x] **F2207** (high, fabrication, easy/small) `molecular_visualization_io/types.rs:576` — DX parser reads data from 'object 2', so APBS files load as zeros
- [x] **F2213** (high, missing-feature, hard/large) `netcdf/functions.rs:15` — No real NetCDF classic binary format; AMBER trajectories unusable externally
- [x] **F2219** (high, fabrication, easy/small) `netcdf/types.rs:326` — Trajectory builder write_cdl writes no data section; frames discarded
- [x] **F2222** (high, bug, medium/small) `numpy/functions.rs:161` — NPY reader ignores fortran_order, returns transposed data
- [x] **F2233** (high, bug, easy/small) `openfoam/functions.rs:1606` — green_gauss_gradient divides by 3V; neighbour indexing panics on standard lists
- [x] **F2234** (high, bug, easy/trivial) `openfoam/types.rs:141` — box_mesh internal y-faces wound so normals point neighbour→owner
- [x] **F2235** (high, bug, medium/small) `openfoam/types.rs:332` — OpenFOAM neighbour file includes -1 boundary entries
- [x] **F2236** (high, bug, medium/small) `openfoam_format.rs:125` — extract_paren_block uses rfind; boundary vectors become internal cells
- [x] **F2246** (high, panic, easy/small) `particle_data_io.rs:317` — Header-declared counts drive Vec::with_capacity: abort on corrupt files
- [x] **F2249** (high, bug, medium/small) `particle_data_io.rs:932` — Delta encoder drops sub-threshold motion; reconstruction drifts
- [x] **F2250** (high, fabrication, medium/small) `particle_data_io.rs:999` — ParticleDiffWriter::reconstruct fabricates mass/radius/time; drops sub-threshold drift
- [x] **F2253** (high, bug, medium/medium) `particle_formats.rs:418` — LAMMPS dump reader zeroes positions for xs/xu columns
- [x] **F2258** (high, bug, medium/medium) `particle_formats.rs:1015` — DCD reader: header underflow panic; unit-cell record read as X
- [x] **F2259** (high, bug, medium/medium) `particle_formats.rs:1061` — DCD reader misreads NAMD unit-cell record as X coordinates
- [x] **F2260** (high, bug, easy/small) `particle_formats.rs:1199` — DCD writer puts DELTA in NAMNF slot; VMD/NAMD misread files
- [x] **F2265** (high, bug, easy/small) `plot3d_format.rs:266` — PLOT3D .q read/written point-interleaved, not variable-major
- [x] **F2268** (high, bug, easy/small) `point_cloud_io/functions.rs:276` — estimate_normals eigen solver returns wrong normals (absolute shift)
- [x] **F2269** (high, bug, medium/small) `point_cloud_io/functions.rs:288` — Normal estimation eigen-solver fails on small-scale or axis-aligned covariances
- [x] **F2279** (high, bug, easy/small) `remote_sensing_io.rs:366` — Dem::flow_accumulation re-adds upstream totals every pass (wrong)
- [x] **F2281** (high, panic, medium/small) `remote_sensing_io.rs:810` — LAS reader takes 1.4 extended count from older files; capacity panic
- [x] **F2282** (high, bug, medium/small) `remote_sensing_io.rs:839` — LAS reader misreads LAS 1.2/1.3 counts and formats 6-10
- [x] **F2283** (high, bug, easy/small) `restart_io.rs:277` — write_binary emits arrays of differing lengths under one count n
- [x] **F2284** (high, bug, easy/small) `restart_io.rs:335` — Restart ASCII/JSON writes 6 decimals; reader zeroes parse errors
- [x] **F2286** (high, bug, medium/medium) `restart_io.rs:466` — JSON restart files cannot be read back; JSON writer is lossy
- [x] **F2287** (high, panic, easy/small) `restart_io.rs:510` — Restart binary reader preallocates from unchecked u64 count
- [x] **F2289** (high, bug, easy/small) `restart_io.rs:768` — save_checkpoint swallows write errors then prunes older good checkpoints
- [x] **F2292** (high, bug, easy/small) `robotics_io/urdf_parser.rs:64` — URDF inertial origin rpy ignored; inertia tensor left in wrong frame
- [x] **F2296** (high, bug, easy/small) `seismic_io.rs:468` — SAC header writes NVHDR at the NZJDAY offset; readers reject files
- [x] **F2299** (high, fabrication, hard/medium) `seismic_io.rs:932` — Butterworth FrequencyFilter::apply is a moving-average difference
- [x] **F2304** (high, bug, easy/trivial) `sensor_io.rs:421` — Madgwick gradient Jacobian has four sign errors; pitch correction diverges
- [x] **F2312** (high, bug, easy/small) `simulation_io.rs:313` — Scene/Checkpoint readers drop orientation, velocities, mass and shape
- [x] **F2313** (high, bug, easy/small) `simulation_io.rs:442` — TRAJ binary format misaligns without velocities; reader trusts header counts
- [x] **F2320** (high, bug, medium/small) `spectroscopy_io/types.rs:1023` — JCAMP-DX reader double-applies YFACTOR; writer output not re-readable
- [x] **F2325** (high, perf, easy/small) `stl/functions_2.rs:174` — STL weld_vertices is O(V^2) linear scan over unique vertices
- [x] **F2326** (high, panic, easy/small) `streaming_io/functions.rs:591` — Chunk decoders overflow size math; ChunkedReader::open can loop forever
- [x] **F2328** (high, bug, easy/trivial) `streaming_io/types.rs:30` — RleEncoder never escapes literal 0xFF bytes; decode corrupts data
- [x] **F2332** (high, bug, easy/small) `streaming_io/types.rs:628` — ChunkedReader::open loops forever on wrapping frame length
- [x] **F2335** (high, bug, easy/trivial) `streaming_io/types.rs:1263` — CompressionStream misreads literal f32 whose low byte is 0xFF
- [x] **F2338** (high, bug, medium/medium) `tecplot_format.rs:280` — Tecplot FE zones written and read without element connectivity
- [x] **F2348** (high, bug, medium/small) `trajectory/functions.rs:134` — Kabsch alignment never rotates: mat3_inv_transpose returns the inverse
- [x] **F2370** (high, bug, medium/medium) `vtu/functions.rs:101` — compute_cell_quality: wrong edges, SJ=0 for every hex/quad
- [x] **F2371** (high, bug, medium/small) `vtu/functions.rs:277` — estimate_gradient underestimates gradients ~3× (inconsistent estimator)
- [x] **F2376** (high, fabrication, medium/large) `vtu/types.rs:814` — VtuReader drops connectivity, returns zero points silently
- [x] **F2380** (high, bug, easy/small) `weather_data_io/functions.rs:165` — GEV L-moment fit: wrong constant, wrong location sign, wrong shape sign
- [x] **F2384** (high, bug, easy/small) `weather_data_io/types.rs:1109` — GribReader::next_message loops forever on total_len=0, overflows on huge
- [x] **F2386** (high, bug, medium/medium) `weather_data_io/types.rs:1208` — GRIB2 Section 3/4/5 fields read at wrong octet offsets
- [x] **F2387** (high, fabrication, medium/medium) `weather_data_io/types.rs:1250` — GRIB2 Section 7 decoded as raw f32, ignoring Section 5 packing
- [x] **F2392** (high, bug, easy/trivial) `weather_io.rs:778` — Bolton LCL uses ln(T−Td) in °C instead of ln(T/Td) in K
- [x] **F2394** (high, bug, easy/trivial) `weather_io.rs:1052` — density_altitude applies 120 ft/°C rule to metre inputs
- [x] **F2401** (high, bug, easy/small) `xarray_io.rs:295` — xarray_to_vtk_structured writes C-order data under x-fastest DIMENSIONS
- [x] **F2407** (high, panic, medium/small) `xdmf/functions.rs:341` — XDMF Mixed topology panics via chunks(0) in three writers
- [x] **F2416** (high, panic, easy/small) `xtc_dcd.rs:233` — Untrusted XTC/DCD header counts drive Vec::with_capacity (abort)
- [x] **F2421** (high, bug, easy/small) `xtc_dcd.rs:775` — xtc_to_dcd skips nm→Å conversion; DCD timestep unit mismatched

### Medium (224)
- [ ] **F1974** (medium, test, easy/medium) `crates/oxiphysics-io/Cargo.toml:19` — No malformed-input property or fuzz tests for any decoder
- [ ] **F1975** (medium, docs, easy/small) `crates/oxiphysics-io/README.md:32` — io/lbm/md/fem README examples use fabricated types and signatures
- [ ] **F1977** (medium, bug, medium/small) `abaqus_format.rs:395` — AbaqusWriter emits invalid deck: undefined sets, *End Part without *Part
- [ ] **F1978** (medium, bug, medium/small) `abaqus_format.rs:481` — Abaqus reader drops 2-D nodes and splits continued element lines
- [ ] **F1985** (medium, perf, easy/small) `amber/types.rs:161` — prmtop residue lookup is O(atoms x residues)
- [ ] **F1986** (medium, fabrication, easy/small) `amber/types.rs:1006` — AmberDcd claims VMD-readable DCD but writes a private format
- [ ] **F1987** (medium, panic, easy/trivial) `ambermd_io.rs:316` — Fixed-width 4-byte string slicing panics on non-ASCII prmtop text
- [ ] **F1991** (medium, bug, easy/trivial) `ambermd_io.rs:1147` — Force-field dihedral barrier not divided by IDIVF
- [ ] **F1992** (medium, bug, easy/small) `amr_io/types.rs:88` — AMR coarsening uses truncating division for negative indices
- [ ] **F1994** (medium, perf, easy/small) `amr_io/types.rs:311` — AMR checkpoint, VTU and VTK writers do unbuffered per-scalar writes
- [ ] **F1996** (medium, fabrication, medium/medium) `amr_io/types.rs:535` — AmrWriter::write_vtk emits a .vthb that no VTK reader accepts
- [ ] **F1997** (medium, fabrication, hard/medium) `amr_io/types.rs:556` — AmrWriter::write_vtk emits XML ParaView cannot load as AMR
- [ ] **F2000** (medium, panic, easy/small) `animation_io.rs:599` — export_bvh panics when a frame has fewer poses than joints
- [ ] **F2001** (medium, stub, medium/medium) `animation_io.rs:1218` — resample_clip ignores the InterpolationMode parameter
- [ ] **F2002** (medium, panic, easy/small) `binary_format.rs:561` — OxiFile deserializers trust untrusted counts; unbounded group recursion
- [ ] **F2004** (medium, panic, easy/small) `binary_format.rs:867` — Untrusted length prefixes drive allocations and arithmetic (OOM/panic)
- [ ] **F2006** (medium, bug, easy/small) `binary_format.rs:1103` — rle_decompress_f64 lets a tiny input allocate unbounded memory
- [ ] **F2007** (medium, bug, easy/small) `binary_formats.rs:215` — BinaryTrajectoryWriter accepts frames whose size differs from header
- [ ] **F2008** (medium, panic, easy/small) `binary_formats.rs:266` — BinaryTrajectoryReader frame-size overflow bypasses length check
- [ ] **F2011** (medium, panic, easy/small) `binary_formats.rs:685` — CheckpointWriter::load_state panics on short or corrupt buffers
- [ ] **F2012** (medium, perf, easy/small) `binary_io.rs:78` — Binary frame read/write issue 3 syscalls per atom on raw File
- [ ] **F2016** (medium, bug, medium/small) `cad_io/types.rs:350` — STEP parser drops or truncates entities that span several lines
- [ ] **F2017** (medium, bug, medium/small) `cad_io/types.rs:361` — StepParser drops multi-line and '#n = ' entities; linear find_by_id
- [ ] **F2019** (medium, bug, easy/small) `calculix_format.rs:110` — C3D20 lines exceed 16 entries; readers ignore continuation lines
- [ ] **F2021** (medium, bug, medium/small) `calculix_format.rs:469` — *NSET/*ELSET GENERATE read as literal ids; element continuation lines
- [ ] **F2022** (medium, bug, easy/small) `calculix_format.rs:480` — CalculiX material ids follow HashMap order; missing *ELASTIC gives E=0
- [ ] **F2024** (medium, bug, easy/small) `cgns_format.rs:356` — CgnsWriter accumulates bases across calls; structured dims and solution lost
- [ ] **F2026** (medium, panic, easy/trivial) `checkpoint_io/types.rs:88` — Checkpoint::from_bytes wrapped bounds check then slice panic
- [ ] **F2028** (medium, panic, easy/small) `checkpoint_io/types.rs:440` — Checkpoint readers trust length fields: capacity and arithmetic overflow
- [ ] **F2029** (medium, bug, medium/small) `checkpoint_io/types.rs:591` — Checkpoint writes non-atomic and footer CRC never verified on read
- [ ] **F2030** (medium, bug, easy/small) `checkpoint_io/types.rs:846` — CheckpointDiff cannot reproduce targets of a different length
- [ ] **F2033** (medium, bug, easy/small) `checkpoint_io/types.rs:999` — CheckpointCompressor silently corrupts offsets above 65535; degenerate Default
- [ ] **F2039** (medium, fabrication, medium/medium) `crystallography_io.rs:380` — structure_factor uses f = Z with no angle dependence; B unused
- [ ] **F2040** (medium, fabrication, medium/medium) `crystallography_io.rs:388` — Structure factor uses f=Z ignoring sinθ/λ and B; unknown elements f=10
- [ ] **F2041** (medium, fabrication, easy/small) `crystallography_io.rs:791` — POSCAR fills missing atoms with zeros; Cartesian conversion assumes an orthogonal cell
- [ ] **F2046** (medium, perf, easy/small) `csv_io/functions.rs:1032` — inner_join_csv is an O(L*R) nested loop cloning keys per pair
- [ ] **F2047** (medium, bug, easy/small) `csv_io/types.rs:268` — CsvFile has no RFC4180 quoting; NaN sort comparator is not a total order
- [ ] **F2049** (medium, bug, easy/trivial) `csv_io/types.rs:322` — write_row_f64 uses 6 decimals, zeroing small physical values
- [ ] **F2051** (medium, bug, medium/small) `database_io/types.rs:391` — Hand-rolled JSON: no escaping, params dropped on deserialize
- [ ] **F2052** (medium, panic, easy/trivial) `ensight_format.rs:380` — EnsightWriter::write_geo panics via chunks(0) for unknown element types
- [ ] **F2053** (medium, bug, medium/small) `ensight_format.rs:391` — EnSight variable writers hardcode part 1; 'given' IDs never written
- [ ] **F2054** (medium, bug, medium/medium) `ensight_format.rs:582` — EnSight read_geo: id-0 underflow, given-ID blocks misread, parts overwrite
- [ ] **F2055** (medium, perf, easy/small) `ensight_format.rs:623` — write_ensight_case writes unbuffered per-value lines; dead buffer clones
- [ ] **F2057** (medium, api, medium/large) `error.rs:10` — Crate Error enum not non_exhaustive and lacks typed variants
- [ ] **F2059** (medium, panic, easy/small) `exodus/functions.rs:155` — Exodus parse_text skips connectivity validation, so downstream quality functions panic
- [ ] **F2060** (medium, bug, easy/trivial) `exodus/functions.rs:578` — tet4_quality scaled Jacobian unnormalised: regular tet scores 0.707
- [ ] **F2061** (medium, bug, easy/small) `exodus/functions.rs:647` — hex8_quality scaled Jacobian evaluated only at corner 0
- [ ] **F2064** (medium, bug, easy/small) `exodus/functions.rs:973` — partition_mesh_round_robin silently drops every block after the first
- [ ] **F2065** (medium, bug, medium/small) `exodus/functions.rs:985` — partition_mesh_round_robin drops every block after the first
- [ ] **F2066** (medium, missing-feature, hard/large) `exodus_format.rs:6` — README advertises Exodus, but output is custom ASCII, not NetCDF
- [ ] **F2067** (medium, missing-feature, hard/large) `exodus_format.rs:224` — No real Exodus II (netCDF) I/O; three incompatible text formats
- [ ] **F2068** (medium, panic, easy/small) `exodus_format.rs:388` — ExodusReader preallocates from untrusted counts; set names with spaces break round-trip
- [ ] **F2069** (medium, bug, easy/trivial) `experimental_data_io.rs:329` — Two inconsistent, mis-normalised von Mises strain formulas
- [ ] **F2073** (medium, fabrication, easy/trivial) `experimental_data_io.rs:1122` — all_metrics masks undefined NRMSE/correlation/R² as 0.0
- [ ] **F2076** (medium, bug, easy/small) `finite_element_io/types.rs:668` — Nastran implied-exponent reals rejected; GRID/FORCE cards silently dropped
- [ ] **F2078** (medium, fabrication, easy/small) `finite_element_io/types.rs:724` — MAT1 parser fabricates nu=0.3 and E/G=0 instead of the Nastran rule
- [ ] **F2081** (medium, bug, medium/medium) `finite_element_io/types.rs:1082` — Abaqus parser: continuation lines become elements; 2-D nodes, GENERATE, BCs lost
- [ ] **F2082** (medium, bug, medium/small) `fluent_format.rs:364` — Fluent writer: ND=3 in 2-D files, all faces in one interior zone
- [ ] **F2086** (medium, bug, medium/small) `foam_io.rs:214` — cell_volumes uses centroid distance instead of normal pyramid height
- [ ] **F2087** (medium, missing-feature, easy/small) `foam_io.rs:703` — OpenFOAM case writers only list paths; module claims polyMesh directory reading
- [ ] **F2089** (medium, bug, easy/small) `geospatial_io.rs:829` — DemRaster::sample documented bilinear but does nearest-cell; panics on empty
- [ ] **F2095** (medium, bug, medium/small) `geospatial_io.rs:1332` — dem_contour_lines: saddle zig-zag, nodata contoured, segments unstitched
- [ ] **F2096** (medium, panic, easy/small) `gltf/functions.rs:130` — parse_obj underflows on index 0, truncates quads, rejects negative indices
- [ ] **F2097** (medium, bug, medium/medium) `gltf/functions.rs:1044` — scene_to_gltf_json emits non-rendering, spec-invalid glTF and discards external buffer
- [ ] **F2099** (medium, fabrication, medium/medium) `gltf/types.rs:474` — export_animation_json drops keyframe values; not glTF animation schema
- [ ] **F2100** (medium, stub, medium/medium) `gltf/types.rs:1197` — GltfScene::to_json references accessors it never writes
- [ ] **F2101** (medium, bug, medium/small) `gltf/types.rs:1446` — GLB writer drops node transforms/hierarchy/materials and writes zero normals
- [ ] **F2103** (medium, duplicate, hard/large) `gromacs/formats.rs:88` — Four DCD and four XTC-named codecs; none reads real XTC
- [ ] **F2104** (medium, bug, easy/trivial) `gromacs/formats.rs:254` — XVG parser takes '@ legend on' options as legends, misses 's0 legend'
- [ ] **F2105** (medium, panic, easy/small) `gromacs/gro_file.rs:56` — GRO reader preallocates from count line; byte slicing panics
- [ ] **F2106** (medium, panic, medium/small) `gromacs/gro_file.rs:82` — GRO reader byte-slices (UTF-8 panic), drops triclinic box; writer overflows columns
- [ ] **F2107** (medium, panic, easy/small) `gromacs/gro_file.rs:88` — Fixed-column &str byte slicing panics on non-ASCII input
- [ ] **F2108** (medium, perf, easy/small) `gromacs/gro_file.rs:315` — GroTrajectory::parse quadratic, CRLF-desyncs, silently truncates on errors
- [ ] **F2110** (medium, bug, easy/small) `gromacs/topology.rs:1131` — [ atomtypes ] writer/parser omit ptype column; incompatible with GROMACS
- [ ] **F2112** (medium, bug, easy/trivial) `hdf5_io/convenience.rs:176` — ChunkIterator loops forever when a chunk dimension is zero
- [ ] **F2115** (medium, missing-feature, medium/large) `hdf5_io/file.rs:106` — Hdf5File has no on-disk save/load; checkpoints lost at exit
- [ ] **F2117** (medium, fabrication, hard/large) `hdf5_io/mod.rs:4` — HDF5 support is an in-memory mock; README and keywords advertise HDF5
- [ ] **F2118** (medium, panic, easy/medium) `hdf5_io/science_io.rs:192` — ~110 assert_eq! inside Hdf5Result-returning public writers
- [ ] **F2119** (medium, panic, easy/medium) `hdf5_io/science_io_md.rs:153` — ~110 assert macros panic inside Hdf5Result-returning public APIs
- [ ] **F2120** (medium, bug, easy/small) `hdf5_io/trajectory.rs:137` — TrajectoryStore drops velocities/forces/box on flush; read_velocities panics
- [ ] **F2121** (medium, bug, easy/small) `hdf5_io/trajectory.rs:217` — Checkpoint truncates step to i32 and silently drops scalars
- [ ] **F2123** (medium, api, easy/small) `hdf5_io/types.rs:14` — Hdf5Error and ComparisonError do not implement std::error::Error
- [ ] **F2124** (medium, api, medium/small) `hdf5_simple/functions.rs:486` — write_xdmf claims SHDF support but always emits Format="HDF"
- [ ] **F2129** (medium, panic, easy/trivial) `hdf5_simple/types.rs:977` — ShdfFile::from_bytes allocates from untrusted header counts (panic/abort)
- [ ] **F2131** (medium, bug, easy/small) `hdf5_simple/types.rs:1025` — SHDF Int64 dtype stored as i32; unchecked reserve
- [ ] **F2133** (medium, panic, easy/trivial) `hpc_io.rs:164` — Hdf5Dataset::read_slice panics when offset exceeds length
- [ ] **F2135** (medium, bug, medium/medium) `hpc_io.rs:276` — hpc_io Hdf5File mock: open() reads nothing, create_group returns root
- [ ] **F2136** (medium, bug, easy/small) `hpc_io.rs:369` — ParallelNetcdf::local_size ignores rank; panics with zero ranks
- [ ] **F2137** (medium, stub, easy/small) `hpc_io.rs:418` — AdiosWriter::perform_puts discards all staged data
- [ ] **F2138** (medium, panic, easy/small) `hpc_io.rs:574` — RestartFile::read_from_buffer panics on 24-31 byte input, drops metadata
- [ ] **F2140** (medium, panic, easy/small) `hpc_io.rs:683` — partition_regular_grid panics for uneven splits and zero ranks
- [ ] **F2142** (medium, bug, medium/small) `hpc_io.rs:905` — ScientificJson::to_json emits invalid JSON and has no reader
- [ ] **F2143** (medium, bug, easy/small) `json_io.rs:60` — JSON writer emits invalid JSON: raw keys, control chars, inf→1e308
- [ ] **F2144** (medium, bug, medium/small) `json_io.rs:499` — Streaming parser splits numbers at chunk boundary, hides syntax errors
- [ ] **F2147** (medium, bug, medium/medium) `json_io.rs:1047` — generate_patch paths incompatible with JsonPatch::apply flat-key logic
- [ ] **F2150** (medium, panic, easy/small) `lammps/types.rs:1150` — Binary dump reader panics on negative n_atoms; mislabelled as LAMMPS
- [ ] **F2152** (medium, bug, medium/small) `lammps_dump.rs:89` — Triclinic dump bounds treated as true box; volume/lengths inflated
- [ ] **F2153** (medium, bug, easy/small) `lammps_dump.rs:223` — Triclinic dump bounding box treated as true xlo/xhi
- [ ] **F2156** (medium, panic, easy/small) `lammps_dump.rs:433` — LammpsDumpWriter::write_frame panics or emits ragged rows on inconsistent frame
- [ ] **F2158** (medium, panic, easy/small) `lammps_dump.rs:712` — LAMMPS custom columns short or zero-filled; extract_forces/write_frame panic
- [ ] **F2159** (medium, bug, easy/small) `lattice_io.rs:84` — Lattice binary: Custom→D3Q19, unchecked alloc, per-f64 syscalls
- [ ] **F2160** (medium, perf, easy/small) `lattice_io.rs:174` — Lattice binary I/O does one syscall per f64; unchecked header allocation
- [ ] **F2161** (medium, bug, easy/small) `lattice_io.rs:207` — Lattice reader allocates nx*ny*nz nodes from header unchecked
- [ ] **F2162** (medium, bug, easy/trivial) `lattice_io.rs:415` — write_poscar emits atoms not grouped by species
- [ ] **F2163** (medium, bug, easy/small) `lattice_io.rs:429` — CIF writer assumes orthogonal cell; HCP written with wrong b and γ
- [ ] **F2165** (medium, api, medium/medium) `lib.rs:46` — Root glob re-exports leak helpers and silently shadow VtkWriter
- [ ] **F2167** (medium, panic, easy/small) `machine_learning_io.rs:89` — Tensor/StateDict from_bytes overflow and huge alloc on crafted input
- [ ] **F2168** (medium, bug, easy/small) `machine_learning_io.rs:190` — DenseLayer::forward uses input.len() as weight-row stride
- [ ] **F2169** (medium, bug, easy/small) `material_db.rs:123` — Water/air defaults put bulk modulus in E and volumetric CTE in α
- [ ] **F2170** (medium, bug, easy/small) `material_db_io.rs:442` — Material JSON writes {:.6}, zeroing small properties
- [ ] **F2177** (medium, fabrication, easy/small) `medical_imaging_io/types.rs:260` — NIfTI scl_slope ignored; unsupported dtypes become zeros
- [ ] **F2178** (medium, bug, medium/medium) `medical_imaging_io/types.rs:288` — NIfTI scl_slope/scl_inter and qform/sform rotation never applied
- [ ] **F2180** (medium, bug, easy/trivial) `medical_imaging_io/types.rs:527` — resample reports spacing inconsistent with its endpoint-aligned grid
- [ ] **F2181** (medium, bug, easy/trivial) `medical_imaging_io/types.rs:630` — Histogram bins inconsistent with centres; equalize panics on 0 bins
- [ ] **F2183** (medium, stub, medium/medium) `medical_io/types.rs:248` — medical_io NiftiReader/DicomReader/NiftiWriter are mocks
- [ ] **F2185** (medium, bug, easy/small) `medical_io/types.rs:988` — DicomDataSet::to_bytes truncates lengths above 65535 via `as u16`
- [ ] **F2186** (medium, bug, easy/small) `medical_io/types.rs:1094` — DVH bins mismatch dose axis; zero-bin/empty panics
- [ ] **F2188** (medium, bug, medium/small) `mesh_export.rs:650` — OBJ negative/invalid face indices silently become vertex 0
- [ ] **F2191** (medium, bug, medium/medium) `mesh_io.rs:141` — NASTRAN field8 overflows 8 cols; no continuation; slice panic
- [ ] **F2192** (medium, bug, medium/small) `mesh_io.rs:466` — glTF: reader takes last FLOAT VEC3 as POSITION; writer drops normals/uvs
- [ ] **F2193** (medium, panic, easy/small) `mesh_io.rs:719` — Gmsh MshReader panics on malformed elements; assumes contiguous tags
- [ ] **F2194** (medium, bug, easy/small) `mesh_io.rs:885` — Abaqus reader drops node ids; element continuation lines misparsed
- [ ] **F2195** (medium, panic, easy/small) `mesh_io.rs:1008` — Nastran parse panics on non-ASCII; field8 overflows 8 columns
- [ ] **F2196** (medium, perf, medium/small) `mesh_io.rs:1131` — MeshConverter::deduplicate_vertices is O(V^2) linear scan
- [ ] **F2200** (medium, fabrication, easy/trivial) `mesh_quality/functions.rs:487` — euler_characteristic assumes E=3F/2 instead of counting edges
- [ ] **F2201** (medium, panic, medium/small) `molecular_docking_io/functions.rs:76` — SMILES tokenizer char/byte index mix panics; formula ignores H
- [ ] **F2202** (medium, stub, hard/large) `molecular_docking_io/types.rs:529` — SMILES parser has no bonds/H; to_smiles fakes canonical output
- [ ] **F2203** (medium, bug, easy/trivial) `molecular_docking_io/types.rs:614` — PDBQT TORSDOF never parsed; 6-byte slice can panic
- [ ] **F2208** (medium, fabrication, easy/small) `molecular_visualization_io/types.rs:597` — DX/CHGCAR/CUBE readers zero-pad truncated volumetric grids
- [ ] **F2209** (medium, bug, easy/small) `molecular_visualization_io/types.rs:762` — POSCAR: Cartesian coords unscaled; negative (volume) scale inverts cell
- [ ] **F2210** (medium, bug, medium/small) `molecular_visualization_io/types.rs:1116` — CUBE reader mishandles negative voxel counts and MO DSET_IDS
- [ ] **F2211** (medium, perf, medium/small) `molecular_viz_io/types.rs:192` — auto_bonds is O(N^2) and allocates two Strings per pair
- [ ] **F2214** (medium, bug, easy/small) `netcdf/functions.rs:224` — CDL parser drops standard global attrs, unlimited size, non-double vars
- [ ] **F2215** (medium, fabrication, easy/trivial) `netcdf/functions.rs:1121` — read_particle_trajectory_cdl fabricates origin atoms for missing data
- [ ] **F2217** (medium, bug, easy/small) `netcdf/types.rs:135` — extend_unlimited can grow a fixed dimension of equal size
- [ ] **F2218** (medium, pure-rust, hard/large) `netcdf/types.rs:193` — No real NetCDF/HDF5 file I/O; Nc4File::to_bytes drops data
- [ ] **F2220** (medium, perf, easy/trivial) `netcdf/types.rs:528` — trajectory_write issues write syscalls per coordinate on unbuffered File
- [ ] **F2223** (medium, bug, easy/small) `numpy/functions.rs:220` — NPY reader ignores fortran_order, shape-product overflow; writer unchecked
- [ ] **F2225** (medium, missing-feature, medium/medium) `numpy/types.rs:534` — NPZ container is not ZIP; incompatible with numpy
- [ ] **F2226** (medium, pure-rust, medium/medium) `numpy/types.rs:544` — .npz is a custom container, not ZIP; numpy cannot read it
- [ ] **F2227** (medium, pure-rust, easy/small) `numpy/types.rs:596` — NpzWriter/NpzArchive emit a private container, not ZIP-based .npz
- [ ] **F2228** (medium, perf, easy/small) `obj/functions.rs:29` — weld_vertices uses an O(n²) linear search
- [ ] **F2229** (medium, panic, easy/small) `obj/types.rs:215` — OBJ face index 0 underflows; negative indices rejected
- [ ] **F2230** (medium, perf, medium/medium) `obj/types.rs:1252` — QEM decimation rescans all faces and edges every collapse
- [ ] **F2232** (medium, bug, medium/small) `openfoam/functions.rs:345` — Nonuniform list parser fails on inline lists; silent Uniform(0.0)
- [ ] **F2237** (medium, bug, easy/small) `openfoam_format.rs:237` — Uniform OpenFOAM fields silently become 1-cell fields
- [ ] **F2238** (medium, fabrication, easy/small) `openfoam_format.rs:246` — detect_n_cells invents a 1-cell field when the count is unknown
- [ ] **F2239** (medium, bug, easy/small) `openfoam_format.rs:335` — write_foam_field omits FoamFile header; dimensions hardcoded
- [ ] **F2242** (medium, bug, easy/small) `parallel_io.rs:553` — DistributedArray::assemble accepts gaps/overlaps; hex decoder can panic
- [ ] **F2243** (medium, bug, easy/small) `parallel_io.rs:559` — DistributedArray::assemble misses gaps/overlaps; text parse can panic
- [ ] **F2244** (medium, bug, easy/trivial) `parallel_io.rs:713` — merge_rank_files silently skips missing rank files
- [ ] **F2247** (medium, fabrication, medium/medium) `particle_data_io.rs:459` — H5part writer is not HDF5; drops ids, fields, attributes
- [ ] **F2248** (medium, fabrication, medium/small) `particle_data_io.rs:508` — ParticleVtkWriter.binary and with_binary() do nothing
- [ ] **F2255** (medium, perf, medium/small) `particle_formats.rs:582` — MSD is a direct O(nf²·np) sum; use an FFT method (oxifft)
- [ ] **F2257** (medium, bug, easy/small) `particle_formats.rs:913` — GRO reader whitespace-splits fixed columns; writer misaligned too
- [ ] **F2261** (medium, bug, easy/small) `pdb.rs:186` — PDB >99,999 atoms breaks round trip; non-ASCII slice panics
- [ ] **F2264** (medium, bug, easy/small) `plot3d_format.rs:127` — PLOT3D dims product overflows and drives huge preallocation
- [ ] **F2270** (medium, fabrication, hard/large) `point_cloud_io/functions.rs:662` — ball_pivoting_mesh is not Ball Pivoting; output always has holes
- [ ] **F2272** (medium, bug, medium/small) `quantum_chemistry_io.rs:672` — ORCA excited-state parser never matches real TDDFT output
- [ ] **F2273** (medium, bug, easy/trivial) `quantum_chemistry_io.rs:686` — ORCA total_energy double-counts dispersion correction
- [ ] **F2275** (medium, missing-feature, medium/large) `remote_sensing_io.rs:21` — GeoTiffReader cannot read GeoTIFF; constructor-fed container only
- [ ] **F2277** (medium, bug, medium/small) `remote_sensing_io.rs:226` — Slope ignores cell size; aspect/hillshade orientation wrong
- [ ] **F2278** (medium, bug, easy/trivial) `remote_sensing_io.rs:334` — D8 flow_direction sends pits and flats east (uphill)
- [ ] **F2280** (medium, perf, easy/trivial) `remote_sensing_io.rs:424` — watershed_delineation scans the whole grid per cell: O((hw)^2)
- [ ] **F2285** (medium, fabrication, medium/medium) `restart_io.rs:417` — Hdf5Like/MessagePack are relabelled binary; JSON restart unreadable
- [ ] **F2288** (medium, duplicate, medium/large) `restart_io.rs:742` — 76 duplicated public type names across overlapping modules
- [ ] **F2290** (medium, perf, easy/small) `restart_io.rs:815` — IncrementalRestart: O(n) dedupe per mark; delta has no indices
- [ ] **F2293** (medium, bug, easy/small) `scientific_formats.rs:112` — Cube reader mishandles negative counts; no data length check
- [ ] **F2294** (medium, fabrication, easy/trivial) `scientific_formats.rs:190` — Molden writer gives every atom atomic number 1 (H)
- [ ] **F2295** (medium, panic, easy/small) `seismic_io.rs:293` — SEG-Y/SAC header decoders panic on short buffers
- [ ] **F2297** (medium, fabrication, medium/medium) `seismic_io.rs:600` — Steim-1 encoder uses inverted codes and one diff per word
- [ ] **F2300** (medium, panic, easy/small) `sensor_data_io.rs:207` — Sensor binary decoders overflow size math and allocate untrusted counts
- [ ] **F2301** (medium, bug, medium/medium) `sensor_data_io.rs:826` — Thermocouple cold-junction compensation adds temperature, not EMF
- [ ] **F2302** (medium, bug, easy/small) `sensor_io.rs:313` — TelemetryProtocol decode never checks the CRC; encode truncates the length
- [ ] **F2303** (medium, bug, medium/small) `sensor_io.rs:395` — Complementary filter clamps roll to ±90° and drags yaw to zero
- [ ] **F2306** (medium, perf, easy/small) `sensor_io.rs:1158` — sensor_fft_spectrum is a naive O(N^2) DFT; oxifft unused
- [ ] **F2307** (medium, panic, easy/small) `sensor_io.rs:1216` — sensor_synchronize panics on NaN, is O(N·M) and is unbounded
- [ ] **F2308** (medium, bug, easy/trivial) `simulation_database.rs:296` — RLE and diff use absolute EPSILON, so small-magnitude fields are lossy
- [ ] **F2310** (medium, bug, easy/small) `simulation_database.rs:1115` — latin_hypercube_sample always yields segment midpoints; not stratified
- [ ] **F2311** (medium, bug, easy/small) `simulation_io.rs:59` — read_xdmf_timesteps cannot parse real XDMF; its test is vacuous
- [ ] **F2315** (medium, bug, easy/small) `simulation_log.rs:500` — SimSnapshot::to_json fixed {:.9} zeroes small values, emits NaN
- [ ] **F2318** (medium, bug, medium/small) `spectroscopy_io/types.rs:276` — ALS baseline uses first-difference penalty, not Eilers' D2
- [ ] **F2319** (medium, panic, easy/trivial) `spectroscopy_io/types.rs:381` — SpectralAnalysis::moving_average panics when half_window > n
- [ ] **F2321** (medium, bug, easy/small) `stl/functions.rs:155` — Absolute 1e-10 area threshold deletes valid micro-scale triangles
- [ ] **F2322** (medium, fabrication, medium/medium) `stl/functions.rs:474` — repair_close_holes adds zero-area caps that add boundary edges
- [ ] **F2324** (medium, bug, easy/small) `stl/functions.rs:564` — obj_to_stl underflows on index 0, fabricates origin, drops quads
- [ ] **F2327** (medium, bug, easy/small) `streaming_io/functions.rs:924` — Streaming XYZ reader rereads the whole file per frame; bad coords become 0
- [ ] **F2329** (medium, bug, medium/small) `streaming_io/types.rs:240` — FrameInterpolator CubicHermite drops tangent terms (smoothstep only)
- [ ] **F2337** (medium, bug, easy/trivial) `su2_format.rs:235` — read_su2 turns I/O errors into empty lines; can loop forever
- [ ] **F2339** (medium, bug, easy/small) `tecplot_format.rs:467` — Tecplot ZONE key lookup by substring; E= matches inside ZONETYPE=
- [ ] **F2340** (medium, pure-rust, easy/small) `time_series_io.rs:459` — Hand-rolled radix-2 FFT, pow2-only; oxifft unused in io crate
- [ ] **F2342** (medium, panic, easy/trivial) `time_series_io.rs:587` — welch_psd underflows when overlap > segment_len; n=1 window is NaN
- [ ] **F2343** (medium, panic, easy/small) `time_series_io.rs:634` — autocorrelation/cross_correlation panic on empty input; NaN when constant
- [ ] **F2344** (medium, bug, medium/small) `time_series_io.rs:716` — detect_peaks prominence uses whole-side minima; O(n^2)
- [ ] **F2346** (medium, panic, easy/small) `time_series_io.rs:938` — channel_from_json panics on truncated input; to_json emits Rust escapes
- [ ] **F2349** (medium, api, easy/small) `trajectory/types.rs:54` — integrate_vacf needs raw VACF but only normalised is public
- [ ] **F2350** (medium, panic, easy/trivial) `trajectory/types.rs:82` — XyzReader preallocates from untrusted atom count: panic or abort
- [ ] **F2355** (medium, panic, easy/small) `trajectory/types.rs:486` — TrajectoryResampler panics on short atom_types and fabricates timestep
- [ ] **F2357** (medium, bug, medium/medium) `visualization_io.rs:904` — glTF physics export maps static to isKinematic; non-Khronos schema
- [ ] **F2359** (medium, fabrication, easy/small) `visualization_io.rs:1374` — Viridis/Plasma colour maps use fabricated polynomial coefficients
- [ ] **F2360** (medium, bug, easy/small) `vtk/functions.rs:530` — Legacy VTK reader allocates from POINTS count, accepts truncation
- [ ] **F2366** (medium, missing-feature, medium/medium) `vtk_writer/types.rs:10` — XML VTU writers are ASCII-only; no binary or zlib-compressed output
- [ ] **F2367** (medium, missing-feature, medium/medium) `vtk_writer/types.rs:81` — VTU writers ASCII-only; binary helper unused; no zlib compression
- [ ] **F2373** (medium, bug, easy/trivial) `vtu/types.rs:279` — VTU 6-component tensor documented in xz/yz-swapped order vs VTK
- [ ] **F2374** (medium, perf, medium/medium) `vtu/types.rs:403` — VtuWriter: per-value String allocation, ASCII-only, whole-file String
- [ ] **F2375** (medium, bug, medium/medium) `vtu/types.rs:776` — VtuReader misparses points and discards connectivity
- [ ] **F2377** (medium, bug, easy/small) `wavefront_extended.rs:138` — OBJ quads/n-gons truncated; negative indices silently dropped
- [ ] **F2378** (medium, bug, medium/medium) `wavefront_extended.rs:252` — Mesh rigid-body mass multiplies area by kg/m³ density
- [ ] **F2379** (medium, panic, easy/small) `weather_data_io/functions.rs:57` — METAR parsers panic on non-ASCII tokens via byte slicing
- [ ] **F2385** (medium, bug, easy/small) `weather_data_io/types.rs:1184` — GRIB reference time uses 365-day years, 30-day months
- [ ] **F2390** (medium, bug, easy/small) `weather_io.rs:320` — NetcdfAtmDataset text round-trip loses packing attrs and precision
- [ ] **F2395** (medium, bug, easy/trivial) `weather_io.rs:1182` — isa_temperature extends +2.8 K/km gradient to 86 km
- [ ] **F2396** (medium, fabrication, medium/small) `weather_io.rs:1386` — Showalter index: discarded work, T_LCL=Td, fixed 6 K/km
- [ ] **F2400** (medium, bug, easy/small) `xarray_io.rs:228` — read_csv_xarray turns malformed rows into index 0 / value 0
- [ ] **F2402** (medium, fabrication, easy/small) `xarray_io.rs:311` — resample_linear claims trilinear interpolation but is nearest-neighbour
- [ ] **F2403** (medium, fabrication, medium/small) `xarray_io.rs:338` — resample_linear documented as trilinear but is nearest-neighbour
- [ ] **F2408** (medium, bug, easy/small) `xdmf/types.rs:184` — HDF5 attribute reference missing ':' separator; no connectivity
- [ ] **F2409** (medium, bug, easy/small) `xdmf/types.rs:718` — XDMF writers emit cell topology without connectivity, wrong count
- [ ] **F2410** (medium, panic, medium/small) `xdmf/types.rs:870` — XDMF Mixed topology panics via chunks(0)
- [ ] **F2413** (medium, missing-feature, hard/large) `xtc_dcd.rs:22` — No real GROMACS XTC codec; simplified LE format only
- [ ] **F2415** (medium, bug, easy/trivial) `xtc_dcd.rs:99` — dcd_to_xtc stamps a fabricated 3 nm box on every frame
- [ ] **F2417** (medium, perf, easy/small) `xtc_dcd.rs:267` — read_frame_at/read_metadata decode the whole trajectory despite docs
- [ ] **F2418** (medium, bug, easy/small) `xtc_dcd.rs:534` — DCD writer emits corrupt records when atom counts mismatch
- [ ] **F2420** (medium, fabrication, easy/trivial) `xtc_dcd.rs:737` — compute_rmsd returns 0.0 (perfect match) for mismatched inputs

### Low (54) — tracked in the external ledger
F1993, F1998, F2010, F2013, F2018, F2023, F2025, F2031, F2056, F2058, F2071, F2084, F2085, F2093, F2094, F2111, F2113, F2114, F2116, F2125, F2126, F2127, F2155, F2164, F2166, F2176, F2179, F2231, F2245, F2251, F2263, F2271, F2298, F2305, F2309, F2316, F2333, F2341, F2352, F2353, F2356, F2361, F2363, F2365, F2368, F2369, F2372, F2381, F2382, F2399, F2404, F2406, F2411, F2412

<!-- AUDIT-2026-10-12 END -->
