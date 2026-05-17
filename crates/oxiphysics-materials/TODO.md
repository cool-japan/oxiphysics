# oxiphysics-materials TODO

Last updated: 2026-05-17 / v0.1.1

## Phase 1: Foundation
- [x] Define core types and traits
- [x] Implement basic error handling
- [x] Add unit tests

## Phase 2: Core Implementation
- [x] Implement primary algorithms (75+ modules, 6,241 public items)
- [x] Add integration tests (4,486 tests total, 0 stubs)
- [x] Elastic, plasticity, hyperelastic, viscoelastic, creep, fatigue, fracture, damage
- [x] Composite / fiber / nanocomposite modules
- [x] Geological, geomechanics, geomaterial_models, porous_media
- [x] Biological, biomaterials, polymer_physics, polymer_mechanics
- [x] Smart materials: shape_memory, shape_memory_alloy, magnetocaloric_materials, metamaterials
- [x] Energy: battery_materials, energy_materials, electrochemistry, thermoelectrics, hydrogen_storage
- [x] Electronic/optical: semiconductor, dielectric, optical, optical_materials, quantum_materials, superconductor
- [x] Thermal, radiation, radiation_shielding, nuclear_materials
- [x] Manufacturing: additive_manufacturing, aerospace, alloy, foam, nano/nanomaterials
- [x] Tribology, tribology_ext, corrosion, multiphysics
- [x] EOS, phase_transform, crystal_plasticity, acoustics
- [x] Contact presets: ContactMaterialPair, FrictionCombineRule, ModulusCombineRule, RestitutionCombineRule
- [x] Performance benchmarks (basic)

## Phase 3: Polish
- [x] Documentation (rustdoc on all public items)
- [x] Extended examples (additional usage examples) — see `materials_bench` module
- [x] Benchmark suite expansion (`materials_bench` — `MatBenchHarness`, elastic/neo-Hookean/viscoplastic/fatigue/thermal batch kernel timing)
- [x] SIMD acceleration paths (`simd_paths` — SoA layout, elastic/neo-Hookean/viscoplastic/fatigue/return-mapping batch kernels)
