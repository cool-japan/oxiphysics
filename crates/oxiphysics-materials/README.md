# oxiphysics-materials

**Status: [Alpha]** — v0.1.1 (2026-05-17)

[![Tests](https://img.shields.io/badge/tests-4486-brightgreen)](https://github.com/cool-japan/oxiphysics)
[![docs.rs](https://img.shields.io/docsrs/oxiphysics-materials)](https://docs.rs/oxiphysics-materials)

Comprehensive material science and constitutive model library for the [OxiPhysics](https://github.com/cool-japan/oxiphysics) engine.
Provides 75+ domain modules covering virtually every material class used in physics simulation.

## Features

- **Structural materials**: elastic, plasticity, hyperelastic, viscoelastic, creep, fatigue, fracture, damage mechanics
- **Composites & fibers**: composite, composite_materials, composite_failure, fiber_composites, nanocomposites
- **Geological & geomechanics**: geological, geomaterial_models, geomechanics, soil, porous_media
- **Biological & biomedical**: biological_materials, biomaterials, polymer_physics, polymer_mechanics
- **Advanced materials**: smart_materials, shape_memory, shape_memory_alloy, magnetocaloric_materials, metamaterials
- **Energy & electrochemistry**: battery_materials, energy_materials, electrochemistry, thermoelectrics, hydrogen_storage
- **Electronic & optical**: semiconductor, dielectric_materials, optical, optical_materials, quantum_materials, superconductor
- **Thermal & radiation**: thermal, radiation, radiation_shielding, nuclear_materials
- **Manufacturing**: additive_manufacturing, aerospace_materials, alloy_materials, foam_materials, nano_materials, nanomaterials
- **Surface & tribology**: tribology, tribology_ext, corrosion, multiphysics
- **Equations of state**: eos, phase_transform, crystal_plasticity, acoustic (acoustics module)
- **Contact presets**: `ContactMaterialPair`, `FrictionCombineRule`, `ModulusCombineRule`, `RestitutionCombineRule`

## Key Exports

```rust
use oxiphysics_materials::{
    ContactMaterialPair, FrictionCombineRule, ModulusCombineRule, RestitutionCombineRule,
};
```

## Notes

- **Standalone**: no dependency on `oxiphysics-core`; can be used in any Rust project
- 6,241 public API items across 75+ modules
- 4,486 tests — 0 stubs

## License

Apache-2.0 — Part of the [OxiPhysics](https://github.com/cool-japan/oxiphysics) project by COOLJAPAN OU (Team Kitasan)
