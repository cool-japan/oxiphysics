"""
Smoke tests for oxiphysics Python bindings.
These verify the compiled cdylib loads and the top-level API is accessible.
Run via `python -m pytest python/tests/` after `maturin develop`.
"""
import sys
import os

# When running under maturin develop, oxiphysics is importable directly.
# When running in CI, the cdylib must be on sys.path or PYTHONPATH.
try:
    import oxiphysics
except ImportError:
    # Skip gracefully if the extension module is not built yet.
    import pytest
    pytest.skip("oxiphysics extension module not built", allow_module_level=True)


def test_version():
    assert hasattr(oxiphysics, "__version__")
    assert isinstance(oxiphysics.__version__, str)


def test_sim_config_constructors():
    cfg = oxiphysics.SimConfig.earth_gravity()
    assert cfg is not None
    cfg2 = oxiphysics.SimConfig.zero_gravity()
    assert cfg2 is not None


def test_rigid_body_config_constructors():
    body = oxiphysics.RigidBodyConfig.dynamic(1.0, 0.0, 5.0, 0.0)
    assert body is not None
    static = oxiphysics.RigidBodyConfig.static_body(0.0, 0.0, 0.0)
    assert static is not None


def test_contact_result_default():
    # ContactResult is not directly constructible from Python in the current API;
    # verify it is importable.
    assert hasattr(oxiphysics, "ContactResult")


def test_physics_world_basic():
    world = oxiphysics.PhysicsWorld.with_earth_gravity()
    # body_count is a method, not a property
    assert world.body_count() == 0
    handle = world.add_sphere_body(1.0, 0.0, 10.0, 0.0, 0.5)
    assert handle >= 0
    world.step(0.016)
    pos = world.get_position(handle)
    # get_position returns a tuple of 3 floats, or None
    assert pos is not None and len(pos) == 3


def test_submodules_exist():
    """Verify that sub-modules are importable (they may be empty until Wave 2 completes)."""
    submodules = [
        "analytics", "constraints", "fem", "geometry",
        "io", "lbm", "materials", "md",
        "rigid", "sph", "vehicle", "viz", "world",
    ]
    for name in submodules:
        assert hasattr(oxiphysics, name), f"oxiphysics.{name} sub-module missing"
