"""Integration tests for oxiphysics Python bindings.

Run after `maturin develop` or `pip install -e .`:
    pytest crates/oxiphysics-python/tests/

These are import and smoke tests only — they verify the wheel loads and
classes instantiate without checking physics correctness.
"""

import pytest

try:
    import oxiphysics
    HAS_OXIPHYSICS = True
except ImportError:
    HAS_OXIPHYSICS = False

pytestmark = pytest.mark.skipif(
    not HAS_OXIPHYSICS,
    reason="oxiphysics wheel not installed — run `maturin develop` first"
)


def test_module_imports():
    """The module must be importable."""
    import oxiphysics  # noqa: F401


def test_physics_world_instantiation():
    from oxiphysics import PhysicsWorld, SimConfig
    cfg = SimConfig()
    world = PhysicsWorld(cfg)
    assert world is not None


def test_sim_config_defaults():
    from oxiphysics import SimConfig
    cfg = SimConfig()
    assert cfg is not None


def test_rigid_body_config():
    from oxiphysics import RigidBodyConfig
    cfg = RigidBodyConfig()
    assert cfg is not None
