"""Tests for core PhysicsWorld lifecycle.

These tests verify the fundamental PhysicsWorld API: construction, stepping,
body addition, and contact result queries.

Run after ``maturin develop --release`` via::

    python -m pytest python/tests/test_world.py -v
"""

import pytest

try:
    import oxiphysics as ox
    HAS_OXIPHYSICS = True
except ImportError:
    HAS_OXIPHYSICS = False

pytestmark = pytest.mark.skipif(not HAS_OXIPHYSICS, reason="oxiphysics not built")


def test_physics_world_creates_with_earth_gravity():
    world = ox.PhysicsWorld.with_earth_gravity()
    assert world is not None


def test_physics_world_creates_with_zero_gravity():
    world = ox.PhysicsWorld(0.0, 0.0, 0.0)
    assert world is not None


def test_physics_world_body_count_starts_zero():
    world = ox.PhysicsWorld.with_earth_gravity()
    assert world.body_count() == 0


def test_physics_world_step_does_not_crash():
    world = ox.PhysicsWorld.with_earth_gravity()
    world.step(1.0 / 60.0)


def test_add_sphere_body():
    world = ox.PhysicsWorld.with_earth_gravity()
    handle = world.add_sphere_body(1.0, 0.0, 10.0, 0.0, 0.5)
    assert handle >= 0


def test_body_count_increments_after_add():
    world = ox.PhysicsWorld.with_earth_gravity()
    world.add_sphere_body(1.0, 0.0, 10.0, 0.0, 0.5)
    assert world.body_count() == 1


def test_get_position_returns_tuple():
    world = ox.PhysicsWorld.with_earth_gravity()
    handle = world.add_sphere_body(1.0, 0.0, 10.0, 0.0, 0.5)
    pos = world.get_position(handle)
    assert pos is not None
    assert len(pos) == 3


def test_step_moves_body_under_gravity():
    world = ox.PhysicsWorld.with_earth_gravity()
    handle = world.add_sphere_body(1.0, 0.0, 10.0, 0.0, 0.5)
    y_before = world.get_position(handle)[1]
    for _ in range(60):
        world.step(1.0 / 60.0)
    y_after = world.get_position(handle)[1]
    assert y_after < y_before, "Body should fall under Earth gravity"


def test_add_rigid_body_via_config():
    world = ox.PhysicsWorld.with_earth_gravity()
    cfg = ox.RigidBodyConfig.dynamic(1.0, 0.0, 5.0, 0.0)
    body_id = world.add_rigid_body(cfg)
    assert body_id >= 0


def test_sim_config_earth_gravity():
    cfg = ox.SimConfig.earth_gravity()
    assert cfg is not None


def test_sim_config_zero_gravity():
    cfg = ox.SimConfig.zero_gravity()
    assert cfg is not None
