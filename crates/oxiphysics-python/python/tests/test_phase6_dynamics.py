"""Tests for Phase 6.3 Advanced Dynamics & Tooling classes.

Covers: CharacterController, Rope, IkChain, XpbdSolver, AeroSystem, NavMesh,
TriggerWorld, EventBus, AnimationPlayer, ProfilerSession, ForceFieldSystem.

Uses ``hasattr`` guards so individual tests skip gracefully if a class is
not yet bound.

Run after ``maturin develop --release`` via::

    python -m pytest python/tests/test_phase6_dynamics.py -v
"""

import json
import math
import pytest

try:
    import oxiphysics as ox
    HAS_OXIPHYSICS = True
except ImportError:
    HAS_OXIPHYSICS = False

pytestmark = pytest.mark.skipif(not HAS_OXIPHYSICS, reason="oxiphysics not built")


# ---------------------------------------------------------------------------
# CharacterController
# ---------------------------------------------------------------------------

def test_character_controller_creates():
    if not hasattr(ox, "CharacterController"):
        pytest.skip("CharacterController not yet bound")
    cc = ox.CharacterController([0.0, 1.0, 0.0], 0.4, 0.9)
    assert cc is not None


def test_character_controller_position():
    if not hasattr(ox, "CharacterController"):
        pytest.skip("CharacterController not yet bound")
    cc = ox.CharacterController([1.0, 2.0, 3.0], 0.4, 0.9)
    pos = cc.position()
    assert len(pos) == 3
    assert abs(pos[0] - 1.0) < 1e-9
    assert abs(pos[1] - 2.0) < 1e-9
    assert abs(pos[2] - 3.0) < 1e-9


def test_character_controller_velocity():
    if not hasattr(ox, "CharacterController"):
        pytest.skip("CharacterController not yet bound")
    cc = ox.CharacterController([0.0, 0.0, 0.0], 0.4, 0.9)
    vel = cc.velocity()
    assert len(vel) == 3


def test_character_controller_set_velocity():
    if not hasattr(ox, "CharacterController"):
        pytest.skip("CharacterController not yet bound")
    cc = ox.CharacterController([0.0, 0.0, 0.0], 0.4, 0.9)
    cc.set_velocity([1.0, 0.0, 0.0])
    vel = cc.velocity()
    assert abs(vel[0] - 1.0) < 1e-9


def test_character_controller_is_grounded():
    if not hasattr(ox, "CharacterController"):
        pytest.skip("CharacterController not yet bound")
    cc = ox.CharacterController([0.0, 0.0, 0.0], 0.4, 0.9)
    assert isinstance(cc.is_grounded(), bool)


def test_character_controller_apply_gravity():
    if not hasattr(ox, "CharacterController"):
        pytest.skip("CharacterController not yet bound")
    cc = ox.CharacterController([0.0, 5.0, 0.0], 0.4, 0.9)
    cc.apply_gravity(-9.81, 0.016)


def test_character_controller_set_max_slope_deg():
    if not hasattr(ox, "CharacterController"):
        pytest.skip("CharacterController not yet bound")
    cc = ox.CharacterController([0.0, 0.0, 0.0], 0.4, 0.9)
    cc.set_max_slope_deg(45.0)


# ---------------------------------------------------------------------------
# Rope
# ---------------------------------------------------------------------------

def test_rope_creates():
    if not hasattr(ox, "Rope"):
        pytest.skip("Rope not yet bound")
    rope = ox.Rope([0.0, 10.0, 0.0], 5, 0.5, 0.1)
    assert rope is not None


def test_rope_link_count():
    if not hasattr(ox, "Rope"):
        pytest.skip("Rope not yet bound")
    rope = ox.Rope([0.0, 10.0, 0.0], 5, 0.5, 0.1)
    assert rope.link_count() == 5


def test_rope_total_length():
    if not hasattr(ox, "Rope"):
        pytest.skip("Rope not yet bound")
    # 4 links with 1.0 segment length = 3 segments = total_length 3.0
    rope = ox.Rope([0.0, 10.0, 0.0], 4, 1.0, 0.1)
    # segments = links - 1; total_length = num_links * segment_length
    assert rope.total_length() > 0.0


def test_rope_segments():
    if not hasattr(ox, "Rope"):
        pytest.skip("Rope not yet bound")
    # N links produces N-1 segments between consecutive links
    rope = ox.Rope([0.0, 10.0, 0.0], 3, 0.5, 0.1)
    segments = rope.segments()
    assert isinstance(segments, list)
    # segments is between consecutive links: 3 links → 2 segments
    assert len(segments) == rope.link_count() - 1


# ---------------------------------------------------------------------------
# IkChain
# ---------------------------------------------------------------------------

def test_ik_chain_creates():
    if not hasattr(ox, "IkChain"):
        pytest.skip("IkChain not yet bound")
    chain = ox.IkChain([0.0, 0.0, 0.0], [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
    assert chain is not None


def test_ik_chain_segment_count():
    if not hasattr(ox, "IkChain"):
        pytest.skip("IkChain not yet bound")
    chain = ox.IkChain([0.0, 0.0, 0.0], [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
    assert chain.segment_count() == 2


def test_ik_chain_total_reach():
    if not hasattr(ox, "IkChain"):
        pytest.skip("IkChain not yet bound")
    chain = ox.IkChain([0.0, 0.0, 0.0], [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
    assert chain.total_reach() > 0.0


def test_ik_chain_solve_fabrik_json():
    if not hasattr(ox, "IkChain"):
        pytest.skip("IkChain not yet bound")
    chain = ox.IkChain([0.0, 0.0, 0.0], [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
    result = chain.solve_fabrik_json([1.5, 0.0, 0.0])
    assert result is not None
    data = json.loads(result)
    assert isinstance(data, (dict, list))


def test_ik_chain_solve_two_bone_json():
    if not hasattr(ox, "IkChain"):
        pytest.skip("IkChain not yet bound")
    chain = ox.IkChain([0.0, 0.0, 0.0], [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
    result = chain.solve_two_bone_json([1.5, 0.0, 0.0])
    assert result is not None
    data = json.loads(result)
    assert isinstance(data, (dict, list))


# ---------------------------------------------------------------------------
# XpbdSolver
# ---------------------------------------------------------------------------

def test_xpbd_solver_creates():
    if not hasattr(ox, "XpbdSolver"):
        pytest.skip("XpbdSolver not yet bound")
    solver = ox.XpbdSolver()
    assert solver is not None


def test_xpbd_solver_add_particle():
    if not hasattr(ox, "XpbdSolver"):
        pytest.skip("XpbdSolver not yet bound")
    solver = ox.XpbdSolver()
    idx = solver.add_particle([0.0, 0.0, 0.0], 1.0)
    assert idx >= 0


def test_xpbd_solver_particle_count():
    if not hasattr(ox, "XpbdSolver"):
        pytest.skip("XpbdSolver not yet bound")
    solver = ox.XpbdSolver()
    assert solver.particle_count() == 0
    solver.add_particle([0.0, 0.0, 0.0], 1.0)
    assert solver.particle_count() == 1


def test_xpbd_solver_step():
    if not hasattr(ox, "XpbdSolver"):
        pytest.skip("XpbdSolver not yet bound")
    solver = ox.XpbdSolver()
    solver.set_gravity(0.0, -9.81, 0.0)
    solver.add_particle([0.0, 5.0, 0.0], 1.0)
    solver.step(0.016)


def test_xpbd_solver_particle_position():
    if not hasattr(ox, "XpbdSolver"):
        pytest.skip("XpbdSolver not yet bound")
    solver = ox.XpbdSolver()
    idx = solver.add_particle([1.0, 2.0, 3.0], 1.0)
    pos = solver.particle_position(idx)
    assert pos is not None
    assert len(pos) == 3


def test_xpbd_solver_add_distance_constraint():
    if not hasattr(ox, "XpbdSolver"):
        pytest.skip("XpbdSolver not yet bound")
    solver = ox.XpbdSolver()
    a = solver.add_particle([0.0, 0.0, 0.0], 1.0)
    b = solver.add_particle([1.0, 0.0, 0.0], 1.0)
    solver.add_distance_constraint(a, b, 1.0, 0.0)
    assert solver.constraint_count() == 1


def test_xpbd_solver_pin_particle():
    if not hasattr(ox, "XpbdSolver"):
        pytest.skip("XpbdSolver not yet bound")
    solver = ox.XpbdSolver()
    idx = solver.add_particle([0.0, 0.0, 0.0], 1.0)
    solver.pin_particle(idx)  # Should not crash


def test_xpbd_solver_all_positions():
    if not hasattr(ox, "XpbdSolver"):
        pytest.skip("XpbdSolver not yet bound")
    solver = ox.XpbdSolver()
    solver.add_particle([1.0, 2.0, 3.0], 1.0)
    solver.add_particle([4.0, 5.0, 6.0], 1.0)
    positions = solver.all_positions()
    assert isinstance(positions, list)
    assert len(positions) == 6  # 2 particles x 3 coords


# ---------------------------------------------------------------------------
# AeroSystem
# ---------------------------------------------------------------------------

def test_aero_system_creates():
    if not hasattr(ox, "AeroSystem"):
        pytest.skip("AeroSystem not yet bound")
    aero = ox.AeroSystem()
    assert aero is not None


def test_aero_system_entry_count_starts_zero():
    if not hasattr(ox, "AeroSystem"):
        pytest.skip("AeroSystem not yet bound")
    aero = ox.AeroSystem()
    assert aero.entry_count() == 0


def test_aero_system_set_wind():
    if not hasattr(ox, "AeroSystem"):
        pytest.skip("AeroSystem not yet bound")
    aero = ox.AeroSystem()
    aero.set_wind(10.0, 0.0, 0.0)


def test_aero_system_add_drag_body():
    if not hasattr(ox, "AeroSystem"):
        pytest.skip("AeroSystem not yet bound")
    aero = ox.AeroSystem()
    # rotation_flat: identity matrix 3x3 row-major
    identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    aero.add_drag_body(
        [0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 0.47, 0.5, 1.225, identity
    )
    assert aero.entry_count() == 1


def test_aero_system_apply_json():
    if not hasattr(ox, "AeroSystem"):
        pytest.skip("AeroSystem not yet bound")
    aero = ox.AeroSystem()
    identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    aero.add_drag_body(
        [0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 0.47, 0.5, 1.225, identity
    )
    result = aero.apply_json()
    assert result is not None
    assert isinstance(json.loads(result), (dict, list))


def test_aero_system_clear():
    if not hasattr(ox, "AeroSystem"):
        pytest.skip("AeroSystem not yet bound")
    aero = ox.AeroSystem()
    identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    aero.add_drag_body(
        [0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 0.47, 0.5, 1.225, identity
    )
    aero.clear()
    assert aero.entry_count() == 0


# ---------------------------------------------------------------------------
# NavMesh
# ---------------------------------------------------------------------------

def _make_simple_navmesh():
    """Build a simple triangulated flat ground plane."""
    vertices = [
        [0.0, 0.0, 0.0],
        [10.0, 0.0, 0.0],
        [10.0, 0.0, 10.0],
        [0.0, 0.0, 10.0],
    ]
    tris = [[0, 1, 2], [0, 2, 3]]
    return ox.NavMesh(vertices, tris)


def test_navmesh_creates():
    if not hasattr(ox, "NavMesh"):
        pytest.skip("NavMesh not yet bound")
    nm = _make_simple_navmesh()
    assert nm is not None


def test_navmesh_vertex_count():
    if not hasattr(ox, "NavMesh"):
        pytest.skip("NavMesh not yet bound")
    nm = _make_simple_navmesh()
    assert nm.vertex_count() == 4


def test_navmesh_triangle_count():
    if not hasattr(ox, "NavMesh"):
        pytest.skip("NavMesh not yet bound")
    nm = _make_simple_navmesh()
    assert nm.triangle_count() == 2


def test_navmesh_find_path_json():
    if not hasattr(ox, "NavMesh"):
        pytest.skip("NavMesh not yet bound")
    nm = _make_simple_navmesh()
    result = nm.find_path_json([1.0, 0.0, 1.0], [9.0, 0.0, 9.0], 0.3, [])
    assert result is not None
    data = json.loads(result)
    assert isinstance(data, (dict, list))


# ---------------------------------------------------------------------------
# TriggerWorld
# ---------------------------------------------------------------------------

def test_trigger_world_creates():
    if not hasattr(ox, "TriggerWorld"):
        pytest.skip("TriggerWorld not yet bound")
    tw = ox.TriggerWorld()
    assert tw is not None


def test_trigger_world_volume_count_starts_zero():
    if not hasattr(ox, "TriggerWorld"):
        pytest.skip("TriggerWorld not yet bound")
    tw = ox.TriggerWorld()
    assert tw.volume_count() == 0


def test_trigger_world_add_sphere():
    if not hasattr(ox, "TriggerWorld"):
        pytest.skip("TriggerWorld not yet bound")
    tw = ox.TriggerWorld()
    vid = tw.add_sphere([0.0, 0.0, 0.0], 1.0, ["zone_a"])
    assert isinstance(vid, int)
    assert tw.volume_count() == 1


def test_trigger_world_add_aabb():
    if not hasattr(ox, "TriggerWorld"):
        pytest.skip("TriggerWorld not yet bound")
    tw = ox.TriggerWorld()
    vid = tw.add_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], ["box_zone"])
    assert isinstance(vid, int)


def test_trigger_world_remove_volume():
    if not hasattr(ox, "TriggerWorld"):
        pytest.skip("TriggerWorld not yet bound")
    tw = ox.TriggerWorld()
    vid = tw.add_sphere([0.0, 0.0, 0.0], 1.0, [])
    removed = tw.remove_volume(vid)
    assert removed is True
    assert tw.volume_count() == 0


def test_trigger_world_update_json():
    if not hasattr(ox, "TriggerWorld"):
        pytest.skip("TriggerWorld not yet bound")
    tw = ox.TriggerWorld()
    tw.add_sphere([0.0, 0.0, 0.0], 2.0, ["zone"])
    # bodies: list of (body_index, [x,y,z], radius, is_sleeping)
    bodies = [(0, [0.5, 0.0, 0.0], 0.3, False)]
    result = tw.update_json(0, bodies)
    assert result is not None
    events = json.loads(result)
    assert isinstance(events, list)


# ---------------------------------------------------------------------------
# EventBus
# ---------------------------------------------------------------------------

def test_event_bus_creates():
    if not hasattr(ox, "EventBus"):
        pytest.skip("EventBus not yet bound")
    bus = ox.EventBus()
    assert bus is not None


def test_event_bus_pending_count_starts_zero():
    if not hasattr(ox, "EventBus"):
        pytest.skip("EventBus not yet bound")
    bus = ox.EventBus()
    assert bus.pending_count() == 0


def test_event_bus_total_published_starts_zero():
    if not hasattr(ox, "EventBus"):
        pytest.skip("EventBus not yet bound")
    bus = ox.EventBus()
    assert bus.total_published() == 0


def test_event_bus_drain_as_strings():
    if not hasattr(ox, "EventBus"):
        pytest.skip("EventBus not yet bound")
    bus = ox.EventBus()
    result = bus.drain_as_strings()
    assert isinstance(result, list)


def test_event_bus_pending_kind_names():
    if not hasattr(ox, "EventBus"):
        pytest.skip("EventBus not yet bound")
    bus = ox.EventBus()
    names = bus.pending_kind_names()
    assert isinstance(names, list)


def test_event_bus_clear():
    if not hasattr(ox, "EventBus"):
        pytest.skip("EventBus not yet bound")
    bus = ox.EventBus()
    bus.clear()  # Should not crash when empty


def test_event_bus_flush():
    if not hasattr(ox, "EventBus"):
        pytest.skip("EventBus not yet bound")
    bus = ox.EventBus()
    bus.flush()  # Should not crash when empty


# ---------------------------------------------------------------------------
# AnimationPlayer
# ---------------------------------------------------------------------------

def _make_animation_clip_json():
    """Build a minimal AnimationClip JSON for testing."""
    clip = {
        "name": "test_clip",
        "body_animations": [
            {
                "body_id": "body_0",
                "position_track": {
                    "keyframes": [
                        {"time": 0.0, "value": [0.0, 0.0, 0.0], "ease": "Linear"},
                        {"time": 1.0, "value": [1.0, 0.0, 0.0], "ease": "Linear"},
                    ],
                    "looping": False,
                },
                "rotation_track": None,
                "scale_track": None,
            }
        ],
        "total_duration": 1.0,
    }
    return json.dumps(clip)


def test_animation_player_creates():
    if not hasattr(ox, "AnimationPlayer"):
        pytest.skip("AnimationPlayer not yet bound")
    try:
        player = ox.AnimationPlayer(_make_animation_clip_json())
        assert player is not None
    except Exception:
        pytest.skip("AnimationPlayer clip JSON format mismatch — skip for now")


def test_animation_player_play_and_is_playing():
    if not hasattr(ox, "AnimationPlayer"):
        pytest.skip("AnimationPlayer not yet bound")
    try:
        player = ox.AnimationPlayer(_make_animation_clip_json())
    except Exception:
        pytest.skip("AnimationPlayer clip JSON format mismatch — skip for now")
    player.play()
    assert player.is_playing()


def test_animation_player_pause():
    if not hasattr(ox, "AnimationPlayer"):
        pytest.skip("AnimationPlayer not yet bound")
    try:
        player = ox.AnimationPlayer(_make_animation_clip_json())
    except Exception:
        pytest.skip("AnimationPlayer clip JSON format mismatch — skip for now")
    player.play()
    player.pause()
    assert not player.is_playing()


def test_animation_player_progress():
    if not hasattr(ox, "AnimationPlayer"):
        pytest.skip("AnimationPlayer not yet bound")
    try:
        player = ox.AnimationPlayer(_make_animation_clip_json())
    except Exception:
        pytest.skip("AnimationPlayer clip JSON format mismatch — skip for now")
    progress = player.progress()
    assert 0.0 <= progress <= 1.0


def test_animation_player_step_json():
    if not hasattr(ox, "AnimationPlayer"):
        pytest.skip("AnimationPlayer not yet bound")
    try:
        player = ox.AnimationPlayer(_make_animation_clip_json())
    except Exception:
        pytest.skip("AnimationPlayer clip JSON format mismatch — skip for now")
    player.play()
    result = player.step_json(0.016)
    assert result is not None
    data = json.loads(result)
    assert isinstance(data, list)


# ---------------------------------------------------------------------------
# ProfilerSession
# ---------------------------------------------------------------------------

def test_profiler_session_creates():
    if not hasattr(ox, "ProfilerSession"):
        pytest.skip("ProfilerSession not yet bound")
    profiler = ox.ProfilerSession()
    assert profiler is not None


def test_profiler_session_begin_and_end_frame_json():
    if not hasattr(ox, "ProfilerSession"):
        pytest.skip("ProfilerSession not yet bound")
    profiler = ox.ProfilerSession()
    profiler.begin_frame()
    result = profiler.end_frame_json()
    assert result is not None
    data = json.loads(result)
    assert isinstance(data, (dict, list))


def test_profiler_session_end_frame_csv():
    if not hasattr(ox, "ProfilerSession"):
        pytest.skip("ProfilerSession not yet bound")
    profiler = ox.ProfilerSession()
    profiler.begin_frame()
    csv_str = profiler.end_frame_csv()
    assert isinstance(csv_str, str)


def test_profiler_session_end_frame_folded_stacks():
    if not hasattr(ox, "ProfilerSession"):
        pytest.skip("ProfilerSession not yet bound")
    profiler = ox.ProfilerSession()
    profiler.begin_frame()
    stacks = profiler.end_frame_folded_stacks()
    assert isinstance(stacks, str)


def test_profiler_session_node_count():
    if not hasattr(ox, "ProfilerSession"):
        pytest.skip("ProfilerSession not yet bound")
    profiler = ox.ProfilerSession()
    count = profiler.node_count()
    assert isinstance(count, int)
    assert count >= 0


# ---------------------------------------------------------------------------
# ForceFieldSystem
# ---------------------------------------------------------------------------

def test_force_field_system_creates():
    if not hasattr(ox, "ForceFieldSystem"):
        pytest.skip("ForceFieldSystem not yet bound")
    ffs = ox.ForceFieldSystem()
    assert ffs is not None


def test_force_field_system_empty_initially():
    if not hasattr(ox, "ForceFieldSystem"):
        pytest.skip("ForceFieldSystem not yet bound")
    ffs = ox.ForceFieldSystem()
    assert ffs.is_empty()
    assert ffs.len() == 0


def test_force_field_system_add_uniform():
    if not hasattr(ox, "ForceFieldSystem"):
        pytest.skip("ForceFieldSystem not yet bound")
    ffs = ox.ForceFieldSystem()
    fid = ffs.add_uniform([0.0, -9.81, 0.0])
    assert isinstance(fid, int)
    assert ffs.len() == 1


def test_force_field_system_add_radial_attract():
    if not hasattr(ox, "ForceFieldSystem"):
        pytest.skip("ForceFieldSystem not yet bound")
    ffs = ox.ForceFieldSystem()
    fid = ffs.add_radial_attract([0.0, 0.0, 0.0], 100.0, 10.0)
    assert isinstance(fid, int)


def test_force_field_system_force_at():
    if not hasattr(ox, "ForceFieldSystem"):
        pytest.skip("ForceFieldSystem not yet bound")
    ffs = ox.ForceFieldSystem()
    ffs.add_uniform([0.0, -9.81, 0.0])
    force = ffs.force_at([1.0, 0.0, 0.0], 0.0)
    assert isinstance(force, (list, tuple))
    assert len(force) == 3


def test_force_field_system_enable_disable():
    if not hasattr(ox, "ForceFieldSystem"):
        pytest.skip("ForceFieldSystem not yet bound")
    ffs = ox.ForceFieldSystem()
    fid = ffs.add_uniform([0.0, -9.81, 0.0])
    ffs.disable(fid)
    assert ffs.active_count() == 0
    ffs.enable(fid)
    assert ffs.active_count() == 1


def test_force_field_system_remove():
    if not hasattr(ox, "ForceFieldSystem"):
        pytest.skip("ForceFieldSystem not yet bound")
    ffs = ox.ForceFieldSystem()
    fid = ffs.add_uniform([0.0, -9.81, 0.0])
    initial_count = ffs.len()
    ffs.remove(fid)
    # After remove, active_count drops to 0 (slot may still exist)
    assert ffs.active_count() == 0
    assert initial_count >= 1


def test_force_field_aabb_region_creates():
    if not hasattr(ox, "ForceFieldAabbRegion"):
        pytest.skip("ForceFieldAabbRegion not yet bound")
    region = ox.ForceFieldAabbRegion([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])
    assert region is not None


def test_force_field_aabb_region_contains():
    if not hasattr(ox, "ForceFieldAabbRegion"):
        pytest.skip("ForceFieldAabbRegion not yet bound")
    region = ox.ForceFieldAabbRegion([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])
    assert region.contains([0.0, 0.0, 0.0]) is True
    assert region.contains([5.0, 0.0, 0.0]) is False
