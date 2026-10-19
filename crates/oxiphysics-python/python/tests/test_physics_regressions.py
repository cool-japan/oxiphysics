"""Regression tests of the physics behind the Python bindings.

Each test pins one audited defect of the bindings against an independent
reference: an analytic solution, a brute-force sum or a conservation law.

Covers: LBM body forcing (D2Q9 and D3Q19 walled channels), the binary VTK
writer, ConvexHull / GJK, the XPBD volume constraint and the molecular dynamics
neighbour list and Nose-Hoover thermostat.

Run after ``maturin develop`` via::

    python -m pytest python/tests/test_physics_regressions.py -v
"""

import math
import random
import struct

import pytest

try:
    import oxiphysics as ox
    HAS_OXIPHYSICS = True
except ImportError:
    HAS_OXIPHYSICS = False

pytestmark = pytest.mark.skipif(not HAS_OXIPHYSICS, reason="oxiphysics not built")


# ---------------------------------------------------------------------------
# LBM: Guo body forcing
# ---------------------------------------------------------------------------

CHANNEL_LAYERS = 21
CENTRE_SPEED = 0.02


def _channel_acceleration(viscosity):
    """Body acceleration that gives the centre speed ``a H^2 / (8 nu) = CENTRE_SPEED``."""
    return 8.0 * viscosity * CENTRE_SPEED / CHANNEL_LAYERS**2


def _poiseuille(accel, viscosity, s):
    """Plane Poiseuille speed at distance ``s`` from a wall of the channel."""
    return accel * s * (CHANNEL_LAYERS - s) / (2.0 * viscosity)


def _steady_steps(viscosity):
    """Ten decay times of the slowest diffusive mode of the channel."""
    return math.ceil(10.0 * CHANNEL_LAYERS**2 / (math.pi**2 * viscosity))


def test_lbm_body_force_adds_its_acceleration_every_step():
    accel = [3.0e-5, -1.2e-5]
    for viscosity in (0.1, 1.0 / 6.0, 0.02):
        cfg = ox.lbm.PyLbmConfig(6, 5, viscosity)
        cfg.body_force = accel
        sim = ox.lbm.PyLbmSimulation(cfg)
        previous = sim.velocity_at(2, 3)
        for _ in range(10):
            sim.step()
            current = sim.velocity_at(2, 3)
            for k in range(2):
                assert current[k] - previous[k] == pytest.approx(accel[k], abs=1e-12)
            previous = current


def test_lbm_2d_channel_reproduces_the_poiseuille_profile():
    viscosity = 1.0 / 6.0
    accel = _channel_acceleration(viscosity)
    cfg = ox.lbm.PyLbmConfig(3, CHANNEL_LAYERS + 2, viscosity)
    cfg.body_force = [accel, 0.0]
    sim = ox.lbm.PyLbmSimulation(cfg)
    sim.add_bottom_wall()
    sim.add_top_wall()
    sim.run(_steady_steps(viscosity))
    # Full-way bounce-back puts the walls half a cell outside the fluid layers.
    for y in range(1, CHANNEL_LAYERS + 1):
        expected = _poiseuille(accel, viscosity, y - 0.5)
        assert sim.velocity_at(1, y)[0] == pytest.approx(expected, abs=0.02 * CENTRE_SPEED)
    centre = sim.velocity_at(1, CHANNEL_LAYERS // 2 + 1)[0]
    assert centre == pytest.approx(CENTRE_SPEED, rel=0.02)
    assert sim.velocity_at(1, 0) == [0.0, 0.0]


def test_lbm_3d_channel_reproduces_the_poiseuille_profile():
    viscosity = 1.0 / 6.0
    accel = _channel_acceleration(viscosity)
    # One cell in x and y: the exact periodic reduction of a laterally uniform channel.
    cfg = ox.lbm.PyLbm3dConfig(1, 1, CHANNEL_LAYERS + 2, viscosity)
    cfg.body_force = [accel, 0.0, 0.0]
    sim = ox.lbm.PyLbm3dSimulation(cfg)
    wall = ox.lbm.LbmBoundary.NoSlipWall
    sim.set_boundary(0, 0, 0, wall)
    sim.set_boundary(0, 0, CHANNEL_LAYERS + 1, wall)
    sim.run(_steady_steps(viscosity))
    for z in range(1, CHANNEL_LAYERS + 1):
        expected = _poiseuille(accel, viscosity, z - 0.5)
        velocity = sim.velocity_at(0, 0, z)
        assert velocity[0] == pytest.approx(expected, abs=0.02 * CENTRE_SPEED)
        assert abs(velocity[1]) < 1e-9 and abs(velocity[2]) < 1e-9


# ---------------------------------------------------------------------------
# VTK writer
# ---------------------------------------------------------------------------

def _vtk_writer(count):
    writer = ox.io.PyVtkWriter("points.vtk")
    for k in range(count):
        writer.add_position(float(k), 0.5 * k, -1.0)
    return writer


def test_vtk_binary_vertex_cell_k_references_point_k():
    count = 5
    data = bytes(_vtk_writer(count).to_vtk_binary())
    marker = b"CELLS 5 10\n"
    start = data.index(marker) + len(marker)
    cells = struct.unpack(">10i", data[start:start + 40])
    assert cells == tuple(v for k in range(count) for v in (1, k))
    assert data[start + 40:start + 41] == b"\n"


def test_vtk_writer_rejects_arrays_with_the_wrong_length(tmp_path):
    writer = _vtk_writer(3)
    writer.add_point_data("pressure", [1.0, 2.0])
    target = tmp_path / "bad.vtk"
    with pytest.raises(ValueError):
        writer.to_vtk_binary()
    with pytest.raises(ValueError):
        writer.write_binary()
    with pytest.raises(ValueError):
        writer.write_to_file(str(target))
    assert not target.exists()


def test_vtk_writer_writes_the_bytes_it_reports(tmp_path):
    writer = _vtk_writer(4)
    writer.add_point_data("pressure", [1.0, 2.0, 3.0, 4.0])
    writer.add_cell_data("id", [0.0, 1.0, 2.0, 3.0])
    target = tmp_path / "ok.vtk"
    writer.write_to_file(str(target))
    assert target.read_bytes() == bytes(writer.to_vtk_binary())
    assert writer.write_binary() == target.stat().st_size


# ---------------------------------------------------------------------------
# Convex hull and GJK
# ---------------------------------------------------------------------------

UNIT_TETRA = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]


def _shifted(points, offset):
    return [[p[k] + offset[k] for k in range(3)] for p in points]


def test_convex_hull_of_a_tetrahedron_is_not_its_bounding_box():
    hull = ox.geometry.PyConvexHull.from_points(UNIT_TETRA)
    assert hull.volume() == pytest.approx(1.0 / 6.0, abs=1e-12)
    assert hull.surface_area() == pytest.approx(1.5 + 0.5 * math.sqrt(3.0), abs=1e-12)
    assert len(hull.get_vertices()) == 4
    assert len(hull.get_faces()) == 12


def test_convex_hull_rejects_degenerate_input():
    with pytest.raises(ValueError):
        ox.geometry.PyConvexHull.from_points(UNIT_TETRA[:3])
    coplanar = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
    with pytest.raises(ValueError):
        ox.geometry.PyConvexHull.from_points(coplanar)
    with pytest.raises(ValueError):
        ox.geometry.PyConvexHull.from_points([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])


def test_convex_hull_gjk_overlap_is_not_a_bounding_sphere_test():
    first = ox.geometry.PyConvexHull.from_points(UNIT_TETRA)
    # Disjoint from the first tetrahedron (x + y + z >= 1.5 against <= 1), yet the
    # bounding spheres of the two overlap.
    far_side = [[0.5, 0.5, 0.5], [1.5, 0.5, 0.5], [0.5, 1.5, 0.5], [0.5, 0.5, 1.5]]
    disjoint = ox.geometry.PyConvexHull.from_points(far_side)
    assert not first.gjk_overlap(disjoint)
    assert not disjoint.gjk_overlap(first)
    # An interpenetrating pair does overlap.
    overlapping = ox.geometry.PyConvexHull.from_points(_shifted(UNIT_TETRA, [0.25] * 3))
    assert first.gjk_overlap(overlapping)
    assert overlapping.gjk_overlap(first)


# ---------------------------------------------------------------------------
# XPBD volume constraint
# ---------------------------------------------------------------------------

def _tetra_volume(positions):
    p = [positions[3 * i:3 * i + 3] for i in range(4)]
    e1, e2, e3 = ([p[i][k] - p[0][k] for k in range(3)] for i in (1, 2, 3))
    cross = [e2[1] * e3[2] - e2[2] * e3[1],
             e2[2] * e3[0] - e2[0] * e3[2],
             e2[0] * e3[1] - e2[1] * e3[0]]
    return sum(e1[k] * cross[k] for k in range(3)) / 6.0


def _falling_tetrahedron(constrained):
    positions = [coordinate for point in UNIT_TETRA for coordinate in point]
    solver = ox.constraints.PyPbdSolver(positions, [0.0, 1.0, 1.0, 1.0])
    if constrained:
        solver.add_volume(0, 1, 2, 3, 0.0)
    for _ in range(20):
        solver.step(1.0 / 60.0)
    return _tetra_volume(solver.positions)


def test_pbd_volume_constraint_holds_the_volume_against_gravity():
    # The apex is pinned and the three base particles fall: without a constraint
    # the volume (1 - d) / 6 shrinks with the fall distance d.
    rest = 1.0 / 6.0
    assert abs(_falling_tetrahedron(False) - rest) / rest > 0.3
    assert abs(_falling_tetrahedron(True) - rest) / rest < 0.05


def test_pbd_add_volume_rejects_bad_particles():
    positions = [coordinate for point in UNIT_TETRA for coordinate in point]
    solver = ox.constraints.PyPbdSolver(positions, [1.0] * 4)
    with pytest.raises(IndexError):
        solver.add_volume(0, 1, 2, 4, 0.0)
    with pytest.raises(ValueError):
        solver.add_volume(0, 1, 1, 3, 0.0)
    solver.add_volume(0, 1, 2, 3, 0.0)


def test_pbd_add_stretch_rejects_bad_particles():
    positions = [coordinate for point in UNIT_TETRA for coordinate in point]
    solver = ox.constraints.PyPbdSolver(positions, [1.0] * 4)
    with pytest.raises(IndexError):
        solver.add_stretch(0, 4, 0.0)
    with pytest.raises(ValueError):
        solver.add_stretch(2, 2, 0.0)
    solver.add_stretch(0, 1, 0.0)


# ---------------------------------------------------------------------------
# Molecular dynamics
# ---------------------------------------------------------------------------

def _jittered_lattice(per_edge, box, jitter, seed):
    rng = random.Random(seed)
    spacing = box / per_edge
    return [
        [(index + 0.5) * spacing + jitter * (2.0 * rng.random() - 1.0)
         for index in (ix, iy, iz)]
        for ix in range(per_edge) for iy in range(per_edge) for iz in range(per_edge)
    ]


def _lj_energy(positions, box, epsilon=1.0, sigma=1.0, cutoff=2.5):
    """Brute-force truncated-shifted Lennard-Jones energy with the minimum image."""
    rc2 = cutoff * cutoff
    s6_cut = (sigma / cutoff) ** 6
    shift = 4.0 * epsilon * s6_cut * (s6_cut - 1.0)
    total = 0.0
    for i, a in enumerate(positions):
        for b in positions[i + 1:]:
            r2 = 0.0
            for k in range(3):
                d = b[k] - a[k]
                d -= box * round(d / box)
                r2 += d * d
            if r2 >= rc2 or r2 < 1e-20:
                continue
            s6 = (sigma * sigma / r2) ** 3
            total += 4.0 * epsilon * s6 * (s6 - 1.0) - shift
    return total


def test_md_neighbour_list_energy_matches_brute_force():
    cfg = ox.md.PyMdConfig()
    positions = _jittered_lattice(7, cfg.box_size, 0.1, 7)
    sim = ox.md.PyMdSimulation(cfg)
    sim.set_thermostat(False)
    for position in positions:
        sim.add_atom(position, 0)
    sim.step(0.0)  # evaluates the forces at the unmoved positions
    expected = _lj_energy(positions, cfg.box_size)
    assert abs(expected) > 1e-3
    assert sim.potential_energy() == pytest.approx(expected, rel=1e-9, abs=1e-9)


def test_md_microcanonical_dynamics_conserve_energy():
    cfg = ox.md.PyMdConfig()
    positions = _jittered_lattice(5, cfg.box_size, 0.1, 21)
    sim = ox.md.PyMdSimulation(cfg)
    sim.set_thermostat(False)
    rng = random.Random(9)
    for i, position in enumerate(positions):
        sim.add_atom(position, 0)
        sim.set_velocity(i, [0.8 * (2.0 * rng.random() - 1.0) for _ in range(3)])
    sim.step(0.0)
    start = sim.total_energy()
    sim.run(0.002, 400)
    assert abs(sim.total_energy() - start) / max(abs(start), 1.0) < 1e-3


def _lj_gas(initial_temperature, target, tau, seed=2024):
    """64 LJ atoms on a jittered 4x4x4 lattice with zero-drift Maxwell-like velocities."""
    cfg = ox.md.PyMdConfig()
    cfg.box_size = 6.0
    cfg.target_temperature = target
    sim = ox.md.PyNvtSimulation(cfg, tau)
    rng = random.Random(seed)
    for position in _jittered_lattice(4, 6.0, 0.05, seed):
        sim.add_atom(position, 0)
    count = sim.atom_count()
    velocities = [[rng.gauss(0.0, 1.0) for _ in range(3)] for _ in range(count)]
    drift = [sum(v[k] for v in velocities) / count for k in range(3)]
    velocities = [[v[k] - drift[k] for k in range(3)] for v in velocities]
    kinetic = sum(0.5 * sum(c * c for c in v) for v in velocities)
    scale = math.sqrt(initial_temperature / (2.0 * kinetic / (3.0 * count)))
    for i, v in enumerate(velocities):
        sim.set_velocity(i, [c * scale for c in v])
    return sim


def test_nose_hoover_samples_the_target_temperature_and_conserves_h():
    target = 1.2
    sim = _lj_gas(2.0, target, 0.1)
    h0 = sim.conserved_energy()
    max_drift = 0.0
    temperatures = []
    for step in range(3000):
        sim.step(0.002)
        max_drift = max(max_drift, abs((sim.conserved_energy() - h0) / h0))
        if step >= 1000:
            temperatures.append(sim.temperature())
    mean = sum(temperatures) / len(temperatures)
    spread = math.sqrt(sum((t - mean) ** 2 for t in temperatures) / len(temperatures)) / mean
    assert abs(mean - target) < 0.05 * target
    # A canonical thermostat lets T fluctuate; a hard velocity rescale would pin it.
    assert spread > 0.02
    assert max_drift < 1e-3
