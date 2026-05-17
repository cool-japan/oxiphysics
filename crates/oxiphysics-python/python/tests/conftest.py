"""Pytest fixtures shared across all oxiphysics test modules."""
import pytest

try:
    import oxiphysics as ox
    HAS_OXIPHYSICS = True
except ImportError:
    HAS_OXIPHYSICS = False

# Skip entire test modules when the native extension is not installed
skip_no_oxiphysics = pytest.mark.skipif(
    not HAS_OXIPHYSICS,
    reason="oxiphysics native extension not built (run: maturin develop --release)",
)


@pytest.fixture(scope="session")
def physics_world():
    """A default PhysicsWorld with Earth gravity for reproducibility."""
    if not HAS_OXIPHYSICS:
        pytest.skip("oxiphysics not built")
    world = ox.PhysicsWorld.with_earth_gravity()
    return world
