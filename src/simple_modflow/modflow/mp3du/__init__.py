"""Public MP3DU helpers for simple_modflow.

The supported particle-tracking surface is intentionally small:

- ``ParticleTrackingInput``
- ``prepare_particle_tracking``
- ``run_particle_tracking``

Legacy FloPy helpers for the native MODFLOW ``PRT`` model live in
``simple_modflow.modflow.mp3du.legacy_prt`` and are not part of the preferred
API.
"""

from .particles import (
    ParticleTrackingInput,
    prepare_particle_tracking,
    run_particle_tracking,
)

__all__ = [
    "ParticleTrackingInput",
    "prepare_particle_tracking",
    "run_particle_tracking",
]
