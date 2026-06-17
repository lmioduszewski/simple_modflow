"""Public MP3DU-engine helpers for myflopy.

The supported particle-tracking surface is intentionally small:

- ``ParticleTrackingInput``
- ``prepare_particle_tracking``
- ``run_particle_tracking``

The preferred native MODFLOW 6 PRT workflow lives in
``myflopy.modflow.mf6.prt``. Legacy experimental PRT package wrappers
remain in ``myflopy.modflow.mp3du.legacy_prt`` for compatibility only.
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
