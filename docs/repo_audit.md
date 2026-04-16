# Repository Audit

## Rearrangement

The repository has been grouped around usage:

- `demos/`: interactive scripts and end-to-end showcases
- `tests/`: verification and diagnostics scripts
- `flight/`: reusable flight dynamics, control, and trajectory helpers
- `math/`: quaternion and rotation utilities
- Existing domain folders remain in place: `motion_planner/`, `radar/`, `terrain/`, `safety/`, `DEM/`

`project/setup_project_paths.m` is the shared bootstrap that adds project folders to the MATLAB path and resets the working directory to the repository root so existing relative DEM paths continue to work after moves.

## Unnecessary Or Likely Unnecessary Code

Strong candidates:

- `placeholder.m`
  - No inbound references in the active repository.
  - Name and size suggest scaffold code rather than production logic.
- `math/R_to_quat.m`
  - No inbound references in the active repository.
  - Safe candidate for removal if you do not plan to add reverse quaternion conversions.

Entrypoint-only files with zero inbound references:

- `demos/demo_dem_radar.m`
- `demos/demo_fmm_bit_star.m`
- `demos/demo_fsm_bit_star_direct.m`
- `demos/demo_mission_1.m`
- `demos/demo_mission_2.m`
- `demos/demo_mission_3.m`
- `demos/demo_radar_fmm_bit.m`
- `demos/demo_terrain_radar.m`
- `demos/final_mission_demo.m`
- `demos/ray_caster_demo.m`
- `tests/test_mission2_traj.m`
- `tests/test_trajectory_smoother.m`
- `tests/test_trajectory_upgrades.m`

These are not dead code by themselves. They are scripts or test entrypoints, so zero inbound references is expected.

## Non-Code Cleanup Candidates

- `final_mission_demo.asv`
  - MATLAB autosave artifact.
- `nul`
  - Stray file with no code role.
- `backup_original/`
  - Archive content rather than active project code.

These were left in place for safety, but they are good cleanup or archive targets if you want a stricter repository surface.
