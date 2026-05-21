# Terrain-Following Low-Observable UAV MATLAB Simulation

A script-driven MATLAB workspace for terrain-following, radar-aware UAV planning and flight simulation. The current active flow focuses on terrain/radar modeling, radar-aware planners, minimum-snap trajectory smoothing, and INDI-based 6-DOF tracking.

## Main Capabilities

| Area | Current scope |
|------|---------------|
| Flight and control | 6-DOF quadrotor dynamics in NED, differential flatness, INDI control, quaternion attitude utilities |
| Trajectory generation | Minimum-snap trajectory smoothing with adaptive segment-time allocation |
| Terrain | Synthetic terrain, GeoTIFF DEM loading, terrain maps/meshes, LOS and ray-casting helpers |
| Radar | Radar-site model, threat maps, terrain masking, SkyMap visualization |
| Planning | Radar-aware RRT*, FMM/FSM-guided BIT*, AIT* direct planner experiments, corridor extraction |
| Safety | CBF safety-filter prototype |
| Demos | Mission, terrain/radar masking, planner benchmark, FMM/FSM/BIT* and AIT* scripts |

## Quick Start

Run from the repository root in MATLAB:

```matlab
setup_project_paths;

% Flight/control demos
run('demos/demo_indi_6dof_ned.m')
run('demos/demo_aerobatics.m')

% Mission demos
run('demos/demo_mission_1.m')
run('demos/demo_mission_2.m')
run('demos/demo_mission_3.m')

% Terrain and radar demos
run('demos/demo_terrain_masking.m')
run('demos/demo_terrain_radar.m')
run('demos/demo_dem_radar.m')
run('demos/ray_caster_demo.m')

% Planner demos and benchmarks
run('demos/demo_radar_binary_masking.m')
run('demos/demo_radar_fmm_bit.m')
run('demos/demo_fmm_bit_star.m')
run('demos/demo_fsm_bit_star_direct.m')
run('demos/demo_ait_star_direct.m')
run('demos/demo_plain_radar_bit_star_performance.m')
run('demos/demo_radar_planner_benchmark.m')

% Integrated mission
run('demos/final_mission_demo.m')
```

## Repository Layout

```text
test_rig/
|-- demos/                  Script entry points and visual demos
|-- tests/                  Lightweight validation scripts
|-- flight/                 Dynamics, INDI, flatness, yaw, aerobatics
|-- math/                   Quaternion and rotation helpers
|-- motion_planner/         RRT*, BIT*/FMM/FSM/AIT* planners and smoothing
|-- radar/                  Radar sites, threat maps, FMM/FSM heuristics, SkyMap
|-- terrain/                DEM loading, terrain maps/meshes, LOS, ray casting
|-- safety/                 CBF safety filter
|-- project/                MATLAB path bootstrap
|-- DEM/                    Local GeoTIFF terrain assets
`-- docs/                   Notes and audits
```

## Active Planner Flow

The actively maintained radar-planning flow is:

```text
terrain/DEM -> threat map -> planner path -> simplify_path
            -> trajectory_smoother -> flight/control demo or visualization
```

`trajectory_smoother` is the current trajectory interface. The older sliding-window trajectory generator has been removed from the active code surface.

## Dependencies

| Dependency | Used for |
|------------|----------|
| Core MATLAB | Scripts, plotting, interpolation, numeric routines |
| UAV Toolbox / Robotics-related trajectory support | `minsnappolytraj` trajectory generation |
| Image Processing Toolbox | `bwdist` in corridor and binary-mask flows |
| Statistics and Machine Learning Toolbox | KD-tree helpers such as `createns`, `knnsearch`, `rangesearch` |
| Parallel Computing Toolbox, optional | `parfor` acceleration where available |

## DEM Assets

`DEM/` contains local GeoTIFF terrain files used by demos. Some scripts fall back to synthetic terrain when the requested DEM is missing.

Known assets in this workspace include:

- `agri.tif`
- `artvin.tif`
- `bosphorus.tif`
- `denizli.tif`
- `konya.tif`
- `oltu.tif`
- `siirt.tif`

Large media files are ignored by `.gitignore`.

## Validation

Use the smallest relevant script first:

```matlab
run('tests/test_trajectory_smoother.m')
run('tests/test_mission2_traj.m')
run('tests/test_ait_star_direct_smoke.m')
run('demos/demo_fsm_bit_star_direct.m')
run('demos/demo_radar_planner_benchmark.m')
```

Full visual demos can be slow because they build terrain, threat grids, and planner trees.

## References

- Tal, E. and Karaman, S., "Accurate Tracking of Aggressive Quadrotor Trajectories Using Incremental Nonlinear Dynamic Inversion and Differential Flatness", IEEE CSL, 2021.
- Mellinger, D. and Kumar, V., "Minimum Snap Trajectory Generation and Control for Quadrotors", IEEE ICRA, 2011.
- Karaman, S. and Frazzoli, E., "Sampling-Based Algorithms for Optimal Motion Planning", IJRR, 2011.
- LaValle, S. M., "Rapidly-Exploring Random Trees: A New Tool for Path Planning", 1998.
- Pelosi, M., Kopp, C. and Brown, M., "Range-Limited UAV Trajectory Using Terrain Masking Under Radar Detection Risk", Applied Artificial Intelligence, 2012.

## Author

4th year Aerospace Engineering student, Istanbul Technical University.
