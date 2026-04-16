# Terrain-Following Low-Observable UAV — MATLAB Simulation

A modular MATLAB simulation framework for **terrain-following, radar-aware UAV flight**.  
The system plans stealthy trajectories through complex terrain using sampling-based and grid-based planners, tracks them with an INDI controller built on differential flatness, and enforces safety constraints via a Control Barrier Function filter.

## Features

| Area | Capability |
|------|-----------|
| **Control** | 6-DOF quadrotor dynamics (NED), INDI outer/inner loop, quaternion PD attitude control |
| **Trajectory** | Minimum-snap polynomial generation, sliding-window online replanning, aerobatic maneuver library (loops, flips, barrel rolls, Immelmann, Split-S, figure-eight) |
| **Terrain** | Synthetic height-map generator (ridge, valley, hills), real GeoTIFF DEM loader (`dem_loader`), fast bilinear interpolation, terrain mesh with spatial indexing |
| **Radar** | Simplified radar equation, SNR-based detection probability, 3-D threat map with precomputed risk grid, terrain masking via LOS sampling |
| **Ray Casting** | Grid-based stepping with binary-search refinement, mesh-based Möller–Trumbore ray-triangle intersection, shadow maps, coverage maps, horizon profiles |
| **Planning — RRT\*** | Radar-aware RRT\* with KD-tree queries, dynamic rewiring radius, informed ellipsoidal sampling, shadow-zone bias, shortcut pruning |
| **Planning — Hybrid FMM+BIT\*** | Fast Marching Method heuristic on a 3-D cost grid → corridor-guided BIT\* with FMM-informed vertex/edge pruning and warm-start path |
| **Planning — Stealth Corridor** | FSM-based Eikonal solve on the visibility grid, gradient backtracking to extract a safe corridor skeleton, distance-transform-based radii |
| **Safety** | CBF quadratic-program filter enforcing terrain clearance, radar exposure limit, and velocity bound |
| **Visualization** | SkyMap app, 3-D flight animation, threat/shadow map overlays, RRT\* live tree growth |

## Quick Start

```matlab
% Bootstrap paths (run once per MATLAB session)
setup_project_paths;

% Core demos
run('demos/demo_indi_6dof_ned.m')        % INDI 6-DOF tracking
run('demos/demo_mission_1.m')            % Waypoint mission
run('demos/demo_mission_2.m')            % Multi-segment mission
run('demos/demo_mission_3.m')            % Complex mission
run('demos/demo_aerobatics.m')           % Aerobatic maneuvers

% Terrain & radar demos
run('demos/demo_terrain_masking.m')      % Terrain masking visualization
run('demos/demo_terrain_radar.m')        % Combined terrain + radar
run('demos/demo_dem_radar.m')            % Real DEM + radar

% Planner demos
run('demos/demo_radar_binary_masking.m') % Binary radar masking planning
run('demos/demo_radar_fmm_bit.m')       % Hybrid FMM + BIT* planner
run('demos/demo_fmm_bit_star.m')        % FMM-BIT* standalone
run('demos/demo_fsm_bit_star_direct.m') % FSM-BIT* direct planner

% Utilities
run('demos/ray_caster_demo.m')           % Ray casting demonstration
run('demos/final_mission_demo.m')        % Integrated end-to-end mission
```

## Repository Layout

```
test_rig/
├── demos/                  Demo and mission scripts (14 scripts)
├── tests/                  Verification scripts
├── flight/                 Flight dynamics, INDI control, flatness, aerobatics
│   ├── flatness.m          Differential flatness mapping
│   ├── quad_params_indi.m  Vehicle parameters
│   ├── aerobatic_maneuvers.m  Maneuver library
│   └── yaw_planner.m      Yaw reference generator
├── math/                   Quaternion and rotation utilities
├── motion_planner/         Path planning and trajectory smoothing
│   ├── rrt_star_radar.m    Radar-aware RRT* (KD-tree, informed sampling)
│   ├── bit_star_planner.m  Corridor-guided BIT*
│   ├── fmm_bit_star_planner.m  Hybrid FMM + BIT*
│   ├── fsm_bit_star_planner_direct.m  FSM-based direct planner
│   ├── compute_stealth_corridor.m  Stealth corridor extraction
│   ├── trajectory_smoother.m  Minimum-snap smoothing
│   ├── generate_trajectory_sliding_window.m  Online replanning
│   └── optimize_time_allocation.m  Segment time optimization
├── radar/                  Radar modeling and threat assessment
│   ├── radar_site.m        Radar site model (radar equation, SNR, Pd)
│   ├── threat_map.m        3-D precomputed threat grid
│   ├── compute_fmm_heuristic.m  FMM heuristic on cost grid
│   ├── compute_fsm_heuristic.m  FSM heuristic variant
│   └── skymap_app.m / skymap_viewer.m  SkyMap visualization
├── terrain/                Terrain modeling and ray casting
│   ├── terrain_generator.m Synthetic terrain (ridge, valley, hills)
│   ├── terrain_map.m       Height map with fast interpolation
│   ├── terrain_mesh.m      Triangle mesh with spatial indexing
│   ├── dem_loader.m        GeoTIFF DEM loader
│   ├── radar_ray_caster.m  Ray casting (grid + Möller–Trumbore)
│   ├── los_checker.m       Line-of-sight checker
│   └── moller_trumbore.m   Ray-triangle intersection kernel
├── safety/
│   └── cbf_safety_filter.m CBF-based safety QP
├── DEM/                    Terrain data assets (GeoTIFF)
├── project/
│   └── setup_project_paths.m  MATLAB path bootstrap
└── docs/
    └── repo_audit.md       Repository cleanup notes
```

## Dependencies

| Toolbox | Required for |
|---------|-------------|
| Core MATLAB | Everything |
| Image Processing Toolbox | `bwdist` in corridor generation and binary masking demos |
| Statistics and Machine Learning Toolbox | `createns`, `knnsearch`, `rangesearch` in RRT\* |
| Parallel Computing Toolbox *(optional)* | `parfor` acceleration in threat map and RRT\* rewiring |

## DEM Assets

Real terrain data under `DEM/`:

| File | Region |
|------|--------|
| `agri.tif` | Ağrı (mountainous) |
| `artvin.tif` | Artvin (steep valleys) |
| `bosphorus.tif` | Bosphorus / İstanbul |

## References

- **Tal, E. and Karaman, S.** — *Accurate Tracking of Aggressive Quadrotor Trajectories Using Incremental Nonlinear Dynamic Inversion and Differential Flatness* (IEEE CSL, 2021)
- **Tal, E. and Karaman, S.** — *Global Incremental Flight Control for Agile Maneuvering of a Tailsitter Flying Wing* (AIAA SciTech, 2022)
- **Mellinger, D. and Kumar, V.** — *Minimum Snap Trajectory Generation and Control for Quadrotors* (IEEE ICRA, 2011)
- **Karaman, S. and Frazzoli, E.** — *Sampling-Based Algorithms for Optimal Motion Planning* (IJRR, 2011)
- **LaValle, S. M.** — *Rapidly-Exploring Random Trees: A New Tool for Path Planning* (1998)
- **Pelosi, M., Kopp, C. and Brown, M.** — *Range-Limited UAV Trajectory Using Terrain Masking Under Radar Detection Risk* (Applied Artificial Intelligence, 2012)

## Author

4th year Aerospace Engineering student @ İstanbul Technical University
