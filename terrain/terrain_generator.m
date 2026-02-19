function terrain = terrain_generator(params)
%TERRAIN_GENERATOR Generate synthetic terrain for simulation
%
% Generates a 2D height map representing terrain elevation.
% Useful for testing terrain-following and radar masking algorithms
% without requiring real DEM data.
%
% Inputs:
%   type   - Terrain type: 'flat', 'hills', 'valley', 'ridge', 'canyon', 'random'
%   params - Struct with parameters:
%            .bounds   - [N_min, N_max, E_min, E_max] in meters (default: [-500,500,-500,500])
%            .resolution - Grid resolution in meters (default: 10)
%            .base_height - Base terrain height in meters (default: 0)
%            .amplitude   - Height variation amplitude (default: 100)
%            .wavelength  - Feature wavelength in meters (default: 200)
%            .seed        - Random seed for reproducibility (default: 42)
%
% Output:
%   terrain - Struct with:
%             .N_grid     - [M x N] North coordinates
%             .E_grid     - [M x N] East coordinates
%             .Z          - [M x N] Terrain height (positive up, NED: -D)
%             .bounds     - [N_min, N_max, E_min, E_max]
%             .resolution - Grid resolution
%             .type       - Terrain type string
%
% Example:
%   terrain = terrain_generator(struct('type','hills', 'amplitude', 50));
%   surf(terrain.N_grid, terrain.E_grid, terrain.Z);
%
% Author: Quadrotor Terrain Following Project

% Default parameters
bounds = get_param(params, 'bounds', [-500, 500, -500, 500]);
resolution = get_param(params, 'resolution', 10);
base_height = get_param(params, 'base_height', 0);
amplitude = get_param(params, 'amplitude', 100);
wavelength = get_param(params, 'wavelength', 200);
seed = get_param(params, 'seed', 42);
type = get_param(params, 'type', 'hills');

% Create grid
N_vec = bounds(1):resolution:bounds(2);
E_vec = bounds(3):resolution:bounds(4);
[N_grid, E_grid] = meshgrid(N_vec, E_vec);

% Set random seed for reproducibility
rng(seed);

% Generate terrain based on type
switch lower(type)
    case 'flat'
        Z = base_height * ones(size(N_grid));

    case 'hills'
        % Rolling hills using superposition of sinusoids
        Z = base_height + ...
            amplitude * 0.5 * sin(2*pi*N_grid/wavelength) .* cos(2*pi*E_grid/wavelength) + ...
            amplitude * 0.3 * sin(2*pi*N_grid/(wavelength*0.7) + 1.2) + ...
            amplitude * 0.2 * cos(2*pi*E_grid/(wavelength*0.5) + 0.8);

    case 'valley'
        % V-shaped valley along North axis
        valley_center_E = (bounds(3) + bounds(4)) / 2;
        valley_width = (bounds(4) - bounds(3)) / 4;
        dist_from_center = abs(E_grid - valley_center_E);
        Z = base_height + amplitude * (dist_from_center / valley_width);
        % Add some noise
        Z = Z + amplitude * 0.1 * randn(size(Z));

    case 'ridge'
        % Ridge running North-South with rugged details
        ridge_center_E = (bounds(3) + bounds(4)) / 2;
        ridge_width = wavelength / 2.5; % Slightly narrower for consistent localized height

        % 1. Main visual mass
        dist_from_ridge = abs(E_grid - ridge_center_E);

        % Asymmetric ridge profile
        Z_main = base_height + amplitude * exp(-(dist_from_ridge.^2) / (2*ridge_width^2));

        % 2. Add structural noise (large features)
        % Variation along the ridge line
        ridge_modulation = 1 + 0.2 * sin(2*pi*N_grid/(wavelength*1.5));
        Z = Z_main .* ridge_modulation;

        % 3. Add rocky detail (high frequency noise)
        % Using multiple octaves for fractal look
        Z = Z + amplitude * 0.15 * sin(2*pi*N_grid/(wavelength*0.5)) .* cos(2*pi*E_grid/(wavelength*0.5));
        Z = Z + amplitude * 0.08 * sin(2*pi*N_grid/(wavelength*0.2)) .* cos(2*pi*E_grid/(wavelength*0.2));

        % 4. Add localized peaks
        peak_period = wavelength * 1.2;
        Z = Z + amplitude * 0.25 * exp(-(dist_from_ridge.^2)/(2*(ridge_width*0.5)^2)) .* ...
            max(0, sin(2*pi*N_grid/peak_period));

    case 'canyon'
        % Canyon with steep walls
        canyon_center_E = (bounds(3) + bounds(4)) / 2;
        canyon_width = wavelength / 3;
        dist_from_center = abs(E_grid - canyon_center_E);
        % Steep walls using tanh
        Z = base_height + amplitude * tanh((dist_from_center - canyon_width) / (canyon_width * 0.3));
        Z = max(Z, base_height);  % Floor at base height

    case 'mountain'
        % Prominent central peak with rugged alpine features
        peak_N = (bounds(1) + bounds(2)) / 2;
        peak_E = (bounds(3) + bounds(4)) / 2;

        % 1. Main Base Shape (Gaussian)
        dist_from_peak = sqrt((N_grid - peak_N).^2 + (E_grid - peak_E).^2);
        Z_base = base_height + amplitude * exp(-(dist_from_peak.^2) / (2*(wavelength*0.8)^2));

        % 2. Radiating Ridges (irregular)
        angle = atan2(E_grid - peak_E, N_grid - peak_N);
        % Use multiple frequencies for ridges to make them less perfect
        ridges = 0.5 * sin(5*angle) + 0.3 * sin(3*angle + 1) + 0.2 * cos(7*angle);
        % Ridges are prominent at mid-slopes, fade at peak and far base
        ridge_mask = exp(-dist_from_peak / (wavelength * 0.6));

        Z = Z_base .* (1 + 0.15 * ridges .* ridge_mask);

        % 3. Vertical Noise (Rock texture)
        noise = 0.05 * amplitude * sin(2*pi*N_grid/50) .* cos(2*pi*E_grid/50);
        Z = Z + noise .* exp(-dist_from_peak / wavelength);

        % 4. Secondary Peak or Shoulder
        shoulder_N = peak_N + wavelength * 0.4;
        shoulder_E = peak_E - wavelength * 0.3;
        dist_shoulder = sqrt((N_grid - shoulder_N).^2 + (E_grid - shoulder_E).^2);
        Z = Z + 0.4 * amplitude * exp(-(dist_shoulder.^2) / (2*(wavelength*0.4)^2));

    case 'random'
        % Fractal-like random terrain using multiple octaves
        Z = base_height * ones(size(N_grid));
        for octave = 1:4
            freq = 2^(octave-1);
            amp = amplitude / freq;
            phase_n = rand * 2 * pi;
            phase_e = rand * 2 * pi;
            Z = Z + amp * sin(2*pi*freq*N_grid/wavelength + phase_n) .* ...
                cos(2*pi*freq*E_grid/wavelength + phase_e);
        end
        % Add some Gaussian bumps
        n_bumps = 5;
        for i = 1:n_bumps
            bump_N = bounds(1) + rand * (bounds(2) - bounds(1));
            bump_E = bounds(3) + rand * (bounds(4) - bounds(3));
            bump_width = wavelength * (0.3 + 0.4 * rand);
            bump_height = amplitude * (0.3 + 0.7 * rand) * (2*rand - 1);
            dist = sqrt((N_grid - bump_N).^2 + (E_grid - bump_E).^2);
            Z = Z + bump_height * exp(-(dist.^2) / (2*bump_width^2));
        end

    case 'coastal'
        % Coastal terrain: flat near water, rising inland
        water_level = base_height;
        shore_E = bounds(3) + (bounds(4) - bounds(3)) * 0.3;
        inland_slope = amplitude / (bounds(4) - shore_E);
        Z = water_level + max(0, E_grid - shore_E) * inland_slope;
        % Add coastal hills
        Z = Z + amplitude * 0.3 * sin(2*pi*N_grid/wavelength) .* ...
            max(0, (E_grid - shore_E) / (bounds(4) - shore_E));

    otherwise
        error('Unknown terrain type: %s', type);
end

% Ensure terrain is positive (above sea level reference)
Z = max(Z, 0);

% Build output structure
terrain.N_grid = N_grid;
terrain.E_grid = E_grid;
terrain.Z = Z;
terrain.N_vec = N_vec;
terrain.E_vec = E_vec;
terrain.bounds = bounds;
terrain.resolution = resolution;
terrain.type = type;
terrain.base_height = base_height;
terrain.amplitude = amplitude;
end

function val = get_param(params, name, default)
if isfield(params, name)
    val = params.(name);
else
    val = default;
end
end
