function terrain_data = dem_loader(tif_filepath, params)
%DEM_LOADER Load GeoTIFF DEM file and convert to terrain_map compatible format
%
% Reads a GeoTIFF (.tif) Digital Elevation Model file using MATLAB's
% readgeoraster function and converts it to the format expected by the
% terrain_map class for use in quadrotor simulation and radar analysis.
%
% Inputs:
%   tif_filepath - Path to .tif GeoTIFF file (string or char)
%   params (optional) - Struct with configuration options:
%     .origin          - [N, E] origin for local NED frame (default: center of DEM)
%     .target_resolution - Resample to this resolution [m] (default: native)
%     .crop_bounds     - [N_min, N_max, E_min, E_max] crop region in local coords
%     .fill_nodata     - Method to fill NaN/nodata: 'nearest', 'linear', 'none'
%     .vertical_datum  - Vertical offset to apply [m] (default: 0)
%
% Outputs:
%   terrain_data - Struct compatible with terrain_map constructor:
%     .N_vec      - North coordinate vector [1 x M]
%     .E_vec      - East coordinate vector [1 x N]
%     .Z          - Height matrix [N x M] (E rows, N cols)
%     .bounds     - [N_min, N_max, E_min, E_max]
%     .resolution - Grid spacing [m]
%     .type       - 'dem_geotiff'
%     .source     - Original filepath
%     .crs_info   - Coordinate reference system information
%
% Example:
%   % Load DEM and create terrain map
%   terrain_data = dem_loader('elevation.tif');
%   tm = terrain_map(terrain_data);
%   tm.plot();
%
%   % Load with custom origin and resolution
%   params.origin = [500000, 4500000];  % UTM coordinates
%   params.target_resolution = 30;      % 30m resolution
%   terrain_data = dem_loader('elevation.tif', params);
%
% Supported formats:
%   - GeoTIFF with UTM projection (meters)
%   - GeoTIFF with geographic coordinates (lat/lon) - will be converted
%
% Requirements:
%   - MATLAB Mapping Toolbox (for readgeoraster)
%
% Author: Quadrotor Terrain Following Project
% See also: terrain_map, terrain_generator, readgeoraster

%% Input validation
if nargin < 1
    error('dem_loader:NoInput', 'Filepath to GeoTIFF is required.');
end

if ~isfile(tif_filepath)
    error('dem_loader:FileNotFound', 'File not found: %s', tif_filepath);
end

% Default parameters
if nargin < 2
    params = struct();
end

origin = get_param(params, 'origin', []);
target_resolution = get_param(params, 'target_resolution', []);
crop_bounds = get_param(params, 'crop_bounds', []);
fill_nodata = get_param(params, 'fill_nodata', 'nearest');
vertical_datum = get_param(params, 'vertical_datum', 0);

%% Read GeoTIFF file
fprintf('Loading GeoTIFF: %s\n', tif_filepath);

try
    [Z_raw, R] = readgeoraster(tif_filepath);
catch ME
    error('dem_loader:ReadError', 'Failed to read GeoTIFF: %s\n%s', ...
        tif_filepath, ME.message);
end

% Convert to double and handle multiple bands
if size(Z_raw, 3) > 1
    warning('dem_loader:MultiBand', 'Multi-band GeoTIFF detected. Using first band.');
    Z_raw = Z_raw(:, :, 1);
end
Z_raw = double(Z_raw);

%% Extract coordinate reference information
crs_info = struct();
crs_info.source_file = tif_filepath;

if isa(R, 'map.rasterref.MapCellsReference') || isa(R, 'map.rasterref.MapPostingsReference')
    % Projected coordinate system (UTM, etc.) - coordinates in meters
    crs_info.type = 'projected';
    crs_info.x_limits = R.XWorldLimits;
    crs_info.y_limits = R.YWorldLimits;
    crs_info.cell_extent_x = R.CellExtentInWorldX;
    crs_info.cell_extent_y = R.CellExtentInWorldY;

    % Create coordinate vectors
    % Note: Y typically increases northward, X increases eastward
    % For NED frame: N = Y (northing), E = X (easting)
    x_vec = R.XWorldLimits(1) + R.CellExtentInWorldX/2 : R.CellExtentInWorldX : R.XWorldLimits(2);
    y_vec = R.YWorldLimits(1) + R.CellExtentInWorldY/2 : R.CellExtentInWorldY : R.YWorldLimits(2);

    native_resolution = mean([R.CellExtentInWorldX, R.CellExtentInWorldY]);

elseif isa(R, 'map.rasterref.GeographicCellsReference') || isa(R, 'map.rasterref.GeographicPostingsReference')
    % Geographic coordinate system (lat/lon) - need to convert to meters
    crs_info.type = 'geographic';
    crs_info.lat_limits = R.LatitudeLimits;
    crs_info.lon_limits = R.LongitudeLimits;

    warning('dem_loader:Geographic', ...
        'Geographic coordinates detected. Converting to approximate UTM meters.');

    % Create lat/lon vectors
    lat_vec = linspace(R.LatitudeLimits(1), R.LatitudeLimits(2), size(Z_raw, 1));
    lon_vec = linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), size(Z_raw, 2));

    % Convert to approximate meters using center latitude
    center_lat = mean(R.LatitudeLimits);
    meters_per_deg_lat = 111320;  % Approximate
    meters_per_deg_lon = 111320 * cosd(center_lat);

    y_vec = (lat_vec - center_lat) * meters_per_deg_lat;
    x_vec = (lon_vec - mean(R.LongitudeLimits)) * meters_per_deg_lon;

    native_resolution = mean([...
        abs(diff(R.LatitudeLimits)) / size(Z_raw, 1) * meters_per_deg_lat, ...
        abs(diff(R.LongitudeLimits)) / size(Z_raw, 2) * meters_per_deg_lon]);
else
    error('dem_loader:UnknownCRS', 'Unknown coordinate reference type: %s', class(R));
end

fprintf('  Native resolution: %.2f m\n', native_resolution);
fprintf('  Grid size: %d x %d\n', size(Z_raw, 2), size(Z_raw, 1));

%% Handle coordinate orientation
% GeoTIFF typically has origin at top-left with Y decreasing downward
% We need Y (North) increasing upward

% Check if Y needs to be flipped (common in GeoTIFFs)
if length(y_vec) > 1 && y_vec(1) > y_vec(end)
    y_vec = flip(y_vec);
    Z_raw = flipud(Z_raw);
end

% Ensure X increases
if length(x_vec) > 1 && x_vec(1) > x_vec(end)
    x_vec = flip(x_vec);
    Z_raw = fliplr(Z_raw);
end

%% Apply origin transformation
% Convert to local NED frame centered at origin
if isempty(origin)
    % Default: center of DEM
    origin = [mean(y_vec), mean(x_vec)];
end

N_vec = y_vec - origin(1);  % Y -> North
E_vec = x_vec - origin(2);  % X -> East

crs_info.origin = origin;
crs_info.origin_description = 'Local NED frame origin in source CRS';

%% Handle nodata/NaN values
nodata_mask = isnan(Z_raw) | Z_raw < -1000 | Z_raw > 10000;
n_nodata = sum(nodata_mask(:));

if n_nodata > 0
    fprintf('  Nodata pixels: %d (%.1f%%)\n', n_nodata, 100*n_nodata/numel(Z_raw));

    switch lower(fill_nodata)
        case 'nearest'
            % Fill using nearest neighbor interpolation
            Z_raw = fillmissing(Z_raw, 'nearest');
            % Handle remaining edge NaNs
            Z_raw(isnan(Z_raw)) = nanmean(Z_raw(:));
        case 'linear'
            % Fill using linear interpolation
            Z_raw = fillmissing(Z_raw, 'linear');
            Z_raw(isnan(Z_raw)) = nanmean(Z_raw(:));
        case 'none'
            % Leave NaN values
        otherwise
            warning('dem_loader:UnknownFill', 'Unknown fill method: %s. Using nearest.', fill_nodata);
            Z_raw = fillmissing(Z_raw, 'nearest');
    end
end

%% Apply vertical datum offset
Z_raw = Z_raw + vertical_datum;

%% Resample if requested
if ~isempty(target_resolution) && target_resolution ~= native_resolution
    fprintf('  Resampling from %.2f m to %.2f m resolution\n', native_resolution, target_resolution);

    % Create new coordinate vectors
    N_vec_new = N_vec(1) : target_resolution : N_vec(end);
    E_vec_new = E_vec(1) : target_resolution : E_vec(end);

    % Interpolate to new grid
    [E_grid_old, N_grid_old] = meshgrid(E_vec, N_vec);
    [E_grid_new, N_grid_new] = meshgrid(E_vec_new, N_vec_new);

    F = griddedInterpolant(E_grid_old', N_grid_old', Z_raw', 'linear', 'nearest');
    Z_new = F(E_grid_new', N_grid_new')';

    N_vec = N_vec_new;
    E_vec = E_vec_new;
    Z_raw = Z_new;
    native_resolution = target_resolution;
end

%% Crop to bounds if specified
if ~isempty(crop_bounds)
    N_min = crop_bounds(1);
    N_max = crop_bounds(2);
    E_min = crop_bounds(3);
    E_max = crop_bounds(4);

    N_idx = N_vec >= N_min & N_vec <= N_max;
    E_idx = E_vec >= E_min & E_vec <= E_max;

    N_vec = N_vec(N_idx);
    E_vec = E_vec(E_idx);
    Z_raw = Z_raw(E_idx, N_idx);

    fprintf('  Cropped to: N=[%.0f, %.0f], E=[%.0f, %.0f]\n', N_min, N_max, E_min, E_max);
end

%% Build output structure
% Note: terrain_map expects Z as [length(E_vec) x length(N_vec)]
% where Z(i,j) is height at (E_vec(i), N_vec(j))

% Transpose Z to match expected format: rows = E, cols = N
Z = Z_raw;  % Already in correct orientation from GeoTIFF read

terrain_data = struct();
terrain_data.N_vec = N_vec(:)';
terrain_data.E_vec = E_vec(:)';
terrain_data.Z = Z;
terrain_data.bounds = [min(N_vec), max(N_vec), min(E_vec), max(E_vec)];
terrain_data.resolution = native_resolution;
terrain_data.type = 'dem_geotiff';
terrain_data.source = tif_filepath;
terrain_data.crs_info = crs_info;

fprintf('  Output bounds: N=[%.0f, %.0f], E=[%.0f, %.0f] m\n', ...
    terrain_data.bounds(1), terrain_data.bounds(2), ...
    terrain_data.bounds(3), terrain_data.bounds(4));
fprintf('  Height range: [%.1f, %.1f] m\n', min(Z(:)), max(Z(:)));
fprintf('DEM loading complete.\n');

end

%% Helper function
function val = get_param(params, name, default)
    if isfield(params, name)
        val = params.(name);
    else
        val = default;
    end
end
