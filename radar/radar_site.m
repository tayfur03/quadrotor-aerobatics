classdef radar_site < handle
%RADAR_SITE Class representing a ground-based radar installation
%
% Models a radar site with position, characteristics, and detection capability.
% Uses simplified radar equation for detection probability calculation.
%
% Properties:
%   position     - [3x1] radar position [N; E; altitude] (altitude positive up)
%   name         - String identifier
%   type         - 'surveillance', 'tracking', 'fire_control'
%   P_t          - Transmit power [W]
%   G            - Antenna gain (linear, not dB)
%   freq         - Operating frequency [Hz]
%   lambda       - Wavelength [m]
%   R_max        - Maximum detection range [m]
%   P_fa         - Probability of false alarm
%   n_pulses     - Number of integrated pulses
%   enabled      - Whether radar is active
%
% Methods:
%   get_detection_probability(target_pos, rcs, los_checker)
%   get_received_power(target_pos, rcs)
%   get_snr(target_pos, rcs)
%
% Example:
%   radar = radar_site([0; 0; 100], 'SAM-1', 'surveillance');
%   P_det = radar.get_detection_probability([500; 200; 50], 0.1, los);
%
% Author: Quadrotor Terrain Following Project

    properties
        position        % [N; E; altitude] in meters
        name            % Radar identifier
        type            % Radar type
        P_t             % Transmit power [W]
        G               % Antenna gain (linear)
        freq            % Frequency [Hz]
        lambda          % Wavelength [m]
        R_max           % Maximum range [m]
        P_fa            % False alarm probability
        n_pulses        % Integrated pulses
        noise_figure    % Receiver noise figure [dB]
        bandwidth       % Receiver bandwidth [Hz]
        enabled         % Is radar active?
        azimuth_limits  % [az_min, az_max] in degrees (empty = 360 deg)
        elevation_limits % [el_min, el_max] in degrees
    end

    properties (Constant)
        k_B = 1.380649e-23;  % Boltzmann constant [J/K]
        T_0 = 290;           % Reference temperature [K]
        c = 299792458;       % Speed of light [m/s]
    end

    methods
        function obj = radar_site(position, name, type, params)
            %RADAR_SITE Constructor
            %   radar = radar_site(pos, name, type)
            %   radar = radar_site(pos, name, type, params)

            obj.position = position(:);
            obj.name = name;

            if nargin < 3
                type = 'surveillance';
            end
            obj.type = type;

            if nargin < 4
                params = struct();
            end

            % Set defaults based on radar type
            switch lower(type)
                case 'surveillance'
                    % Long-range, lower resolution
                    default_Pt = 1e6;        % 1 MW
                    default_G = 1000;        % ~30 dB
                    default_freq = 3e9;      % S-band (3 GHz)
                    default_Rmax = 200e3;    % 200 km
                    default_n_pulses = 20;

                case 'tracking'
                    % Medium range, high accuracy
                    default_Pt = 100e3;      % 100 kW
                    default_G = 5000;        % ~37 dB
                    default_freq = 10e9;     % X-band (10 GHz)
                    default_Rmax = 100e3;    % 100 km
                    default_n_pulses = 50;

                case 'fire_control'
                    % Short range, very accurate
                    default_Pt = 50e3;       % 50 kW
                    default_G = 10000;       % ~40 dB
                    default_freq = 15e9;     % Ku-band (15 GHz)
                    default_Rmax = 50e3;     % 50 km
                    default_n_pulses = 100;

                case 'manpads'
                    % Man-portable, short range
                    default_Pt = 1e3;        % 1 kW
                    default_G = 100;         % ~20 dB
                    default_freq = 35e9;     % Ka-band
                    default_Rmax = 10e3;     % 10 km
                    default_n_pulses = 10;

                otherwise
                    % Generic medium radar
                    default_Pt = 100e3;
                    default_G = 1000;
                    default_freq = 5e9;
                    default_Rmax = 100e3;
                    default_n_pulses = 30;
            end

            % Apply parameters
            obj.P_t = get_param(params, 'P_t', default_Pt);
            obj.G = get_param(params, 'G', default_G);
            obj.freq = get_param(params, 'freq', default_freq);
            obj.R_max = get_param(params, 'R_max', default_Rmax);
            obj.n_pulses = get_param(params, 'n_pulses', default_n_pulses);
            obj.P_fa = get_param(params, 'P_fa', 1e-6);
            obj.noise_figure = get_param(params, 'noise_figure', 5);  % 5 dB
            obj.bandwidth = get_param(params, 'bandwidth', 1e6);      % 1 MHz
            obj.enabled = get_param(params, 'enabled', true);
            obj.azimuth_limits = get_param(params, 'azimuth_limits', []);
            obj.elevation_limits = get_param(params, 'elevation_limits', [-5, 85]);

            % Compute wavelength
            obj.lambda = obj.c / obj.freq;
        end

        function P_r = get_received_power(obj, target_pos, rcs)
            %GET_RECEIVED_POWER Calculate received power from target
            %   P_r = get_received_power(target_pos, rcs)
            %
            %   Uses radar range equation:
            %   P_r = (P_t * G^2 * lambda^2 * sigma) / ((4*pi)^3 * R^4)

            target_pos = target_pos(:);
            R = norm(target_pos - obj.position);

            if R < 1  % Avoid singularity
                R = 1;
            end

            % Radar range equation
            P_r = (obj.P_t * obj.G^2 * obj.lambda^2 * rcs) / ...
                  ((4*pi)^3 * R^4);
        end

        function snr = get_snr(obj, target_pos, rcs)
            %GET_SNR Calculate signal-to-noise ratio
            %   snr = get_snr(target_pos, rcs)
            %
            %   Returns SNR in linear scale (not dB)

            P_r = obj.get_received_power(target_pos, rcs);

            % Noise power
            F = 10^(obj.noise_figure / 10);  % Convert NF from dB
            P_n = obj.k_B * obj.T_0 * obj.bandwidth * F;

            % SNR with pulse integration (coherent assumed)
            snr = (P_r / P_n) * obj.n_pulses;
        end

        function P_det = get_detection_probability(obj, target_pos, rcs, los_obj)
            %GET_DETECTION_PROBABILITY Calculate probability of detection
            %   P_det = get_detection_probability(target_pos, rcs)
            %   P_det = get_detection_probability(target_pos, rcs, los_checker)
            %
            %   If los_checker is provided, terrain masking is considered.
            %   Uses Swerling Case 1 detection model (fluctuating target).

            if ~obj.enabled
                P_det = 0;
                return;
            end

            target_pos = target_pos(:);

            % Check range
            R = norm(target_pos - obj.position);
            if R > obj.R_max
                P_det = 0;
                return;
            end

            % Check angular coverage
            if ~obj.is_in_coverage(target_pos)
                P_det = 0;
                return;
            end

            % Check line of sight (terrain masking)
            if nargin >= 4 && ~isempty(los_obj)
                % Convert positions to LOS checker format [N; E; alt]
                radar_los_pos = obj.position;
                target_los_pos = target_pos;

                if ~los_obj.has_los(radar_los_pos, target_los_pos)
                    P_det = 0;  % Terrain blocked
                    return;
                end
            end

            % Calculate SNR
            snr = obj.get_snr(target_pos, rcs);

            % Detection probability using Swerling 1 model
            % Approximation: P_d ≈ P_fa^(1/(1+snr)) for Swerling 1
            % More accurate: use Marcum Q function approximation

            if snr <= 0
                P_det = 0;
            else
                % Swerling 1 detection probability approximation
                % Based on Albersheim's equation for detection threshold
                snr_db = 10 * log10(snr);

                % Required SNR for P_d = 0.9, P_fa = 1e-6 is about 13.2 dB
                % Adjust for actual P_fa
                snr_threshold_db = 10.8 - 10*log10(obj.n_pulses) + ...
                                   5*log10(-log10(obj.P_fa));

                % Sigmoid approximation for P_d vs SNR
                snr_margin = snr_db - snr_threshold_db;
                P_det = 1 / (1 + exp(-0.8 * snr_margin));

                % Ensure physical bounds
                P_det = max(0, min(1, P_det));
            end
        end

        function in_coverage = is_in_coverage(obj, target_pos)
            %IS_IN_COVERAGE Check if target is within radar angular coverage

            target_pos = target_pos(:);
            rel_pos = target_pos - obj.position;

            % Calculate azimuth and elevation
            azimuth = atan2d(rel_pos(2), rel_pos(1));  % degrees from North
            range_horiz = norm(rel_pos(1:2));
            elevation = atan2d(rel_pos(3), range_horiz);

            in_coverage = true;

            % Check azimuth limits
            if ~isempty(obj.azimuth_limits)
                az_min = obj.azimuth_limits(1);
                az_max = obj.azimuth_limits(2);

                % Handle wraparound
                if az_min <= az_max
                    in_coverage = in_coverage && (azimuth >= az_min && azimuth <= az_max);
                else
                    in_coverage = in_coverage && (azimuth >= az_min || azimuth <= az_max);
                end
            end

            % Check elevation limits
            if ~isempty(obj.elevation_limits)
                in_coverage = in_coverage && ...
                    (elevation >= obj.elevation_limits(1) && ...
                     elevation <= obj.elevation_limits(2));
            end
        end

        function info = get_info(obj)
            %GET_INFO Return radar parameters as string
            info = sprintf('%s (%s): P_t=%.0f kW, G=%.0f dB, f=%.1f GHz, R_max=%.0f km', ...
                obj.name, obj.type, obj.P_t/1e3, 10*log10(obj.G), ...
                obj.freq/1e9, obj.R_max/1e3);
        end

        function plot_coverage(obj, ax, max_range)
            %PLOT_COVERAGE Visualize radar coverage area
            if nargin < 2 || isempty(ax)
                figure;
                ax = gca;
            end
            if nargin < 3
                max_range = obj.R_max;
            end

            % Draw range ring
            theta = linspace(0, 2*pi, 100);
            x = obj.position(1) + max_range * cos(theta);
            y = obj.position(2) + max_range * sin(theta);

            plot(ax, x, y, 'r--', 'LineWidth', 1.5);
            hold(ax, 'on');

            % Draw radar position
            plot(ax, obj.position(1), obj.position(2), 'r^', ...
                 'MarkerSize', 12, 'MarkerFaceColor', 'r');

            % Label
            text(ax, obj.position(1), obj.position(2) + max_range*0.1, ...
                 obj.name, 'FontWeight', 'bold', 'Color', 'r');
        end
    end
end

function val = get_param(params, name, default)
    if isfield(params, name)
        val = params.(name);
    else
        val = default;
    end
end
