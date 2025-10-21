%% Pointing induced spot size vs range
% INPUT: set exactly one of the two below
point_err_deg  = 0.035;   % total pointing error in degrees, leave [] if using mrad
point_err_mrad = [];     % total pointing error in milliradians, leave [] if using deg

% Convert to radians
if ~isempty(point_err_deg) && isempty(point_err_mrad)
    pe_rad = deg2rad(point_err_deg);
elseif isempty(point_err_deg) && ~isempty(point_err_mrad)
    pe_rad = point_err_mrad * 1e-3;
else
    error('Set exactly one of point_err_deg or point_err_mrad');
end

% Pretty print
fprintf('Total pointing error = %.4f deg  (%.0f urad)\n', pe_rad*180/pi, pe_rad*1e6);

% Range vector in meters
R = linspace(1e3, 650e3, 1000);   % 1 km to 200 km

% Linear footprint at target due to pointing error
spot_size = pe_rad .* R;

% Plot
figure;
plot(R/1e3, spot_size, 'LineWidth', 1.5);
xlabel('Range to target [km]');
ylabel('Spot size at target [m]');
title(sprintf('Sensor Pointing Footprint (Spot Size) vs Range [Error = %.3f°]', pe_rad*180/pi));

grid on;

% Reference markers and labels
hold on;
ref_ranges = [10e3, 100e3, 250e3 430e3];  % 10, 30, 50, 100 km
for rr = ref_ranges
    sz = pe_rad * rr;
    % hollow red circles
    plot(rr/1e3, sz, 'ro', 'MarkerFaceColor', 'none');
    % text placement like your original version
    text(rr/1e3, sz, sprintf('  %.2f m @ %.0f km', sz, rr/1e3));
end
hold off;
