%% Compute spot size (pixel footprint) as a function of target range
% Sensor / optics specs
nx = 1280;
ny = 1024;
pixel_pitch = 8e-6;  % meters
f = 0.018;  % focal length in meters

% Exact per-pixel IFOV (rad)
ifov = 2 * atan((pixel_pitch/2)/f);

% Exact full-sensor FOVs (rad)
fov_x = 2*atan((nx*pixel_pitch)/(2*f));
fov_y = 2*atan((ny*pixel_pitch)/(2*f));

% Print nicely
fprintf("IFOV = %.3f deg (%.0f µrad)\n", ifov*(180/(2*pi)), ifov*1e6);
fprintf("FOVx = %.2f deg\n", fov_x*180/pi);
fprintf("FOVy = %.2f deg\n", fov_y*180/pi);

% Define a range vector (target distances) over which to compute footprint
R = linspace(1e3, 200e3, 1000);  % from 1 km to 100 km, say (in meters)

% Compute linear footprint size = IFOV × R
spot_size = ifov * R;

% Plot
figure;
plot(R/1e3, spot_size);  % x axis in km
xlabel('Range to target (km)');
ylabel('Spot size (m) — linear footprint of one pixel');
title(sprintf('Pixel footprint vs Range (FOV_x = %.1f, FOVy = %.1f)', fov_x*180/pi, fov_y*180/pi));
grid on;

% Optionally, you can overlay lines or reference spot sizes of interest
hold on;
ref_ranges = [10e3, 30e3, 50e3, 100e3];  % 10 km, 30 km, 50 km
for rr = ref_ranges
    sz = ifov * rr;
    plot(rr/1e3, sz, 'ro');
    text(rr/1e3, sz, sprintf('  %.1f m @ %.0f km', sz, rr/1e3));
end
hold off;
