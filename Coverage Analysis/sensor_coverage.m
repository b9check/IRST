%% Earth with balloon at 20 km, boresight to horizon, shaded FOV inverse cone
clear; close all; clc;

%% Parameters
Re      = 6371e3;           % Earth radius [m]
h       = 20e3;             % Balloon altitude [m]
FOVdeg  = 32;               % Full field of view [deg]
halfAng = deg2rad(FOVdeg/2);
coneLen = 150e3;            % Visualized cone length [m]

%% Sensor placement at equator for clarity
r_hat = [1 0 0];                  % Local radial up
P0    = (Re + h) * r_hat;         % Sensor position
t_hat = [0 1 0];                  % Local horizon boresight direction

% Orthonormal frame that spans the cone cross section
u_hat = cross(t_hat, r_hat);  u_hat = u_hat / norm(u_hat);
v_hat = cross(t_hat, u_hat);  v_hat = v_hat / norm(v_hat);

%% Build shaded cone surface
nAz   = 240;                      % azimuth samples around the cone
nAx   = 80;                       % samples along the cone axis
phi   = linspace(0, 2*pi, nAz);
s     = linspace(0, coneLen, nAx);           % distance along boresight
[Phi, S] = meshgrid(phi, s);

R = S * tan(halfAng);                          % cone radius at distance S

% Proper broadcasting: expand u_hat and v_hat to 3rd dimension
DirCross = bsxfun(@times, u_hat, cos(Phi)) + bsxfun(@times, v_hat, sin(Phi));

XS = P0(1) + S .* t_hat(1) + R .* DirCross(:,:,1);
YS = P0(2) + S .* t_hat(2) + R .* DirCross(:,:,2);
ZS = P0(3) + S .* t_hat(3) + R .* DirCross(:,:,3);

%% Draw Earth
figure('Color','w'); hold on; axis equal;
[xe, ye, ze] = sphere(240);
surf(Re*xe, Re*ye, Re*ze, 'FaceColor', [0.7 0.85 1.0], 'EdgeColor', 'none');
camlight headlight; lighting gouraud;

%% Draw sensor and boresight
plot3(P0(1), P0(2), P0(3), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
text(P0(1), P0(2), P0(3) + 0.02*Re, 'Balloon', 'HorizontalAlignment','center', 'FontWeight','bold');
qLen = 60e3;
quiver3(P0(1), P0(2), P0(3), qLen*t_hat(1), qLen*t_hat(2), qLen*t_hat(3), ...
    0, 'LineWidth', 2, 'Color', [0.1 0.5 0.1]);
text(P0(1)+qLen*t_hat(1), P0(2)+qLen*t_hat(2), P0(3)+qLen*t_hat(3), ...
    'Boresight', 'Color', [0.1 0.5 0.1]);

%% Shaded inverse cone
surf(XS, YS, ZS, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'FaceColor', [0 0 1]);

%% Labels and view
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title(sprintf('Balloon at 20 km, 32° FOV inverse cone pointed at horizon'));
grid on; view(38, 22);
