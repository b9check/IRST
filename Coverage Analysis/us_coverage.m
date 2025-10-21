%% Balloon Coverage Mesh over the Continental United States
% Requires Mapping Toolbox
close all; clear;

%% CONFIG
% Geographic domain for CONUS
latlim = [24 50];              % degrees north
lonlim = [-125 -66];           % degrees east

% Earth and altitudes
re   = 6371e3;                 % Earth radius [m]
hB   = 20e3;                   % balloon altitude [m]
hM   = 30e3;                   % target altitude [m]

% Coverage radius choice
useHorizon      = true;        % true uses horizon based coverage
manualRange_km  = 500;         % used only if useHorizon equals false

% Mesh controls
sep_km   = 200;                % desired neighbor spacing [km] for hex mesh
landOnly = true;               % keep nodes on land inside CONUS bbox

% Paint grid resolution for coverage counts
degRes = 0.25;                 % degrees

% Plot look
baseMapName = "satellite";

%% COVERAGE RADIUS
if useHorizon
    phiB = acos(re/(re+hB));
    phiM = acos(re/(re+hM));
    theta = phiB + phiM;
    r_km = theta*(re+hM)/1e3;  % projected ground radius at target altitude
else
    r_km = manualRange_km;
end

%% BUILD HEX MESH OVER CONUS BBOX
% Hex geometry in km
row_step_km = sep_km*sqrt(3)/2;

% Convert row step to degrees latitude
deg_per_km_lat = 1/111.32;
dlat = row_step_km * deg_per_km_lat;

% Generate row latitudes
lat_rows = (latlim(1)+dlat/2):dlat:(latlim(2)-dlat/2);

clat = []; clon = [];
for j = 1:numel(lat_rows)
    lat_j = lat_rows(j);

    % Longitude step depends on latitude
    deg_per_km_lon = 1./(111.32*cosd(lat_j)+eps);
    dlon = sep_km * deg_per_km_lon;

    if isinf(dlon) || dlon <= 0
        continue
    end

    % Odd rows are offset by half a step
    offset = 0;
    if mod(j,2) == 0
        offset = dlon/2;
    end

    lon_row = (lonlim(1)+offset):dlon:(lonlim(2)-dlon/2);
    lat_row = repmat(lat_j, size(lon_row));

    clat = [clat, lat_row]; %#ok<AGROW>
    clon = [clon, lon_row]; %#ok<AGROW>
end

%% FILTER TO LAND INSIDE CONUS BBOX
if landOnly
    try
        S = shaperead('landareas', 'UseGeoCoords', true);
        idxUSA = find(strcmp({S.Name}, 'United States'), 1);
        if ~isempty(idxUSA)
            polyLat = S(idxUSA).Lat;
            polyLon = S(idxUSA).Lon;

            % Limit polygon to CONUS bbox to avoid Alaska parts
            in_bbox = clat >= latlim(1) & clat <= latlim(2) & clon >= lonlim(1) & clon <= lonlim(2);
            in_poly = inpolygon(clat, clon, polyLat, polyLon);
            keep    = in_bbox & in_poly;

            clat = clat(keep);
            clon = clon(keep);
        else
            warning('United States polygon not found in landareas. Using bbox only.');
            in_bbox = clat >= latlim(1) & clat <= latlim(2) & clon >= lonlim(1) & clon <= lonlim(2);
            clat = clat(in_bbox); clon = clon(in_bbox);
        end
    catch
        warning('Could not read landareas shapefile. Using bbox only.');
        in_bbox = clat >= latlim(1) & clat <= latlim(2) & clon >= lonlim(1) & clon <= lonlim(2);
        clat = clat(in_bbox); clon = clon(in_bbox);
    end
else
    in_bbox = clat >= latlim(1) & clat <= latlim(2) & clon >= lonlim(1) & clon <= lonlim(2);
    clat = clat(in_bbox); clon = clon(in_bbox);
end

numNodes = numel(clat);

%% COVERAGE COUNTS GRID
latVec = latlim(1):degRes:latlim(2);
lonVec = lonlim(1):degRes:lonlim(2);
[LonGrid, LatGrid] = meshgrid(lonVec, latVec);
counts = zeros(size(LatGrid), 'uint16');

% Great circle distances for each node to all grid cells
for i = 1:numNodes
    dAng = distance("gc", clat(i), clon(i), LatGrid, LonGrid);
    dKm  = deg2rad(dAng)*re/1e3;
    counts = counts + uint16(dKm <= r_km);
end

%% SIMPLE COVERAGE METRICS
areaCell_km2 = (111.32*degRes) * (111.32*degRes) .* cosd(LatGrid);  % approximate cell area
areaCell_km2(areaCell_km2<0) = 0;

A_total   = nansum(areaCell_km2(:));
A_cov1    = nansum(areaCell_km2(counts>=1));
A_cov2    = nansum(areaCell_km2(counts>=2));
A_cov3    = nansum(areaCell_km2(counts>=3));

fprintf('Nodes: %d, sep approx %g km, range approx %.0f km\n', numNodes, sep_km, r_km);
fprintf('Area any coverage: %.1f percent of domain\n', 100*A_cov1/A_total);
fprintf('Area two plus:      %.1f percent of domain\n', 100*A_cov2/A_total);
fprintf('Area three plus:    %.1f percent of domain\n', 100*A_cov3/A_total);

%% PLOT
figure('Color','w');
gx = geoaxes('Basemap', baseMapName); hold(gx,'on');
geolimits(gx, latlim, lonlim);

% Paint coverage
idx1 = counts==1;
idx2 = counts>=2;

if any(idx1(:))
    geoscatter(gx, LatGrid(idx1), LonGrid(idx1), degRes*200, 'b', 'filled', 'MarkerFaceAlpha', 1.0);
end
if any(idx2(:))
    geoscatter(gx, LatGrid(idx2), LonGrid(idx2), degRes*200, 'r', 'filled', 'MarkerFaceAlpha', 1.0);
end

% Plot nodes
geoscatter(gx, clat, clon, 24, 'k', 'filled');

% Optional state outlines if available
try
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    for k = 1:numel(states)
        geoplot(gx, states(k).Lat, states(k).Lon, 'w', 'LineWidth', 0.5);
    end
end

legend({'Single Balloon Coverage','Two Plus Balloons Coverage','Balloons'}, 'Location','bestoutside');

if useHorizon
    modeStr = sprintf('Horizon based range about %.0f km', r_km);
else
    modeStr = sprintf('Manual range equals %.0f km', r_km);
end
title(gx, sprintf('CONUS hex mesh — %d nodes, sep equals %g km, %s', numNodes, sep_km, modeStr));
hold(gx,'off');
