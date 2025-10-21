%% Balloon Coverage with Overlap Highlighting — concentric midpointed rings
close all; clear;

% --- CONFIG ---
lat0 = 13.4443; lon0 = 144.7937;   % Guam
re   = 6371e3;                     % Earth radius [m]

% Altitudes
hB = 20e3;   % balloon altitude [m]
hM = 30e3;   % target altitude [m]

% Coverage radius choice
useHorizon      = true;    % true uses horizon based coverage, false uses manualRange_km
manualRange_km  = 430;

% Swarm layout controls
numB          = 100;        % total balloons
sep_km        = 400;       % desired neighbor spacing in kilometers
includeCenter = false;     % place a balloon at Guam if true

% Map resolution for overlap paint
degRes = 0.1;

% --- COVERAGE RADIUS ---
if useHorizon
    phiB = acos(re/(re+hB));
    phiM = acos(re/(re+hM));
    theta = phiB + phiM;
    r_km = theta*(re+hM)/1e3;       % projected ground radius at target altitude
else
    r_km = manualRange_km;
end

% --- POSITIONS: concentric midpointed rings on the sphere ---
% Use hex like radial step so across ring spacing resembles sep
dr_km = sep_km*sqrt(3)/2;

clat = []; clon = [];
if includeCenter
    clat = lat0; clon = lon0;
end

remaining = numB - numel(clat);
k = 1;
while remaining > 0
    r_km_ring = k * dr_km;

    % Balloons on this ring so that arc spacing around the ring is about sep
    circ_km = 2*pi*r_km_ring;
    m_k = max(6, round(circ_km/sep_km));   % at least a hexagon
    take = min(m_k, remaining);

    % Bearings with midpointing
    % Odd rings have zero offset
    % Even rings shift by half of their own step so they fall between neighbors on odd rings
    step_deg = 360/m_k;
    offset_deg = 0;
    if mod(k,2) == 0
        offset_deg = step_deg/2;
    end
    bearings = (0:take-1)*step_deg + offset_deg;

    % Great circle offset from Guam by central angle for this ring
    ang_deg = rad2deg((r_km_ring*1e3)/re);
    [lat_ring, lon_ring] = reckon(lat0, lon0, ang_deg, bearings, "degrees");

    clat = [clat, lat_ring]; %#ok<AGROW>
    clon = [clon, lon_ring]; %#ok<AGROW>

    remaining = numB - numel(clat);
    k = k + 1;
end

% --- Reference horizon at target altitude (optional) ---
theta0G   = acos(re/(re+hM));
horizG_deg= rad2deg(theta0G)
[latrG, lonrG] = scircle1(lat0, lon0, horizG_deg, "degrees");

% --- GRID for overlap paint ---
centerAng_deg = distance("gc", lat0, lon0, clat, clon);
centerDist_km  = deg2rad(centerAng_deg)*re/1e3;
maxCenterDist_km = max([0, centerDist_km]);
maxReach_km = maxCenterDist_km + r_km;

latPad = maxReach_km / 111.32;
lonPad = maxReach_km / (111.32*cosd(lat0)+eps);

latlim = [lat0 - latPad, lat0 + latPad];
lonlim = [lon0 - lonPad, lon0 + lonPad];

latVec = latlim(1):degRes:latlim(2);
lonVec = lonlim(1):degRes:lonlim(2);
[LonGrid, LatGrid] = meshgrid(lonVec, latVec);
counts = zeros(size(LatGrid));

% --- PAINT COVERAGE COUNTS ---
for i = 1:numel(clat)
    dAng = distance("gc", clat(i), clon(i), LatGrid, LonGrid);
    dKm  = deg2rad(dAng)*re/1e3;
    mask = (dKm <= r_km);
    counts(mask) = counts(mask) + 1;
end

% --- PLOT ---
figure("Color","w");
gx = geoaxes("Basemap","satellite"); hold(gx,"on");

idx1 = counts==1;
idx2 = counts>=2;
if any(idx1(:))
    h1 = geoscatter(gx, LatGrid(idx1), LonGrid(idx1), degRes*200, "b", "filled", "MarkerFaceAlpha", 1.0);
else
    h1 = geoscatter(gx, NaN, NaN);
end
if any(idx2(:))
    h2 = geoscatter(gx, LatGrid(idx2), LonGrid(idx2), degRes*200, "r", "filled", "MarkerFaceAlpha", 1.0);
else
    h2 = geoscatter(gx, NaN, NaN);
end

pgG = geopolyshape(latrG, wrapTo360(lonrG));
hG  = geoplot(gx, pgG, "FaceColor", [0.6 0.3 0], "FaceAlpha", 1.0, "EdgeColor", "none");

hC = geoscatter(gx, lat0, lon0, 100, "y", "filled");
hBalloons = geoscatter(gx, clat, clon, 50, "k", "filled");

legend([h1, h2, hBalloons, hG, hC], ...
    {"Single Balloon Coverage","Two Plus Balloons Coverage","Balloons","Ground Radar Coverage","Guam"}, ...
    "Location","bestoutside");

if useHorizon
    modeStr = sprintf("Horizon based range about %.0f km", r_km);
else
    modeStr = sprintf("Manual range equals %.0f km", r_km);
end
% title(gx, sprintf("Concentric midpointed rings — %d nodes, sep equals %g km, %s", numel(clat), sep_km, modeStr));

hold(gx,"off");
