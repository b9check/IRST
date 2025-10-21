% Acceptable reaction time
t_react = 300; % seconds

% Drag parameters
Cd = 0.08;      % drag coefficient (–)
A  = 0.8;      % cross‐sectional area [m^2]
m  = 2500;     % mass [kg]

% Missile parameters (inputs)
v0s = linspace(341*5, 341*15, 1000) ; % m/s
hMs = linspace(30e3, 100e3, 1000); % m
dist0M = 2500e3; % Missile initial distance

% Balloon parameters
d_b_off = 1000e3; % m
hB = 20e3; % m

% Constants
re    = 6371e3;

ground_defendable = false(numel(v0s), numel(hMs));
air_defendable    = false(numel(v0s), numel(hMs));
ground_react_times = zeros(numel(v0s), numel(hMs));
air_react_times = zeros(numel(v0s), numel(hMs));
dist0B_mat        = zeros(numel(v0s), numel(hMs));
dist0G_mat        = zeros(numel(v0s), numel(hMs));

for i = 1:length(v0s)
    v0 = v0s(i);
    for j = 1:length(hMs)
        hM = hMs(j);
        rho = 1.225 * exp(-hM * 1.41e-4);
        k = 0.5 * Cd * A * rho / m;   % [1/m]

        % AIR RADAR CALCS
        phi1 = d_b_off/re;
        phi2 = acos(re/(re+hB));
        phi3 = acos(re/(re+hM));
        theta0B = phi1+phi2+phi3;  
        dist0B  = theta0B*(re + hM);  
        dist0B_mat(i,j) = dist0B;
        if dist0B < dist0M
            % closed‐form time‐to‐impact (seconds):
            v0B = v0 * exp(-k*(dist0M - dist0B));
            time_to_impactB = (exp(k*dist0B) - 1) / (k * v0B);
            air_react_times(i, j) = time_to_impactB/60;
        else
            time_to_impactB = (exp(k*dist0M) - 1) / (k * v0);
            air_react_times(i, j) = time_to_impactB/60;
        end


        if time_to_impactB > t_react || dist0B > dist0M
            air_defendable(i, j) = 1;
        end

        
        % GROUND RADAR CALCS
        % Calculate time to impact
        theta0G = acos(re/(re + hM));     
        dist0G  = theta0G*(re + hM);  
        dist0G_mat(i,j) = dist0G;
        v0G = v0 * exp(-k*(dist0M - dist0G));
        
        if dist0G < dist0M
        % closed‐form time‐to‐impact (seconds):
            time_to_impactG = (exp(k*dist0G) - 1) / (k * v0G);
            ground_react_times(i, j) = time_to_impactG/60;
        else
            time_to_impactG = (exp(k*dist0M) - 1) / (k * v0);
            ground_react_times(i, j) = time_to_impactG/60;
        end



        if time_to_impactG > t_react ||  dist0G > dist0M
            ground_defendable(i, j) = 1;
        end
    end
end

newC = zeros(size(air_defendable));
newC( air_defendable & ~ground_defendable ) = 1;
newC( ground_defendable )                   = 2;

% new 3-entry colormap: gray, orange, green
map3 = [ ...
  0.8 0.8 0.8;   % 0 = neither
  1   0.5 0;     % 1 = air only
  0   1   0;     % 2 = ground (& air)
];

figure(1);
imagesc(hMs/1e3, v0s/341, newC);
set(gca,'YDir','normal');
colormap(map3);
caxis([0 2]);
xlabel('Missile Cruising Altitude (km)');
ylabel('Missile Initial Speed (Mach)');
title('"Defendable" (5+ Minutes of Tracking before Impact) Regions for an Incoming HGV');

% legend for just the three states
hold on;
h = gobjects(3,1);
names = {'Not Defendable','Defendable with Balloons + Ground','Defendable by Just Ground'};
for idx = 0:2
  h(idx+1) = plot(NaN,NaN,'s', ...
    'MarkerFaceColor',map3(idx+1,:), ...
    'MarkerEdgeColor','k');
end
legend(h(2:3), names(2:3),'Location','northeast');
hold off;



% …[all the prior code above]…

figure(2);
% ------------------------
% Ground‐radar subplot
% ------------------------
subplot(2,1,1)
imagesc(hMs/1e3, v0s/341, ground_react_times);
set(gca,'YDir','normal');
caxis([0 12]);      % saturate at 0–11 min
cb1 = colorbar;
ylabel(cb1, 'Reaction Time (min)');
ylabel('Missile Initial Speed (Mach)');
title('Ground Radar Reaction Time');

% annotate Mach 10 @ 40 km
x_pt = 40;   % km
y_pt = 7;   % Mach
[~, i_idx] = min(abs(v0s/341 - y_pt));
[~, j_idx] = min(abs(hMs/1e3 - x_pt));
tG = ground_react_times(i_idx, j_idx);
hold on;
plot(x_pt, y_pt, 'ko', 'MarkerSize', 8, 'LineWidth', 2);
text(x_pt + 2, y_pt + 0.5, sprintf('%.1f minutes to react', tG), ...
     'Color','k', 'FontWeight','bold');
hold off;

% ------------------------
% Balloon subplot
% ------------------------
figure(2)
subplot(2,1,2)
imagesc(hMs/1e3, v0s/341, air_react_times);
set(gca,'YDir','normal');
caxis([0 12]);
cb2 = colorbar;
ylabel(cb2, 'Reaction Time (min)');
xlabel('Missile Cruising Altitude (km)');
ylabel('Missile Initial Speed (Mach)');
title('Balloon Reaction Time');

% annotate Mach 10 @ 40 km
tB = air_react_times(i_idx, j_idx);
hold on;
plot(x_pt, y_pt, 'ko', 'MarkerSize', 8, 'LineWidth', 2);
text(x_pt + 2, y_pt + 0.5, sprintf('%.1f minutes to react', tB), ...
     'Color','k', 'FontWeight','bold');
hold off;