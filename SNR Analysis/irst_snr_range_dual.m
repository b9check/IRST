function irst_snr_range_dual(range_km, isHumid)
% Plot HyperWatch & SBIR SNR vs HyperWatch slant-range (10–300 km etc.).
%
% Inputs
%   range_km : [Rmin Rmax] 2-vector  → auto-fills 200 points
%            : OR explicit vector of slant ranges to evaluate.
%   isHumid  : true  → use Humid atmosphere table
%            : false → use Standard atmosphere table
%
% Example
%   irst_snr_range_dual([10 300], true);   % Humid, sweep 10–300 km
%
% Requires irst_params_default.m and irst_calc_point.m on the MATLAB path.

% ---------- build vector of ranges ----------
if numel(range_km)==2                 % treat as [min max]
    range_km = linspace(range_km(1), range_km(2), 200);
end
scenario = tern(isHumid,'Humid','Std');

% ---------- base parameter struct ----------
p              = irst_params_default();
Re = p.Re_km;                 % Earth radius (km)
rs = Re + p.z_sensor_km;      % HW sensor radius
rt = Re + p.z_target_km;      % target radius
rsat = Re + p.z_sbir_km;      % satellite radius
speed = p.Mach_tgt;
alt = p.z_target_km;
time_int = p.fixed_tint_s;

% ---------- allocate result arrays ----------
n        = numel(range_km);
SNR_HW   = zeros(1,n);
SNR_SB   = zeros(1,n);
R_SB_km  = zeros(1,n);

% ---------- sweep ----------
for i = 1:n
    R_hw = range_km(i);

    % invert curved-Earth slant-range → ground arc ψ
    cospsi = (rs^2 + rt^2 - R_hw^2) / (2*rs*rt);
    if cospsi < -1 || cospsi > 1
        SNR_HW(i)=NaN; SNR_SB(i)=NaN; R_SB_km(i)=NaN;
        continue
    end
    psi = acos(cospsi);                  % rad
    p.ground_sep_km = psi * Re;          % km

    res = irst_calc_point(p);
    SNR_HW(i)  = res.SNR_HW;
    SNR_SB(i)  = res.SNR_SBIR;
    R_SB_km(i) = res.R_SBIR_m/1e3;
end

% ---------- single semilog-y plot ----------
figure; hold on; grid on;
semilogy(range_km, SNR_HW,  '-', 'LineWidth',1.6);
semilogy(range_km, SNR_SB, '--', 'LineWidth',1.6);
yline(10, 'k--', 'LineWidth', 1.4, 'Label', 'Detection Threshold (SNR=10)', ...
    'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','top');
set(gca,'YScale','log')   % <-- keep axis in log space
xlabel('Slant range (km)');
ylabel('SNR per Integration');
% title(sprintf('SNR vs Slant Range — Mach %.0f Missile at %.0f km Altitude (Neutrino SX8, t_{int}=%.2f ms)', speed, alt, time_int*1000));
legend({'HyperWatch','LEO Satellite'},'Location','best');
end

% ---------- tiny ternary helper ----------
function y = tern(cond,a,b); if cond, y=a; else, y=b; end; end
