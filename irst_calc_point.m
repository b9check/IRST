function res = irst_calc_point(p)
% Core IRST radiometry/SNR engine.
% Inputs:
%   p = parameter struct from irst_params_default (fields documented there).
% Output:
%   res = struct with geometry, intermediate rates, electrons, SNR, NEI.

%% ---------- Derived geometry ----------
psi = p.ground_sep_km / p.Re_km;  % rad central angle
rs   = p.Re_km + p.z_sensor_km;
rt   = p.Re_km + p.z_target_km;
rsat = p.Re_km + p.z_sbir_km;

R_HW_km   = sqrt(rs^2   + rt^2   - 2*rs*rt*cos(psi));
R_SBIR_km = sqrt(rsat^2 + rt^2   - 2*rsat*rt*cos(psi));

R_HW_m   = R_HW_km   * 1e3;
R_SBIR_m = R_SBIR_km * 1e3;

%% ---------- Atmosphere: vertical τ tables ----------
[z_tab, tau_tab] = tau_table_embedded(p.ATM_SCENARIO);
if p.use_external_tau
    [z_tab, tau_tab] = tau_table_external(p.tau_table_path);
end
tau2space = @(zkm) interp1(z_tab, tau_tab, zkm, 'linear','extrap');

tau_SBIR_vert = tau2space(0);                   % surface->space
tau_HW_20     = tau2space(p.z_sensor_km);       % 20->space
tau_HW_40     = tau2space(p.z_target_km);       % 40->space
tau_HW_vert_20to40 = tau_HW_20 / max(tau_HW_40,1e-6);

% slant airmass scaling for SBIR (secant approx; clamp)
airmass_SB = 1 / max(cos(psi), 0.1);
tau_SBIR_slant = tau_SBIR_vert^airmass_SB;

% HW angle dependence negligible over 20–40 km shell
tau_HW_slant = tau_HW_vert_20to40;

%% ---------- Background radiances ----------
L_bg_HW_W = blackbodyRadianceBand(p.T_bg_sky_HW_K, p.band_lo_um, p.band_hi_um, 1.0);

L_bg_earth_W = blackbodyRadianceBand(p.T_bg_earth_K, p.band_lo_um, p.band_hi_um, 1.0);
switch lower(p.ATM_SCENARIO)
    case 'std',   path_frac = 0.10;
    case 'dry',   path_frac = 0.05;
    case 'humid', path_frac = 0.20;
    otherwise,    path_frac = 0.10;
end
L_path_SBIR_W = path_frac * L_bg_earth_W;
L_bg_SBIR_W   = tau_SBIR_slant * L_bg_earth_W + L_path_SBIR_W;

%% ---------- Target temperature ----------
T_amb_K = ussaTemp(p.z_target_km);
eta_surface = 0.4;   % fraction of stagnation T rise (20% typical)
T_tgt_K = T_amb_K + eta_surface * (T_amb_K * ((1 + p.recovery*(p.gamma_air-1)/2 * p.Mach_tgt^2) - 1));
disp(T_tgt_K)
L_tgt_W = blackbodyRadianceBand(T_tgt_K, p.band_lo_um, p.band_hi_um, p.eps_tgt);

% differential radiance (clip >=0)
dL_HW_W   = max(L_tgt_W * tau_HW_slant   - L_bg_HW_W,   0);
dL_SBIR_W = max(L_tgt_W * tau_SBIR_slant - L_bg_SBIR_W, 0);

%% ---------- Sensor factors ----------
lambda_ave_m = p.lambda_ave_um * 1e-6;
f_m          = p.Fnum * p.D_ap_m;
Omega_pix    = (p.dcc_m / f_m)^2; % sr approx

% dark current
A_pix_m2   = p.dcc_m^2;
A_pix_cm2  = A_pix_m2 * 1e4;
q_e        = 1.602176634e-19;
Idark_A    = p.dark_A_per_cm2 * A_pix_cm2;
ndark_per_s= Idark_A / q_e;

% photon conversion
h = 6.62607015e-34;
c = 2.99792458e8;
Ephoton_J  = h*c / lambda_ave_m;
phot_per_W = 1 / Ephoton_J;

% PVF
F_lambda_over_d = p.Fnum * p.lambda_ave_um / (p.dcc_m*1e6);
PVF = 0.0335 + 0.9665*exp(-0.887*F_lambda_over_d);

%% ---------- Signal electron rates ----------
G_HW   = (pi*p.D_ap_m^2/4) / (R_HW_m^2);
G_SBIR = (pi*p.D_ap_m^2/4) / (R_SBIR_m^2);
Atgt_cm2 = p.Atgt_m2 * 1e4;

P_sig_HW_W   = dL_HW_W   * Atgt_cm2 * G_HW;
P_sig_SBIR_W = dL_SBIR_W * Atgt_cm2 * G_SBIR;

Phi_sig_HW   = P_sig_HW_W   * phot_per_W;
Phi_sig_SBIR = P_sig_SBIR_W * phot_per_W;
Nsig_HW_per_s   = Phi_sig_HW   * p.tau_opt * PVF * p.QE;
Nsig_SBIR_per_s = Phi_sig_SBIR * p.tau_opt * PVF * p.QE;

%% ---------- Background electron rates ----------
Aap_m2 = pi*(p.D_ap_m^2)/4;
Lp_HW_Wm2sr   = L_bg_HW_W   * 1e4;
Lp_SBIR_Wm2sr = L_bg_SBIR_W * 1e4;
Pbg_HW_W   = Lp_HW_Wm2sr   * Omega_pix * Aap_m2 * p.tau_opt;  % no τ
Pbg_SBIR_W = Lp_SBIR_Wm2sr * Omega_pix * Aap_m2 * p.tau_opt;  % already slant
Phi_bg_HW   = Pbg_HW_W   * phot_per_W;
Phi_bg_SBIR = Pbg_SBIR_W * phot_per_W;
Nbg_HW_per_s   = Phi_bg_HW   * p.QE;
Nbg_SBIR_per_s = Phi_bg_SBIR * p.QE;

%% ---------- Lens emission ----------
if isempty(p.lens_T_K)
    Nlens_per_s = 0;
else
    L_lens_W = blackbodyRadianceBand(p.lens_T_K, p.band_lo_um, p.band_hi_um, p.lens_emiss);
    L_lens_Wm2sr = L_lens_W * 1e4;
    P_lens_W = L_lens_Wm2sr * Omega_pix * Aap_m2 * p.tau_opt;
    Phi_lens = P_lens_W * phot_per_W;
    Nlens_per_s = Phi_lens * p.QE;
end

%% ---------- Integration times ----------
switch lower(p.OP_MODE)
    case 'auto_both'
        denom_HW   = max(Nbg_HW_per_s + Nlens_per_s + ndark_per_s, 1);
        denom_SBIR = max(Nbg_SBIR_per_s + Nlens_per_s + ndark_per_s, 1);
        tint_HW_s   = min(p.max_tint_s, (p.wcap_e*p.design_fill) / denom_HW);
        tint_SBIR_s = min(p.max_tint_s, (p.wcap_e*p.design_fill) / denom_SBIR);

    case 'fixed'
        tint_HW_s   = min(p.fixed_tint_s, p.max_tint_s);
        tint_SBIR_s = tint_HW_s;

    case 'hybrid'
        % HW integrates to cap (long integration freedom in cold sky)
        tint_HW_s   = p.max_tint_s;
        % SBIR background-limited
        denom_SBIR  = max(Nbg_SBIR_per_s + Nlens_per_s + ndark_per_s, 1);
        tint_SBIR_s = min(p.max_tint_s, (p.wcap_e*p.design_fill) / denom_SBIR);

    otherwise
        error('Unknown OP_MODE.');
end

%% ---------- Electrons / integration ----------
Nsig_HW   = Nsig_HW_per_s   * tint_HW_s;
Nbg_HW    = Nbg_HW_per_s    * tint_HW_s;
Nlens_HW  = Nlens_per_s     * tint_HW_s;
Ndark_HW  = ndark_per_s     * tint_HW_s;

Nsig_SBIR = Nsig_SBIR_per_s * tint_SBIR_s;
Nbg_SBIR  = Nbg_SBIR_per_s  * tint_SBIR_s;
Nlens_SB  = Nlens_per_s     * tint_SBIR_s;
Ndark_SB  = ndark_per_s     * tint_SBIR_s;

%% ---------- Noise ----------
Nnoise_HW   = sqrt( Nbg_HW + Nlens_HW + Ndark_HW + p.nread_rms_e^2 );
Nnoise_SBIR = sqrt( Nbg_SBIR + Nlens_SB + Ndark_SB + p.nread_rms_e^2 );

%% ---------- SNR ----------
SNR_HW   = Nsig_HW   / Nnoise_HW;
SNR_SBIR = Nsig_SBIR / Nnoise_SBIR;

%% ---------- NEI scaling ----------
hc = 6.62607015e-34 * 2.99792458e8;
K_common  = (pi*(p.D_ap_m^2)/4) * p.tau_opt * PVF * p.QE * (lambda_ave_m/hc);
NEI_HW_ph  = (Nnoise_HW)   / (K_common * tint_HW_s);
NEI_SB_ph  = (Nnoise_SBIR) / (K_common * tint_SBIR_s);

%% ---------- Fill fractions ----------
fill_HW = (Nbg_HW + Nlens_HW + Ndark_HW)/p.wcap_e;
fill_SB = (Nbg_SBIR + Nlens_SB + Ndark_SB)/p.wcap_e;

%% ---------- Package ----------
res = struct();
res.p              = p;
res.psi_rad        = psi;
res.R_HW_m         = R_HW_m;
res.R_SBIR_m       = R_SBIR_m;
res.tau_SBIR_vert  = tau_SBIR_vert;
res.tau_SBIR_slant = tau_SBIR_slant;
res.tau_HW_slant   = tau_HW_slant;
res.L_bg_HW_W      = L_bg_HW_W;
res.L_bg_SBIR_W    = L_bg_SBIR_W;
res.T_tgt_K        = T_tgt_K;
res.dL_HW_W        = dL_HW_W;
res.dL_SBIR_W      = dL_SBIR_W;
res.Nsig_HW        = Nsig_HW;
res.Nbg_HW         = Nbg_HW;
res.Nlens_HW       = Nlens_HW;
res.Ndark_HW       = Ndark_HW;
res.Nnoise_HW      = Nnoise_HW;
res.Nsig_SBIR      = Nsig_SBIR;
res.Nbg_SBIR       = Nbg_SBIR;
res.Nlens_SB       = Nlens_SB;
res.Ndark_SB       = Ndark_SB;
res.Nnoise_SBIR    = Nnoise_SBIR;
res.SNR_HW         = SNR_HW;
res.SNR_SBIR       = SNR_SBIR;
res.NEI_HW_ph      = NEI_HW_ph;
res.NEI_SB_ph      = NEI_SB_ph;
res.tint_HW_s      = tint_HW_s;
res.tint_SBIR_s    = tint_SBIR_s;
res.fill_HW        = fill_HW;
res.fill_SB        = fill_SB;
res.Nbg_HW_per_s   = Nbg_HW_per_s;
res.Nbg_SBIR_per_s = Nbg_SBIR_per_s;
res.PVF            = PVF;
res.F_lambda_over_d= F_lambda_over_d;
end

%% --------- Helper: built-in τ(z→space) tables --------------------------
function [z_tau_km, tau_tab] = tau_table_embedded(scenario)
% Engineering LWIR 8–10.1 µm band transmittance to space (mid-lat clear)
z_tau_km = [0 5 10 15 20 25 30 35 40 45 50]';
tau_Std   = [0.80 0.88 0.93 0.96 0.97 0.98 0.985 0.990 0.994 0.996 0.998]';
tau_Dry   = [0.88 0.93 0.96 0.98 0.99 0.995 0.997 0.998 0.999 1.000 1.000]';
tau_Humid = [0.68 0.77 0.84 0.89 0.92 0.94 0.96 0.97 0.98 0.985 0.990]';
switch lower(scenario)
    case 'std',   tau_tab = tau_Std;
    case 'dry',   tau_tab = tau_Dry;
    case 'humid', tau_tab = tau_Humid;
    otherwise,    error('Bad ATM_SCENARIO string.');
end
end

%% --------- Helper: load & monotonic-clean external τ table -------------
function [z, tau] = tau_table_external(path)
tbl = readmatrix(path);
z   = tbl(:,1);
tau = tbl(:,2);
tau = cummax(tau);  % enforce monotonic increasing
end


%% --------- Helper: cumulative max (monotonic enforce) -------------------
function out = cummax(v)
out = v;
for i = 2:numel(v)
    if out(i) < out(i-1), out(i) = out(i-1); end
end
end

%% --------- Helper: US Standard Atmos temperature vs altitude ------------
function T = ussaTemp(zkm)
% Piece-wise approximate US Std Atmos 1976 temperature (K) up to 86 km.
z = zkm*1000;  % m
if z<11000
    T = 288.15 - 6.5e-3*z;
elseif z<20000
    T = 216.65;
elseif z<32000
    T = 216.65 + (z-20000)*1.0e-3;
elseif z<47000
    T = 228.65 + (z-32000)*2.8e-3;
elseif z<51000
    T = 270.65;
elseif z<71000
    T = 270.65 - (z-51000)*2.8e-3;
else
    T = 214.65 - (z-71000)*2.0e-3;
end
end

%% --------- Helper: Band-integrated black-body radiance -------------------
function Lband = blackbodyRadianceBand(T_K,lam_lo_um,lam_hi_um,epsi)
% Returns band-integrated radiance [W/cm²·sr] from lam_lo to lam_hi.
% epsi = emissivity (1 for black-body).
if nargin<4, epsi = 1; end
N = 2000;
lam = linspace(lam_lo_um,lam_hi_um,N);
Lspec = epsi * planck_radiance(lam,T_K);   % W/cm²·sr·µm
Lband = trapz(lam,Lspec);                  % integrate over µm
end

%% --------- Helper: Planck spectral radiance (W/cm²·sr·µm) ---------------
function L = planck_radiance(lambda_um,T)
c1 = 3.741771852e8;   % W·µm⁴/m² (2πhc² *1e24 /1e4)
c2 = 1.438776877e4;   % µm·K   (hc/k)
lam = lambda_um;
expo = exp(c2./(lam*T)) - 1;
M_W_m2_um = c1 ./ (lam.^5 .* expo);    % spectral exitance [W/m²·µm]
L = (M_W_m2_um/pi) / 1e4;              % convert to radiance [W/cm²·sr·µm]
end
