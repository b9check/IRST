function p = irst_params_default()
% Default parameter struct for IRST HyperWatch vs SBIR comparisons.
% Values from Driggers Ch17 Table 17.1 (LWIR 8–10.1 µm) + engineering assumptions.

p = struct();

% Atmos scenario
p.ATM_SCENARIO = 'Std';   % 'Std'|'Dry'|'Humid'

% Operating mode (integration policy)
% 'hybrid'    : HW integrates to cap; SBIR auto-fill.
% 'auto_both' : both auto-fill to design_fill.
% 'fixed'     : both use p.fixed_tint_s (capped by max_tint_s).
p.OP_MODE = 'fixed';
p.fixed_tint_s = 0.05e-3;   % used only if OP_MODE='fixed'

% Geometry
p.Re_km        = 6371;
p.z_sensor_km  = 20;       % HW
p.z_target_km  = 50;       % missile
p.z_sbir_km    = 750;      % LEO
p.ground_sep_km= 0;        % ground arc (km) between platform & target

% Target aero heating
p.Mach_tgt     = 8;
p.gamma_air    = 1.4;
p.recovery     = 0.9;
p.eps_tgt      = 0.9;
d = 1.5;
p.Atgt_m2      = (pi/4)*d^2;

% Sensor: Teledyne FLIR Neutrino SX8 (MWIR)
p.band_lo_um     = 3.4;      % MWIR band (was 8.0)
p.band_hi_um     = 5.1;      % MWIR band (was 10.1)
p.lambda_ave_um  = (p.band_lo_um + p.band_hi_um)/2;  % ~4.25 µm
p.Fnum = 2.5;
p.f_m  = 0.03;         % 18 mm lens
p.D_ap_m = p.f_m / p.Fnum;   % 7.2 mm aperture
p.Fnum           = 2.5;      % f/2.5 option (was 1.5) - or use 3.0 or 4.0
p.dcc_m          = 8e-6;     % 8 µm pitch (was 15 µm) - SX8 specific
p.tau_opt        = 0.8;     % 70% - reasonable ✓
p.QE             = 0.90;     % 75% - reasonable for HOT FPA ✓
p.wcap_e         = 2.6e6;    % 2.6 million electrons (was 12e6) - from datasheet
p.design_fill    = 0.65;     % 65% - reasonable ✓
p.nread_rms_e    = 100;      % ~100 e- (was 1000) - typical for modern HOT FPA
p.dark_A_per_cm2 = 1e-9;     % ~1 nA/cm² (was 6.1e-7) - typical for HOT MWIR
p.frame_rate_Hz  = 60;       % 60 Hz ✓
p.max_tint_s     = 16e-3;    % 16 ms max (was 1.7ms) - from datasheet

% Lens emission
p.lens_T_K   = [];
p.lens_emiss = 0.9;

% Backgrounds
p.T_bg_sky_HW_K  = 3;
p.T_bg_earth_K   = 290;

% External tau file (optional)
p.use_external_tau = false;
p.tau_table_path   = 'tau_Lbg_vs_alt_clean.txt';

%% ========== CALCULATE AND STORE NEdT ==========
p.NETD_mK = calculate_NETD(p);

end

%% ========== NEdT CALCULATION FUNCTION ==========
function NETD_mK = calculate_NETD(p)
% Calculate NEdT with current parameters
% Assumes maximum integration time for best-case NEdT

fprintf('\n╔════════════════════════════════════════════════════════╗\n');
fprintf('║              NEdT Calculation (Neutrino SX8)         ║\n');
fprintf('╚════════════════════════════════════════════════════════╝\n\n');

%% Test conditions
T_test = 303;  % 30°C blackbody (standard NEdT test condition)
tint = p.max_tint_s;  % Use maximum integration time

fprintf('Test Conditions:\n');
fprintf('  Blackbody temp:     %.0f K (%.0f°C)\n', T_test, T_test-273.15);
fprintf('  Integration time:   %.2f ms\n', tint * 1000);
fprintf('  Aperture:           %.2f cm\n', p.D_ap_m * 100);
fprintf('  F-number:           f/%.1f\n\n', p.Fnum);

%% Calculate Planck radiance and derivative
L_test = blackbodyRadianceBand(T_test, p.band_lo_um, p.band_hi_um, 1.0);
L_test_plus = blackbodyRadianceBand(T_test + 1, p.band_lo_um, p.band_hi_um, 1.0);
dL_dT = L_test_plus - L_test;  % W/cm²·sr·K

%% Geometric factors
f_m = p.Fnum * p.D_ap_m;
F_lambda_over_d = p.Fnum * p.lambda_ave_um / (p.dcc_m * 1e6);
PVF = 0.0335 + 0.9665 * exp(-0.887 * F_lambda_over_d);
Omega_pix = (p.dcc_m / f_m)^2;  % steradians
Aap_m2 = pi * (p.D_ap_m^2) / 4;

fprintf('Optical Parameters:\n');
fprintf('  F*λ/d:              %.3f\n', F_lambda_over_d);
fprintf('  PVF:                %.3f\n\n', PVF);

%% Photon conversion
h = 6.62607015e-34;
c = 2.99792458e8;
lambda_ave_m = p.lambda_ave_um * 1e-6;
phot_per_W = lambda_ave_m / (h * c);

%% Background electrons per second
P_bg_W = (L_test * 1e4) * Omega_pix * Aap_m2 * p.tau_opt;
Nbg_per_s = P_bg_W * phot_per_W * p.QE;

%% Dark current electrons per second
A_pix_cm2 = (p.dcc_m^2) * 1e4;
q_e = 1.602176634e-19;
Idark_A = p.dark_A_per_cm2 * A_pix_cm2;
Ndark_per_s = Idark_A / q_e;

%% Electrons during integration
Nbg = Nbg_per_s * tint;
Ndark = Ndark_per_s * tint;
N_noise = sqrt(Nbg + Ndark + p.nread_rms_e^2);

%% Signal per Kelvin
dP_dT_W = (dL_dT * 1e4) * Omega_pix * Aap_m2 * p.tau_opt;
dN_dT = dP_dT_W * phot_per_W * p.QE * PVF * tint;

%% NEdT
NETD_K = N_noise / dN_dT;
NETD_mK = NETD_K * 1000;

fprintf('Noise Budget:\n');
fprintf('  Background:         %.0f e⁻\n', Nbg);
fprintf('  Dark current:       %.0f e⁻\n', Ndark);
fprintf('  Read noise:         %.0f e⁻\n', p.nread_rms_e);
fprintf('  Total noise:        %.0f e⁻ rms\n\n', N_noise);

fprintf('Signal:\n');
fprintf('  dN/dT:              %.0f e⁻/K\n\n', dN_dT);

fprintf('╔════════════════════════════════════════════════════════╗\n');
fprintf('║  Calculated NEdT:  %6.1f mK                          ║\n', NETD_mK);
fprintf('║  (Datasheet spec: <38 mK)                            ║\n');
fprintf('╚════════════════════════════════════════════════════════╝\n\n');

if NETD_mK <= 38
    fprintf('✓ Meets datasheet specification!\n\n');
else
    fprintf('⚠️  Does not meet spec. Need to adjust:\n');
    fprintf('   - Increase τ_opt × QE (current: %.3f)\n', p.tau_opt * p.QE);
    fprintf('   - Or reduce read noise (current: %.0f e⁻)\n\n', p.nread_rms_e);
end

end

%% ========== HELPER: Band-integrated Planck radiance ==========
function Lband = blackbodyRadianceBand(T_K, lam_lo_um, lam_hi_um, epsi)
% Returns band-integrated radiance [W/cm²·sr] from lam_lo to lam_hi
if nargin < 4, epsi = 1; end
N = 2000;
lam = linspace(lam_lo_um, lam_hi_um, N);
Lspec = epsi * planck_radiance(lam, T_K);
Lband = trapz(lam, Lspec);
end

%% ========== HELPER: Planck spectral radiance ==========
function L = planck_radiance(lambda_um, T)
% Returns spectral radiance [W/cm²·sr·µm]
c1 = 3.741771852e8;   % W·µm⁴/m²
c2 = 1.438776877e4;   % µm·K
lam = lambda_um;
expo = exp(c2./(lam*T)) - 1;
M_W_m2_um = c1 ./ (lam.^5 .* expo);
L = (M_W_m2_um/pi) / 1e4;  % Convert to W/cm²·sr·µm
end