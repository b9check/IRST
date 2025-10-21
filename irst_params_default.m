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
p.OP_MODE = 'hybrid';
p.fixed_tint_s = 100e-6;   % used only if OP_MODE='fixed'

% Geometry
p.Re_km        = 6371;
p.z_sensor_km  = 20;       % HW
p.z_target_km  = 40;       % missile
p.z_sbir_km    = 750;      % LEO
p.ground_sep_km= 0;        % ground arc (km) between platform & target

% Target aero heating
p.Mach_tgt     = 8;
p.gamma_air    = 1.4;
p.recovery     = 1.0;
p.eps_tgt      = 1.0;
p.Atgt_m2      = 1.0;

% Sensor: Teledyne FLIR Neutrino SX8 (MWIR)
p.band_lo_um     = 3.4;      % MWIR band (was 8.0)
p.band_hi_um     = 5.1;      % MWIR band (was 10.1)
p.lambda_ave_um  = (p.band_lo_um + p.band_hi_um)/2;  % ~4.25 µm
p.D_ap_m         = 2.83e-2;  % 2.83 cm (reasonable for f/2.5 or custom lens)
p.Fnum           = 2.5;      % f/2.5 option (was 1.5) - or use 3.0 or 4.0
p.dcc_m          = 8e-6;     % 8 µm pitch (was 15 µm) - SX8 specific
p.tau_opt        = 0.70;     % 70% - reasonable ✓
p.QE             = 0.75;     % 75% - reasonable for HOT FPA ✓
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
p.T_bg_sky_HW_K  = 220;
p.T_bg_earth_K   = 290;

% External tau file (optional)
p.use_external_tau = false;
p.tau_table_path   = 'tau_Lbg_vs_alt_clean.txt';
end

