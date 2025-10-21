%% --------- Helper: built-in τ(z→space) tables (MWIR 3.4–5.1 µm) ----------
function [z_tau_km, tau_tab] = tau_table_embedded(scenario)
% Engineering MWIR 3.4–5.1 µm band transmittance to space (mid-latitude clear)
% Values roughly derived from MODTRAN5 MWIR band averages (0 deg zenith).

z_tau_km = [0 2 5 10 15 20 25 30 35 40 45 50]';

% --- Transmittance from altitude → space ---
% MWIR is nearly transparent above 10–15 km; humidity dominates near surface.
tau_Std =  [0.85 0.90 0.94 0.975 0.985 0.990 0.994 0.996 0.997 0.998 0.999 1.000]';
tau_Dry =  [0.90 0.94 0.97 0.985 0.992 0.996 0.998 0.999 0.9995 1.000 1.000 1.000]';
tau_Humid=[0.75 0.83 0.90 0.955 0.975 0.985 0.990 0.993 0.995 0.997 0.998 0.999]';

switch lower(scenario)
    case 'std',   tau_tab = tau_Std;
    case 'dry',   tau_tab = tau_Dry;
    case 'humid', tau_tab = tau_Humid;
    otherwise,    error('Bad ATM_SCENARIO string.');
end
end
