function irst_compare_hw_sbir
% Wrapper: load defaults, run calc, print results.

p = irst_params_default();

% >>> EDIT headline comparison here <<<
% Example: choose 98 km ground arc for ~100 km HW slant:
% p.ground_sep_km = 98;     % uncomment to change
% p.OP_MODE = 'hybrid';     % 'hybrid'|'auto_both'|'fixed'
% p.ATM_SCENARIO = 'Std';   % 'Std'|'Dry'|'Humid'

res = irst_calc_point(p);
irst_print_result(res);
end

function irst_print_result(res)
p = res.p;
fprintf('---- SENSOR SUMMARY (Driggers Table 17.1 8–10.1 µm) ----\n');
fprintf('PVF(avg) = %.3f (Fλ/d = %.3f)\n',res.PVF,res.F_lambda_over_d);
fprintf('Design fill = %.0f%% of %.2f Me- well\n',p.design_fill*100,p.wcap_e/1e6);
fprintf('Mode = %s; Atmos = %s\n',p.OP_MODE,p.ATM_SCENARIO);
fprintf('Max tint cap = %.3g s; Frame period = %.3g s (%.1f Hz)\n',...
    p.max_tint_s,1/p.frame_rate_Hz,p.frame_rate_Hz);

fprintf('\n-- Geometry --\n');
fprintf('HW Range = %.0f km; SBIR Range = %.0f km  (ground sep=%.0f km)\n',...
    res.R_HW_m/1e3, res.R_SBIR_m/1e3, p.ground_sep_km);

fprintf('\n-- Atmosphere --\n');
fprintf('tau_SBIR(vertical surface->space)=%.3f; slant=%.3f\n',...
    res.tau_SBIR_vert,res.tau_SBIR_slant);
fprintf('tau_HW(20->40 km slant)=%.3f\n',res.tau_HW_slant);

fprintf('\n-- Target Thermal --\n');
fprintf('T_tgt = %.1f K\n',res.T_tgt_K);

fprintf('\n[DEBUG] Background rates: HW=%.3g e-/s, SBIR=%.3g e-/s (SBIR/HW=%.2f)\n',...
    res.Nbg_HW_per_s,res.Nbg_SBIR_per_s, res.Nbg_SBIR_per_s/max(res.Nbg_HW_per_s,eps));
fprintf('[DEBUG] Fill fractions: HW=%.2f, SBIR=%.2f (design=%.2f)\n',...
    res.fill_HW,res.fill_SB,p.design_fill);

fprintf('\n---- HYPERWATCH ----\n');
fprintf('tint = %.3g s\n',res.tint_HW_s);
fprintf('Signal e- = %.3g, Background e- = %.3g, Lens e- = %.3g, Dark e- = %.3g\n',...
    res.Nsig_HW,res.Nbg_HW,res.Nlens_HW,res.Ndark_HW);
fprintf('SNR = %.3g\n',res.SNR_HW);
fprintf('NEI ~ %.3g photons/cm^2\n',res.NEI_HW_ph);

fprintf('\n---- SBIR ----\n');
fprintf('tint = %.3g s\n',res.tint_SBIR_s);
fprintf('Signal e- = %.3g, Background e- = %.3g, Lens e- = %.3g, Dark e- = %.3g\n',...
    res.Nsig_SBIR,res.Nbg_SBIR,res.Nlens_SB,res.Ndark_SB);
fprintf('SNR = %.3g\n',res.SNR_SBIR);
fprintf('NEI ~ %.3g photons/cm^2\n',res.NEI_SB_ph);
end
