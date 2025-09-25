function eps_eff = wetsnow_permittivity_tinga73_extended(frequency, temperature, density, liquid_water, ...
    ice_perm_func, water_perm_func, params)

% WETSNOW_PERMITTIVITY_TINGA73_EXTENDED
% Extended Tinga model with density scaling and Bruggeman high-LWC blend

% Defaults
if nargin < 5 || isempty(ice_perm_func)
    ice_perm_func = @(f, T) ice_permittivity_maetzler06(f, T);
end
if nargin < 6 || isempty(water_perm_func)
    water_perm_func = @(f, T) water_permittivity_maetzler87(f, T);
end
if nargin < 7
    params = struct( ...
    'eps_s', 88, ...
    'eps_inf', 4.9, ...
    'tau', 2e-10, ...
    'alpha', 0.2, ...
    'sigma', 0.001);
    % depol_bruggeman = 0.2262;%1/3;
end
if isa(water_perm_func,'function_handle')
    eps_w = water_perm_func(frequency, params);

else
    eps_w = water_perm_func;
end

% Constants
rho_ice = 917;
rho_water = 1000;
eps_air = 1;

% Power law thresholds
% Thresholds for blend zone (density-dependent)
rho_c = 1000; alpha = 2;
A_pf = 0.085; A_fs = 0.3;
LWC_pf = A_pf * (1 - (density / rho_c).^alpha);
LWC_fs = A_fs * (1 - (density / rho_c).^alpha);


% Convert LWC mass fraction to volume fraction
% vw = liquid_water .* rho_water ./ (liquid_water .* rho_water + (1 - liquid_water) .* rho_ice);
% Assume we are inputing Volume Fractions!
vw = liquid_water;

% Get permittivities
% params = struct( ...
%     'eps_s', 88, ...
%     'eps_inf', 4.9, ...
%     'tau', 1e-9, ...
%     'alpha', 0.0896, ...
%     'sigma', 0.001);
% eps_w = water_perm_func(frequency, params);
% eps_w = water_perm_func(frequency, temperature);
eps_i = ice_perm_func(frequency, temperature);

% Volume ratios
Vw_i = 1 + rho_ice / rho_water * vw ./ (1 - vw);
Va_i = (rho_ice ./ density) .* (1 + vw ./ (1 - vw));

% Depolarization (fixed)
n = 0.2262;
prefactor = (1 - n) / n;

% Tinga formulation
alpha_val = 2 * eps_w + eps_i;
diff_wi = eps_w - eps_i;
diff_wa = eps_w - eps_air;

numer = Vw_i ./ Va_i .* diff_wa .* alpha_val - (1 ./ Va_i) .* diff_wi .* (2 * eps_w + eps_air);
denom = (2 * eps_air + eps_w) .* alpha_val ...
      - 2 ./ Vw_i .* diff_wa .* diff_wi ...
      - (Vw_i ./ Va_i) .* diff_wa .* alpha_val ...
      + (1 ./ Va_i) .* diff_wi .* (2 * eps_w + eps_air);

eps_tinga = eps_air .* (1 + prefactor .* numer ./ denom);
eps_bruggeman = zeros(1,length(liquid_water));
for ii = 1:length(liquid_water)
% Bruggeman model (symmetric, 3-phase)
[frac_volume, fi, fw] = compute_frac_volumes(density, liquid_water(ii));
frac_volume = min(frac_volume, 0.999);  % Clamp total inclusion
fa = max(1 - frac_volume, eps);  % Clamp air
% fa = 1 - frac_volume;
deps = [eps_i, eps_w, eps_air];
f = [fi, fw, fa];
% nvals = repmat(depol_bruggeman, 1, 3);
nvals = [1/3,0.05,0.1];
% Newton Solver
eps_bruggeman(ii) = bruggeman_newton_solver(f, deps, nvals);
% Genreal fsolve
% eps_bruggeman(ii) = bruggeman_symmetric_solver(f, deps, nvals);
end
% Join Bruggeman and Tinga Consistently
[~,scaleIx] = min(abs(liquid_water-A_pf));
if abs(eps_bruggeman(scaleIx)) < 1e-3
    epsScale = 1;  % avoid blowup
else
    epsScale = eps_tinga(scaleIx) ./ eps_bruggeman(scaleIx);
end
% epsScale = eps_tinga(scaleIx)./eps_bruggeman(scaleIx);
eps_bruggeman = eps_bruggeman.*epsScale;
% Smooth blend between Tinga and Bruggeman
w = ones(size(vw));
blend_zone = vw > LWC_pf & vw < LWC_fs;
w(vw >= LWC_fs) = 0;
w(blend_zone) = cos(pi/2 * (vw(blend_zone) - LWC_pf) / (LWC_fs - LWC_pf)).^2;

eps_eff = w .* eps_tinga + (1 - w) .* eps_bruggeman;

% Enforce positive imaginary
eps_eff = real(eps_eff) + 1i * abs(imag(eps_eff));

end