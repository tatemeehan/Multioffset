function Es = wetsnow_permittivity_tinga73(frequency, temperature, density, liquid_water, n, ice_permittivity_model, water_permittivity_model)
    % Computes the effective permittivity based on Tinga et al. (1973) for three-component mixing.
    % Parameters:
    % - frequency: frequency in Hz
    % - temperature: temperature in K
    % - density: density of snow (kg/m^3)
    % - liquid_water: liquid water fraction
    % - n: depolarization factor 0 - Needle-Like; 1/3 - Spherical; 1 - Disk
    % - ice_permittivity_model: function handle for ice permittivity model
    % - water_permittivity_model: function handle for water permittivity model
    
    % Constants
    FREEZING_POINT = 273.15; % Freezing point in Kelvin
    DENSITY_OF_ICE = 917; % Density of ice (kg/m^3)
    DENSITY_OF_WATER = 1000; % Density of water (kg/m^3)

    % Validate inputs
    if (temperature < FREEZING_POINT) && any(liquid_water > 0)
        error('Liquid water is positive, but temperature is below freezing. This is incompatible.');
    end

    % Calculate wetness (W)
    W = liquid_water .* DENSITY_OF_WATER ./ (liquid_water .* DENSITY_OF_WATER + (1 - liquid_water) .* DENSITY_OF_ICE);

    % Calculate volume ratios
    Vw_i = 1 + DENSITY_OF_ICE ./ DENSITY_OF_WATER .* W ./ (1 - W);
    Va_i = (DENSITY_OF_ICE ./ density) .* (1 + W ./ (1 - W));

    % Default permittivity models
    if nargin < 5
        n = 0.2262; % Depolarization Factor More Needle-Like
    end
    if nargin < 6 || isempty(water_permittivity_model)
        water_permittivity_model = @(f, T) water_permittivity_tiuri80(f, T); % Define this function separately
    end
    if nargin < 7 || isempty(ice_permittivity_model)
        ice_permittivity_model = @(f, T) ice_permittivity_tiuri84(f, T); % Define this function separately
    end

    % SMRT Fellers
    % Permittivity of components
    eps_a = 1; % Background permittivity
    eps_w = water_permittivity_model(frequency, FREEZING_POINT);
    eps_i = ice_permittivity_model(frequency, temperature);

    % Compute alpha and other terms
    alpha = 2.*eps_w+eps_i;
    diff_wi = eps_w - eps_i;
    diff_wa = eps_w - eps_a;
    
    % Depolarization Parameterization
    depolarization = (1-n)./n; % Default n = 0.2262

    denominator = (2 * eps_a + eps_w) .* alpha ...
                - 2 ./ Vw_i .* diff_wa .* diff_wi ...
                - (Vw_i ./ Va_i) .* diff_wa .* alpha ...
                + (1 ./ Va_i) .* diff_wi .* (2 * eps_w + eps_a);

    Es = eps_a .* (1 + depolarization .* ((Vw_i ./ Va_i .* diff_wa .* alpha ...
         - 1 ./ Va_i .* diff_wi .* (2 * eps_w + eps_a)) ./ denominator));

end
