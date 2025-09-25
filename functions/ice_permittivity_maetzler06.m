function eps_ice = ice_permittivity_maetzler06(frequency, temperature)
% Calculates the complex ice dielectric constant depending on the frequency and temperature
% Based on Mätzler, C. (2006). Thermal Microwave Radiation: Applications for Remote Sensing, p456-461.
% 
% Inputs:
%   - frequency: Frequency in Hz
%   - temperature: Temperature in K
% 
% Outputs:
%   - eps_ice: Complex permittivity of pure ice
% 
% Raises:
%   - Error if temperature <= 0 K (unphysical input)

    % Define constants
    FREEZING_POINT = 273.15; % Freezing point of water in Kelvin

    % Check for valid temperature
    if temperature <= 0
        error('Temperature must be greater than 0 K.');
    end

    % Convert frequency to GHz
    freqGHz = frequency / 1e9;

    % Real part of the permittivity
    Ereal = 3.1884 + 9.1e-4 * (temperature - FREEZING_POINT);

    % Dimensionless temperature parameter
    theta = 300.0 / temperature - 1.0;

    % Calculate alpha
    alpha = (0.00504 + 0.0062 * theta) * exp(-22.1 * theta);

    % Parameters for beta calculation
    B1 = 0.0207;
    B2 = 1.16e-11;
    b = 335.0;

    % Delta beta
    deltabeta = exp(-9.963 + 0.0372 * (temperature - FREEZING_POINT));

    % Betam calculation
    betam = (B1 / temperature) * (exp(b / temperature) / (exp(b / temperature) - 1)^2) + B2 * freqGHz^2;

    % Beta calculation
    beta = betam + deltabeta;

    % Imaginary part of the permittivity
    Eimag = alpha / freqGHz + beta * freqGHz;

    % Complex permittivity
    eps_ice = Ereal + 1i * Eimag;
end
