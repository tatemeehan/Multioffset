function eps_ice = ice_permittivity_tiuri84(frequency, temperature)
    % Calculates the complex ice dielectric constant based on Tiuri et al. (1984)
    % The Complex Dielectric Constant of Snow at Microwave Frequencies.
    % IEEE Journal of Oceanic Engineering, vol. 9, no. 5., pp. 377-382.
    %
    % Parameters:
    % frequency: Frequency in Hz
    % temperature: Temperature in Kelvin
    %
    % Returns:
    % eps_ice: Complex permittivity of pure ice

    % Constants
    FREEZING_POINT = 273.15; % Freezing point of water in Kelvin
    DENSITY_OF_ICE = 917; % Density of ice in kg/m^3 (approximate value)

    % Convert temperature to Celsius
    tempC = temperature - FREEZING_POINT;

    % Validate temperature (must be <= freezing point)
    if any(tempC > 0)
        error('The ice temperature must be lower or equal to the freezing point (273.15 K).');
    end

    % Units conversion: Density in g/cm^3
    density_gm3 = DENSITY_OF_ICE * 1e-3;

    % Eq (1) - Real part
    Ereal = 1 + 1.7 * density_gm3 + 0.7 * density_gm3^2;

    % Eq (6) - Imaginary part
    Eimag = 1.59e6 * ...
        (0.52 * density_gm3 + 0.62 * density_gm3^2) * ...
        (frequency^-1 + 1.23e-14 * sqrt(frequency)) * exp(0.036 * tempC);

    % Complex permittivity
    eps_ice = Ereal + 1i * Eimag;
end
