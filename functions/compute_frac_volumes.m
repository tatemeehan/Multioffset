function [frac_volume, fi, fw] = compute_frac_volumes(density, liquid_water)
    % Computes the fractional volume of ice+water, the fractional volume of ice,
    % and the fractional volume of water based on the snow density and the
    % liquid water volume fraction (relative to ice + water).
    %
    % Inputs:
    %   density     - Density of the snow, including the ice and water phases (kg/m^3).
    %   liquid_water - Fractional volume of liquid water with respect to the ice+water mixture (0 <= liquid_water <= 1).
    %
    % Outputs:
    %   frac_volume - Total fractional volume of ice and water (ice + water).
    %   fi          - Fractional volume of ice.
    %   fw          - Fractional volume of water.

    % Constants
    DENSITY_OF_ICE = 917;  % Density of ice in kg/m^3
    DENSITY_OF_WATER = 1000;  % Density of liquid water in kg/m^3

    % Calculate the density of the mixture
    density_melange = DENSITY_OF_ICE * (1 - liquid_water) + DENSITY_OF_WATER * liquid_water;

    % Calculate the total fractional volume
    frac_volume = density ./ density_melange;

    % Calculate the fractional volume of ice
    fi = frac_volume .* (1 - liquid_water);

    % Calculate the fractional volume of water
    fw = frac_volume .* liquid_water;

    % Return the fractional volumes
    return
end
