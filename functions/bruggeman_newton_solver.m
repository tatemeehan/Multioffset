function epsilon_eff = bruggeman_newton_solver(f, epsilon_i, nvals)
% BRUGGEMAN_NEWTON_SOLVER
% Fast Newton–Raphson solver for Bruggeman symmetric mixing
%   Inputs:
%       f        - volume fractions [f_ice, f_water, f_air]
%       epsilon_i - permittivities [eps_ice, eps_water, eps_air]
%       nvals    - depolarization factors [n_ice, n_water, n_air]
%   Output:
%       epsilon_eff - effective permittivity (complex)

    % Initial guess: weighted average
    epsilon_eff = sum(f .* epsilon_i);

    % Newton–Raphson iteration
    for iter = 1:5  % fixed 5 iterations
        denom = epsilon_eff + nvals .* (epsilon_i - epsilon_eff);
        numer = f .* (epsilon_i - epsilon_eff) ./ denom;
        
        f_val = sum(numer);
        
        % derivative terms
        d_denom = (1 - nvals);
        d_numer = -f .* (1 ./ denom.^2) .* (epsilon_i - epsilon_eff).^2 .* (1 + nvals);
        
        df_val = sum(d_numer);
        
        % Newton update
        epsilon_eff = epsilon_eff - f_val / df_val;
    end

    % enforce positive imaginary part
    epsilon_eff = real(epsilon_eff) + 1i * abs(imag(epsilon_eff));
end