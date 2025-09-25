function d_enh = fce_filter(data, nx, nt, nk, dx, dt)
% Fast Coherent Enhancement (FCE) Filter for GPR data
% INPUTS:
%   data - [nt x nx] matrix, GPR data (time x trace)
%   nx   - number of traces used for stacking (odd number recommended)
%   nt   - number of samples in the coherence time window
%   nk   - number of dip scan angles
%   dx   - trace spacing (m)
%   dt   - time step (s)
% OUTPUT:
%   d_enh - enhanced data of same size as input

[NT, NX] = size(data);
half_ap = floor(nx / 2);
dip_angles = linspace(-10, 10, nk); % ns/m

% Preallocate
stacked = zeros(NT, nk);     % For a single trace
C = zeros(NT, nk);           % Coherence scores
d_enh = zeros(size(data));   % Output

for x0 = 1+half_ap : NX-half_ap
    d0 = data(:, x0);  % current trace

    % For each dip angle
    for k = 1:nk
        slope = dip_angles(k); % ns/m
        
        % Stack over nx traces at this dip
        temp = zeros(NT, nx);
        for i = -half_ap : half_ap
            xi = x0 + i;
            time_shift = round(slope * i * dx / dt); % convert ns/m to index
            shifted = circshift(data(:, xi), -time_shift);
            temp(:, i + half_ap + 1) = shifted;
        end
        stacked(:, k) = mean(temp, 2);
        C(:, k) = d0 .* stacked(:, k); % pointwise product
    end

    % Recursive coherence window
    C_sum = movsum(C, [floor(nt/2) floor(nt/2)], 1);
    [C_max, best_k] = max(C_sum, [], 2);

    % Normalize weights
    weights = C_max ./ max(C_max);
    d_enh(:, x0) = data(:, x0) .* weights;
end
end
