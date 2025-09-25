function Ei = blendCosineVertical(Ei, bottomIx, N)
    % Applies a vertical cosine blend across bottomIx +/- N
    % Inputs:
    %   Ei        - [nz x nx] matrix with extrapolated values already applied
    %   bottomIx  - row index at transition
    %   N         - number of rows above/below to blend
    % Output:
    %   Ei        - smoothed across transition zone

    [nz, nx] = size(Ei);
    z_blend = (bottomIx - N):(bottomIx + N);

    % Make sure rows are in valid range
    z_blend = z_blend(z_blend >= 1 & z_blend <= nz);
    w = cos(pi * ((z_blend - (bottomIx - N)) / (2*N))).^2;

    % Reference row (midpoint of transition)
    refRow = bottomIx;

    for i = 1:length(z_blend)
        row = z_blend(i);
        alpha = w(i);              % taper weight
        Ei(row, :) = alpha * Ei(row, :) + (1 - alpha) * Ei(refRow, :);
    end
end
