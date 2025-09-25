function field = generateSeparableGRF(x, z, sill, rangeX, rangeZ)
% Generates a separable 2D Gaussian random field
%   Cov(x1,z1;x2,z2) = Cx(x1,x2) * Cz(z1,z2)
%
% Inputs:
%   x       - [1 x nx] vector of uniformly spaced x-coordinates
%   z       - [nz x 1] vector of nonuniform z-coordinates
%   sill    - variance (scalar)
%   rangeX  - correlation range in x (same units as x)
%   rangeZ  - correlation range in z (same units as z)
%
% Output:
%   field   - [nz x nx] Gaussian random field

    nx = length(x);
    nz = length(z);

    %% --- Construct horizontal covariance (Cx) ---
    dx = x(2) - x(1); % must be uniform!
    x_shifted = ifftshift((-floor(nx/2):ceil(nx/2)-1) * dx);
    Cx = sill * exp(-abs(x_shifted) / rangeX);
    Sx = real(fft(Cx)); % power spectral density
    Sx(Sx < 0) = 0;      % numerical fix

    %% --- Construct vertical covariance (Cz) ---
    Dz = abs(z - z');   % [nz x nz] pairwise distances
    Cz = sill * exp(-Dz / rangeZ);

    % Cholesky decomposition in z
    Lz = chol(Cz + 1e-10*eye(nz), 'lower'); % stabilize

    %% --- Simulate separable field ---
    % White noise: [nz x nx]
    W = randn(nz, nx);

    % Apply vertical correlation (Cholesky)
    Zcorr = Lz * W;

    % Apply horizontal correlation via FFT
    Zcorr_fft = fft(Zcorr, [], 2);             % FFT along x
    Zcorr_fft = Zcorr_fft .* sqrt(Sx);         % Scale
    field = real(ifft(Zcorr_fft, [], 2));      % Back to spatial domain

    % Normalize (optional)
    field = field - mean(field(:));
    field = field * sqrt(sill) / std(field(:));
end
