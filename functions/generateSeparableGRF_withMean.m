function field = generateSeparableGRF_withMean(x, z, meanField, sill, rangeX, rangeZ)
% Generates a separable 2D Gaussian random field with a mean trend
%   Cov(x1,z1;x2,z2) = Cx(x1,x2) * Cz(z1,z2)
%   Uses Gaussian (squared exponential) kernel in both directions.
%
% Inputs:
%   x         - [1 x nx] vector of uniform x-coordinates
%   z         - [nz x 1] vector of nonuniform z-coordinates
%   meanField - [nz x nx] matrix of mean values
%   sill      - scalar variance of the field
%   rangeX    - scalar correlation range in x
%   rangeZ    - scalar correlation range in z
%
% Output:
%   field     - [nz x nx] Gaussian random field with mean added

    nx = length(x);
    nz = length(z);

    %% --- Horizontal (x) Gaussian kernel ---
    dx = x(2) - x(1); % must be uniform
    x_shifted = ifftshift((-floor(nx/2):ceil(nx/2)-1) * dx);
    Cx = sill * exp(-(x_shifted.^2) / (2 * rangeX^2));
    Sx = real(fft(Cx));
    Sx(Sx < 0) = 0; % stability fix

    %% --- Vertical (z) Gaussian kernel ---
    Dz = abs(z - z');
    Cz = sill * exp(-(Dz.^2) / (2 * rangeZ^2));
    Lz = chol(Cz + 1e-10*eye(nz), 'lower'); % stabilize

    %% --- Correlated noise generation ---
    W = randn(nz, nx);          % white noise
    Zcorr = Lz * W;             % vertical correlation
    Zcorr_fft = fft(Zcorr, [], 2);         % horizontal FFT
    Zcorr_fft = Zcorr_fft .* sqrt(Sx);     % scale by spectrum
    noise = real(ifft(Zcorr_fft, [], 2));  % inverse FFT

    %% --- Add mean trend ---
    if ~isequal(size(meanField), [nz, nx])
        error('meanField must be of size [nz x nx].');
    end

    field = meanField + noise;

    % Optional: rescale to desired sill (field variance)
    field = field - mean(field(:));
    field = field * sqrt(sill) / std(field(:));
    field = field + mean(meanField(:));
end
