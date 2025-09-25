function field = generateGRF_Nystrom(Xgrid, Zgrid, meanField, sill, rangeX, rangeZ, nLandmarks)
% Generate 2D GRF using Nyström approximation (for large grids)
%
% Inputs:
%   Xgrid, Zgrid - [nz x nx] coordinate matrices
%   meanField    - [nz x nx] mean trend field
%   sill         - scalar variance
%   range        - correlation range (Gaussian kernel)
%   nLandmarks   - number of Nyström sample points (e.g., 500–1000)
%
% Output:
%   field        - [nz x nx] random field with mean

    [nz, nx] = size(Xgrid);
    N = nz * nx;

    % Flatten coordinates
    Xvec = Xgrid(:);
    Zvec = Zgrid(:);
    coords = [Xvec, Zvec];

    %% Step 1: Randomly select landmark points
    rng(1); % for reproducibility
    perm = randperm(N);
    idx_landmarks = perm(1:nLandmarks);
    landmarks = coords(idx_landmarks, :);

    %% Step 2: Compute small kernel K (m x m) and cross-covariance W (N x m)
    % D_K = pdist2(landmarks, landmarks);    % small kernel
    % K = sill * exp(-D_K.^2 / (2 * range^2)) + 1e-8*eye(nLandmarks); % stabilize
    Dx_K = pdist2(landmarks(:,1), landmarks(:,1));
    Dz_K = pdist2(landmarks(:,2), landmarks(:,2));
    D_K = sqrt((Dx_K./rangeX).^2 + (Dz_K./rangeZ).^2);

    % D_W = pdist2(coords, landmarks);       % cross-distance
    % W = sill * exp(-D_W.^2 / (2 * range^2));
    Dx_W = pdist2(coords(:,1), landmarks(:,1));
    Dz_W = pdist2(coords(:,2), landmarks(:,2));
    D_W = sqrt((Dx_W./rangeX).^2 + (Dz_W./rangeZ).^2);

    K = sill * exp(-0.5 * D_K.^2) + 1e-8*eye(nLandmarks);
    W = sill * exp(-0.5 * D_W.^2);

    %% Step 3: Eigendecompose K
    [Uk, Sk] = eig(K);
    Lk = diag(Sk);
    Lk(Lk < 0) = 0; % numerical fix
    sqrtInvK = Uk * diag(1 ./ sqrt(Lk)) * Uk';

    %% Step 4: Simulate
    Z = randn(nLandmarks, 1);
    sample = W * sqrtInvK * Z;

    %% Step 5: Reshape and add mean
    noise = reshape(sample, nz, nx);
    field = meanField + noise;

    % Optional: Normalize
    % field = field - mean(field(:));
    % field = field * sqrt(sill) / std(field(:));
    % field = field + mean(meanField(:));
end
