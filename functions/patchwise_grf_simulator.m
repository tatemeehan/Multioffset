% PATCHWISE_GRF_SIMULATOR
% Simulates a large 2D Gaussian Random Field (GRF) in memory-efficient patches
% with overlap and blending using anisotropic Gaussian kernels.
% Author: ChatGPT Advanced Version

function field = patchwise_grf_simulator(X, Z, meanField, sill, rangeX, rangeZ, patchSize, overlap)
    % Inputs:
    %   X, Z         - [nz x nx] coordinate matrices
    %   meanField    - mean trend field [nz x nx]
    %   sill         - variance of GRF
    %   rangeX, rangeZ - correlation lengths in x and z
    %   patchSize    - [pz, px] size of each patch (rows, cols)
    %   overlap      - [oz, ox] number of overlapping rows/cols between patches
    %
    % Output:
    %   field        - full simulated GRF [nz x nx]

    [nz, nx] = size(X);
    pz = patchSize(1); px = patchSize(2);
    oz = overlap(1); ox = overlap(2);

    field = zeros(nz, nx);
    weight = zeros(nz, nx);

    % Loop over patches
    for i = 1:(px - ox):(nx - 1)
        for j = 1:(pz - oz):(nz - 1)
            % Patch index limits (including overlap)
            i1 = i;
            i2 = min(i + px - 1, nx);
            j1 = j;
            j2 = min(j + pz - 1, nz);

            Xp = X(j1:j2, i1:i2);
            Zp = Z(j1:j2, i1:i2);
            Mp = meanField(j1:j2, i1:i2);
            [nzi, nxi] = size(Xp);

            % Flatten coordinates for local covariance
            coords = [Xp(:), Zp(:)];
            Dx = pdist2(coords(:,1), coords(:,1));
            Dz = pdist2(coords(:,2), coords(:,2));
            D = sqrt((Dx./rangeX).^2 + (Dz./rangeZ).^2);
            C = sill * exp(-0.5 * D.^2);

            % Cholesky and simulate
            try
                L = chol(C + 1e-10*eye(size(C)), 'lower');
            catch
                L = chol(C + 1e-6*eye(size(C)), 'lower');
            end
            z = randn(nxi*nzi, 1);
            patch = reshape(L*z, nzi, nxi);

            % Create blending weights (2D cosine window)
            wx = cos(linspace(-pi/2, pi/2, nxi)).^2;
            wz = cos(linspace(-pi/2, pi/2, nzi)).^2;
            W = wz' * wx;

            % Update output with weighted average
            field(j1:j2, i1:i2) = field(j1:j2, i1:i2) + W .* (patch + Mp);
            weight(j1:j2, i1:i2) = weight(j1:j2, i1:i2) + W;
        end
    end

    % Normalize weighted output
    field = field ./ max(weight, 1e-12);
end
