% extrapolateImaginaryPermittivity
% Performs column-wise exponential fitting to a mean permittivity profile,
% smooths the resulting (a, b, c) parameters using localized weighted smoothing,
% and reconstructs a 2D extrapolated permittivity field.
% Now includes interface-level rescaling to match real data before noise addition,
% with horizontal smoothing of the real interface to avoid vertical striping.

function Ei = extrapolateImaginaryPermittivity(Ei, bottomIx, extendN, Dvec, thresholdR, GPR, ii)

% Smooth Permittivity Data for Model Fitting
eiMean = zeros(numel(500:bottomIx), numel(Dvec));
sigma = thresholdR;

for kk = 1:numel(Dvec)
    dist = abs(Dvec - Dvec(kk));
    w = exp(-(dist.^2) / (2*sigma^2)); w = w / sum(w);

    % Weighted average profile (smooth horizontal blend)
    eiMean(:, kk) = Ei(500:bottomIx, :) * w';
end

[nz, nx] = size(Ei);
baseZ = (1:(bottomIx - 500 + 1))';
zExtrap = (1:extendN)';

% Step 1: Exponential fitting
a_params = nan(1, nx);
b_params = nan(1, nx);
c_params = nan(1, nx);
for kk = 1:nx
    y = Ei(500:bottomIx, kk);
    z = baseZ;
    if any(isnan(y)) || all(y <= 0), continue; end
    try
        fitObj = fit(z, y, fittype('a*exp(-b*x)+c'), 'StartPoint', [y(1)-y(end), 0.01, y(end)]);
        a_params(kk) = fitObj.a;
        b_params(kk) = fitObj.b;
        c_params(kk) = fitObj.c;
    catch, end
end

% Step 2: Mask and smooth
% Define validity masks
a_mask = isfinite(a_params) & (a_params) < 0.1 & (a_params) > -0.1;   % limit magnitude
b_mask = isfinite(b_params) & b_params < 2.5e-2;  % limit decay rate
c_mask = isfinite(c_params) & (c_params) < 0.1 & (c_params) > -0.1;   % clamp based on your scale

% Combined mask (use only fully valid sets)
valid_mask = a_mask & b_mask & c_mask;

% Replace invalid entries with NaN for smoothing
a_params(~valid_mask) = NaN;
b_params(~valid_mask) = NaN;
c_params(~valid_mask) = NaN;

% Smooth using only valid values
a_smooth = smoothParam(a_params, Dvec, sigma);
b_smooth = smoothParam(b_params, Dvec, sigma);
c_smooth = smoothParam(c_params, Dvec, sigma);

% Step 3: Build extrapolated region
for kk = 1:nx
    if any(isnan([a_smooth(kk), b_smooth(kk), c_smooth(kk)])), continue; end
    Ei(bottomIx+1:bottomIx+extendN, kk) = a_smooth(kk) * exp(-b_smooth(kk) * zExtrap) + c_smooth(kk);
end

% Step 4: Horizontal interface rescaling
real_rows = bottomIx-5:bottomIx;
extrap_rows = bottomIx+1:bottomIx+5;
avg_real = mean(movmean(Ei(real_rows, :), 151, 2), 1);
avg_extrap = mean(Ei(extrap_rows, :), 1);
scale = avg_real ./ avg_extrap; scale(~isfinite(scale)) = 1;
for kk = 1:nx
    Ei(bottomIx+1:end, kk) = Ei(bottomIx+1:end, kk) * scale(kk);
end

% Step 5: Add vertically correlated Gaussian process noise
sill = 2.5; range = 0.25; % 25 cm Vertical Correlation Length
rangeIx = ceil(range / GPR.D.dx{ii}) * 10; % 10 m Spatial Correlation Length
distanceMat = GPR.Geolocation.Distance{ii} .* ones(size(Ei(bottomIx+1:end, :)));
depth = GPR.D.Depth{ii}(bottomIx+1:end, :);
stdEi = std(Ei(bottomIx-50:bottomIx, :), 0, 'all').* 0.85;

for kk = 1:rangeIx:(nx - rangeIx)
    dSub = depth(:, kk:kk+rangeIx-1);
    xSub = distanceMat(:, kk:kk+rangeIx-1);
    D = pdist2([xSub, dSub], [xSub, dSub]);
    C = sill - (sill * (1 - exp(-D ./ range)));
    W = chol(C + 1e-8 * eye(size(C)));

    randEi = normrnd(0, stdEi, size(W, 1), 1);
    tmp = W * randEi;
    tmp = tmp - mean(tmp);
    tmp = repmat(tmp, 1, rangeIx);
    tmp = tmp .* exp(-((1:size(tmp,1))') / 300);

    meanEi = Ei(bottomIx+1:end, kk:kk+rangeIx-1);
    EiOut =  tmp + meanEi;
    EiOut = max(6e-5, EiOut);
    Ei(bottomIx+1:end, kk:kk+rangeIx-1) = EiOut;
end

% Final smoothing
Ei(bottomIx-5:end, :) = imgaussfilt(Ei(bottomIx-5:end, :), 2);
end

function paramSmooth = smoothParam(paramVec, Dvec, sigma)
nx = length(paramVec);
paramSmooth = nan(1, nx);
for kk = 1:nx
    dist = abs(Dvec - Dvec(kk));
    w = exp(-(dist.^2) / (2 * sigma^2));
    w(isnan(paramVec)) = 0; w = w / sum(w);
    paramSmooth(kk) = nansum(paramVec .* w);
end
end
