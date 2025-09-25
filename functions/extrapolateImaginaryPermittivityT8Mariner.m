% extrapolateImaginaryPermittivityT8Mariner.m
% -------------------------------------------------------------------------
% Created by T8 & GPT Mariner
%
% Extrapolates complex dielectric permittivity below the deepest GPR return
% using exponential curve fitting and physically grounded Gaussian noise.
%
% This function models the sub-signal permittivity structure by:
%   • Fitting exponential decay curves to the observed real-valued data
%   • Smoothing model parameters laterally for spatial coherence
%   • Seamlessly extending into unmeasured depths
%   • Adding spatially correlated, depth-decaying Gaussian noise
%
% Developed for forward and inverse modeling of wet snow and firn
% over the Juneau Ice Field.
%
% This implementation underpinned the final retrieval of ε′ + jε″ over the JIF
% using multi-offset GPR wave speed and frequency-time analyses.
% -------------------------------------------------------------------------

function Ei = extrapolateImaginaryPermittivityT8Mariner(Ei, bottomIx, extendN, Dvec, thresholdR, GPR, ii)

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

    % Step 1: Subsampled exponential fitting with parfor
    p = 0.25;
    % fitIx = 1:nx;
    % sampleCols = datasample(fitIx,round(p.*numel(fitIx)),'replace',false);
    % fitSpacing = round(thresholdR / GPR.D.dx{ii});
    % sampleCols = 1:fitSpacing:nx;
    % Set buoy coverage percentage
    p = 0.25;
    nBuoys = round(p * nx);
    fitSpacing = round(nx / nBuoys);
    sampleCols = 1:fitSpacing:nx;  % Regularly spaced buoys based on percentage
    a_params = nan(1, nx);
    b_params = nan(1, nx);
    c_params = nan(1, nx);

    a_sub = nan(size(sampleCols));
    b_sub = nan(size(sampleCols));
    c_sub = nan(size(sampleCols));

    parfor idx = 1:length(sampleCols)
        kk = sampleCols(idx);
        y = Ei(500:bottomIx, kk);
        z = baseZ;
        if any(isnan(y)) || all(y <= 0), continue; end
        try
            fitObj = fit(z, y, fittype('a*exp(-b*x)+c'), 'StartPoint', [y(1)-y(end), 0.01, y(end)]);
            a_sub(idx) = fitObj.a;
            b_sub(idx) = fitObj.b;
            c_sub(idx) = fitObj.c;
        catch, end
    end
    a_params = a_sub; b_params = b_sub; c_params = c_sub;
        validIdx = isfinite(a_sub) & (a_sub) < 0.2 & (a_sub) > -0.2 & ...
        isfinite(b_sub) & b_sub < 5e-2 & ...
        isfinite(c_sub) & (c_sub) < 0.15 & (c_sub) > -0.15;
    a_params(~validIdx) = NaN; b_params(~validIdx) = NaN; c_params(~validIdx) = NaN;

    % Interpolate fit parameters back to full domain
    % validIdx = isfinite(a_sub) & isfinite(b_sub) & isfinite(c_sub);
    a_params = interp1(Dvec(sampleCols(validIdx)), a_params(validIdx), Dvec, 'linear', 'extrap');
    b_params = interp1(Dvec(sampleCols(validIdx)), b_params(validIdx), Dvec, 'linear', 'extrap');
    c_params = interp1(Dvec(sampleCols(validIdx)), c_params(validIdx), Dvec, 'linear', 'extrap');

    % Fill any large gaps from interpolation
    a_params = fillmissing(a_params, 'nearest');
    b_params = fillmissing(b_params, 'nearest');
    c_params = fillmissing(c_params, 'nearest');

    % Mask and smooth
    % T8's trusted tightrope bounds
    % valid = isfinite(a_params) & (a_params) < 0.2 & (a_params) > 0 & ...
    %     isfinite(b_params) & b_params > 5e-2 & ...
    %     isfinite(c_params) & (c_params) < 0.15 & (c_params) > 0;

    sigma = thresholdR / 2;
    a_smooth = smoothParam(a_params, Dvec, sigma);
    b_smooth = smoothParam(b_params, Dvec, sigma);
    c_smooth = smoothParam(c_params, Dvec, sigma);

    % Reconstruct extrapolated region
    for kk = 1:nx
        if any(isnan([a_smooth(kk), b_smooth(kk), c_smooth(kk)])), continue; end
        Ei(bottomIx+1:bottomIx+extendN, kk) = a_smooth(kk) * exp(-b_smooth(kk) * zExtrap) + c_smooth(kk);
    end

    % Rescale interface
    real_rows = bottomIx-5:bottomIx;
    extrap_rows = bottomIx+1:bottomIx+5;
    avg_real = mean(movmean(Ei(real_rows, :), 25, 2), 1);
    avg_extrap = mean(Ei(extrap_rows, :), 1);
    scale = avg_real ./ avg_extrap; scale(~isfinite(scale)) = 1;
    for kk = 1:nx
        Ei(bottomIx+1:end, kk) = Ei(bottomIx+1:end, kk) * scale(kk);
    end

    % Smoothing Of Extrpolation here
    Ei(bottomIx+1:end, :) = imgaussfilt(Ei(bottomIx+1:end, :), [5 51]);

    % Add depth-correlated Gaussian process noise
    sill = 2.5; range = 0.25;
    rangeIx = ceil(range / GPR.D.dx{ii}) * 10;
    distanceMat = GPR.Geolocation.Distance{ii} .* ones(size(Ei(bottomIx+1:end, :)));
    depth = GPR.D.Depth{ii}(bottomIx+1:end, :);
    stdEi = std(Ei(bottomIx-50:bottomIx, :), 0, 'all') * 0.95;

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
        if ~isequal(size(tmp), size(meanEi))
            tmp = tmp(:, 1:size(meanEi,2));
        end
        EiOut = max(6e-5, tmp + meanEi);
        Ei(bottomIx+1:end, kk:kk+rangeIx-1) = EiOut;
    end
    
    % Continuity
    for kk = 1:nx
    Ei(bottomIx-5:bottomIx+1, kk) = interp1([[bottomIx-11:bottomIx-6],[bottomIx+1:bottomIx+6]],[Ei(bottomIx-11:bottomIx-6,kk);Ei(bottomIx+2:bottomIx+7,kk)],bottomIx-5:bottomIx+1,'pchip');
    end
     % Final smoothing
    Ei(bottomIx-5:end, :) = imgaussfilt(Ei(bottomIx-5:end, :), 2);
    Ei(bottomIx-15:bottomIx-5, :) = imgaussfilt(Ei(bottomIx-15:bottomIx-5, :), 1);
    Ei(bottomIx-25:bottomIx-15, :) = imgaussfilt(Ei(bottomIx-25:bottomIx-15, :), 0.5);
    Ei(bottomIx-35:bottomIx-25, :) = imgaussfilt(Ei(bottomIx-35:bottomIx-25, :), 0.25);
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
