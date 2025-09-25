% EXPONENTIALEXTRAPSMOOTH2D
% Performs column-wise exponential fitting to a mean permittivity profile,
% smooths the resulting (a, b, c) parameters using localized weighted smoothing,
% and reconstructs a 2D extrapolated permittivity field.

function Ei = exponentialExtrapSmooth2D(Ei, bottomIx, extendN, Dvec, method, thresholdR)
% Inputs:
%   Ei        - [nz x nx] permittivity field (bottomIx+1:end rows are extrapolated)
%   bottomIx  - row index where extrapolation begins
%   extendN   - number of rows to extrapolate downward
%   Dvec      - [1 x nx] horizontal distances
%   method    - 'exp' or similar (currently only 'exp' supported)
%   thresholdR - horizontal distance for parameter smoothing
eiMean = zeros(numel(500:bottomIx), numel(Dvec));
sigma = thresholdR;

for kk = 1:numel(Dvec)
    dist = abs(Dvec - Dvec(kk));
    w = exp(-(dist.^2) / (2*sigma^2)); w = w / sum(w);

    % Weighted average profile (smooth horizontal blend)
    eiMean(:, kk) = Ei(500:bottomIx, :) * w';
end

[nz, nx] = size(Ei);
baseZ = (1:(bottomIx - 500 + 1))'; % depth indices of known portion
zExtrap = (1:extendN)';

a_params = nan(1, nx);
b_params = nan(1, nx);
c_params = nan(1, nx);

% Step 1: Fit exponential to each column
for kk = 1:nx
    % y = Ei(500:bottomIx, kk);
    y = eiMean(:,kk);
    z = baseZ;

    % Basic check for valid data
    if any(isnan(y)) || all(y <= 0)
        continue;
    end

    % Fit: y = a * exp(-b*z) + c
    expModel = fittype('a*exp(-b*x)+c', 'independent', 'x', 'coefficients', {'a','b','c'});
    try
        fitObj = fit(z, y, expModel, 'StartPoint', [y(1)-y(end), 0.01, y(end)]);
        a_params(kk) = fitObj.a;
        b_params(kk) = fitObj.b;
        c_params(kk) = fitObj.c;
    catch
        % Leave NaN if fit fails
    end
end
% figure();
% subplot(3,1,1); plot(a_params); title('a'); ylabel('a'); grid on;
% subplot(3,1,2); plot(b_params); title('b'); ylabel('b'); grid on;
% subplot(3,1,3); plot(c_params); title('c'); ylabel('c'); grid on;

% Step 2: Smooth parameters with localized Gaussian weights
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
% 
% figure();
% subplot(3,1,1); plot(a_smooth); title('a'); ylabel('a'); grid on;
% subplot(3,1,2); plot(b_smooth); title('b'); ylabel('b'); grid on;
% subplot(3,1,3); plot(c_smooth); title('c'); ylabel('c'); grid on;
% % a_smooth = smoothParam(a_params, Dvec, sigma);
% % b_smooth = smoothParam(b_params, Dvec, sigma);
% % c_smooth = smoothParam(c_params, Dvec, sigma);
% figure();
% plot([a_smooth(:), b_smooth(:), c_smooth(:)]);
% legend('a','b','c'); title('Smoothed Parameters');

% Step 3: Reconstruct extrapolated values from smoothed parameters
for kk = 1:nx
    a = a_smooth(kk);
    b = b_smooth(kk);
    c = c_smooth(kk);

    if any(isnan([a, b, c]))
        continue;
    end

    Ei(bottomIx+1:bottomIx+extendN, kk) = a * exp(-b * zExtrap) + c;
    Ei(bottomIx+1:bottomIx+extendN, kk) = max(Ei(bottomIx+1:bottomIx+extendN, kk),6e-5);
end

    % % Step 5: Interface-level rescaling before noise
    % real_rows = bottomIx-5:bottomIx;
    % extrap_rows = bottomIx+1:bottomIx+5;
    % Ei_real_smoothed = movmean(Ei(real_rows, :), 25, 2); % smooth across x with window of 25 columns
    % avg_real = mean(Ei_real_smoothed, 1);
    % avg_extrap = mean(Ei(extrap_rows, :), 1);
    % % 
    % % avg_real = mean(Ei(real_rows, :), 1);
    % % avg_extrap = mean(Ei(extrap_rows, :), 1);
    % 
    % scale = avg_real ./ avg_extrap;
    % scale(~isfinite(scale)) = 1;
    % 
    % for kk = 1:nx
    %     Ei(bottomIx+1:end, kk) = Ei(bottomIx+1:end, kk) * scale(kk);
    % end
end

function paramSmooth = smoothParam(paramVec, Dvec, sigma)
nx = length(paramVec);
paramSmooth = nan(1, nx);
for kk = 1:nx
    dist = abs(Dvec - Dvec(kk));
    w = exp(-(dist.^2) / (2*sigma^2));
    w(isnan(paramVec)) = 0; % ignore missing fits
    w = w / sum(w);
    paramSmooth(kk) = nansum(paramVec .* w);
end
end
