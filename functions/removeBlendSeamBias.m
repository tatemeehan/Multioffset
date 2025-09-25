function blendedStack = removeBlendSeamBias(blendedStack, twtOverlap, dt, window)
% removeBlendSeamBias - Corrects lateral bias artifact at blend transition
%
% Inputs:
%   blendedStack - [nt x nx] radargram after blending
%   twtOverlap   - Blend start time (ns)
%   dt           - Time step (ns)
%   window       - Half-width of sample window around seam (default: 10)
%
% Output:
%   blendedStack - Bias-corrected stack

if nargin < 4
    window = 10;
end

ix0 = round(twtOverlap / dt);
[nt, nx] = size(blendedStack);

% Ensure we don't exceed bounds
ix1 = max(ix0 - window, 1);
ix2 = min(ix0 + window, nt);

% Compute vertical profile difference just above and below seam
above = mean(blendedStack(ix1:ix0-1,:), 1);
below = mean(blendedStack(ix0+1:ix2,:), 1);

% Estimate column-wise seam bias
bias = below - above;

% Subtract seam bias from lower half
for j = 1:nx
    blendedStack(ix0+1:end, j) = blendedStack(ix0+1:end, j) - bias(j);
end

end