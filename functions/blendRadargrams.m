function [blendedStack, nearW, farW] = blendRadargrams(nearStack, farStack, twtOverlap, dt, varargin)
% blendRadargrams - Smoothly blends two radargram stacks (e.g., near- and far-offset)
%
% Inputs:
%   nearStack   - [nt x nx] matrix, radargram from near offsets
%   farStack    - [nt x nx] matrix, radargram from far offsets
%   twtOverlap  - Overlap depth in nanoseconds (ns) for start of blending
%   dt          - Sample interval in nanoseconds
%
% Optional:
%   'WindowLength' - Number of samples used for transition blending (default: 200)
%   'CorrectBias'  - Logical flag, apply mean correction between stacks (default: true)
%
% Outputs:
%   blendedStack - Final blended radargram [nt x nx]
%   nearW        - Weighting function applied to nearStack
%   farW         - Weighting function applied to farStack

% Parse inputs
p = inputParser;
addParameter(p, 'WindowLength', 200, @isscalar);
addParameter(p, 'CorrectBias', true, @islogical);
parse(p, varargin{:});
nOverlap = p.Results.WindowLength;
correctBias = p.Results.CorrectBias;

[nt, nx] = size(nearStack);

% Bias correction (optional)
if correctBias
    trans_ix = round(twtOverlap/dt) + (-floor(nOverlap/4):floor(nOverlap/4));
    trans_ix = trans_ix(trans_ix > 0 & trans_ix <= nt); % bound indices
    nearMean = mean(mean((nearStack(trans_ix,:)), 2));
    farMean  = mean(mean((farStack(trans_ix,:)), 2));
    ampCorr  = mean(farMean - nearMean);
    farStack = farStack - ampCorr;
end

% Construct Tukey window
blendWin = tukeywin(nOverlap, 0.8);
nearWin = [ones(round(twtOverlap/dt),1); blendWin(floor(nOverlap/2)+1:end)];
farWin  = [zeros(round(twtOverlap/dt),1); blendWin(1:floor(nOverlap/2))];

% Ensure length
nearW = zeros(nt,1);
farW  = ones(nt,1);
nearW(1:numel(nearWin)) = nearWin;
farW(1:numel(farWin))   = farWin;

% Normalize weights to prevent gain
weightSum = nearW + farW;
% nearW = nearW./weightSum;
% farW = farW./weightSum;
% weightSum = nearW + farW;

weightSum(weightSum == 0) = 1; % prevent divide by zero

% Apply blending
blendedStack = (nearStack .* nearW + farStack .* farW) ./ weightSum;

end
