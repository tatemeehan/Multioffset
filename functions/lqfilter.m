function [filtData,badix] = lqfilter(data,X,Y,r,q,qglobal)
% lqfilter.m is a local spatial quantile filter. The desired outlier
% quantiles q = [qmin, qmax] determines a local threshold for data
% rejection within a search radius of r. The locally rejected data is
% gathered into a global bin for which qthresh is applied to determine if
% the local outliers are also global outliers. The globaly rejected
% are then interpolated between by linear piecewise functions with linear
% extrapolation.
%
% Inputs - data vector
%        - X x-coordinate vector
%        - Y y-coordinate vector
%        - r search radius (determined from X and Y)
%        - q local quantiles of outliers default q = [0.05, 0.95];
%        - qglobal global quantile for outlier rejection default 0.95
%
% Outputs - filtData the filtered Data vector
%         - badix the indicies of the filtered data points
% Tate Meehan SnowEx 2021

% Check inputs
if nargin < 6
    qglobal = 0.95;
end
if nargin < 5
    q = [0.05,0.95];
end
% Ensure Vectors
data = data(:); X = X(:); Y = Y(:);
n = length(data);
% Check if Parallel Pool is active
isPar = ~isempty(gcp('nocreate'));
% Allocation
badix = [];
countBadIx = zeros(1,n);
% Local Quantile Filter
if isPar
    parfor kk = 1:n
        dist = sqrt((X-X(kk)).^2+(Y-Y(kk)).^2);
        binIx = find(dist<r);
        pickbin = data(binIx);
        qt = quantile(pickbin,q);
        tmpbadix = find(pickbin<qt(1) | pickbin>qt(2));
        tmpbadix = binIx(tmpbadix);
        badix = [tmpbadix;badix];
    end
    parfor kk = 1:length(data)
        countBadIx(kk) = sum(badix == kk);
    end
else
    for kk = 1:n
        dist = sqrt((X-X(kk)).^2+(Y-Y(kk)).^2);
        binIx = find(dist<r);
        pickbin = data(binIx);
        qt = quantile(pickbin,q);
        tmpbadix = find(pickbin<qt(1) | pickbin>qt(2));
        tmpbadix = binIx(tmpbadix);
        badix = [tmpbadix;badix];
    end
    for kk = 1:n
        countBadIx(kk) = sum(badix == kk);
    end
end
% Apply Global Outlier Threshold
threshold = quantile(countBadIx,qglobal);
rmvIx = find(countBadIx>threshold);
goodix = find(~ismember(1:n,rmvIx));
% Filter Outliers with Linear Interpolation
filtData = interp1(goodix,data(goodix),1:n,'pchip','extrap');
end

