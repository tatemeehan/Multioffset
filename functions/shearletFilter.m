function [data] = shearletFilter(data,t0,dt,numScales,pow)
% shearletFilter.m performs the shearlet transform on the radargram and
% filters the coeffiecient matricies to remove the background, direct wave,
% and random noise. The background is removed by zeroing the lowpass filter
% of the shearlet coefficients. Vertical Noise Spikes are removed by
% zeroing the vertically polarized coefficients and for lower frequencies
% also the adjacent coefficients. The direct wave is removed by combination
% of median subtraction applied to the horizontally polarized coefficients
% and with the use of an exponentially tapered window. The window length is
% determined by the time-zero plus 50 samples. This approach is similar to
% that of Wang and Liu (2016), but perhaps a bit more sophisticated.
%
% Inputs - data - a matrix (should be approximately square)
%        - t0 - a scalar indicating the time of the airwave arrival
%        - dt  - the sample interval
%        - numScales - The number of Shearlet Scales - Default is 5
%        - pow - The power of the exponential window taper - Default is 5
%
% Outputs - data - the filtered data
%
% Tate Meehan, SnowEx 2021

% Get Image Size
[m,n] = size(data);
 % Set Exponential Power of Window Ramp
 if nargin < 5
     pow = 5;
 end
% Set Number of Shearlet Scales
 if nargin < 4
     defNumScales = floor(log2(min([m,n]))-3);
     if defNumScales < 5
         numScales = defNumScales;
     else
         numScales = 5;
     end
 end
 if numScales > 5
     warning(['Supplied parameter numScales is larger than algorithm support. '...
         'Setting numScales to Default.'])
     numScales = 5;
 end
 if nargin < 3
     winIx = 50;
     warning(['Time-zero information not supplied. Using Default window taper.'])
 else
     % Set Taper Window
     winIx = t0./dt+50;
 end

 % Build Shearlet System
 sls = shearletSystem('ImageSize',[m n],'PreserveEnergy',true,'NumScales',numScales);
 % Perform the 2D Shearlet Transform
 c = sheart2(sls,data);
 % Filter the Low Pass Coefficient Matrix
 c(:,:,1) = zeros(m,n);
 % Filter the Vertically and Horizontally Polarized Coefficient Matricies
 v = [3,10,20,31,46]; v = v(1:numScales);   % Vertical   Coefficients
 h = [6,15,25,38,55]; h = h(1:numScales);   % Horizontal Coefficients
 % Set Vertical Coefficient Matricies to Zero
 for kk = 1:numScales
     c(:,:,v(kk)) = zeros(m,n);
     if kk < 4 & kk > 1
         % Zero Adjacent Coefficients
         c(:,:,v(kk)+1) = zeros(m,n);
         c(:,:,v(kk)-1) = zeros(m,n);
     end
 end
 % Subtract Rowwise Median from Horizontal Coefficient Matricies
 for jj = 1:numScales
     M = zeros(m,1);
     for kk = 1:m
         M(kk) = median(c(kk,:,h(jj)));
     end
     % Form the Outter Product
     M = M*ones(1,n);
     % Subtract Median
     c(:,:,h(jj)) = c(:,:,h(jj))-M;
 end
 % Filter All Coefficient Matricies with Tapered Window
 window = linspace(0,1,winIx).^pow';
 for kk = 1:size(c,3)
     c(1:winIx,:,kk) = window.*c(1:winIx,:,kk).*ones(winIx,n);
 end
 % Reconstruct the Data by the Inverse Shearlet Transform
 data = isheart2(sls,c);
end

