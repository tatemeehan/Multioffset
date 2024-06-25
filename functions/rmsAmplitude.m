function [rmsRad] = rmsAmplitude(Rad,L)
% rmsAmplitude.m  Computes the Root-Mean-Square Amplitude of the Radargram

% Tate Meehan - Juneau Ice Field Research Program, June 2024
[m,n] = size(Rad);
rmsRad = zeros(m,n);
if nargin < 2
    L = 50;
end
L = floor(L./2);
w = hamming(2.*L+1);
for kk = 1:n
    trc = Rad(:,kk).^2;
    rms = sqrt(conv(trc,w,'same')./sum(w));
    rmsRad(:,kk) = rms;
end