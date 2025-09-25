function blendedStack = blendRadargrams_fft(nearStack, farStack, sigma)
% blendRadargrams_fft - Frequency-domain blending of near and far offset radargrams
%
% Inputs:
%   nearStack - [nt x nx] radargram with high resolution (near offsets)
%   farStack  - [nt x nx] radargram with deep penetration (far offsets)
%   sigma     - (optional) Gaussian rolloff width in frequency space (default: 0.12)
%
% Output:
%   blendedStack - [nt x nx] blended radargram

if nargin < 3
    sigma = 0.12;
end

[nt, nx] = size(nearStack);

% 2D FFT of inputs
FFT_near = fft2(nearStack);
FFT_far  = fft2(farStack);

% Create frequency grids
f_t = (-floor(nt/2):ceil(nt/2)-1)' / nt;
f_x = (-floor(nx/2):ceil(nx/2)-1) / nx;
[F_t, F_x] = meshgrid(f_x, f_t);
R = sqrt(F_t.^2 + F_x.^2); % radial frequency

% Construct Gaussian weighting masks
W_low  = exp(-(R / sigma).^2);      % far (low-pass)
W_high = 1 - W_low;                 % near (high-pass)

% Shift FFTs for blending
FFT_near_shifted = fftshift(FFT_near);
FFT_far_shifted  = fftshift(FFT_far);

% Blend in frequency domain
FFT_blend = W_high .* FFT_near_shifted + W_low .* FFT_far_shifted;

% Inverse FFT to time-space
blendedStack = real(ifft2(ifftshift(FFT_blend)));

end
