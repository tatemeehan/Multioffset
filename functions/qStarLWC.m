function [GPR] = qStarLWC(GPR,f0,fmethod,qmethod,chan,isPlot)
% qStar.m computes the instantaneous frequency of the GPR image and 
% estimates the frequency downshift from f0, the instantaneous frequency of a
% reference airwave signal. The reference and measured frequency is supplied
% into Bradford's 2007 method for estimation of Q* attentuation. Q* is used
% to estimate the imaginary component of snow dielectric permittivity,
% supplied the real part. Snow/firn liquid water is estimated from the
% imaginary component of the permittivity following Shivola & Tiuri (1984).

% This function requires Stacking and Interval Velocity Models.
%
% Inputs: GPR - the GPR data Structure
%          f0 - apriori air wave instantaneous frequency
%     fmethod - Method for picking instant freq of ground ('avg' or default: 'max')
%     qmethod - Q* Method ('Bradford07' or default: 'Bradford09')
%        chan - Channel or Offset number used in LWC estimation (default 1)
%      isPlot - Plot the Instantaneous Frequency and Downshift
%
% Outputs: GPR    - The GPR data Structure
%          Output Variables:
%           GPR.D.stackingImagPermittivity - Imaginary Permittivity
%           relative to Stacking Velocity
%           GPR.D.intervalImagPermittivity - Imaginary Permittivity
%           relative to Interval Velocity
%           GPR.D.qLWC - Snow/Firn LWC Estimated from Q* Attenuation
%           * Interval Velocity is used to estimate LWC
%
% Tate Meehan - CRREL, Juneau IceField Research Project, 3/31/2025

%% Default Parameters
if nargin < 3
    fmethod = 'avg';
    qmethod = 'Bradford07';
    isPlot = 0;
end
if nargin < 4
    qmethod = 'Bradford07';
    isPlot = 0;
end
if nargin < 5
    % Select Channel (near offset)
    chan = 1;
    isPlot = 0;
end
if nargin < 6;
        isPlot = 0;
end

%% Forward Model
% Tau For Bruggman polarizations n = [1/3,0.05,0.1];
% 1.3 GHz
% tauParams = [2.70604364852145e-10, -8.95922225296886e-13, 1.10036955414893e-15 ,-4.79169239741482e-19];
% tauFun = @(p,rho) p(1) + p(2).*rho + p(3).*rho.^2 + p(4).*rho.^3;
% 50 MHz
% tauParams = [1.12500000001000e-10, -1.23881816659395e-11, 0.01, 275, ...
%     -8.28312725292080e-12, 0.01, 475, -2.02719661105318e-13, 1.12011837199496e-16];
% 0.5 GHz
tauParams = [1.12500000000915e-10, -1.93658735175413e-11, 0.01, 275,...
    -6.38731614738757e-12, 0.01, 475, -1.27852360690599e-13, 4.32973804133215e-17];
tauFun = @(p, rho) ...
    p(1) + ...
    p(2) * tanh(p(3)*(rho - p(4))) + ...
    p(5) * tanh(p(6)*(rho - p(7))) + ...
    p(8)*rho + p(9)*rho.^2;

% Model Parameters
frequency = (f0(chan)./1.13).*10.^9;  % 1.3 GHz
temperature = 273.15; % K

% Define density and LWC ranges
rho_range = linspace(100, 917, 1000);   % Density from 200 to 550 kg/m^3
lwc_range = linspace(0, 0.35, numel(rho_range));    % LWC from 0 to 15% volume fraction
[RHO, LWC] = meshgrid(rho_range, lwc_range);

% Initialize error grid
eps_model = zeros(size(RHO));

% Compute Permittivity Forward Model
% for idx = 1:numel(RHO)
parfor kk = 1:length(rho_range)
    tau = tauFun(tauParams,rho_range(kk));      % Relaxation Time

    params = struct( ...
        'eps_s', 87.5, ...
        'eps_inf', 4.9, ...
        'tau', tau, ...
        'alpha', 0, ...
        'sigma', 0.001);

    eps_w = water_permittivity_colecole(frequency, params);

    eps_model(:,kk) = [wetsnow_permittivity_tinga73_extended(frequency, temperature, rho_range(kk), lwc_range, ...
        @ice_permittivity_maetzler06, eps_w)].';
    % eps_model(:,kk) = wetsnow_permittivity_tinga73(frequency, temperature, rho_range(kk), lwc_range, 0.2262, @ice_permittivity_maetzler06, eps_w)
end

%% Interpolation Look-up Approach for Inversion
% Reshape inputs
eps_real = real(eps_model(:));
eps_imag = imag(eps_model(:));
rho_vals = RHO(:);
lwc_vals = LWC(:);

% Ensure uniqueness by grouping duplicate ε*, and averaging ρ/LWC
[eps_unique, ~, idxu] = unique([eps_real, eps_imag], 'rows');
rho_mean = accumarray(idxu, rho_vals, [], @mean);
lwc_mean = accumarray(idxu, lwc_vals, [], @mean);

% Build scattered interpolants
rhoInvInterp = scatteredInterpolant(eps_unique(:,1), eps_unique(:,2), rho_mean, 'natural', 'none');
lwcInvInterp = scatteredInterpolant(eps_unique(:,1), eps_unique(:,2), lwc_mean, 'natural', 'none');

%% Compute Instantaneous Frequency
for ii = 1:GPR.MD.nFiles
% HH Polarization Signal
% Raw Signal
% sHH = GPR.D.MxRadar{1}(:,1:2:end);
% Processed Signal
signal = GPR.D.RadarNMO{ii}{chan};

% Sample Interval
dt = GPR.D.dt{1};
dx = GPR.D.dx{1};
% Number of samples for Smoothing
n = 1;

% Instantaneous Frequency - Barnes Eqn 3
% method ... 1 Unwrap the instantaneous phase and differentiate (Barnes eqn 2)
%            2 Differentiate real and imaginary parts (Barnes Eqn 3)
%            3 Claerbout's method (Barnes eqn 14)
%            4 Barnes first approximation (Barnes eqn 9)
%            5 Barnes second approximation (Barnes eqn 12)
method = 2;

% Absolute Value of Frequency - 1
% flag ... 1 Calculate the absolute value of the instantaneous frequency
%      ... 2 Calculate the signed instantaneous frequency
flag = 1;
% Median Smoothing Window
filtm = round(1./dt); % 1 ns
filtn = round(5./dx); % 5 m
if mod(filtm,2)==0
    filtm = filtm+1;
end
if mod(filtn,2)==0
    filtn = filtn+1;
end
instantFreq = medfilt2(ins_freq(signal,dt,n,method,flag),[filtm,filtn],'symmetric');
freq = instantFreq;
fMean = mean(freq(6:10,:));
freq(1:5,:) = ones(5,1)*fMean;
%% Estimate Q*
% Scalar Peak Frequency For Ricker Wavelet
f0 = f0(chan);
f0rick = f0./1.13;
freq = freq./1.13;

% Angular Frequency
w0 = 2.*pi.*f0rick; %w0 = w0(:);
wt = 2.*pi.*(freq); %wt = wt(:);
df = w0-wt;

% Bradford Retreival Method
if strcmp(qmethod,'Bradford07')
    % Bradford 2007
    qWeight = ones(size(freq)).*0.25; qWeight(20:end,:) = 4;
    H = hamming(21);
    qWeight = conv2(qWeight,H./sum(H),"same");
    invQstar = (qWeight./(GPR.D.TimeAxis{ii}(:))).*((w0.^2-wt.^2)./(w0.^2.*wt));
elseif strcmp(qmethod,'Bradford09')
    % Bradford 2009
    qWeight = ones(size(freq)); qWeight(10:end,:) = 2;
    H = hamming(21);
    qWeight = conv2(qWeight,H./sum(H),"same");
    invQstar = (qWeight./(pi.*GPR.D.TimeAxis{ii}(:))).*((w0.^2-wt.^2)./(w0.^2.*wt));
else
    qWeight = ones(size(freq)); qWeight(10:end,:) = 3;
    H = hamming(21);
    qWeight = conv2(qWeight,H./sum(H),"same");
    invQstar = (qWeight./(GPR.D.TimeAxis{ii}(:))).*((w0.^2-wt.^2)./(w0.^2.*wt));
end

% Median Filter
invQstar = nanmedfilt2(invQstar,[filtm, filtn]);
% invQstar(1:25,:) = imgaussfilt(invQstar(1:25,:),1);
% Q*
Qstar = abs(medfilt2(1./invQstar,[filtm, filtn],'symmetric'));
Qstar = imgaussfilt(Qstar,1);
% Qstar = abs(nanmedfilt2(1./invQstar,[filtm, filtn]));
% % Q* Mean and Standard Deviation
% QstarStd = movstd(Qstar,L);
% % SG
% % Qstar = movmean(Qstar,L./2);
% O = 5; % polynomial order
% Qstar= sgolayfilt(Qstar,O,L); % m

%% Density Liquid Water Content Estimation
% Store Imaginary Component Estimate
GPR.D.stackingImagPermittivity{ii} = GPR.D.stackingPermittivity{ii}./(2.*Qstar);
GPR.D.intervalImagPermittivity{ii} = GPR.D.intervalPermittivity{ii}./(2.*Qstar);
NaNix = find(GPR.D.intervalImagPermittivity{ii} < 0.001 |...
    GPR.D.intervalImagPermittivity{ii} > max(eps_imag(:)));
GPR.D.intervalImagPermittivity{ii}(NaNix) = NaN;
bottomIx = 1200;
extendN = length(GPR.D.TimeAxis{ii})-bottomIx;
GPR.D.intervalImagPermittivity{ii}(1200:end,:) = NaN;
% DryIx = find(GPR.D.intervalImagPermittivity{ii} < 0.001);
% GPR.D.intervalImagPermittivity{ii}(DryIx) = 0.001;
% WetIx = find(GPR.D.intervalImagPermittivity{ii} > max(eps_imag(:)));
% GPR.D.intervalImagPermittivity{ii}(WetIx) = max(eps_imag(:));

% Density and LWC Inversion Using Interpolation Look-up
% Invert For Retrieved Complex Permittivity
GPR.D.intervalImagPermittivity{ii} = inpaint_nans(GPR.D.intervalImagPermittivity{ii},5);
WetIx = find(GPR.D.intervalImagPermittivity{ii} > max(eps_imag(:)));
GPR.D.intervalImagPermittivity{ii}(WetIx) = max(eps_imag(:));
DryIx = find(GPR.D.intervalImagPermittivity{ii} < min(eps_imag(:)));
GPR.D.intervalImagPermittivity{ii}(DryIx) = min(eps_imag(:));
Er = GPR.D.intervalPermittivity{ii};
Ei = GPR.D.intervalImagPermittivity{ii};

% Extrapolate Imaginary Permittivity
thresholdR = 250;
Dvec = GPR.Geolocation.Distance{ii};  % 1 x N
% Ei = extrapolateImaginaryPermittivity(Ei, bottomIx, extendN, Dvec, thresholdR, GPR, ii);
Ei = extrapolateImaginaryPermittivityT8Mariner(Ei, bottomIx, extendN, Dvec, thresholdR, GPR, ii);

% Density and LWC Inversion Tinga & Bruggeman Method
rho = rhoInvInterp(Er, Ei);
lwc = lwcInvInterp(Er, Ei);
rho = inpaint_nans(rho,5);
lwc = inpaint_nans(lwc,5);
lwc(lwc<0) = 0;
GPR.D.qDensity{ii} = rho;
GPR.D.qLWC{ii} = lwc;
% Add Water Weight to Snow Density
GPR.D.qDensity{ii} = GPR.D.qDensity{ii}+GPR.D.qLWC{ii}.*1000;
GPR.D.qDensity{ii}(GPR.D.qDensity{ii}>917) = 917; % Literally 1 data point

% % LWC inversion Tiuri Method
% Tk = 273.15;
% kwimag = imag(water_permittivity_maetzler87(GPR.D.f0GHz{ii}.*10.^9, Tk));
% % kimag = GPR.D.stackingImagPermittivity{ii};
% % Use Instantaneous Velocity for Instantaneous Frquency
% kimag = GPR.D.intervalImagPermittivity{ii};
% % Minimization Scheme
% % Estimate LWC
% testlwc = [0:0.001:0.25];
% % Shivola & Tiuri (1984) Inversion
% lwc = zeros(size(kimag));
% for kk = 1:numel(kimag)
%     testLWC = [(0.1.*testlwc+0.8.*testlwc.^2).*kwimag - kimag(kk)]';
%     [~,lwcIx] = min(abs(testLWC));
%     lwc(kk) = testlwc(lwcIx);
% end
% 
% % Set LWC Thresholds
% lwc(lwc>=0.25) = NaN;
% lwc(lwc<0) = NaN;
% % lwc(1:10,:) = NaN;
% lwc(1200:end,:) = NaN;
% lwc = inpaint_nans(lwc,5);
% % H = hamming(11);
% % lwcsmooth = conv2(lwc,H./sum(H),"same");
% % GPR.D.qLWC{ii} = lwc;

%% Plot Instantaneous Frequency Downshift
if isPlot
figure();
subplot(2,1,1)
imagesc(GPR.Geolocation.Distance{ii},GPR.D.TimeAxis{ii},instantFreq);
cb = colorbar;colormap(bone)
cb.Label.String = 'Instantaneous Frequency (GHz)';cb.Ticks = [.25,.5,0.75];
cb.Location = "northoutside";
cb.FontSize = 14;
clim([0.25 .75])
set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
ylabel('Trave-time (ns)');
% title('Instantaneous Frequency')
subplot(2,1,2)
imagesc(GPR.Geolocation.Distance{ii},GPR.D.TimeAxis{ii},f0-freq);
% grid on; grid minor;
% xlim([min(GPR.Geolocation.Distance{ii}),max(GPR.Geolocation.Distance{ii})])
% title('Frequency Shift'); 
cb = colorbar;colormap(bone)
cb.Label.String = '\Delta f (GHz)';cb.Ticks = [0,.25,.5];
cb.Location = "northoutside";
cb.FontSize = 14;
xlabel('Distance (m)');ylabel('Trave-time (ns)');clim([0 0.5])
set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
end
end
end