function [GPR] = qStar(GPR,f0,fmethod,qmethod,chan,isPlot)
% qStar.m computes the instantaneous frequency of the GPR image and
% extracts the frequency content of the ground reflection. The frequency
% downshift is estimated from f0, the instantaneous frequency of a
% reference airwave signal. The reference and measured frequency is supplied
% into Bradford's 2009 method for estimation of Q* attentuation. Q* is used
% to estimate the imaginary component of snow dielectric permittivity,
% supplied the real part, which then can estimate snow liquid water
% content.
% This function requires that the two-way travel-time of the ground
% reflection be known and stored in the GPR structure.
%
% Inputs: GPR - the GPR data Structure
%          f0 - apriori air wave instantaneous frequency
%     fmethod - Method for picking instant freq of ground ('avg' or default: 'max')
%     qmethod - Q* Method ('Bradford07' or default: 'Bradford09')
%
% Outputs: Qstar    - Frequency Attenuation Factor Q* = eps'./2.*eps"
%          QstarStd - Standard Deviation of Q*
% Tate Meehan - 12/12/24 while he was not at AGU2024

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
n = 10;

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

%% Estimate Q*
% Scalar Peak Frequency For Ricker Wavelet
f0 = f0(chan);
f0 = f0./1.13;
freq = freq./1.13;

% Angular Frequency
w0 = 2.*pi.*f0; %w0 = w0(:);
wt = 2.*pi.*(freq); %wt = wt(:);
df = w0-wt;

% Bradford Retreival Method
if strcmp(qmethod,'Bradford07')
    % Bradford 2007
    invQstar = (4./(GPR.D.TimeAxis{ii}(:))).*((w0.^2-wt.^2)./(w0.^2.*wt));
elseif strcmp(qmethod,'Bradford09')
    % Bradford 2009
    invQstar = (2./(pi.*GPR.D.TimeAxis{ii}(:))).*((w0.^2-wt.^2)./(w0.^2.*wt));
else
    invQstar = (3./(GPR.D.TimeAxis{ii}(:))).*((w0.^2-wt.^2)./(w0.^2.*wt));
end

% Median Filter
invQstar = medfilt2(invQstar,[filtm, filtn]);
% Q*
Qstar = 1./invQstar;
% % Q* Mean and Standard Deviation
% QstarStd = movstd(Qstar,L);
% % SG
% % Qstar = movmean(Qstar,L./2);
% O = 5; % polynomial order
% Qstar= sgolayfilt(Qstar,O,L); % m

%% Liquid Water Content Estimation
% Store Imaginary Component Estimate
GPR.D.stackingImagPermittivity{ii} = GPR.D.stackingPermittivity{ii}./(2.*Qstar);
GPR.D.intervalImagPermittivity{ii} = GPR.D.intervalPermittivity{ii}./(2.*Qstar);

% LWC inversion
Tk = 273.15;
kwimag = imag(water_permittivity_maetzler87(GPR.D.f0GHz{ii}.*10.^9, Tk));
% kimag = GPR.D.stackingImagPermittivity{ii};
% Use Instantaneous Velocity for Instantaneous Frquency
kimag = GPR.D.intervalImagPermittivity{ii};
% Minimization Scheme
% Estimate LWC
testlwc = [0:0.001:0.5];
% Shivola & Tiuri (1984) Inversion
lwc = zeros(size(kimag));
for kk = 1:numel(kimag)
    testLWC = [(0.1.*testlwc+0.8.*testlwc.^2).*kwimag - kimag(kk)]';
    [~,lwcIx] = min(abs(testLWC));
    lwc(kk) = testlwc(lwcIx);
end
GPR.D.qLWC{ii} = lwc;

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