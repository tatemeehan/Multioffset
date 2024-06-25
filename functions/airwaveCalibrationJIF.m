function [ timeZeroData, twt, to ] = airwaveCalibrationJIF( data, t, dt, chanNo, offset, path,file)
% airwaveCalibrationJIF.m applies predetermined travel-time perturbations
% to the common offset gather radargrams. Airwave calibrations are picked 
% prior using pickAirwaveCalibraiton.m. The data files are lengthened as
% to not truncate data at the far offsets.
%
%   Input:      data   - Uncorrected Radargram Data Matrix
%               t      - Two-way Travel Time Window   
%               dt     - Sample interval [ns]
%               chanNo - Channel Number of Radargram
%               offset - Tx - Rx Offset of Radargram
%               path   - The Directory containing the calibration file
%               file   - The calibration .mat file

%   Output:     timeZeroData - Time Zero Correced Data
%               twt          - New (lengthened) twt axis
%               to           - Intercept Time (Removed Time)

% Written by Tate Meehan, CRREL, Juneau Icefield Research Project 2024

% Using Prior AirWave Calibration Picks Adjust the RadarGrams
% Give the Calibration Picks a touch of overhead Travel-time

% Size of Radargram
[m,n] = size(data);
% Airwave Velocity
c = 0.3;
% Load Airwave Calibraiton File
load(fullfile(path,file));
% Calcuate Median of Airwave Calibrations
perturbation = round(median(cat(1,AirwaveCalibration.ixShift{:})));
% Determine Length of New Gathers (Add some Samples to Lengthen Window)
nsamp = numel(t) + max(perturbation);
% Round to Next 50 ns
tmp = (nsamp-1).*dt;
tmp = round(tmp/50)*50;
nsamp = (tmp./dt)+1;
% Configure New TWT Axis
twt = 0:dt:(nsamp-1).*dt;
% Allocate Radargram
timeZeroData = zeros(nsamp,n);
% Truncate Radargram at Airwave
t0ix = round(median(cat(1,AirwaveCalibration.ixPick{:})));
t0ix = t0ix(chanNo)-1;
t0calIx = round((offset./c)./dt);
zeroData = mean(mean(data(1:t0ix,:)));
% Append Padded Data to New Radargram
datIx = m-t0ix;
padIx = nsamp - (datIx);
% Time-Zero Correction
timeZeroData(1:t0calIx-1,:) = ones(t0calIx-1,n).*zeroData;
timeZeroData(t0calIx:datIx+t0calIx,:) = data(t0ix:m,:);
timeZeroData(datIx+1:nsamp,:) = zeroData.*ones(padIx,n);
% The Theroretical Time-Zero
to = t0calIx.*dt;
end


