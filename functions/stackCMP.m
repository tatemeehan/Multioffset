function [GPR] = stackCMP(GPR,ds,winsize)
% stackCMP.m stacks adjacent CMP gathers by rolling average.
% Inputs
% GPR - the data strcuture containing 2D matricies of CMP gathers
% ds - the linear sample interval distance (meters)
% winsize - the averaging window size meters

% Outputs
% GPR.D.CMP - Stacked CMP Gathers (OverWrites Original GPR.D.CMP)
% GPR.Geolocation.cmpDistance - Distance Axis of Stacked CMP Gathers
% GPR.Geolocation.cmpX - X Coordinate of Stacked CMP Gathers
% GPR.Geolocation.cmpY - Y Coordinate of Stacked CMP Gathers
% GPR.Geolocation.cmpZ - Z Coordinate of Stacked CMP Gathers

% Written By: Tate Meehan, CRREL, Juneau Icefield Research Project 2024

% Set Default Stacking Window Size to 5 m
if nargin < 2
    ds = 5;
    winsize = 5;
elseif nargin < 3
    winsize = 5;
end
for ii = 1 : GPR.MD.nFiles
    S = GPR.Geolocation.Distance{ii};
    C = GPR.D.CMP{ii};
    newS = min(S):ds:max(S);
    newC = cell(1, length(newS));
    [ntrc,nchan] = size(GPR.D.CMP{ii}{1});
    for kk = 1:length(newS)
%     parfor (kk = 1:length(newS),GPR.MD.nWorkers)
        % Calulate Distance
        dist = abs(S-newS(kk));
        % Indicies within insize
        distIx = find(dist<winsize);
        ndist = length(distIx);
        % Create a Temp Matrix of the many C's
        tmpC = zeros(ntrc,nchan,ndist);
        for jj = 1:ndist
            tmpC(:,:,jj) = C{distIx(jj)};
        end
        % Mean CMP Stacking
        newC{kk} = mean(tmpC,3);
    end
    % Store Stacked CMP gathers
    GPR.D.CMP{ii} = newC;
    % Store New Disctance Axis
    GPR.Geolocation.cmpDistance{ii} = newS;
    % Interpolate CMP X,Y,Z Coordinates
    newX = interp1(GPR.Geolocation.Distance{ii},GPR.Geolocation.X{ii},newS,'pchip','extrap');
    newY = interp1(GPR.Geolocation.Distance{ii},GPR.Geolocation.Y{ii},newS,'pchip','extrap');
    newZ = interp1(GPR.Geolocation.Distance{ii},GPR.Geolocation.Z{ii},newS,'pchip','extrap');
    % Store CMP Coordinates
    GPR.Geolocation.cmpX{ii} = newX; GPR.Geolocation.cmpY{ii} = newY; GPR.Geolocation.cmpZ{ii} = newZ;
end
% End of Function
end