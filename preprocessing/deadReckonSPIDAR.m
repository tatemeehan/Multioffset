function [trhd] = deadReckonSPIDAR( trhd, chan, delta )
% Extract Multiplexed Traces
trcIx = find(trhd(23,:)==chan);
warning('off','MATLAB:mir_warning_changing_try_catch')
% Initialize Delta if Not Supplied
try delta;
    % Get Delta for the Channel
    if size(delta,1)>1
        delta = delta(chan,:);
    end
catch
    % Init 0 perturbations
    delta = [0,0,0];
end

% Dead Reckoning
theta = (mod(trhd(19,trcIx),90));
for kk = 1:length(trcIx)
    % Convert Heading Relative to Left Axis of Quadrant
    tmp = mod(trhd(19,trcIx(kk)),360);
    if tmp >= 0 && tmp < 90 %Q1
        dlx = delta(1).*cosd(theta(kk));
        dly = delta(1).*sind(theta(kk));
        dx = (delta(2)).*sind(theta(kk))+dlx;
        dy = (delta(2)).*cosd(theta(kk))-dly;
        dz = ([delta(2)].*sind(trhd(17,trcIx(kk))))+delta(3);
        % Overwrite Antenna Midpoint Positions in Trace Header
        trhd(13:15,trcIx(kk)) = trhd(13:15,trcIx(kk)) + ([-dx;-dy;-dz]);
    elseif tmp >= 90 && tmp < 180%Q4
        dlx = delta(1).*sind(theta(kk));
        d1y = delta(1).*cosd(theta(kk));
        dx = (delta(2)).*cosd(theta(kk))-dlx;
        dy = delta(2).*sind(theta(kk))+d1y;
        dz = ([delta(2)].*sind(trhd(17,trcIx(kk))))+delta(3);
        % Overwrite Antenna Midpoint Positions in Trace Header
        trhd(13:15,trcIx(kk)) = trhd(13:15,trcIx(kk)) + ([-dx;dy;-dz]);
    elseif tmp >= 180 && tmp < 270%Q3
        dlx = delta(1).*cosd(theta(kk));
        dly = delta(1).*sind(theta(kk));
        dx = (delta(2)).*sind(theta(kk))+dlx;
        dy = (delta(2)).*cosd(theta(kk))-dly;
        dz = ([delta(2)].*sind(trhd(17,trcIx(kk))))+delta(3);
        % Overwrite Antenna Midpoint Positions in Trace Header
        trhd(13:15,trcIx(kk)) = trhd(13:15,trcIx(kk)) + ([dx;dy;-dz]);
    else%Q2
        dlx = delta(1).*sind(theta(kk));
        d1y = delta(1).*cosd(theta(kk));
        dx = (delta(2)).*cosd(theta(kk))-dlx;
        dy = delta(2).*sind(theta(kk))+d1y;
        dz = ([delta(2)].*sind(trhd(17,trcIx(kk))))+delta(3);
        % Overwrite Antenna Midpoint Positions in Trace Header
        trhd(13:15,trcIx(kk)) = trhd(13:15,trcIx(kk)) + ([dx;-dy;-dz]);
    end
end

% Distance Array Corrected to MidPoint Centers
trhd(16,trcIx) = trhd(16,trcIx) + trhd(21,trcIx);
% Unwrap Heading to 0 - 360
trhd(19,trcIx) = mod(trhd(19,trcIx),360);
% Ensure that no points equal 360 exactly
ix360 = trhd(19,trcIx)==360;
trhd(19,trcIx(ix360))=0;
end