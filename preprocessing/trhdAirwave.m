function [trhd,Rad] = trhdAirwave(trhd,Rad,GPS,Geometry)
% trhdAirwave is a basic trace referencing procedure. Where airwave files
% are stationary.
% Inputs
% trhd = the trace header
% GPS = [trcno,unixTime,longitude,latitude,elevation,geiodN,heading,speed,pdop,hdop,vdop];
% Geometry - srtucture with gpr antenna array & gps geometry
% Outputs
% trhd - the georeferrenced trace header (40 row format)
% Rad -  the preprocessed GPR multiplexed traces (static traces removed)

% Tate Meehan 6/23/24 Juneau Icefield Research Project 2024

%% Check Inputs
if isfield(Geometry,'delta')
else
    delta = zeros(numel(unique(trhd(23,:))),3);
end
if isfield(Geometry,'latency')
else
    latency = zeros(numel(unique(trhd(23,:))),1);
end
% Default Calculation of Heading
if nargin<4
    error('Insufficient Inputs');
end
%% Start of Processing
% Determine number of Channels
nChan = numel(unique(trhd(23,:)));
% Single Channel
if nChan == 1
    trhd(23,:) = 1;
end
multiplexNtrcs = size(trhd,2);
% Remove any unbinned Traces ..
xtraTrcs = mod(multiplexNtrcs,nChan);
if xtraTrcs > 0
    trhd((multiplexNtrcs-xtraTrcs)+1:end) = [];
    Rad((multiplexNtrcs-xtraTrcs)+1:end) = [];
end
multiplexNtrcs = multiplexNtrcs-xtraTrcs;
nchanTrcs = multiplexNtrcs./nChan;
% Move Shot Number
trhd(30,:) = trhd(1,:);
% Correct Trace Number
trhd(1,:) = 1:multiplexNtrcs;

%% Configure GPS Information (may not Exist)
% Remove Non Unique Values
[~,unIx] = unique(GPS(:,2)); % Unique Times
GPS = GPS(unIx,:);
trcno = GPS(:,1);time = GPS(:,2);
% Remove Non Unique Positions (Rare!)
tmp = [GPS(:,3),GPS(:,4)];
[~,unIx] = unique(tmp,'rows','stable');
GPS = GPS(unIx,:);
% Remove NaN Values
GPS = GPS(~isnan(GPS(:,3)),:);
% Remove Zero Values
GPS = GPS(~~(GPS(:,3)),:);
% Define Variables
lon = GPS(:,3); lat = GPS(:,4);z = GPS(:,5);
geoid = GPS(:,6);
if isempty(GPS)
    warning('No GPS information in Airwave Trace Header!')
    GEOID = zeros(1,multiplexNtrcs);
    X = zeros(1,multiplexNtrcs);
    Y = zeros(1,multiplexNtrcs);
    ZONE = zeros(1,multiplexNtrcs);
    % Include Elevations
    Z = zeros(1,multiplexNtrcs);
    try
        time = interp1(trcno,time,multiplexNtrcs,'linear');
    catch
        time = mean(time).*ones(1,multiplexNtrcs);
    end
else
    %% Convert to UTM
    [E,N,utmzone] = deg2utm(lat,lon);
    % Convert UTM zone to interger format
    % Positive Norhern Hemisphere, Negative Southern Hemisphere
    % Possible UTM zones
    zonekey = {'C','D','E','F','G','H','I','J','L','M','N','P','Q','R','S','T','U','V','W'};
    zone = zeros(length(E),1);
    for kk = 1:length(E)
        tmpzone = strsplit(utmzone(kk,:));
        zoneix = find(contains(zonekey,tmpzone{2}));
        if zoneix > 10
            % Northern Hemisphere
            zone(kk) = str2num(tmpzone{1});
        else
            % Southern Hemisphere
            zone(kk) = -str2num(tmpzone{1});
        end
    end
    clear('tmpzone')
    %% GPS Interpolation (OutPerforms MLS)
    GEOID = ones(1,multiplexNtrcs).*mean(geoid);
    X = ones(1,multiplexNtrcs).*mean(E);
    Y = ones(1,multiplexNtrcs).*mean(N);
    ZONE = ones(1,multiplexNtrcs).*mean(zone);
    % Include Elevations
    Z = ones(1,multiplexNtrcs).*mean(z);
    try
        time = interp1(trcno,time,multiplexNtrcs,'linear');
    catch
        time = mean(time).*ones(1,multiplexNtrcs);
    end
end
distance = zeros(1,multiplexNtrcs);
% Append Distance to Trace Header
trhd(2,:) = distance;
% Append GPS Data to Trace Header
trhd(4,:) = time;
% Antenna Midpoint
trhd(13,:) = X;
trhd(14,:) = Y;
trhd(15,:) = Z;
Slope = zeros(1,multiplexNtrcs);
Velocity = zeros(1,multiplexNtrcs);
Heading = zeros(1,multiplexNtrcs);
% Append GPS to Trace Header
trhd(17:19,:) = [Slope;Velocity;Heading];

clear('GPSix','GPSixEdges','dS','dX','dY','dZ''dT','ds','dx','dy','dz''dt',...
    'Distance','Slope','Velocity','Heading','Tailing','Time','tmpTime','GPS',...
    'interpIx','heading','dheading','wrapIx','distance');
% Configure Array Geometry
geometryArray = Geometry.midpointy;  % Midpoint Locations Relative to 0
trhd(21,:) = geometryArray([trhd(23,:)]);
% Append the Offset Array to Trace Header
offsetArray = Geometry.offset;
channelArray = trhd(23,:);
offsetAppend = offsetArray(channelArray);
trhd(22,:) = offsetAppend;
%% Dead Reckoning
% Average Positions for Bin Centers (The GPS Location and Distance)
for jj = 1:nChan:size(trhd,2)
    % Store Bin Centers [m]
    trhd(10:12,jj:jj+nChan-1) = mean(trhd(13:15,jj:jj+nChan-1),2)*ones(1,nChan);
    % Overwrite Distance with Bin Center Position [m]
    tmp = mean(trhd(16,jj:jj+nChan-1),2)*ones(1,nChan);
    trhd(2,jj:jj+nChan-1) = tmp;
    % Overwrite Tailing with Average Bin Center Heading
    trhd(20,jj:jj+nChan-1) = mean(trhd(19,jj:jj+nChan-1),2)*ones(1,nChan);
end

% Zero Starting Distance
tmp = trhd(2,1);
trhd(2,:) = trhd(2,:) - tmp;
trhd(16,:) = trhd(16,:)-tmp;

clear('tmp')

%% Reorganize Trace Header
% GPS time
tmp4 = trhd(4,:);
% frequency MHz
tmp9 = trhd(9,:);
% Move GPS time to Geolocation block
trhd(4,:) = tmp9;
% Move Frequency
trhd(9,:) = tmp4;
clear('tmp4')
% Eliminate Irrelavant array
trhd(24,:) = zeros(1,length(trhd(24,:)));
% channel number
tmp23 = trhd(23,:);
% midpoint geometery
tmp21 = trhd(21,:);
% offset
tmp22 = trhd(22,:);
% shot gather bin center distance
tmp2 = trhd(2,:);
% Antenna Midpoint Distance
tmp16 = trhd(16,:);
% Move Channel Number
trhd(2,:) = tmp23;
clear('tmp23')
% nsamps
tmp3 = trhd(3,:);
% Move Frequency
trhd(3,:) = tmp9;
clear('tmp9')
% dt
tmp7 = trhd(7,:);
% Move dt
trhd(4,:) = tmp7;%./1000;
clear('tmp7')
% Move nsamps
trhd(5,:) = tmp3;
clear('tmp3')
% Move midpoint Geometery
trhd(6,:) = tmp21;
clear('tmp21')
% Move Offset array
trhd(7,:) = tmp22;
clear('tmp22')
% Move Bin Center Geolocation
trhd(8:11,:) = trhd(9:12,:);
% Move Bin Center Distance
trhd(12,:) = tmp2;
clear('tmp2')
% Move Midpoint Geoloctaion
trhd(21:23,:) = trhd(13:15,:);
% Move Antenna Midpoint Distance
trhd(24,:) = tmp16;
clear('tmp16')
% Move Antenna Center Heading
trhd(25,:) = trhd(19,:);
% Move Velocity
trhd(13,:) = trhd(18,:);
% Move Slope
trhd(15,:) = trhd(17,:);
% Move Heading
trhd(14,:) = trhd(20,:);
% Move Midpoint Geolocations
trhd(16:20,:) = trhd(21:25,:);
% Remove Extra Rows
% trhd(21:25,:) = [];
% Use psn2ll and Convert EGM08 to WGS84 bc dead reckoning
[trhd(22,:),trhd(21,:)] = utm2ll(trhd(9,:),trhd(10,:),ZONE);
% Geoid Height: EGM2008
trhd(23,:) = trhd(11,:)+GEOID;
% Midpoint Locations
[trhd(25,:),trhd(24,:)] = utm2ll(trhd(16,:),trhd(17,:),ZONE);
trhd(26,:) = trhd(18,:)+GEOID;
% Add Extra Row for Shot Number
trhd(3:size(trhd,1)+1,:) = trhd(2:end,:);
% Move Shot Number
trhd(2,:) = trhd(31,:);
% Move other GPR meta data
tmp78 = trhd(7:8,:);
trhd(6:8,:) = trhd(4:6,:);
% Replace offset and midpoint
trhd(4:5,:) = flipud(tmp78);
clear('tmp78')
% Add 5 Rows for Array Geomertry
trhd(11:size(trhd,1)+5,:) = trhd(6:end,:);
% Append Geometry
% Move y midpoint
trhd(6,:) = trhd(5,:);
% Add X Midpoint
geometryArray = Geometry.midpointx;  % Midpoint Locations Relative to 0
trhd(5,:) = geometryArray([trhd(3,:)]);
% Add Transmitter X
geometryArray = Geometry.tx;
trhd(7,:) = geometryArray([trhd(3,:)]);
% Add Transmitter Y
geometryArray = Geometry.ty;
trhd(8,:) = geometryArray([trhd(3,:)]);
% Add Reveiver X
geometryArray = Geometry.rx;
trhd(9,:) = geometryArray([trhd(3,:)]);
% Add Reveiver Y
geometryArray = Geometry.ry;
trhd(10,:) = geometryArray([trhd(3,:)]);
% Geoid EGM96
trhd(33,:) = GEOID';
% UTM Zone
trhd(34,:) = ZONE';
% GPS Geometry
trhd(35:37,:) = repmat(Geometry.gps(:),1,size(trhd,2));
% GPS delta
geometryArray = Geometry.delta;
trhd(38:40,:) = geometryArray([trhd(3,:)],:)';
% Correct Midpoint Distance (without Dead Reckoning)
trhd(25,:) = trhd(18,:)+trhd(6,:);

end

