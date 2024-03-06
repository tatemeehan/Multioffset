function [trhd,Rad] = trhdGeoref(trhd,Rad,GPS,Geometry,velocityThreshold,calcVelocityHeading,isInterpolateTrhd,dxGrid,isrmvStaticTrc)
% trhdGeoref performs the preprocessing and interpolation of GPS data.
% Stationary positions and traces are removed. The positions are migrated
% to the midpoints of the of the GPR array. 
% Inputs
% trhd = the trace header
% GPS = [trcno,unixTime,longitude,latitude,elevation,geiodN,heading,speed,pdop,hdop,vdop];
% Geometry - srtucture with gpr antenna array & gps geometry
% velocityThreshold - float (optional - default 0.35 m/s)
% calcVelocityHeading - calculate velocity and heading (if $GPVTG optional - default 1)
% Outputs
% trhd - the georeferrenced trace header (40 row format)
% Rad -  the preprocessed GPR multiplexed traces (static traces removed)

% Tate Meehan 1/9/21 SnowEx 2021

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
if nargin<6
calcVelocityHeading = 1;
end
% Default Velocity Threshold for Static Positions
if nargin<5
velocityThreshold = 0.35;%[m/s]
end
if nargin < 9
    isrmvStaticTrc = 1;
end
%% Start of Processing
% Determine number of Channels
nChan = numel(unique(trhd(23,:)));
% Single Channel
if nChan == 1
    trhd(23,:) = 1;
end
% Remove Non Unique Values
[~,unIx] = unique(GPS(:,2)); % Unique Times
GPS = GPS(unIx,:);
% Remove Non Unique Positions (Rare!)
tmp = [GPS(:,3),GPS(:,4)];
[~,unIx] = unique(tmp,'rows','stable');
GPS = GPS(unIx,:);
% Remove NaN Values
GPS = GPS(~isnan(GPS(:,3)),:);
% Remove Zero Values
GPS = GPS(~~(GPS(:,3)),:);
% Define Variables
lon = GPS(:,3); lat = GPS(:,4);trcno = GPS(:,1);sec = GPS(:,2);z = GPS(:,5);
geoid = GPS(:,6);
multiplexNtrcs = size(trhd,2);
% Move Shot Number
trhd(30,:) = trhd(1,:);
% Correct Trace Number
trhd(1,:) = 1:multiplexNtrcs;

%% Correct for Pauses in Time that will effect interpolation
% Check for Trace Pauses Either Free Run or Odometer..
% Time keeps ticking, but trc number does not
tthresh = 15; % Minimum Pause Time in Seconds
dt = diff(sec'); dtix1 = find(dt>tthresh);dtix2 = dtix1+1;
if ~isempty(dtix1)
    dtix1 = [dtix1,length(trcno)];dtix2 = [1,dtix2];
    segIx = [dtix2;dtix1];
    % We could either adjust the time, or the trace number..
    % We will temporaraily adjust the time
    t = zeros(length(dtix1)-1,1);
    for kk = 1:length(dtix1)-1
        % Calculate Time Difference + one sample interval
        t(kk) = sec(segIx(2,kk)) - sec(segIx(1,kk+1))+mean(diff(sec(segIx(1,kk):segIx(2,kk))));
        % Correct Segment of Time for Pause
        sec(segIx(1,kk+1):segIx(2,kk+1)) = sec(segIx(1,kk+1):segIx(2,kk+1))+t(kk);
    end
    
    % Linear Time Sample Interpolation
    s = sec(:); % Copy of shifted raw times
    sec = interp1(trcno,sec,[1:multiplexNtrcs]','linear','extrap');
    
    % Shift back to original Times
    tmpIx = segIx(:);
    outSeg = zeros(numel(segIx),1);
    for kk = 1:numel(segIx)
        [~,outSeg(kk)] = min(abs(trcno(tmpIx(kk))-[1:multiplexNtrcs]));
    end
    % Ensure Neigboring MidPoints
    for kk = 2:2:numel(segIx)-1
        outSeg(kk) = floor(mean([outSeg(kk),outSeg(kk+1)]));
        outSeg(kk+1) = outSeg(kk)+1;
    end
    % Ensure Correct Ending Indicies
    outSeg(1) = 1; outSeg(numel(segIx)) = multiplexNtrcs;
    % Reshape 
    outSeg = reshape(outSeg,size(segIx));
    segIx = outSeg;
    secOut = sec; % Copy
        % adjust the time
    for kk = 1:size(segIx,2)-1
        % Correct Interpolated Segment of Time for Pause
        secOut(segIx(1,kk+1):segIx(2,kk+1)) = secOut(segIx(1,kk+1):segIx(2,kk+1))-t(kk);
    end
        % Put Interpolated Time in Out Matrix
        timeout = secOut;
        clear('secOut','segIx','outSeg','tmpIx','tthresh')
else
    % Linear Time Sample Interpolation
    s = sec(:); % Copy of raw times
    % Interpolate Trace Times
    sec = interp1(trcno,sec,[1:multiplexNtrcs]','linear','extrap');
    % Put Interpolated Time in Out Matrix
    timeout = sec;
end
%% Convert to Time of Day
isConvertTime = 0;
if isConvertTime
    % Check if Seconds or unix time
    if all(timeout>86400) %unixTime
       time = datetime(timeout,'ConvertFrom','posixtime')
    else
        decsec = timeout./(24*3600);
        time = datestr(decsec,'HH:MM:SS.FFF');
          % Convert to num
%         tmp = datestr(decsec,'HH:MM:SS.FFF');
%         tmp(:,[3,6]) = [];
%         time = str2num(tmp);
    end
else
    time = timeout;
end

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
%% PCHIP Interpolation (OutPerforms MLS)
% LON = pchip(s,lon,sec);
% LAT = pchip(s,lat,sec);
GEOID = pchip(s,geoid,sec);
X = pchip(s,E,sec);
Y = pchip(s,N,sec);
ZONE = interp1(s,zone,sec,'nearest','extrap');
% Include Elevations
Z = pchip(s,z,sec);
% ZWGS84 = Z+GEOID;

% Smooth X,Y,Z
% Filter Smoothing Distance
R = 51;
X = movmean(X,R);
Y = movmean(Y,R);
Z = movmean(Z,R);

% Append GPS Data to Trace Header
trhd(4,:) = time;
% Antenna Midpoint
trhd(13,:) = X;
trhd(14,:) = Y;
trhd(15,:) = Z;
% Put in GPSout

%% Compute Distance, Heading, Velocity, & Slope
X = X(:); Y = Y(:); Z = Z(:);
% Interpolated Sites
dX = [0;diff(X)];
% dX = medfilt1(dX,R);
dX = movmean(dX,R);
dY = [0;diff(Y)];
% dY = medfilt1(dY,R);
dY = movmean(dY,R);
dZ = diff(Z);
dZ = [dZ(1);dZ];
% dZ = medfilt1(dZ,R);
dZ = movmean(dZ,R);
% Raw Sites
dx = [0;diff(E)];
dy = [0;diff(N)];
dz = diff(z);
dz = [dz(1);dz];

% Compute Distance
if all(isnan(Z))
    % If elevation is unavailable
    dS = sqrt(dX.^2+dY.^2);
    ds = sqrt(dx.^2+dy.^2);
    Ds2d = dS;
    ds2d = ds;
else
    dS = sqrt(dX.^2+dY.^2+dZ.^2);
    ds = sqrt(dx.^2+dy.^2+dz.^2);
    dS2d = sqrt(dX.^2+dY.^2);
    ds2d = sqrt(dx.^2+dy.^2);
end
    Distance = cumsum(dS); %[m]
    distance = cumsum(ds);

% Smooth Delta Distance
% dS = medfilt1(dS,R);
dS = movmean(dS,R);
dS2d = movmean(dS2d,R);
trhd(2,:) = Distance; % Configure Distance [m]

% Check if GPVTG was provided
if all(isnan(GPS(:,7)))
    calcVelocityHeading = 1;
end
% Interpolate or Calculate Velocity and Heading
if ~calcVelocityHeading
    % Interpolate Velocity and Heading
    Velocity = pchip(s,GPS(:,8),sec);% pchip(s,lon,sec);
    % Unwrap Heading then Interpolate
    Heading = zeros(size(Velocity));
    dheading = [0;diff(GPS(:,7))];
    %  heading degree threshold
    headThresh = 300;
    wrapIx = find(abs(dheading)>headThresh);
    % Find Static Piositions
    dumIx = find(GPS(:,8)<velocityThreshold);
    rmvIx = find(ismember(wrapIx,dumIx));
    % Remover Static Positions
    wrapIx(rmvIx) = [];
    wrapIx = wrapIx - 1;
    if ~isempty(wrapIx)
        if wrapIx(1) == 1
            GPS(1,7) = GPS(2,7);
            wrapIx(1) = [];
        end
        if any(diff(wrapIx)==1)
            % Locate Bad Bin
            diffIx = find(diff(wrapIx)==1);
            % Get GPS index
            headIx = wrapIx(diffIx)+1;
            % Remove Bad Bin
            wrapIx(diffIx) = [];
            % Correct Heading Data using nearest neighbors
            for kk = 1:length(headIx)
                if headIx~=1
                    GPS(headIx(kk),7) = mean([GPS(headIx(kk)-1,7),GPS(headIx(kk)+1,7)]);
                elseif headIx~=length(GPS(:,7))
                    GPS(headIx(kk),7) = GPS(headIx(kk)+1,7);
                else
                    GPS(headIx(kk),7) = GPS(headIx(kk)-1,7);
                end
            end
        end
        % Interpolate UnWrapped Segments
        for kk = 1:length(wrapIx)+1
            if kk == 1
                [~,interpIx] = min(abs(Distance-distance(wrapIx(kk))));
                interpIx = 1:interpIx;
                Heading(interpIx) = pchip(distance(1:wrapIx(kk)),GPS(1:wrapIx(kk),7),Distance(interpIx));
            elseif kk == length(wrapIx)+1
                %                 [~,interpIx] = min(abs(Si-GPS(wrapIx(kk-1)+1,4)));
                interpIx = interpIx(end)+1;
                interpIx = interpIx:length(Distance);
                Heading(interpIx) =  pchip(distance(wrapIx(kk-1)+1:end),GPS(wrapIx(kk-1)+1:end,7),Distance(interpIx));
            else
                %                 [~,interpIx1] = min(abs(Si-GPS(wrapIx(kk-1)+1,4)));
                [~,interpIx2] = min(abs(Distance-distance(wrapIx(kk))));
                interpIx = interpIx(end)+1:interpIx2;
                %                 interpIx = interpIx1:interpIx2;
                Heading(interpIx) = pchip(distance(wrapIx(kk-1)+1:wrapIx(kk)),GPS(wrapIx(kk-1)+1:wrapIx(kk),7),Distance(interpIx));
            end
        end
    else
        % PCIHP
        Heading = pchip(distance,GPS(:,7),Distance);
    end
    % Ensure that no points equal 360 exactly
    ix360 = Heading==360;
    Heading(ix360)=0;
else
    % Compute Velocity and Heading
    Time = sec;
    dT = [mean(diff(Time));diff(Time)];
    dt = [mean(diff(s));diff(s)];
    Velocity2d = dS2d./dT;
    velocity2d = ds2d./dt;
    Velocity = dS./dT; % [m/s]
    velocity = ds./dt;
    Velocity = movmean(Velocity,R);
    Velocity2d = movmean(Velocity2d,R);
    
    % Find Static Positions
    dumIx = find(velocity2d<velocityThreshold);
    
    % Heading
    dH = [0,0;diff([X,Y])];
    dh = [0,0;diff([E,N])];
    Heading = mod(atan2(dH(:,1), dH(:,2))*180/pi, 360);
    heading = mod(atan2(dh(:,1), dh(:,2))*180/pi, 360);
    dheading = [0;diff(heading)];
    %  heading degree threshold
    headThresh = 300;
    wrapIx = find(abs(dheading)>headThresh);
    rmvIx = find(ismember(wrapIx,dumIx));
    % Remover Static Positions
    wrapIx(rmvIx) = [];
    wrapIx = wrapIx - 1;
    for kk = 1:length(wrapIx)
        [~,interpIx(kk)] = min(abs(Time-s(wrapIx(kk))));
        Heading(interpIx(kk)+1:end) = Heading(interpIx(kk)+1:end)-360;
    end
    % Ensure that no points equal 360 exactly
    ix360 = Heading==360;
    Heading(ix360)=0;
end
% Slope
Slope = atand(dZ./dS);
Slope = movmean(Slope,R);

% Append GPS to Trace Header
trhd(17:19,:) = [Slope';Velocity';Heading'];

%             GPS = [GPS,Distance,Slope,Velocity,Heading];
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

%% Static/Skipped Trace Removal

% Determine Data Acquisition  Method
if all(trhd(5,:) == 0)
    isFreeRun = 1;
    isOdometer = 0;
else
    isOdometer = 1;
    isFreeRun = 0;
end
% Remove Static Traces if is Free Run Acquisition
if isFreeRun
    if isrmvStaticTrc
    % Remove Static Trace Headers from Multiplexed Record
    dupIx = removeStaticPositions(trhd,velocityThreshold,Velocity2d);
    else
        dupIx = [];
    end
    % Check if GP2 was Corrupted
%     w = warning('query','last');
%     if strcmp(w.identifier,'GP2:Corrupted')
%         % Three the hard way
%         rmvIx = 60725;
%         % Ensure that we remove Channel 1 onward
%         rmvChan = trhd(23,rmvIx);
%         rmvIx = rmvIx-(rmvChan-1);
%         rmvIx = [dupIx;[rmvIx:max(dupIx)]'];
%         rmvIx = unique(rmvIx);
%     end
    trhd(:,dupIx) = [];
    % Configure Trace Indicies
    trhd(1,:) = 1:length(trhd);
    % Configure Shot Number
    trhd(30,:) = trhd(30,:)-trhd(30,1)+1;
    % re-Incorporate GPS here
    X = trhd(13,:)';Y = trhd(14,:)';Z = trhd(15,:)';
    % Linear Approximation
    isLinear = 0;
    isArclength = 1;
    if isLinear
    dX = [0;diff(X)];
%     dX = medfilt1(dX,R);
    dX = movmean(dX,R);
    dY = [0;diff(Y)];
%     dY = medfilt1(dY,R);
    dY = movmean(dY,R);
    dZ = diff(Z);
    dZ = [dZ(1);dZ];
%     dZ = medfilt1(dZ,R);
    dZ = movmean(dZ,R);
    % Compute Distance
    if quantile(Z,0.5 < 0)
        % If elevation is unavailable
        dS = sqrt(dX.^2+dY.^2);
    else
        dS = sqrt(dX.^2+dY.^2+dZ.^2);
    end
%     dS = medfilt1(dS,R);
    dS = movmean(dS,R);
    Distance = cumsum(dS); %[m]
%     if ii>1
%         Distance = endDist+cumsum(dS);
%     else
%         Distance = cumsum(dS); %[m]
%     end
%     endDist = Distance(end);
    end
    if isArclength
        % Use Arclength Calculation along smooth curve
        [~,seglen] = arclength(X,Y,Z,'pchip');
        Distance = [0;cumsum(seglen)]';
    end
    trhd(2,:) = Distance; % Configure Distance
    % Configure Distance with midpoint geometry (in Direction of Travel y)
    trhd(16,:) = geometryArray([trhd(23,:)])+trhd(2,:); % Append Midpoint Locations
    % Remove Static Positions from GPSout
    GEOID(dupIx) = [];
    ZONE(dupIx) = [];
    % Remove Static Traces from Multiplexed Data
%     if strcmp(w.identifier,'GP2:Corrupted')
%             Rad(:,rmvIx) = [];
%     else
            Rad(:,dupIx) = [];
%     end
    %             xArray = trhd(2,:);         % Define Configured Distance xArray    
end

%% Fix the Fucked Up File
% for jj = 1:nChan
%     ix = find(trhd(23,:) == jj);
%     rmvtrhd(jj,:) = ix(3523:3523+505);
%     Radar{jj} = Rad(:,ix);
%     if trcnoBad(jj)
% %         rmvtrc(jj,:) = ix(4481:4981);
% %         rmvtrc(jj,:) = ix(4478:4981);
%         Radar{jj}(:,4478:4981) = [];
%         % Remove nans?
%         Radar{jj}(:,[3475,3523]) = [];
% 
%     else
% %         rmvtrc(jj,:) = ix(3523:3523+500);
% %         rmvtrc(jj,:) = ix(3523:3523+503);
% Radar{jj}(:,3523:3523+503) = [];
%         Radar{jj}(:,[3475,3523]) = [];
% 
% 
%     end   
% end
% Multiplex the Radar grams again fml
% Radder = zeros(size(Radar{1},1),size(Radar{1},2).*9);
% for jj = 1:nChan
%     Radder(:,jj:nChan:end) = Radar{jj};
% end
% trhd(:,rmvtrhd(:)) = [];
% GEOID(rmvtrhd(:)) = [];
% ZONE(rmvtrhd(:)) = [];
% Rad = Radder;
% tmpRad(:,sort(rmvtrc(:))) = [];

%% Remove Skip Traces if Wheel Odometer Acquisition
if isOdometer
    disp('Somethings Wrong,Because Data is Free Run')
    dupIx = find(~(diff(trhd(23,:)))); % Find Skipped Traces
    
    for jj = dupIx
        trhd(1,jj:end) = trhd(1,jj:end) - 1; % ReConfigure Trace Indicies
    end
    
    trhd(:,dupIx) = []; % Remove Skipped Traces from Trace Header
    Rad(:,dupIx) = []; % Remove Skipped Traces from Multiplexed Data
%     xArray = trhd(2,:); % Define ReConfigured Distance as xArray
% Remove Static Positions from GPSout
GEOID(dupIx) = [];
ZONE(dupIx) = [];
end


%% Dead Reckoning
% Calculate True Antenna Midpoint Positions using Array Geometry
% delta is the correction from GPS antenna center to reflection midpoint 
tmptrhd = trhd;
for kk = 1:nChan
    [trhd] = deadReckonSPIDAR(trhd,kk,Geometry.delta(kk,:));
end
% Smooth Positions post-Reckoning
% R = 5;
% for kk = 13:15
%     for jj = 1:nChan
%         chanIx = find(trhd(23,:)==jj);
%     trhd(kk,chanIx) = movmean(trhd(kk,chanIx),R);
%     end
% end

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
[trhd(22,:),trhd(21,:)] = utm2ll(trhd(9,:),trhd(10,:),ZONE');
% Geoid Height: EGM2008
trhd(23,:) = trhd(11,:)+GEOID';
% Midpoint Locations
[trhd(25,:),trhd(24,:)] = utm2ll(trhd(16,:),trhd(17,:),ZONE');
trhd(26,:) = trhd(18,:)+GEOID';
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

% If Data Size does not fit Trace Header Size - Interpolate Trace Header
% Allocation
[mtrhd,ntrhd] = size(trhd);ntrcs = size(Rad,2);
if ntrhd~= ntrcs
    keyboard
    for kk = 1:nChan
        tmptrhd = trhd(:,kk:nChan:ntrhd);
        tmpntrcs = size(kk:nChan:ntrcs,2);
        intrptrhd = zeros(mtrhd,tmpntrcs);
        for ll = 1:mtrhd
            intrptrhd(ll,:) = interp1(tmptrhd(1,:),tmptrhd(ll,:),kk:nChan:ntrcs,'nearest','extrap');
        end
        newtrhd{kk} = intrptrhd;
    end
end
% Default Interpolation
if nargin < 7
    isInterpolateTrhd = 1;
    dxGrid = 0.1;
end
% Default dx
if nargin < 8
    dxGrid = 0.1;
end
% Interpolate Distance to Equal Trace Interval
if isInterpolateTrhd
% Calculate ArcLength Along GPR Transect
totalLength = trhd(18,end);
% totalLength = floor(totalLength.*(1./dxGrid))./(1./dxGrid); %Rounding Error
% % Interpolate Coordinates Along Path
% pt = interparc(0:(dxGrid/totalLength):1,trhd(15,1:nChan:end),trhd(16,1:nChan:end),trhd(17,1:nChan:end),'pchip');
% % Recalculate Distance Axis for Header interpolation
% [~,seglen] = arclength(pt(:,1),pt(:,2),pt(:,3),'pchip');
% seglen = round(seglen.*(1./dxGrid))./(1./dxGrid);% Tamp Rounding Error
% Distance = round(cumsum([0;seglen]).*(1./dxGrid))./(1./dxGrid);% Rounding
Distance = [0:dxGrid:totalLength];
oldDistance = trhd(18,1:nChan:end);
intrpTrhd = zeros(size(trhd,1),length(Distance).*nChan);
% FID
intrpTrhd(1,:) = 1:length(Distance).*nChan;
% Shot No.
tmpShotNo = repmat(1:length(Distance),nChan,1);
intrpTrhd(2,:) = [tmpShotNo(:)]';
% Channel Number
intrpTrhd(3,:) = repmat(1:nChan,1,length(Distance));
% Offset
intrpTrhd(4,:) = offsetArray(intrpTrhd(3,:));
% Midpoint X
intrpTrhd(5,:) = Geometry.midpointx(intrpTrhd(3,:));
% Midpoint Y
intrpTrhd(6,:) = Geometry.midpointy(intrpTrhd(3,:));
% TX RX Geomoetery
intrpTrhd(7:10,:) = [Geometry.tx(intrpTrhd(3,:));Geometry.ty(intrpTrhd(3,:));...
    Geometry.rx(intrpTrhd(3,:));Geometry.ry(intrpTrhd(3,:))];
% GPR Frequency Parameters
intrpTrhd(11:13,:) = trhd(11:13,intrpTrhd(3,:));
% Unix Time
for jj = 1:nChan
intrpTrhd(14,jj:nChan:length(Distance).*nChan) = interp1(oldDistance,trhd(14,jj:nChan:end),Distance,'linear');
end
% Midpoint Coordinates etc.
for kk = 15:33
for jj = 1:nChan
intrpTrhd(kk,jj:nChan:length(Distance).*nChan) = interp1(oldDistance,trhd(kk,jj:nChan:end),Distance,'pchip');
end
end
% UTM Zone
for jj = 1:nChan
intrpTrhd(34,jj:nChan:length(Distance).*nChan) = interp1(oldDistance,trhd(34,jj:nChan:end),Distance,'nearest');
end
% GPR Antenna Position
intrpTrhd(35:37,:) = [repmat(Geometry.gps(1),1,length(Distance).*nChan);...
    repmat(Geometry.gps(2),1,length(Distance).*nChan);...
    repmat(Geometry.gps(3),1,length(Distance).*nChan)];
% Dead Reckon Lever Arm
intrpTrhd(38:40,:) = [Geometry.delta(intrpTrhd(3,:),1)';...
    Geometry.delta(intrpTrhd(3,:),2)';Geometry.delta(intrpTrhd(3,:),3)'];
% Interpolate the Multiplexed Image
newRad = zeros(size(Rad,1),length(Distance).*nChan);
for kk = 1:size(Rad,1)
    for jj = 1:nChan
    newRad(kk,jj:nChan:length(Distance).*nChan) = interp1(oldDistance,Rad(kk,jj:nChan:end),Distance,'pchip');
    end
end

% Assign Trhd and Rad
trhd = intrpTrhd;
Rad = newRad;
end

end

