function [GPS, delta, latency] = readGP2(filename,OGdate)
%readGP2 extracts GPS data from NMEA string embedded in S&S .GP2 files.
%   Input :
%       filename = /directory/filename.GP2 - the fullfile sting 
%   Output:
%       GPS = [trcno,unixTime,longitude,latitude,elevation,geiodN,heading,speed,pdop,hdop,vdop];
%           trcno - Trace Number
%           unixTime - The UTC date and time Zulu in unix time format
%           longitude [DD]
%           latitude  [DD]
%           elevation [masl] EGM96
%           geiodN [N] Geiod Height EGM96
%           heading [DD] 
%           speed [m/s]
%           pdop - position dilution of precision
%           hdop - horizontal dilution of precision
%           vdop - vertical dilution of precision
%       delta - [X,Y,Z] relative position of GPS antenna [m]
%       latency - the GPS position to system latency time [s] - unsused 

% Eligible NMEA Strings
% GPGGA
% GPVTG
% GPZDA
% GPRMC
% GPGSA

% Tate Meehan 12/21/2020
% Read Date
OGdate = strsplit(OGdate,'-');
% Open the .GP2 file
fid=fopen(filename);
% Read the whole lines of .GP2 file into a cell array
rawGPS=textscan(fid,'%s','Delimiter','\n');
% Close the file
fclose(fid);
% Reduce Dimension of Cell Array
rawGPS = rawGPS{1};
% Extract Relative Antenna Position (delta perturbations in deadReckon)
tmpdelta = strsplit(rawGPS{5},'=');
delta = str2num(tmpdelta{2});
% Extract GPS Latency (Currently Unused)
tmplatency = strsplit(rawGPS{6},'=');
latency = str2num(tmplatency{2});
% Remove Unwanted Header Rows
rawGPS(1:8) = [];
m = size(rawGPS,1);
% Get Indicies of NMEA Strings
for kk = 1:m
    tmp1 = strsplit(rawGPS{kk},'"');
    tmp2 = strsplit(tmp1{2},',');
    NMEA(kk,:) = tmp2{1};
    if any(strcmp(NMEA(kk,:),["$GPGGA","$GPVTG","$GPZDA","$GPRMC","$GPGSA"])) && ~exist("NMEA1","var")
        NMEA1 = NMEA(kk,:);
    end
    if exist("NMEA1","var")
    iterNMEA(kk) = strcmp(NMEA(kk,:),NMEA1);
    end
    rmvNMEA(kk) = ~any(strcmp(NMEA(kk,:),["$GPGGA","$GPVTG","$GPZDA","$GPRMC","$GPGSA"]));
end
% Remove Undigestible NMEA Strings
NMEA(rmvNMEA,:) = [];
rawGPS(rmvNMEA) = [];
iterNMEA(rmvNMEA) = [];
m = size(rawGPS,1);
% Check that all NMEA Strings were Read by SPIDAR
% Roughly Determine number of NMEA Strings
% nNMEA = quantile(diff(find(iterNMEA)),0.99);
nNMEA = size(unique(NMEA,'rows'),1);

% Find Index of missing NMEA sequence
noNMEA = find(diff(find(iterNMEA))<nNMEA);
tmpIx = find(iterNMEA);
% Change iteration index to false
iterNMEA(tmpIx(noNMEA)) = false;
iter = 1;
trcno = zeros(numel(find(iterNMEA)),1);
% Initialize String Condiditons
isGGA = 0;
isZDA = 0;
isVTG = 0;
isRMC = 0;
isGSA = 0;
for kk = 1:m
    tmprawGPS = erase(rawGPS{kk},'"');
    gpsCell = strsplit(tmprawGPS,',','CollapseDelimiters', false);
    % Reconstruct Trace Number
    if iterNMEA(kk)
        trc1 = str2num(gpsCell{1});
        trcno(kk) = trc1;
    else
        trcno(kk) = trc1;
    end
    
    if kk > 1
       dtrcno = trcno(kk)-trcno(kk-1);
       if dtrcno > 0
           iter = iter+1;
           isGGA = 0;
           isZDA = 0;
           isVTG = 0;
           isRMC = 0;
           isGSA = 0;
       end
    end
    % Extract GPS Data
    % GPGGA
    % Lon, Lat, Elevation, Time
    if strcmp(gpsCell{5},'$GPGGA')
        isGGA = 1;
%         isGPS = str2num(gpsCell{6})~=0; % GPS fix
        isGPS = ~isempty(gpsCell{8}); % (str2num(gpsCell{8}) == 1 | str2num(gpsCell{8}) == 2) |%length(gpsCell)>10; % GPS fix
        if isGPS
            % Longitude
            lon = str2num(gpsCell{9}(1:3));
            lonm = str2num(gpsCell{9}(4:end));
            % Convert to Decimal Degees
            lonm = lonm./60;
            lon = lon+lonm;
            % Check Sign
            if (gpsCell{10}) == 'W'
                lon = -lon;
            end
            % Latittude
            % Get the Latitude DDMM.MMMMM
            lat = str2double(gpsCell{7}(1:2));
            latm = str2double(gpsCell{7}(3:end));
            % Convert to Decimal Degees
            latm = latm./60;
            lat = lat+latm;
            % Check Sign
            if (gpsCell{8}) == 'S'
                lat = -lat;
            end
            % Elevation
            z = str2double(gpsCell{14});
            N = str2double(gpsCell{16});
            if ~isZDA
                % Time        % Time of Day
                t = gpsCell{6};
                H = str2double(t(1:2));
                M = str2double(t(3:4));
                S = str2double(t(5:end));
                s = H.*3600+M.*60+S; % Seconds of Day
                if iter == 1
                    s1 = s;
                end
                d = str2double(OGdate{3});
                % Check if Next Day
                if s < s1
                    d = d+1;
                end
                mnth = month(datetime(OGdate{2},'InputFormat','MMM'));
                y = str2double(OGdate{1});
                dt = datetime([y mnth d H M S]);
                unixTime(iter) = posixtime(dt); % Convert to unixTime

%             six = find(s<s(1));
%             s(six) = s(six)+86400;
%             clear('t','H','M','S');
%             unixTime(iter) = s;
            end
        else
            lon = NaN;
            lat = NaN;
            z = NaN;
            N = NaN;
            s = NaN;
        end
        longitude(iter) = lon;
        latitude(iter) = lat;
        elevation(iter) = z;
        geoidN(iter) = N;
        clear('lon','lat','z','N','s')
    end
    % GPVTG
    % True Heading, Speed km/h
    if strcmp(gpsCell{5},'$GPVTG')
        isVTG = 1;
        isGPS = ~strcmp(gpsCell{6},'T');
        if isGPS
            heading(iter) = str2double(gpsCell{6});
            speedk = str2double(gpsCell{12}); % km/h
            speed(iter) = (speedk.*1000)./(3600); %km/h 2 m/s
            clear('speedk')
        else
            heading(iter) = NaN;
            speed(iter) = NaN;
        end
    end
    % GPZDA
    % Date and Time
    if strcmp(gpsCell{5},'$GPZDA')
        isZDA = 1;
        % Time of Day
        t = gpsCell{6};
        H = str2double(t(1:2));
        M = str2double(t(3:4));
        S = str2double(t(5:end));
        d = str2double(gpsCell{7});
        mnth = str2double(gpsCell{8});
        y = str2double(gpsCell{9});
        dt = datetime([y mnth d H M S]);
        unixTime(iter) = posixtime(dt); % Convert to unixTime
        clear('t','H','M','S','d','mnth','y','dt');
    end
    % GPRMC
    % True Heading, Speed km/h, Time & Date
    if strcmp(gpsCell{5},'$GPRMC')
        isRMC = 1;
        isGPS = ~strcmp(gpsCell{6},'V');
        if ~isGGA % Get Lon, Lat (Elevation Unavailable)
            if isGPS
                % Longitude
                lon = str2num(gpsCell{10}(1:3));
                lonm = str2num(gpsCell{10}(4:end));
                % Convert to Decimal Degees
                lonm = lonm./60;
                lon = lon+lonm;
                % Check Sign
                if (gpsCell{11}) == 'W'
                    lon = -lon;
                end
                % Latittude
                % Get the Latitude DDMM.MMMMM
                lat = str2double(gpsCell{8}(1:2));
                latm = str2double(gpsCell{8}(3:end));
                % Convert to Decimal Degees
                latm = latm./60;
                lat = lat+latm;
                % Check Sign
                if (gpsCell{9}) == 'S'
                    lat = -lat;
                end
                % Elevation
                z = NaN;
                N = NaN;
            else
                lon = NaN;
                lat = NaN;
                z = NaN;
                N = NaN;
            end
                longitude(iter) = lon;
                latitude(iter) = lat;
                elevation(iter) = z;
                geoidN(iter) = N;
                clear('lon','lat','z','N')
        end
        if ~isZDA % Get Date & Time
            if isGPS
                % Time of Day
                t = gpsCell{6};
                H = str2double(t(1:2));
                M = str2double(t(3:4));
                S = str2double(t(5:end));
                date = gpsCell{14};
                d = str2double(date(1:2));
                mnth = str2double(date(3:4));
                y = str2double(date(5:6))+2000;
                dt = datetime([y mnth d H M S]);
                unixTime(iter) = posixtime(dt); % Convert to unixTime
                clear('t','H','M','S','d','mnth','y','dt','date');
            else
                unixTime(iter) = NaN;
            end
        end
        if ~isVTG % Get Heading and Speed
            if isGPS
                heading(iter) = str2double(gpsCell{13});
                speedkn = str2double(gpsCell{12}); % knots
                speed(iter) = speedkn./1.944;   % knots 2 m/s
                clear('speedkn')
            else
                heading(iter) = NaN;
                speed(iter) = NaN;
            end
        end
    end  
    % GPGSA
    if strcmp(gpsCell{5},'$GPGSA')
        isGSA = 1;
        isGPS = gpsCell{7}~=1;
        if isGPS
            % Get Dilution of Precision
            pdop(iter) = str2double(gpsCell{length(gpsCell)-2});
            hdop(iter) = str2double(gpsCell{length(gpsCell)-1});
            tmpvdop = split(gpsCell{length(gpsCell)},'*');
            vdop(iter) = str2double(tmpvdop{1});
            clear('tmpvdop')
        else
            pdop(iter) = NaN;
            hdop(iter) = NaN;
            vdop(iter) = NaN;
        end
    end
% Completed Read Iteration    
end
%% Write GPS Matrix
% Format:
% [trcno,unixTime,longitude,latitude,elevation,geiodN,heading,speed,pdop,hdop,vdop];

% Get Trace Numbers
% [~,sortIx] = unique(trcno);
% Patch for Corrupt or Bad Files
% badIx = find(diff(sortIx)<1,1);
badIx = find(diff(trcno)<0,1)+1;

if ~isempty(badIx)
    errorstr = [filename,' is corrupted..!'];
    warning(errorstr)
%     trcno(sortIx(badIx):end) = [];
    trcno(badIx:end) = [];
    [trcno,trcIx] = unique(trcno);
    trcno(end) = []; % pop off the end for padding just bc?
% end
else
[trcno,trcIx] = unique(trcno);
end
ntrcs = length(trcno);
if length(trcno) ~= length(unixTime)
    warning('GP2:Corrupted','Something is not right')
%     keyboard
end
% Allocate GPS Matrix
GPS = nan(numel(trcno),11);
% Put trcNo
GPS(:,1) = trcno(:);
% Check that we have all the variables
if exist('unixTime','var')
    % Put unixTime
    %     GPS(:,2) = unixTime(:);
    % Error Checking
    if numel(unixTime)==numel(trcno)
        GPS(:,2) = unixTime(1:ntrcs)';
    else
        GPS(:,2) = interp1(1:numel(unixTime),unixTime,1:numel(trcno),'linear','extrap');
    end
end
if exist('longitude','var')
    % Put longitude
    %     GPS(:,3) = longitude(:);
    %     GPS(:,3) = longitude(1:ntrcs)';
    % Error Checking
    if numel(longitude)==numel(trcno)
        GPS(:,3) = longitude(1:ntrcs)';
    else
        GPS(:,3) = interp1(1:numel(longitude),longitude,1:numel(trcno),'linear','extrap');
    end
end
if exist('latitude','var')
    % Put latitude
%     GPS(:,4) = latitude(:);
%     GPS(:,4) = latitude(1:ntrcs)';
    % Error Checking
    if numel(latitude)==numel(trcno)
        GPS(:,4) = latitude(1:ntrcs)';
    else
        GPS(:,4) = interp1(1:numel(latitude),latitude,1:numel(trcno),'linear','extrap');
    end
end
if exist('elevation','var')
    % Put elevtion
%     GPS(:,5) = elevation(:);
%     GPS(:,5) = elevation(1:ntrcs)';
    % Error Checking
    if numel(elevation)==numel(trcno)
        GPS(:,5) = elevation(1:ntrcs)';
    else
        GPS(:,5) = interp1(1:numel(elevation),elevation,1:numel(trcno),'linear','extrap');
    end
end
if exist('geoidN','var')
    % Put geiodN
%     GPS(:,6) = geoidN(:);
%     GPS(:,6) = geoidN(1:ntrcs)';
    % Error Checking
    if numel(geoidN)==numel(trcno)
        GPS(:,6) = geoidN(1:ntrcs)';
    else
        GPS(:,6) = interp1(1:numel(geoidN),geoidN,1:numel(trcno),'linear','extrap');
    end
end
if exist('heading','var')
    % Put heading
%     GPS(:,7) = heading(:);
%     GPS(:,7) = heading(1:ntrcs)';
    % Error Checking
    if numel(heading)==numel(trcno)
        GPS(:,7) = heading(1:ntrcs)';
    else
        GPS(:,7) = interp1(1:numel(heading),heading,1:numel(trcno),'linear','extrap');
    end
end
if exist('speed','var')
    % Put speed
%     GPS(:,8) = speed(:);
%     GPS(:,8) = speed(1:ntrcs)';
    % Error Checking
    if numel(speed)==numel(trcno)
        GPS(:,8) = speed(1:ntrcs)';
    else
        GPS(:,8) = interp1(1:numel(speed),speed,1:numel(trcno),'linear','extrap');
    end
end
if exist('pdop','var')
    % Put pdop
%     GPS(:,9) = pdop(:);
%     GPS(:,9) = pdop(1:ntrcs)';
    % Error Checking
    if numel(pdop)==numel(trcno)
        GPS(:,9) = pdop(1:ntrcs)';
    else
        GPS(:,9) = interp1(1:numel(pdop),pdop,1:numel(trcno),'linear','extrap');
    end
end
if exist('hdop','var')
    % Put hdop
%     GPS(:,10) = hdop(:);
%     GPS(:,10) = hdop(1:ntrcs)';
    % Error Checking
    if numel(hdop)==numel(trcno)
        GPS(:,10) = hdop(1:ntrcs)';
    else
        GPS(:,10) = interp1(1:numel(hdop),hdop,1:numel(trcno),'linear','extrap');
    end
end
if exist('vdop','var')
    % Put vdop
%     GPS(:,11) = vdop(:);
%     GPS(:,11) = vdop(1:ntrcs)';
    % Error Checking
    if numel(vdop)==numel(trcno)
        GPS(:,11) = vdop(1:ntrcs)';
    else
        GPS(:,11) = interp1(1:numel(vdop),vdop,1:numel(trcno),'linear','extrap');
    end
end

% End of Function
end

