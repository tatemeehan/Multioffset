function [GPS,delta,latency] = readGPS(filename, OGdate)
%readGP2 extracts GPS data from NMEA string embedded in S&S .GP2 files.
%   Input :
%       filename = /directory/filename.GPS - the fullfile sting 
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
% Open the .GPS file
fid=fopen(filename);
% Read the whole lines of .GP2 file into a cell array
rawGPS=textscan(fid,'%s','Delimiter','\n');
% Close the file
fclose(fid);
% Reduce Dimension of Cell Array
rawGPS = rawGPS{1};
% remove any empty cells
rmvix = ~cellfun(@isempty,rawGPS);
rawGPS = rawGPS(rmvix);
% Get Indicies of NMEA Strings
nmeaStrings = ["$GPGGA","$GPVTG","$GPZDA","$GPRMC","$GPGSA"];
whichNMEA = zeros(length(rawGPS),numel(nmeaStrings));
for kk = 1:length(rawGPS)
    ix(kk) = strcmp(rawGPS{kk}(1),'$');
    GSAix(kk) = strcmp(rawGPS{kk}(1:6),'$GPGSA');
    % Supported NMEA Strings "$GPGGA","$GPVTG","$GPZDA","$GPRMC","$GPGSA"
    nmeas =(strcmp(rawGPS{kk}(1:6),[nmeaStrings]));
    rmvNMEA(kk) = ~any(nmeas);
    % Create Logical Matrix of NMEA
    whichNMEA(kk,:) = nmeas;
end
% Trace Number Rows
trcix = find(~ix); 
trcstr = rawGPS(trcix);
m = length(trcix);
trcno = zeros(m,1);
% Remove nonRecognized NMEAs
rawGPS(rmvNMEA,:) = [];
whichNMEA(rmvNMEA,:) = [];
% determine number of unique NMEAs & Which
NMEAs = nmeaStrings(sum(whichNMEA,1)>0);
nNMEA = numel(NMEAs);
% Sometimes $GPGSA is not recorded at 1 Hz we need to infill this
nmea1 = find((whichNMEA(1,:)));
nmea1Ix = find(whichNMEA(:,nmea1));
nnmea1 = numel(nmea1Ix);
if median(diff(nmea1Ix))~=nNMEA
    %%% Assumes $GPGSA is the only String needing infilled. %%%
    synthGPS = cell(nnmea1.*nNMEA,1);
    if any(strcmp(NMEAs,"$GPGGA"))
        gpggaIx = find(whichNMEA(:,1));
        gpgga1Ix = gpggaIx(1);
        synthGPGGAix = gpgga1Ix:nNMEA:nnmea1.*nNMEA;
        for ftp = 1:nnmea1
            synthGPS{synthGPGGAix(ftp)} = rawGPS{gpggaIx(ftp)};
        end
    end
    if any(strcmp(NMEAs,"$GPVTG"))
        gpvtgIx = find(whichNMEA(:,2));
        gpvtg1Ix = gpvtgIx(1);
        synthGPVTGix = gpvtg1Ix:nNMEA:nnmea1.*nNMEA;
        for ftp = 1:nnmea1
        synthGPS{synthGPVTGix(ftp)} = rawGPS{gpvtgIx(ftp)};
        end
    end
    if any(strcmp(NMEAs,"$GPZDA"))
        gpzdaIx = find(whichNMEA(:,3));
        gpzda1Ix = gpzdaIx(1);
        synthGPZDAix = gpzda1Ix:nNMEA:nnmea1.*nNMEA;
        for ftp = 1:nnmea1
        synthGPS{synthGPZDAix(ftp)} = rawGPS{gpzdaIx(ftp)};
        end
    end
    if any(strcmp(NMEAs,"$GPRMC"))
        gprmcIx = find(whichNMEA(:,4));
        gprmc1Ix = gprmcIx(1);
        synthGPRMCix = gprmc1Ix:nNMEA:nnmea1.*nNMEA;
        for ftp = 1:nnmea1
        synthGPS{synthGPRMCix(ftp)} = rawGPS{gprmcIx(ftp)};
        end
    end
    if any(strcmp(NMEAs,"$GPGSA"))
        gpgsaIx = find(whichNMEA(:,5));
        gpgsa1Ix = gpgsaIx(1);
        synthGPGSAix = gpgsa1Ix:nNMEA:nnmea1.*nNMEA;
        for ftp = 1:nnmea1
            [~,minIx] = min(abs(gpgsaIx - synthGPGSAix(ftp)));
            % Nearest Neighbor Infill
            synthGPS{synthGPGSAix(ftp)} = rawGPS{gpgsaIx(minIx)};
        end
    end
    rawGPS = synthGPS;
    clear('synthGPS');
end

delta = [0,0,0];
latency = [0,0,0];
ix = 0;
% Initialize String Condiditons
isGGA = 0;
isZDA = 0;
isVTG = 0;
isRMC = 0;
isGSA = 0;
for kk = 1:m
    % Extract Trace Number
    tmptrcstr = strsplit(trcstr{kk});
    trcno(kk) = str2num(tmptrcstr{2}(2:end));
    iter = kk;
    isGGA = 0;
    isZDA = 0;
    isVTG = 0;
    isRMC = 0;
    isGSA = 0;
for jj = 1:nNMEA
    ix = ix+1;
    % Extract GPS Data
    tmprawGPS = rawGPS{ix};
    % Add Dummy Cols to mimic .gp2
    tmprawGPS = ['1,1,1,1,',tmprawGPS];
    gpsCell = strsplit(tmprawGPS,',','CollapseDelimiters', false);
    % GPGGA
    % Lon, Lat, Elevation, Time
    if strcmp(gpsCell{5},'$GPGGA')
        isGGA = 1;
        isGPS = str2num(gpsCell{6})~=0; % GPS fix
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
            % Time of Day
            t = gpsCell{6};
            H = str2double(t(1:2));
            M = str2double(t(3:4));
            S = str2double(t(5:end));
            s = H.*3600+M.*60+S; % Seconds of Day
            % Initial Second of Day
                if iter == 1
                    s1 = s;
                end
                d = str2double(OGdate{3});
                % Check if Next Day
                if s < s1
                    d = d+1;
                end
                mnth = month(datetime(OGdate{2},'InputFormat','MM'));
                y = str2double(OGdate{1});
                dt = datetime([y mnth d H M S]);
                % Convert to unixTime
                unixTime(iter) = posixtime(dt); 
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
        geiodN(iter) = N;
        unixTime(iter) = s;
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
                % Latitude
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
                geiodN(iter) = N;
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
end
% Completed Read Iteration    
end
%% Write GPS Matrix
% Format:
% [trcno,unixTime,longitude,latitude,elevation,geiodN,heading,speed,pdop,hdop,vdop];

% Get Trace Numbers
trcno = unique(trcno);
% Allocate GPS Matrix
GPS = nan(numel(trcno),11);
% Put trcNo
GPS(:,1) = trcno(:);
% Check that we have all the variables
if exist('unixTime','var')
    % Put unixTime
    GPS(:,2) = unixTime(:);
end
if exist('longitude','var')
    % Put longitude
    GPS(:,3) = longitude(:);
end
if exist('latitude','var')
    % Put latitude
    GPS(:,4) = latitude(:);
end
if exist('elevation','var')
    % Put elevtion
    GPS(:,5) = elevation(:);
end
if exist('geiodN','var')
    % Put geiodN
    GPS(:,6) = geiodN(:);
end
if exist('heading','var')
    % Put heading
    GPS(:,7) = heading(:);
end
if exist('speed','var')
    % Put speed
    GPS(:,8) = speed(:);
end
if exist('pdop','var')
    % Put pdop
    GPS(:,9) = pdop(:);
end
if exist('hdop','var')
    % Put hdop
    GPS(:,10) = hdop(:);
end
if exist('vdop','var')
    % Put vdop
    GPS(:,11) = vdop(:);
end

% End of Function
end

