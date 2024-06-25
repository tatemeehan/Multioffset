function [dat,hdr1,trhd,dt,f0,xArray,dx,date,GPS,delta,latency] = rSS(filename)
%
% USAGE: [dat,hdr1,trhd,dt,f0,xArray] = readSensorsSoftwareData(filename)
%
% Description: Read a sensors and software GPR file, both hdr and dt1
% components. The data matrix is returned as well as the header
% information.
%
% INPUT:
%   filename = full filename including the path if not in current directory
%             * Exclude the file Extension
%   example: \my_data\LINE00
% OUTPUT:
%   dat    = data matrix ( row=time, column=space)
%   hdr1   = survey geometry headers
%   trhd   = trace headers
%   dt     = sample interval [ns]
%   f0     = antenna central frequency
%   xArray = position vector corresponding to columns in dat
%   dx     = Not Based on GPS
%   date   = Not Based on GPS
%   GPS    = [trcno,unixTime,longitude,latitude,elevation,geiodN,heading,...
%           speed,pdop,hdop,vdop];
%   delta  = [X,Y,Z] relative position of GPS antenna [m]
%   latency= the GPS position to system latency time [s] - unsused 
%
% Written by: Dylan Mikesell (dylanmikesell@boisestate.edu)
% Last modified: 2 November 2016
% By: Tate Meehan - dt calculation was edited
% Original: dt  = t / ( ns - 1 );
% Modified: dt  = t / ( ns );
% Based on code from John Bradford
%
% V1: wrote code to read 500 MHz demultiplexed data
% V2: rewrote code to handle non-demultiplexed 100 MHz data as well. This
% is a more robust function now but assumes that headers are the same
% within the two data types.

% Ensure that filename excludes an extesion
tmpfilename = split(filename,'.');
filename = tmpfilename{1};
clear('tmpfilename')
%--------------------------------------------------------------------------
% Load the HEADER file
hdfile = [filename '.HD'];

if ~exist(hdfile,'file')
    disp(hdfile)
    error('MATLAB/readSensorsSoftwareData: HEADER file does not exist.');
end

fid = fopen(hdfile,'r');
hdr1 = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
hdr1 = hdr1{1}; % a cell array of strings, one per line in the file.

% Find all the blank lines using cellfun and remove
hdr1( cellfun( @isempty, hdr1 ) ) = [];
% Grab Date Stamp
date = hdr1{3};
% remove the first 3 header lines
hdr1(1:3) = [];

for ii = 1 : numel(hdr1) 
    C{ii} = strtrim( strsplit(hdr1{ii},'=') ); % separate variables and trim whitespace
end

% Find the parts we need
ntr = round( str2double( C{1}{1,2} ) ); % number of traces 
ns  = round( str2double( C{2}{1,2} ) ); % number of sample per trace 
t   = str2double( C{4}{1,2} ); % [nanosecond] t_max 
dt  = t / ( ns ); % [nanosecond] time sample interval

xstart = str2double( C{5}{1,2} ); % [m] start of survey
xstop  = str2double( C{6}{1,2} ); % [m] end of survey
dx     = str2double( C{7}{1,2} ); % [m] sample interval in space
xArray = xstart : dx : xstop; % [m] position vector

f0 = round( str2double( C{9}{1,2} ) ); % [MHz] frequency 

if strcmp(C{15}{1}(1:9),'ELEVATION')
    % DATA TYPE field D.N.E. in O.G. DVL .hd
    try % GPS fix
        datatype = C{18}{1,2}; % Data Byte Type
    catch
        datatype = 'I*2';
    end
else
    try % No GPS
        datatype = C{16}{1,2}; % Data Byte Type
    catch
        datatype = 'I*2';
    end
end

%--------------------------------------------------------------------------
% Load the DATA file
dtfile = [filename '.DT1'];

if ~exist(dtfile,'file')
    disp(dtfile)
    error('MATLAB/readSensorsSoftwareData: DATA file does not exist.');
end

% allocate
dat  = zeros( ns, ntr );
trhd = zeros( 25, ntr );

fid=fopen(dtfile,'r','n');
k = 0;
while (k < ntr)
    trhd(1:22,k+1) = fread( fid, 22, 'real*4' ); % read trace header
    fread( fid, 1, 'integer*1' ); % read blank 
    trhd(23,k+1)    = fread(fid,1,'integer*1'); % read trace header
    fread(fid,2,'integer*1'); % read blank
    trhd(24:25,k+1) = fread(fid,2,'real*4'); % read trace header
    comment         = fread(fid,28); % read trace 'comment'
    if strcmp(datatype,'I*2')
        % Integer Data
        dat(:,k+1)      = fread(fid,ns,'int16'); % read data
    elseif strcmp(datatype,'F*4')
        % Floating Point Data
        dat(:,k+1)      = fread(fid,ns,'real*4'); % read data
    else
        error('unknown data type!')
    end
    k = k + 1; % upate counter
end
fclose(fid);
% Ensure Frequency in Header [MHz]
trhd(9,:) = f0;
% Ensure dt in Header [ns]
trhd(7,:) = dt;
% Ensure nsamps in Header
trhd(3,:) = ns;
%--------------------------------------------------------------------------
% Read GPS Data

% Check File Type
isGPS = exist([filename,'.gps'],'file') == 2;
isGP2 = exist([filename,'.gp2'],'file') == 2;
if ~isGPS && ~isGP2
    warning('No GPS data!')
%     error('GPS data is corrupt.')
else
    if isGPS % Read .GPS file
        [GPS, delta, latency] = readGPS([filename,'.gps'],date);
    end
    if isGP2 % Read .GP2 file
        [GPS, delta, latency] = readGP2([filename,'.gp2'],date);
    end
end


fprintf('Done reading %s.\n',filename);

return