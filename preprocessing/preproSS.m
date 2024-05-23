%% Preprocess Sensors&Software
clear; close all; clc;
%% Establish Directories and Files
workingDirectory = pwd;
% Enter Data Directory
directories = {'D:\alicia\GPR Data - Share w T8'};
% Enter Line Numbers
Lines = {[1]};
% Controls
isWrite = 1; % Write netCDF file
isFID = 1;   % Create Unique Trace ID for each file
isKillChan = 0; % Remove Bad Channels
%% Preprocess
for ff = 1:length(directories)
dataDir = directories{ff};
writeDir = dataDir;
fileNames = dir([dataDir,'\*.DT1']);
lineNo = Lines{ff};
nFiles = length(lineNo);
MxData = cell(nFiles,1);
MxTrhd = MxData;
tmpLineNo = zeros(length(fileNames),1);
% Read All Files in Directory
% getfile index properly
for kk = 1:length(fileNames)
    tmpname = fileNames(kk).name;
    splitname = strsplit(tmpname,'.');
    splitname = strsplit(splitname{1},'-');
    if length(splitname)>=1
    tmpLineNo(kk) = str2num(splitname{1}(5:end));
    else
        error('File Cannot be Loaded, Check the Naming Convention.')
    end
end

for ii = 1:length(lineNo)
fileIx = find(tmpLineNo == lineNo(ii));
lineNames = fileNames(fileIx);
for kk=1:length(lineNames)
    tmpname = fileNames(fileIx(kk)).name;
    splitname = strsplit(tmpname,'.');
    splitname = strsplit(splitname{1},'-');
    if length(splitname)>=2
    ChanNo{ii}(kk) = str2num(splitname{2}(3:end));
    else
        ChanNo{ii}(kk) = 1;
        SqNo{ii}(kk) = 1;
    end
    if length(splitname)>=3
        SqNo{ii}(kk) = str2num(splitname{3}(3:end));
    else
        SqNo{ii}(kk) = 1;
    end
end
nChan(ii) = length(unique(ChanNo{ii})).*length(unique(SqNo{ii}));
for jj = 1:length(fileIx)
    filename = fileNames(fileIx(jj)).name;
    tmpfilename = [fileNames(jj).folder,'\',filename(1:end-4)];
        [dat{jj,ii},hdr1{jj,ii},trhd{jj,ii},~,~,~,~,~,GPS{jj,ii},delta(jj,:),latency(jj)] = rSS(tmpfilename);
    % ntraces
    ntrc(jj) = size(trhd{jj,ii},2);
    % Tx Rx Geometry
    Geometry(ii).tx(jj) = trhd{jj,ii}(18,1);
    Geometry(ii).ty(jj) = trhd{jj,ii}(19,1);
    Geometry(ii).rx(jj) = trhd{jj,ii}(15,1);
    Geometry(ii).ry(jj) = trhd{jj,ii}(16,1);
end
mintrc = min(ntrc);
clear('ntrc')
nsamp = size(dat{1},1);
% Determine Offsets
Geometry(ii).offset = sqrt((Geometry(ii).tx-Geometry(ii).rx).^2+(Geometry(ii).ty-Geometry(ii).ry).^2);
% Sort Channels by Offset
[Geometry(ii).offset,Geometry(ii).sortIx] = sort(Geometry(ii).offset);
Geometry(ii).tx = Geometry(ii).tx(Geometry(ii).sortIx);
Geometry(ii).rx = Geometry(ii).rx(Geometry(ii).sortIx);
Geometry(ii).ty = Geometry(ii).ty(Geometry(ii).sortIx);
Geometry(ii).ry = Geometry(ii).ry(Geometry(ii).sortIx);
%Midpoint
Geometry(ii).midpointx = mean([Geometry(ii).tx;Geometry(ii).rx],1);  % Midpoint Locations Relative to 0
Geometry(ii).midpointy = mean([Geometry(ii).ty;Geometry(ii).ry],1);  % Midpoint Locations Relative to 0
if exist('GPS','var')
% Sort Trace Headers, Data, and GPS by offset
for kk = 1:length(fileIx)
tmptrhd{kk,1} = trhd{Geometry(ii).sortIx(kk),ii};
tmpdat{kk,1} = dat{Geometry(ii).sortIx(kk),ii};
tmpGPS{kk,1} = GPS{Geometry(ii).sortIx(kk),ii};
% Correct GPS trace numbers
tmpIx = find(tmpGPS{kk,1}(:,1)>mintrc);
% Remove Extraneous Trace Numbers
tmpGPS{kk,1}(tmpIx,:) = [];
% Renumber the Traces of each channel
tmpIx = find(tmpGPS{kk,1}(:,1)>=1);
tmpGPS{kk,1}(tmpIx,1) = tmpGPS{kk,1}(tmpIx,1).*nChan(ii);
tmpGPS{kk,1}(:,1) = (tmpGPS{kk,1}(:,1)+(kk-1));
end
for kk = 1:length(fileIx)
trhd{kk,ii} = tmptrhd{kk,1};
dat{kk,ii} = tmpdat{kk,1};
GPS{kk,ii} = tmpGPS{kk,1};
end
clear('tmpGPS','tmptrhd','tmpdat','tmpIx')
if size(delta,2)==3 % .gp2 file loaded
    % Change Sign of delta
    Geometry(ii).delta = delta(Geometry(ii).sortIx,:);
    Geometry(ii).latency = latency(Geometry(ii).sortIx);
    % GPS Location :eyeroll
    Geometry(ii).gps = [median(Geometry(ii).delta(:,1)+Geometry(ii).midpointx(:)),...
        median(Geometry(ii).delta(:,2)+Geometry(ii).midpointy(:)),median(Geometry(ii).delta(:,3))];
    % Recalculate Delta :eyeroll
    Geometry(ii).delta = Geometry(ii).gps ...
        - [Geometry(ii).midpointx(:),Geometry(ii).midpointy(:),zeros(nChan(ii),1)];
else % .gps file loaded
    Geometry(ii).latency = zeros(1,nChan(ii));
    Geometry(ii).gps = [0,0,0]; % Enter GPS Antenna Position
    % Calculate Delta
    Geometry(ii).delta = Geometry(ii).gps ...
        - [Geometry(ii).midpointx(:),Geometry(ii).midpointy(:),zeros(nChan(ii),1)];
    warning(['GPS Antenna Position Unknown. Assuming ', num2str(Geometry(ii).gps)])
end

for jj = 1:length(fileIx)
    trcnoBad(jj) = any(diff(trhd{jj,ii}(1,:))~=1);
    trcNo(jj) = length((trhd{jj,ii}(1,:)));
    badNos{jj} = find(diff(trhd{jj,ii}(1,:))~=1)+1;
end
if any(trcnoBad)
    warning('Trace Numbers Automatically Repaired.')
    for jj = 1:length(fileIx)
        % Remove Corrupt Traces
        trhd{jj,ii}(:,badNos{jj}) = [];
        dat{jj,ii}(:,badNos{jj}) = [];
        % ntraces
        ntrc(jj) = size(trhd{jj,ii},2);
        % Pull out Bad Traces & Create New Trace Numbers
        trhd{jj,ii}(1,:) = 1:length(trhd{jj,ii}(1,:));
    end
        mintrc = min(ntrc);
end

for jj = 1:length(fileIx)
    % Add Channel Number to Header
    trhd{jj,ii}(23,:) = jj;
    % Multiplex The Files :eyeroll
    MxTrhd{ii}(:,jj:nChan(ii):mintrc.*nChan(ii)) = trhd{jj,ii}(:,1:mintrc);
    MxData{ii}(:,jj:nChan(ii):mintrc.*nChan(ii)) = dat{jj,ii}(:,1:mintrc);
end

% Clear Undefined Loop Variables
clear('mintrc','nsamp','ntrc','delta','latency')
%% GeoReference
velocityThreshold = 1;
calcVelocityHeading = 1;
InterpolateTrhd = 1;
dx = 0.1;
rmvStaticTrc = 1;
[MxTrhd{ii},MxData{ii}] = trhdGeoref(MxTrhd{ii},MxData{ii},GPS{1,ii},Geometry(ii),velocityThreshold,calcVelocityHeading,InterpolateTrhd,dx,rmvStaticTrc);%,trcnoBad);

%% Create Unique File ID
% For all files in Queue Create unique trace numbers, shot numbers, and/or
% distances
if isFID
    if ii > 1
        % Extract End of Previous File
        maxtrc = max(MxTrhd{ii-1}(1,:));
        maxshot = max(MxTrhd{ii-1}(2,:));
        maxshotdist = max(MxTrhd{ii-1}(18,:));
        maxmiddist = max(MxTrhd{ii-1}(25,:));       
        % Increment File by last EoF
        MxTrhd{ii}(1,:) = MxTrhd{ii}(1,:)+maxtrc;
        MxTrhd{ii}(2,:) = MxTrhd{ii}(2,:)+maxshot;
        MxTrhd{ii}(18,:) = MxTrhd{ii}(18,:)+maxshotdist...
            + mean(diff(MxTrhd{ii}(18,1:nChan(ii):end)));
        MxTrhd{ii}(25,:) = MxTrhd{ii}(25,:)+maxmiddist ...
            + mean(diff(MxTrhd{ii}(25,:)));
    end
    if ii == nFiles
        % Extract End of Previous File
        maxtrc = max(MxTrhd{ii}(1,:));
        maxshot = max(MxTrhd{ii}(2,:));
        maxshotdist = max(MxTrhd{ii}(18,:));
        maxmiddist = max(MxTrhd{ii}(25,:));   
    end
    % 
    if ff > 1 & ii == 1
        % Increment File by last EoF
        MxTrhd{ii}(1,:) = MxTrhd{ii}(1,:)+maxtrc;
        MxTrhd{ii}(2,:) = MxTrhd{ii}(2,:)+maxshot;
        MxTrhd{ii}(18,:) = MxTrhd{ii}(18,:)+maxshotdist...
            + mean(diff(MxTrhd{ii}(18,1:nChan(ii):end)));
        MxTrhd{ii}(25,:) = MxTrhd{ii}(25,:)+maxmiddist ...
            + mean(diff(MxTrhd{ii}(25,:)));
    end
end
%% Remove Bad Channels
if isKillChan
    badChan = []; % Enter Channel Numbers
    badIx = ismember(MxTrhd{ii}(3,:),badChan);
    MxTrhd{ii}(:,badIx) = [];
    MxData{ii}(:,badIx) = [];
end
%% Write netCDF
if isWrite
    ncname = strsplit(filename,'-');
    if length(ncname) == 1
        % Catch This
        ncname = strsplit(filename,'.');
    end
    ncname = ncname{1};
    writeNC(MxTrhd{ii},MxData{ii},ncname,writeDir,workingDirectory);
end
end
end
% clear('MxTrhd','MxData','GPS','Geometry')
clear('ChanNo','SqNo','nChan','dat','trhd','hdr1')
end