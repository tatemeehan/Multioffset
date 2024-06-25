%% Quick Pick Airwave Calibration
%clear; close all; clc;
addpath(genpath('C:\Users\RDCRLTGM\Desktop\git-repository\Multioffset'))
%% Read .nc file
% Establish Directories and Files
workingDirectory = pwd;
% Enter Data Directory
% directories = {'F:\500MHz_mo_2024\JIF_MO\MT_divide2_2024\raw\Lineset'};
directories = {['F:\500MHz_mo_2024\JIF_MO\C18 _2024\raw\airwave'],['F:\500MHz_mo_2024\JIF_MO\DEM1_2024\raw\airwave']...
    ['F:\500MHz_mo_2024\JIF_MO\MT_divide_2024\raw\airwave'],['F:\500MHz_mo_2024\JIF_MO\MT_divide2_2024\raw\airwave']...
    ['F:\500MHz_mo_2024\JIF_MO\NWBN_2024\raw\airwave'],['F:\500MHz_mo_2024\JIF_MO\NWBN2_2024\raw\airwave']};
% Enter Line Numbers
Lines = {[1,3],[1,10],[1],[2],[1],[1]};
nChan = 16;
chans = 1:nChan;
% Allocate Radargram
RAD = [];
distanceAxe = [];
isPlot = 0;
iter = 1;
for ff = 3:numel(directories)
% Extract Common Offset from .nc
for ii = 1:numel(Lines{ff})
    filename = ['Line',num2str(Lines{1}(ii)),'.nc'];
    filepath = fullfile(directories{1},filename);
    % Read netCDF data file
    disp(' ')
    fprintf('Reading .nc File \n')
    tic
    ncRad = ncread(filepath,'DATA');
    MxTrhd{ii} = ncRad(1:40,:);
    MxData{ii} = ncRad(41:end,:);
        % Radar Parameters
    TWT = [0:MxTrhd{ii}(12,1):(MxTrhd{ii}(13,1)-1).*MxTrhd{ii}(12,1)];
    dt = MxTrhd{ii}(12,1);
    f0 = MxTrhd{ii}(11,1);
    dx = median(diff(distanceAxe));
    for jj = 1:nChan
% Channel
channel = chans(jj);
    offset = MxTrhd{ii}(4,channel);

    % Extract Common Offset Data
    RAD = [MxData{ii}(:,channel:nChan:end)];
    % Quick Process
    [ RAD ] = processCommonOffset( RAD, f0, dt, TWT, offset, distanceAxe, dx, channel );
    % Airwave Picking
    [pickX,pickT] = polarPickerT8( {RAD} );
    % Layman's Deconvolution
    prompt = ['*',blanks(2),'Perform Layman`s Deconvolution?',blanks(2),'*'];
    titleTxt = ['Layman`s Deconvolution'];
    laymanQuest = questdlg(prompt, titleTxt,'Yes','No','Yes');
    if strcmp(laymanQuest,'Yes')
        laymanDecon = inputdlg('  Number of Samples to Subtract from Pick  :',...
            'Layman`s Deconvolution', 1,{'10'});
        laymanDecon = str2num(cell2mat(laymanDecon));
    else
        laymanDecon = 0;
    end
    FBpick(jj) = median(pickT{1}-laymanDecon);
%     % Parameters
%     R = [100,100,250,250,250,250,250,250,250,250,250,250,250,250,250,250];
%     merR = R(jj);   %Rank of MER window [ns]
%     powMER = 3; % Power of MER operation
%     for kk = 1:length(RAD(1,:))
%         % First Break Picking
%         tmp = RAD(:,kk);
%         [CommonFBix, ~, ~] = wong_mer(tmp, merR, dt, powMER);
%         CommonFBPick(kk) = TWT(CommonFBix);
%     end
%     
%     % Travel Time Correction is the Mean Static Adjustment
%     % FBpick = mean(CommonFBPick);
%     FBpick(jj) = median(CommonFBPick);
    tpick(jj) = FBpick(jj).*dt;
    c = 0.3; % Airwave Velocity
    tcal = offset./c;
    tshift(jj) = tcal-tpick(jj);
    shiftIx(jj) = tshift(jj)./dt;
    offsetArray(jj) = offset;
    if isPlot
        figure();imagesc(RAD);colormap(bone);title(num2str(offset))%(length(directories{ii})+2:end))
    hold on; plot(1:length(RAD(1,:)),ones(length(RAD(1,:)),1).*FBpick(jj),'m')
    end
    end
    
    tmpFile = strsplit(directories{ff},'\');
    AirwaveCalibration.project{iter} = tmpFile{4};
    AirwaveCalibration.line{iter} = ['Line',num2str(Lines{ff}(ii))];
    AirwaveCalibration.channel{iter} = 1:16;
    AirwaveCalibration.offset{iter} = offsetArray;
    AirwaveCalibration.tPick{iter} = tpick;
    AirwaveCalibration.ixPick{iter} = FBpick;
    AirwaveCalibration.tShift{iter} = tshift;
    AirwaveCalibration.ixShift{iter} = shiftIx;
    AirwaveCalibration.laymanDecon{iter} = laymanDecon;
    toc
    iter = iter+1;
    save("F:\500MHz_mo_2024\JIF_MO\AirwaveCalibration.mat","AirwaveCalibration")
end
end
