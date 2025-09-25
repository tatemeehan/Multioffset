%% Airwave Instantaneous Frequency Calibration
clear; close all; clc;
addpath(genpath('C:\Users\RDCRLTGM\Desktop\git-repository\Multioffset'))
%% Read .nc file
% Establish Directories and Files
workingDirectory = pwd;
% Enter Data Directory
% directories = {'F:\500MHz_mo_2024\JIF_MO\MT_divide2_2024\raw\Lineset'};
directories = {['E:\JIF24\JIF_MO\C18 _2024\raw\airwave'],['E:\JIF24\JIF_MO\DEM1_2024\raw\airwave']...
    ['E:\JIF24\JIF_MO\MT_divide_2024\raw\airwave'],['E:\JIF24\JIF_MO\MT_divide2_2024\raw\airwave']...
    ['E:\JIF24\JIF_MO\NWBN_2024\raw\airwave'],['E:\JIF24\JIF_MO\NWBN2_2024\raw\airwave']};
load('E:\JIF24\JIF_MO\AirwaveCalibration.mat');
% Enter Line Numbers
Lines = {[1,3],[1,10],[1],[2],[1],[1]};
nChan = 16;
chans = 1:nChan;
% Allocate Radargram
RAD = [];
distanceAxe = [];
isPlot = 1;
iter = 1;
for ff = 1:numel(directories)
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
            instantFreq = ins_freq(RAD,dt,10,2,1);
            filtInstantFreq = medfilt2(instantFreq,[10,50]);
            initFreq = mean(filtInstantFreq(AirwaveCalibration.ixPick{ii}(jj):AirwaveCalibration.ixPick{ii}(jj)+30,:));
            initFreq2 = max(filtInstantFreq(AirwaveCalibration.ixPick{ii}(jj):AirwaveCalibration.ixPick{ii}(jj)+30,:));
            initFreq3 = min(filtInstantFreq(AirwaveCalibration.ixPick{ii}(jj):AirwaveCalibration.ixPick{ii}(jj)+30,:));
            initFreq4 = median((filtInstantFreq(AirwaveCalibration.ixPick{ii}(jj):AirwaveCalibration.ixPick{ii}(jj)+30,:)));
            % Smooth Freq
            L = 11;
            H = hamming(L);
            % SMooth intiFreq
            initFreq = conv(initFreq,H./sum(H),"same");
            initFreq(1:L) = mean(initFreq); initFreq(numel(initFreq)-L:end)=mean(initFreq);
            initFreq = mean(initFreq);
            fmean(jj) = initFreq;
            % SMooth intiFreq
            initFreq2 = conv(initFreq2,H./sum(H),"same");
            initFreq2(1:L) = mean(initFreq2); initFreq2(numel(initFreq2)-L:end)=mean(initFreq2);
            initFreq2 = mean(initFreq2);
            fmax(jj) = initFreq2;
            % SMooth intiFreq
            initFreq3 = conv(initFreq3,H./sum(H),"same");
            initFreq3(1:L) = mean(initFreq3); initFreq3(numel(initFreq3)-L:end)=mean(initFreq3);
            initFreq3 = mean(initFreq3);
            fmin(jj) = initFreq3;
            % SMooth intiFreq
            initFreq4 = conv(initFreq4,H./sum(H),"same");
            initFreq4(1:L) = mean(initFreq4); initFreq4(numel(initFreq4)-L:end)=mean(initFreq4);
            initFreq4 = mean(initFreq4);
            fmedian(jj) = initFreq4;
            % Quick Process
            %     [ RAD ] = processCommonOffset( RAD, f0, dt, TWT, offset, distanceAxe, dx, channel );

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
            %     tpick(jj) = FBpick(jj).*dt;
            %     c = 0.3; % Airwave Velocity
            %     tcal = offset./c;
            %     tshift(jj) = tcal-tpick(jj);
            %     shiftIx(jj) = tshift(jj)./dt;
                offsetArray(jj) = offset;
            if isPlot
                figure();imagesc(filtInstantFreq);colormap(bone);colorbar;title(num2str(offset))%(length(directories{ii})+2:end))
                hold on;plot(ones(size(RAD,2),1).*AirwaveCalibration.ixPick{ii}(jj),'m')            
            end
        end

        tmpFile = strsplit(directories{ff},'\');
        AirwaveFrequency.project{iter} = tmpFile{4};
        AirwaveFrequency.line{iter} = ['Line',num2str(Lines{ff}(ii))];
        AirwaveFrequency.channel{iter} = 1:16;
        AirwaveFrequency.offset{iter} = offsetArray;
        AirwaveFrequency.fmean{iter} = fmean;
        AirwaveFrequency.fmedian{iter} = fmedian;
        AirwaveFrequency.fmin{iter} = fmin;
        AirwaveFrequency.fmax{iter} = fmax;
        toc
        iter = iter+1;
        save("E:\JIF24\JIF_MO\AirwaveFrequency.mat","AirwaveFrequency")
    end
end
