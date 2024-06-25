%% Read MxData
clear; %close all; clc;
addpath(genpath('C:\Users\RDCRLTGM\Desktop\git-repository\Multioffset'))
%% Read .nc file
% Establish Directories and Files
workingDirectory = pwd;
% Enter Data Directory
% directories = {'F:\500MHz_mo_2024\JIF_MO\MT_divide2_2024\raw\Lineset'};
directories = {'F:\500MHz_mo_2024\JIF_MO\NWBN2_2024\raw\Lineset'};
% Enter Line Numbers
Lines = {[2,4,5,6,7,8]};
% Channel
channel = 12;
% Allocate Radargram
RAD = [];
distanceAxe = [];
% Extract Common Offset from .nc
for ii = 1:numel(Lines{1})
filename = ['Line',num2str(Lines{1}(ii)),'.nc'];
    filepath = fullfile(directories{1},filename);
    % Read netCDF data file
    disp(' ')
    fprintf('Reading .nc File \n')
    tic
    ncRad = ncread(filepath,'DATA');
    MxTrhd{ii} = ncRad(1:40,:);
    MxData{ii} = ncRad(41:end,:);
    % Concatenate Data
    RAD = [RAD,MxData{ii}(:,channel:16:end)];
    distanceAxe = [distanceAxe,MxTrhd{ii}(18,1:16:end)];
    toc
end
%% Make some Radargrams

% Radar Parameters
offset = MxTrhd{1}(4,channel);
TWT = [0:MxTrhd{1}(12,1):(MxTrhd{1}(13,1)-1).*MxTrhd{1}(12,1)];
dt = MxTrhd{1}(12,1);
f0 = MxTrhd{1}(11,1);
dx = median(diff(distanceAxe));
% Quick Process
[ RAD ] = processCommonOffset( RAD, f0, dt, TWT, offset, distanceAxe, dx, channel );

% Colormap
cmap = csvread('C:\Users\RDCRLTGM\Desktop\git-repository\Multioffset\colormaps\radarlove.csv');

% Radagram
figure();imagesc(distanceAxe./1000,TWT,RAD);colormap(cmap);
xlabel('Distance (km)')
ylabel('Two-Way Travel-time (ns)')
title(['Offset: ',num2str(offset), ' m']);
