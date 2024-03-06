% 
clear;close all; clc;
%% Directories
directories = [{'E:\SnowEx23AK\GPR\UKT030823'},{'E:\SnowEx23AK\GPR\UKT030923'},...
    {'E:\SnowEx23AK\GPR\UKT031023'},{'E:\SnowEx23AK\GPR\UKT031123'},...
    {'E:\SnowEx23AK\GPR\UKT031223'},{'E:\SnowEx23AK\GPR\UKT031323'},...
    {'E:\SnowEx23AK\GPR\UKT031423'},{'E:\SnowEx23AK\GPR\UKT031523'}];
GPRcsv = [];
GPRll = [];
%% Process Traveltime Data
for hh = 1:size(directories,2)
    if hh == 1
    try GPR;
    catch
        GPR.MD.dataDir = directories{hh};
        GPR.MD.workDir = pwd;
        GPR.MD.Dir = dir([GPR.MD.dataDir,'/','*.nc']);
    end
    else
        GPR.MD.dataDir = directories{hh};
        GPR.MD.workDir = pwd;
        GPR.MD.Dir = dir([GPR.MD.dataDir,'/','*.nc']);
    end
% if ~isLoadMat
for ii = 1:length(GPR.MD.Dir)
GPR.MD.dataDir = directories{hh};
GPR.MD.workDir = pwd;
% GPR.MD.Dir = dir([GPR.MD.dataDir,'/','*.nc']);
GPR.MD.fileNames = GPR.MD.Dir(ii).name;
GPR.MD.lineNo = ii;%1:length(GPR.MD.fileNames);                   % Array of data "LINE" numbers
GPR.MD.nFiles = length(GPR.MD.fileNames);        % Number of Files
% GPR.Geometry.badChan = [1,2,3]; % Temporary Fix for Bad Channels


%% Read GPR Data
[GPR] = readGPRnc(GPR);
% Store X,Y,Z
GPRcsv = [GPRcsv;[GPR.Geolocation.X{1}(:),GPR.Geolocation.Y{1}(:),GPR.Geolocation.Z{1}(:)]];
GPRll = [GPRll;[GPR.Geolocation.Latitude{1}(:),GPR.Geolocation.Longitude{1}(:)]];
if ii == length(GPR.MD.Dir)
    maxDist(hh) = max(GPR.Geolocation.Distance{1});
end
end

end

%% Figure
figure();plot(GPRcsv(1:10:end,1)./1000,GPRcsv(1:10:end,2)./1000,'ok','markersize',1,'markerfacecolor','k')
title('Toolik MultiPolarization GPR')
xlabel('Easting (km)')
ylabel('Northing (km)')
daspect([1,1,1])
set(gca,'fontname','serif','fontsize',12,'fontweight','bold')
grid on; grid minor;

%% Export .csv
T = table(GPRll(:,1),GPRll(:,2),'VariableNames',{'Latitude','Longitude'});
writetable(T,'E:\SnowEx23AK\GPR\UKTmultipolarizationGPRcoords.csv')
