%% Multi-offset GPR Snow Analysis
clear; close all; clc;
%% Establish Directories and Files
% Data Directory
directories = {'E:\Keegan\NM-02-25-24\queue'};
% Paths
addpath(genpath('C:\Users\RDCRLTGM\Desktop\git-repository\Multioffset'))
%% Processing WorkFlow Controls
% Parallel Computing Enabled
isParallel = 1;
% Read Data
isTrimTWT = 0;          % Truncate Recorded Data
isReduceData = 0;       % Thin Traces
isEditPicks = 1;        % Edit Picks

% Write SWE Data
isWriteCSVhva = 0; % .csv output
isSaveMat = 0;
isSavePlots = 0;
% Load Processed Data
isLoadMat = 0;

% Load Color Maps
yetBlack = load('yetBlack.txt');

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
if ~isLoadMat
if isempty(GPR.MD.Dir)
    warning("Data Queue is Empty!")
end
for ii = 1:length(GPR.MD.Dir)
GPR.MD.dataDir = directories{hh};
GPR.MD.workDir = pwd;
GPR.MD.fileNames = GPR.MD.Dir(ii).name;
GPR.MD.lineNo = ii;% Array of data "LINE" numbers
GPR.MD.nFiles = 1; % Number of Files (Only handles one at a time)

%% Control Parallelization
if isParallel
    % Wake Parallel Computing
    if isempty(gcp('nocreate'))     
        nWorkers = 4;
        p = parpool(nWorkers);
    else nWorkers = 4;
    end
else
    nWorkers = 1;
end
GPR.MD.isParallel = isParallel;
GPR.MD.nWorkers = nWorkers;

%% Read GPR Data
[GPR] = readGPRnc(GPR);

%% Process GPR Data
display('Processing GPR Data')
rmNtrc = 0; padding = 0; isKillChan = 0; killArray = [];
[GPR] = processSWEPR(GPR,isTrimTWT,isReduceData,rmNtrc,padding,isKillChan,killArray);
%% Automatic Travel-time Picking
isAutoPickSurface = 0;
if isAutoPickSurface
    display('Auto Picking Snow Surface')
    [GPR] = autoPickSurface2(GPR,1);
    if GPR.Geometry.nChan{1} > 1
        [GPR] = autoEstimateSurface(GPR);
    end
end
%% Interperet Surface
isInterpretSurface = 0;
if isInterpretSurface
    display('Interpret the Surface Horizon')
    [GPR] = interpretSurface(GPR);
end
%% GPR on Surface
if isAutoPickSurface == 0 && isInterpretSurface == 0
    isOnSurface = 1;
    [GPR] = onSurface(GPR);
else
    isOnSurface = 0;
end
%% Flatten Surface Structure
isSOF = 0;
if isSOF
display('Flattening Surface')
[GPR] = flattenSurface(GPR);
%% Structure Oriented Filter
display('Performing Structure Oriented Filter')
isFXdeconvolution = 0;
isKuwahara = 1;
if isFXdeconvolution
[GPR] = fxDeconvolution(GPR);
end
if isKuwahara
[GPR] = kuwaharaFilter(GPR);
end
end
%% Auto Pick Ground
isAutoPickGround = 0;
if isAutoPickGround
    isQuickLook = 1;
        display('Auto Picking Ground Surface')
        [GPR] = autoPickGround3(GPR,isQuickLook);
end
%% Interpret Ground
isInterpretGround = 1;
if isInterpretGround
    display('Interpret the Ground Horizon')
    [GPR] = interpretGround(GPR);
end
%% QC Picks
if isEditPicks
    display('QC Picks')
    isLoadPicks = 1;
    [GPR] = pickEditorMx(GPR,isLoadPicks);
end
%% Horizon Velocity Analysis
isHVA = 1;
if isHVA
    display('Horizon Velocity Analysis')
[GPR] = HVA(GPR);
isRayTrace = 0;
if isRayTrace
% Raytracing Constrained Inversion
density = 300; % Approximate Average Snow Density
[GPR] = raytraceConstraint(GPR,density);
% Recompute HVA
[GPR] = HVA(GPR);
end
%% Estimate Snow Properties
display('Calculating Snow Properties')
[GPR] = estimateSnow(GPR,isOnSurface);
end
%% Save GPR.mat
if isSaveMat
    display('Writing .mat File')

    matname = strsplit(GPR.MD.Dir(ii).folder,'\');
    linename = strsplit(GPR.MD.fileNames,'.');
    linename = regexp(linename{1},'\d*','Match');
    matname = [matname{length(matname)-1},'-',linename{1}];

    cd(GPR.MD.dataDir)
    save([matname,'.mat'],'GPR','-v7.3');
    cd(GPR.MD.workDir)
end
%% Write .csv
% Snow Properties from HVA
if isWriteCSVhva
        display('Writing .csv File')
T = table(GPR.Geolocation.Longitude{1}(:),GPR.Geolocation.Latitude{1}(:),...
    GPR.Geolocation.Elevation{1}(:),GPR.Snow.Depth{1}(:),GPR.Snow.Density{1}(:),...
    GPR.Snow.Porosity{1}(:),GPR.Snow.SWE{1}(:));
T.Properties.VariableNames{'Var1'} = 'Longitude';
T.Properties.VariableNames{'Var2'} = 'Latitude';
T.Properties.VariableNames{'Var3'} = 'Elevation';
T.Properties.VariableNames{'Var4'} = 'Depth';
T.Properties.VariableNames{'Var5'} = 'Density';
T.Properties.VariableNames{'Var6'} = 'Porosity';
T.Properties.VariableNames{'Var7'} = 'SWE';
csvname = strsplit(GPR.MD.Dir(ii).folder,'\');
linename = strsplit(GPR.MD.fileNames,'.');
linename = regexp(linename{1},'\d*','Match');
csvname = [csvname{length(csvname)-1},'-',linename{1}];
cd(GPR.MD.dataDir)
writetable(T,[csvname,'.csv']);
cd(GPR.MD.workDir)
end
end
else
    % Load .mat
    cd(GPR.MD.dataDir)
    load('')
    cd(GPR.MD.workDir)
end
%% Create Figures
% Check if pcolor needed
isResampled = abs((GPR.Geolocation.Distance{1}(2)-GPR.Geolocation.Distance{1}(1))-...
        (GPR.Geolocation.Distance{1}(end)-GPR.Geolocation.Distance{1}(end-1))) <= 0.1; 
if length(GPR.Geolocation.Distance{1})<1000
    isResampled = 0; % use pcolor for smaller transects
end
titstr = strsplit(GPR.MD.fileNames,'.');
projstr = strsplit(GPR.MD.Dir(ii).folder,'\');
% Calculate RMS Distance
distRMS = (sqrt((GPR.Geolocation.X{1}(1:end-1)-GPR.Geolocation.X{1}(2:end)).^2+(GPR.Geolocation.Y{1}(1:end-1)-GPR.Geolocation.Y{1}(2:end)).^2));
% Calculate RMS Position
posRMS = sqrt(var(GPR.Geolocation.X{1})+var(GPR.Geolocation.Y{1}));
if posRMS < 5
    isStationary = 1;
    isResampled = 0;
end
% Make Figures
display('Displaying Figures')

% GPR TWT at 3 offsets    
    figure();
    subplot(3,1,1)
    if isResampled
        imagesc(GPR.Geolocation.Distance{1},GPR.D.TimeAxis{1},GPR.D.Radar{1}{1});
        colormap(bone); caxis([quantile(GPR.D.Radar{1}{1}(:),[0.005,0.995])]);
    else
    pcolor(GPR.Geolocation.Distance{1},GPR.D.TimeAxis{1},GPR.D.Radar{1}{1});
    shading interp;axis ij; colormap(bone);caxis([quantile(GPR.D.Radar{1}{1}(:),[0.005,0.995])]);
    end
    hold on; plot(GPR.Geolocation.Distance{1},GPR.D.surfaceTWT{1}(1,:),'m','linewidth',2)
    plot(GPR.Geolocation.Distance{1},GPR.D.groundTWT{1}(1,:),'m','linewidth',2);
    xlabel('Distance (m)'); ylabel('Travel-time (ns)');title(['Offset ',num2str(round(GPR.Geometry.offset{1}(1),2)),' (m)'])
    set(gca,'fontsize',12,'fontweight','bold','fontname','serif');

    subplot(3,1,2)
    if isResampled
        imagesc(GPR.Geolocation.Distance{1},GPR.D.TimeAxis{1},GPR.D.Radar{1}{2});
        colormap(bone); caxis([quantile(GPR.D.Radar{1}{2}(:),[0.005,0.995])]);
    else
    pcolor(GPR.Geolocation.Distance{1},GPR.D.TimeAxis{1},GPR.D.Radar{1}{2});
    shading interp;axis ij; colormap(bone);caxis([quantile(GPR.D.Radar{1}{2}(:),[0.005,0.995])]);
    end
    hold on; plot(GPR.Geolocation.Distance{1},GPR.D.surfaceTWT{1}(2,:),'m','linewidth',2)
    plot(GPR.Geolocation.Distance{1},GPR.D.groundTWT{1}(2,:),'m','linewidth',2);
    xlabel('Distance (m)'); ylabel('Travel-time (ns)');title(['Offset ',num2str(round(GPR.Geometry.offset{1}(2),2)),' (m)'])
    set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
    
    subplot(3,1,3)
    if isResampled
        imagesc(GPR.Geolocation.Distance{1},GPR.D.TimeAxis{1},GPR.D.Radar{1}{3});
        colormap(bone);caxis([quantile(GPR.D.Radar{1}{3}(:),[0.005,0.995])]);
    else
    pcolor(GPR.Geolocation.Distance{1},GPR.D.TimeAxis{1},GPR.D.Radar{1}{3});
    shading interp;axis ij; colormap(bone);caxis([quantile(GPR.D.Radar{1}{3}(:),[0.005,0.995])]);
    end
    hold on; plot(GPR.Geolocation.Distance{1},GPR.D.surfaceTWT{1}(2,:),'m','linewidth',2)
    plot(GPR.Geolocation.Distance{1},GPR.D.groundTWT{1}(3,:),'m','linewidth',2);
    xlabel('Distance (m)'); ylabel('Travel-time (ns)');title(['Offset ',num2str(round(GPR.Geometry.offset{1}(3),2)),' (m)'])
    set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
    % make sup title
    sgtitle([projstr{end-1},' ',titstr{1}],'fontsize',14,'fontweight','bold','fontname','serif')
    
    % Enlarge figure to full screen.
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
    
    % Save Figure
    if isSavePlots
        display('Saving Figure')
        saveas(gcf,[GPR.MD.dataDir,'\',[projstr{end-1},titstr{1}],'.png'])
    end

is2Dmap = 1; % 2D Map View
is3Dmap = 0; % X,Y,Z, map view

figure();
if is3Dmap
hmap = scatter3(GPR.Geolocation.X{1}./1000,GPR.Geolocation.Y{1}./1000,GPR.Geolocation.Z{1},5,GPR.D.groundT0{ii}-GPR.D.surfaceTWT{ii}(1,:),'filled');
end
if is2Dmap
hmap = scatter(GPR.Geolocation.X{1}./1000,GPR.Geolocation.Y{1}./1000,5,GPR.D.groundT0{ii}-GPR.D.surfaceTWT{ii}(1,:),'filled');
end
title([projstr{end-1},' ',titstr{1}])
xlabel('Easting (km)'); ylabel('Northing (km)');colormap(bone);hc = colorbar;
zlabel('Elevation (masl)');ylabel(hc,'Travel-Time (ns)');
caxis([quantile(GPR.D.groundT0{1}-GPR.D.surfaceTWT{1}(1,:),[0.01,0.99])])
set(hc,'fontsize',12,'fontweight','bold','fontname','serif');
set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
ax = ancestor(hmap, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.2f')
ax.YAxis.Exponent = 0;
ytickformat('%.2f')
daspect([1,1,1])
grid on; grid minor;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Save Map View TWT
if isSavePlots
    display('Saving Figure')
saveas(gcf,[GPR.MD.dataDir,'\',[projstr{end-1},titstr{1}],'map','.png'])
end

if isHVA
% Load Colormap
load('yetBlack.txt')
% Density
figure();
hsnow = scatter(GPR.Geolocation.X{ii}./1000,GPR.Geolocation.Y{ii}./1000,5,GPR.Snow.Density{ii});colormap(yetBlack)
colorbar;
caxis([quantile(GPR.Snow.Density{ii},[0.05,0.95])])
daspect([1,1,1])
title('Density (kg/m^3)')
xlabel('Easting (km)')
ylabel('Northing (km)')
set(gca,'fontsize',12,'fontweight','bold','fontname','serif')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
grid on; grid minor;
ax = ancestor(hsnow, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.2f')
ax.YAxis.Exponent = 0;
ytickformat('%.2f')
    if isSavePlots
        display('Saving Figure')
        saveas(gcf,[GPR.MD.dataDir,'\',[projstr{end-1},titstr{1}],'density.png'])
    end
% Porosity
figure();
hsnow = scatter(GPR.Geolocation.X{ii}./1000,GPR.Geolocation.Y{ii}./1000,5,GPR.Snow.Porosity{ii});colormap(yetBlack)
colorbar;
caxis([quantile(GPR.Snow.Porosity{ii},[0.05,0.95])])
daspect([1,1,1])
hold on;grid on; grid minor;
title('Porosity')
xlabel('Easting (km)')
ylabel('Northing (km)')
set(gca,'fontsize',12,'fontweight','bold','fontname','serif')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
ax = ancestor(hsnow, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.2f')
ax.YAxis.Exponent = 0;
ytickformat('%.2f')
    if isSavePlots
        display('Saving Figure')
        saveas(gcf,[GPR.MD.dataDir,'\',[projstr{end-1},titstr{1}],'porosity.png'])
    end

% Depth
figure();scatter(GPR.Geolocation.X{ii}./1000,GPR.Geolocation.Y{ii}./1000,5,GPR.Snow.Depth{ii});colormap(yetBlack)
colorbar;
caxis([quantile(GPR.Snow.Depth{ii},[0.05,0.95])])
daspect([1,1,1])
grid on; grid minor;
title('Depth (cm)')
xlabel('Easting (km)')
ylabel('Northing (km)')
set(gca,'fontsize',12,'fontweight','bold','fontname','serif')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf,'units','normalized','outerposition',[0 0 1 1])
ax = ancestor(hsnow, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.2f')
ax.YAxis.Exponent = 0;
ytickformat('%.2f')
    if isSavePlots
        display('Saving Figure')
        saveas(gcf,[GPR.MD.dataDir,'\',[projstr{end-1},titstr{1}],'depth.png'])
    end
% SWE
figure();scatter(GPR.Geolocation.X{ii}./1000,GPR.Geolocation.Y{ii}./1000,5,GPR.Snow.SWE{ii});colormap(yetBlack)
colorbar;
caxis([quantile(GPR.Snow.SWE{ii},[0.05,0.95])])
daspect([1,1,1])
grid on; grid minor;
title('SWE (mm)')
xlabel('Easting (km)')
ylabel('Northing (km)')
set(gca,'fontsize',12,'fontweight','bold','fontname','serif')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
ax = ancestor(hsnow, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.2f')
ax.YAxis.Exponent = 0;
ytickformat('%.2f')
    if isSavePlots
        display('Saving Figure')
        saveas(gcf,[GPR.MD.dataDir,'\',[projstr{end-1},titstr{1}],'swe.png'])
    end
end
end
