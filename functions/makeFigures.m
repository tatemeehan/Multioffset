function [GPR] = makeFigures(GPR)
tmp = strsplit(GPR.MD.dataDir,'\');
tmp(length(tmp)-1:end) = [];
figDir = strjoin(tmp,'\');
figDir = [figDir,'\figures'];
if ~exist(figDir, 'dir')
    mkdir(figDir)
end
% Load Color Maps
yetBlack = load('yetBlack.txt');
radarlove = csvread('radarlove.csv');
RdYlBu = csvread('RdYlBu.csv');
lateNite = csvread('LateNite.csv');

lines = 1:GPR.MD.nFiles-1;
% Stacking Velocity
Vstack = [];
Distance = [];
for ii = 1 : GPR.MD.nFiles
% for ii = lines
    Distance = [Distance,GPR.Geolocation.Distance{ii}];
    Vstack = [Vstack,GPR.D.stackingVelocity{ii}];
end

figure();
imagesc(Distance./1000,GPR.D.TimeAxis{1},Vstack);
colormap(yetBlack);
cb = colorbar;
cb.Label.String = 'Stacking Velocity (m/ns)';
xlabel('Distance (km)')
ylabel('Two-way Travel-time (ns)')
ylim([0,250])
daspect([1,10,1])
set(gca,'fontsize',14,'fontweight','bold','fontname','serif')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
f = gcf;
exportgraphics(f,[figDir,'\','stackingVelocity.png'],'Resolution',300)

% Interval Velocity
Vint = [];
for ii = 1 : GPR.MD.nFiles
% for ii = lines
    Vint = [Vint,GPR.D.intervalVelocity{ii}];
end
Perm = (0.3./Vint).^2;

figure();
imagesc(Distance./1000,GPR.D.TimeAxis{1},Vint);
colormap(yetBlack);
cb = colorbar;
cb.Label.String = 'Interval Velocity (m/ns)';
xlabel('Distance (km)')
ylabel('Two-way Travel-time (ns)')
ylim([0,250])
daspect([1,10,1])
set(gca,'fontsize',14,'fontweight','bold','fontname','serif')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
f = gcf;
exportgraphics(f,[figDir,'\','intervalVelocity.png'],'Resolution',300)

% Dry Firn Density
D = [];
for ii = 1 : GPR.MD.nFiles
% for ii = lines
    D = [D,GPR.D.Density{ii}];
end

figure();
imagesc(Distance./1000,GPR.D.DepthAxis{1},D);
colormap(yetBlack);
cb = colorbar;
cb.Label.String = 'Density (kg/m^3)';
xlabel('Distance (km)')
ylabel('Depth (m)')
ylim([0,25])
daspect([1,2,1])
set(gca,'fontsize',14,'fontweight','bold','fontname','serif')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
f = gcf;
exportgraphics(f,[figDir,'\','Density.png'],'Resolution',300)

% Wet Firn Density
D = [];
Z = [];
try GPR.D.qLWC;
    LWC = [];
    qLWC = 1;
catch
    qLWC = 0;
end
for ii = 1 : GPR.MD.nFiles
% for ii = lines
    % D = [D,GPR.D.wDensity{ii}];
    D = [D,GPR.D.qDensity{ii}];
    Z = [Z,GPR.Geolocation.Z{ii}];
    if qLWC
        LWC = [LWC,GPR.D.qLWC{ii}];
    end
end
if ~qLWC
% Regress LWC against Elevation
minLWC = 0.005; maxLWC = 0.03;
w = ((minLWC-maxLWC)./(max(Z)-min(Z))).*Z;
LWC = w+(minLWC-min(w)).*ones(size(D));
iceIx = find(D>.650);
LWC(iceIx) = 0;
% Convolution Smoothing
[nr,nc] = size(LWC);
R=250;
kernel = hamming(2.*R +1);
% Pad Data with End Traces prior to Convolution
% padl = ones(nr,R).*LWC(:,1);
% padr = ones(nr,R).*mean(LWC(:,nc-2.*R+1:nc),2);
convdata = LWC;%padl,LWC,padr];%
% Row Wise Convolution
padu = ones(R,1).*(convdata(1,:));
padd = ones(R,1).*convdata(nr,:);
convdata = [padu;convdata;padd];
LWC = conv2(convdata(:,:),kernel./sum(kernel),'valid');
clear('convdata','padd','padu','kernel');
end
% Remove Mean Bias of Various Files
% if qLWC
%     meanLWC = mean(LWC(100:1000,:));
%     % Curve Fitting
%     G = [ones(numel(meanLWC),1),Distance(:).^2,Distance(:),movmean(Z(:),11).^2,movmean(Z(:),11)];
%     d = meanLWC(:);
%     m = G\d;
%     dcal = G*m;
%     % Mean Matching
%     tmpLWC = LWC-meanLWC+dcal';
%     tmpLWC(tmpLWC<0) = 0;
%     LWC = tmpLWC;
%     clear("tmpLWC");
% end
% % Compute Wet Firn Density
% isCrim = 0;
% isTinga = 1;
% if isCrim
% wD = WetCrim(V,LWC);
% end
% if isTinga
%       % Tinga73
%       tingaBias = 0.15;
%         % Invert For Density - Line Search
%         testDensity = 100:0.1:950;
%         tmpPerm = Perm;
%         tmpDensity = zeros(numel(tmpPerm),1);
%         kwreal = real(water_permittivity_maetzler87(1.*10.^9, 273.15));
%         frequency = 1e9;
%         temperature = 273.15;
%         liquid_water = LWC;
%         parfor kk = 1:numel(tmpPerm)
%             calculatedK = wetsnow_permittivity_tinga73(frequency, ...
%                 temperature, testDensity, liquid_water(kk))+tingaBias;
%             [~, densityIx] = min(abs(real(calculatedK)-tmpPerm(kk)));
%             tmpDensity(kk) = testDensity(densityIx);
%         end
%         wD = reshape(tmpDensity,size(Perm));
% end
if qLWC
    % Depth Conversion
    nTrcs = size(Vstack,2);
    dwD = zeros(length(GPR.D.DepthAxis{1}),nTrcs);
    dLWC = dwD;
    Zinterp = Vstack.*GPR.D.TimeAxis{1}./2;
    Zaxe = GPR.D.DepthAxis{1};
    for jj = 1:nTrcs
        %parfor (jj = 1:nTrcs, GPR.MD.nWorkers) % parfor is slower
        % Depth Conversion by Simple Stretch
        dwD(:,jj) = interp1(Zinterp(:,jj),D(:,jj),Zaxe,'linear','extrap');
        dLWC(:,jj) = interp1(Zinterp(:,jj),LWC(:,jj),Zaxe,'linear','extrap');
    end
    wD = dwD;
    LWC = dLWC;
    clear("dwD","dLWC");
end
% Overwrite Model in Data Structure!!!
sumncols = 0;
for ii = 1 : GPR.MD.nFiles
    ncols = size(GPR.D.Density{ii},2);
    sumncols = sumncols+ncols;
    if ii == 1
        GPR.D.wDensity{ii} = wD(:,1:ncols);
        GPR.D.LWC{ii} = LWC(:,1:ncols);
    else
        GPR.D.wDensity{ii} = wD(:,prevNcols+1:sumncols);
        GPR.D.LWC{ii} = LWC(:,prevNcols+1:sumncols);
    end
        prevNcols = sumncols;
end
% Density Image
figure();
imagesc(Distance./1000,GPR.D.DepthAxis{1},wD);
colormap(yetBlack);
cb = colorbar;
cb.Label.String = ['Density (kg/m^3)'];%['Density (kg/m^3): LWC ',num2str(GPR.D.LWC{1}.*100),' %'];
xlabel('Distance (km)')
ylabel('Depth (m)')
ylim([0,25])
clim([300,900])
daspect([1,2,1])
set(gca,'fontsize',14,'fontweight','bold','fontname','serif')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
f = gcf;
exportgraphics(f,[figDir,'\','qDensityTinga.png'],'Resolution',300)

% LWC Image
figure();
% imagesc(Distance./1000,GPR.D.DepthAxis{1},LWC.*100);
imagesc(Distance./1000,GPR.D.DepthAxis{1},LWC.*100);
colormap(yetBlack);
cb = colorbar;
cb.Label.String = ['LWC (%)'];%['Density (kg/m^3): LWC ',num2str(GPR.D.LWC{1}.*100),' %'];
xlabel('Distance (km)')
ylabel('Depth (m)')
ylim([0,25])
clim([0,25])
daspect([1,2,1])
set(gca,'fontsize',14,'fontweight','bold','fontname','serif')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
f = gcf;
exportgraphics(f,[figDir,'\','qlwcTinga5.png'],'Resolution',300)

% Depth Image
Rad = [];
for ii = 1 : GPR.MD.nFiles
    Rad = [Rad,GPR.D.RadarDepth{ii}];
end
% Rad = imadjust(Rad);
isImadjust = 1;
figure();imagesc(Distance./1000,GPR.D.DepthAxis{1},Rad);colormap(cmapAdapt(Rad,flipud(lateNite)))
xlabel('Distance (km)')
ylabel('Depth (m)')
ylim([0,25])
daspect([1,2,1])
set(gca,'fontsize',14,'fontweight','bold','fontname','serif')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
f = gcf;
% Adjust Contrast
if isImadjust
    keyboard
%     imcontrast % Run imcontrast command when finished adjusting, exit
%     debugging mode
end
exportgraphics(f,[figDir,'\','radargramfinal25-2.png'],'Resolution',300)

% Time Image
Rad = [];
for ii = 1 : GPR.MD.nFiles
    Rad = [Rad,GPR.D.RadarStack{ii}];
    % Rad = [Rad,GPR.D.RadarNMO{ii}{2}];

end
% Rad = rmsAmplitude(Rad,5);
% Rad = spiking(Rad,250,100);
figure();imagesc(Distance./1000,GPR.D.TimeAxis{1},Rad);colormap(cmapAdapt(Rad,flipud(lateNite)))
xlabel('Distance (km)')
ylabel('Travel-Time (ns)')
ylim([0,250])
daspect([1,20,1])
set(gca,'fontsize',14,'fontweight','bold','fontname','serif')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
f = gcf;
exportgraphics(f,[figDir,'\','radargramtime.png'],'Resolution',300)

% Elevation Plot
figure();plot(Distance./1000,Z,'k','linewidth',2)
xlabel('Distance (km)')
ylabel('Elevation (masl)')
set(gca,'fontsize',14,'fontweight','bold','fontname','serif')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf,'Position',[0    0.4    1    .5]);
ylim([1000 2000])
xlim([0,max(Distance)./1000])
daspect([1,100,1])
grid on; grid minor
f = gcf;
exportgraphics(f,[figDir,'\','elevation.png'],'Resolution',300)

end
