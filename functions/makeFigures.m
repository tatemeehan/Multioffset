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

% Stacking Velocity
V = [];
Distance = [];
for ii = 1 : GPR.MD.nFiles
    Distance = [Distance,GPR.Geolocation.Distance{ii}];
    V = [V,GPR.D.stackingVelocity{ii}];
end

figure();
imagesc(Distance./1000,GPR.D.TimeAxis{1},V);
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
V = [];
for ii = 1 : GPR.MD.nFiles
    V = [V,GPR.D.intervalVelocity{ii}];
end

figure();
imagesc(Distance./1000,GPR.D.TimeAxis{1},V);
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
for ii = 1 : GPR.MD.nFiles
    D = [D,GPR.D.wDensity{ii}];
    Z = [Z,GPR.Geolocation.Z{ii}];
end
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
% Compute Wet Firn Density
wD = WetCrim(V,LWC);
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
figure();
imagesc(Distance./1000,GPR.D.DepthAxis{1},wD);
colormap(yetBlack);
cb = colorbar;
cb.Label.String = ['LWC Density (kg/m^3)'];%['Density (kg/m^3): LWC ',num2str(GPR.D.LWC{1}.*100),' %'];
xlabel('Distance (km)')
ylabel('Depth (m)')
ylim([0,30])
daspect([1,2,1])
set(gca,'fontsize',14,'fontweight','bold','fontname','serif')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
f = gcf;
exportgraphics(f,[figDir,'\','lwcDensity.png'],'Resolution',300)

% LWC Image
figure();
imagesc(Distance./1000,GPR.D.DepthAxis{1},LWC.*100);
colormap(yetBlack);
cb = colorbar;
cb.Label.String = ['LWC (%)'];%['Density (kg/m^3): LWC ',num2str(GPR.D.LWC{1}.*100),' %'];
xlabel('Distance (km)')
ylabel('Depth (m)')
ylim([0,30])
daspect([1,2,1])
set(gca,'fontsize',14,'fontweight','bold','fontname','serif')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
f = gcf;
exportgraphics(f,[figDir,'\','lwc.png'],'Resolution',300)

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
exportgraphics(f,[figDir,'\','radargram.png'],'Resolution',300)

% Time Image
Rad = [];
for ii = 1 : GPR.MD.nFiles
%     Rad = [Rad,GPR.D.RadarStack{ii}];
    Rad = [Rad,GPR.D.RadarNMO{ii}{2}];

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
