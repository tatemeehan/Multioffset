vmin = 0.15;
vmax = 0.25;
nv = 100;
v = linspace(vmin,vmax,nv);
for kk = 1:numel(GPR.D.velocityCoherence{1})
    % Pick Maximum Coherence for each
[~,vix] = maxk(GPR.D.velocityCoherence{1}{kk}(:),25);
[vr,vc] = ind2sub(size(GPR.D.velocityCoherence{1}{kk}),vix);
vPick(kk) = mean(v(vc));
vPickStd(kk) = std(v(vc));
tPick(kk) = mean(GPR.D.stackingTimeAxis{1}(vr)) - 10.*GPR.D.dt{1}; 
tPickStd(kk) = std(GPR.D.stackingTimeAxis{1}(vr)); 

end
% T at Offset
tmpT = sqrt( tPick.^2 + (GPR.Geometry.offset{1}(1)./vPick).^2 );    
depth = (tPick.*vPick)./2;
density = DryCrim(vPick);
depth = movmean(depth,101);
density = movmean(density,101);
% swe = movmean(depth,50).*movmean(density,50);
swe = depth.*density;

% figure();plot(GPR.Geolocation.Distance{1},movmean(depth,50));
% figure();plot(GPR.Geolocation.Distance{1},movmean(density,50));


figure();
subplot(1,3,1)
scatter(GPR.Geolocation.X{1},GPR.Geolocation.Y{1},5,depth,'filled');colormap(yetBlack)
clim(quantile(depth(:),[0.05,0.95]))
colorbar
daspect([1,1,1])
title('Depth (m)')
subplot(1,3,2)
scatter(GPR.Geolocation.X{1},GPR.Geolocation.Y{1},5,density,'filled');colormap(yetBlack)
clim(quantile(density(:),[0.05,0.95]))
colorbar
daspect([1,1,1])
title('Density (kg./m^3)')
subplot(1,3,3)
scatter(GPR.Geolocation.X{1},GPR.Geolocation.Y{1},5,swe,'filled');colormap(yetBlack)
clim(quantile(swe(:),[0.05,0.95]))
colorbar
daspect([1,1,1])
title('SWE (m)')

figure();imagesc(GPR.Geolocation.Distance{1},GPR.D.TimeAxis{1},GPR.D.Radar{1}{1});colormap(bone)
hold on; plot(GPR.Geolocation.Distance{1},movmean(tmpT,101),'.m')
ylim([0 10])
xlabel('Distance (m)')
ylabel('TWT (ns)')

figure();imagesc(GPR.Geolocation.Distance{1},GPR.D.TimeAxis{1},GPR.D.Radar{1}{3});colormap(bone)
hold on; plot(GPR.Geolocation.Distance{1},movmean(tmpT,101),'.m')
ylim([0 10])
xlabel('Distance (m)')
ylabel('TWT (ns)')