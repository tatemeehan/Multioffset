function [GPR] = estimateSnow(GPR,isOnSurface)
for ii = 1:GPR.MD.nFiles
vrms = GPR.D.groundVRMS{ii};
vrmsMC = GPR.D.MCvrms;
if isOnSurface
%     GPR.Snow.Depth{ii} =  GPR.D.groundH{ii}(:).*100;
%     [GPR.Snow.Density{ii},GPR.Snow.Porosity{ii}] = DryCrim(vrms(:));
%     GPR.Snow.Density{ii} = GPR.Snow.Density{ii}.*1000;
%     GPR.Snow.SWE{ii} = GPR.Snow.Depth{ii}.*GPR.Snow.Density{ii}./100;

    % Extract Mean and Standard Deviation
    for kk = 1:numel(vrms)
        GPR.Snow.Depth{ii}(kk) = mean(GPR.D.MCh{ii}(1,kk,:)).*100;
        GPR.Snow.DepthStd{ii}(kk) = std(GPR.D.MCh{ii}(1,kk,:)).*100;
        [tmp1,tmp2] = (DryCrim(vrmsMC{ii}(1,kk,:)));
        GPR.Snow.Density{ii}(kk) = mean(tmp1).*1000;
        GPR.Snow.Porosity{ii}(kk) = mean(tmp2);
        GPR.Snow.DensityStd{ii}(kk) = std(tmp1).*1000;
        GPR.Snow.PorosityStd{ii}(kk) = std(tmp2);
        GPR.Snow.SWE{ii}(kk) = GPR.Snow.Depth{ii}(kk).*GPR.Snow.Density{ii}(kk)./100;
        GPR.Snow.SWEstd{ii}(kk) = GPR.Snow.SWE{ii}(kk).*sqrt((GPR.Snow.DensityStd{ii}(kk)...
            ./GPR.Snow.Density{ii}(kk)).^2+(GPR.Snow.DepthStd{ii}(kk)./GPR.Snow.Depth{ii}(kk)).^2);
    end
else
v = [0.3.*ones(1,length(vrms));vrms];
z = [GPR.D.surfaceH{ii};GPR.D.groundH{ii}];
% z = [mean(GPR.D.surfaceH{ii}).*ones(1,length(vrms));GPR.D.groundH{ii}];
% z = [movmean(GPR.D.surfaceH{ii},51);GPR.D.groundH{ii}];
% z = [ones(1,length(vrms));GPR.D.groundH{ii}];

vInterval = zeros(size(v));
% Dix Inversion for Interval Velocity
for kk = 1:length(vrms)
[vInterval(:,kk)] = intervalVelocity(v(:,kk),z(:,kk));
end
GPR.D.vInterval{ii} = vInterval;
GPR.Snow.Depth{ii} =  (GPR.D.groundH{ii}(:) - (GPR.D.surfaceH{ii}(:))).*100;
%(((GPR.D.groundT0{ii}(:)-GPR.D.surfaceT0{ii}(:)).*vInterval(2,:)')./2).*100;%nonParametricSmooth( Distance,vInterval(2,:),Distance,25))./2;
[GPR.Snow.Density{ii},GPR.Snow.Porosity{ii}] = DryCrim(vInterval(2,:)');
% [GPR.Snow.Density{ii},GPR.Snow.Porosity{ii}] = DryCrim(GPR.D.groundVRMS{ii}(:));
% [GPR.Snow.Density{ii},GPR.Snow.Porosity{ii}] = WetCrim( vInterval(2,:)',0.005 );
GPR.Snow.Density{ii} = GPR.Snow.Density{ii}.*1000;
GPR.Snow.SWE{ii} = GPR.Snow.Depth{ii}.*GPR.Snow.Density{ii}./100;
end
end
end

