function [GPR] = onSurface(GPR)  
for ii = 1:GPR.MD.nFiles
    for jj = 1:GPR.Geometry.nChan{ii}
        [m,n] = size(GPR.D.Radar{ii}{jj});
        GPR.D.surfaceTWT{ii}(jj,:) = zeros(1,n);%GPR.D.t0{ii}(jj).*ones(1,n);
        GPR.D.surfaceT0{ii}(jj,:) = zeros(1,n);%GPR.D.t0{ii}(jj).*ones(1,n);
        GPR.D.surfaceH{ii}(jj,:) = zeros(1,n);
        GPR.D.groundTWT{ii}(ii,:) = zeros(1,n);
    end
end
end