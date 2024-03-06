function [GPR] = raytraceConstraint(GPR,ps)
for ii = 1:GPR.MD.nFiles
    if nargin<2
        % Default Approximate Snow Density
        ps = 295;
    end
% Desired VRMS
va = 0.3; vs = DryCrimVRMS(ps); ha = mean(GPR.D.surfaceH{ii}); hg = mean(GPR.D.groundH{ii});
dvrms = stackingVelocity([va;vs],[ha;hg]);
perturbVRMS = GPR.D.groundVRMS{ii}-dvrms;
% Update VRMS with Measured Variability
VRMS = (GPR.D.groundVRMS{ii}-perturbVRMS)+perturbVRMS./10;
groundH = GPR.D.groundT0{ii}.*VRMS./2;
% Raytracing Step
synthTT = zeros(size(GPR.D.groundTWT{ii}));
% Raytrace Parameters
xcap = 0.025;
pfan = -1;
itermax = 20;
optflag = 1;
pflag = 0;
dflag = 0;
for kk = 1:length(GPR.D.groundH{ii})
synthTT(:,kk) = traceray_pp([VRMS(kk);VRMS(kk)],...
    [0;groundH(kk)],0,0,groundH(kk),GPR.Geometry.offset{ii},xcap,pfan,itermax,optflag,pflag,dflag);
end
% Update Travel Times
 GPR.D.groundTWT{ii} = synthTT;
end
end

