function [GPR] = DepthConversion(GPR)
% DepthConversion.m Performs Time-to-Depth Conversion on the Stacked Radar
% Section by linear interpolation to a linear depth axis. The 2-D velocity
% and depth space is generated in expVelocityModeling.m. Normal moveout
% corection and image processing must occur before depth conversion!
%
% Created by Tate Meehan, CRREL, Juneau IceField Research Project, 2024

for ii = 1:GPR.MD.nFiles
    nTrcs = size(GPR.D.RadarStack{ii},2);
    depthRadar = zeros(length(GPR.D.DepthAxis{ii}),nTrcs);
    tmpRadar = GPR.D.RadarStack{ii};
    Z = GPR.D.Depth{ii};
    Zaxe = GPR.D.DepthAxis{ii};
    for jj = 1:nTrcs
        %parfor (jj = 1:nTrcs, GPR.MD.nWorkers) % parfor is slower
        % Depth Conversion by Simple Stretch
        depthRadar(:,jj) = interp1(Z(:,jj),tmpRadar(:,jj),Zaxe,'linear','extrap');
    end
    GPR.D.RadarDepth{ii} = depthRadar;
end
% End of Function
end