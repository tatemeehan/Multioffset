function [GPR] = gatherCMP(GPR)
% Sort CMP Gathers by Midpoint GPS distance
GPR.D.CMP = cell(1,GPR.MD.nFiles);
for ii = 1 : GPR.MD.nFiles
    chan1Ix = find(GPR.D.trhd{ii}(3,:)==1);
    tmpRadar = cell(1,GPR.Geometry.nChan{ii});
    ntrc = numel(find(GPR.D.trhd{ii}(3,:)==1));
for kk = 1:ntrc
    GPR.D.CMP{ii}{kk}(:,1) = GPR.D.Radar{ii}{1}(:,kk);
    tmpRadar{1}(:,kk) = GPR.D.Radar{ii}{1}(:,kk);
    chan1Dist = GPR.D.trhd{ii}(25,chan1Ix(kk));
    for tt = 2:GPR.Geometry.nChan{ii}
        chanIx = find(GPR.D.trhd{ii}(3,:)==tt);
        chanDist = GPR.D.trhd{ii}(25,chanIx);
        [~,trcIx] = min(abs(chan1Dist-chanDist));
        % Sort CMP Gather
        GPR.D.CMP{ii}{kk}(:,tt) = GPR.D.Radar{ii}{tt}(:,trcIx);
        % Re-Sort Common Offsets by Midpoints
        tmpRadar{tt}(:,kk) = GPR.D.Radar{ii}{tt}(:,trcIx);
    end
end
% Overwrite Radargrams with CMP Sorties
GPR.D.Radar{ii} = tmpRadar;
end
% End of Function
end