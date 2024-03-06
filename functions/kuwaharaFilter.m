function [GPR] = kuwaharaFilter(GPR)
for ii = 1:GPR.MD.nFiles
    nChan = GPR.Geometry.nChan{ii};
    Radar = GPR.D.Radar{ii};
    for jj = 1:nChan
        Radar{jj} = Kuwahara(Radar{jj},5);
    end
    GPR.D.Radar{ii} = Radar;
end
end

