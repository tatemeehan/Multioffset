function [GPR] = flattenSurface(GPR)
for ii = 1:GPR.MD.nFiles
    nChan = GPR.Geometry.nChan{ii};
    taxe = GPR.D.TimeAxis{ii};
    for jj = 1:nChan
        % Extract GPR Data
        rad = GPR.D.Radar{ii}{jj};
        % Get Size
        [m,n] = size(rad);
        % Allocate Matrix
        flatRad = zeros(m,n);
        % Extract Surface Travel-times
        surfaceTWT = GPR.D.surfaceTWT{ii}(jj,:);
        % Create Datum for Trace Interpolation
        surfaceDatum = round(mean(surfaceTWT),1);
        % Create Grid of Perturbations
        surfacePerturb = surfaceDatum - surfaceTWT;
        % Matrix of Final Traveltimes
        tStak = repmat(taxe,1,n);
        % Matrix of Starting Traveltimes
        tFlat = tStak + surfacePerturb;
        for kk = 1:n
        % Flatten Surface with Linear Interpolation
        flatRad(:,kk) = interp1(tFlat(:,kk),rad(:,kk),tStak(:,kk),'linear','extrap');
        end
        % Save Flattened Image
        GPR.D.Radar{ii}{jj} = flatRad;
        % Correct Travel-Times
        GPR.D.surfaceTWT{ii}(jj,:) = surfaceDatum.*ones(1,n);
    end
    
    % Correct Hieght of Antennas
    tmpH = GPR.D.surfaceH{ii};
    tmpDatumH = mean(tmpH);
    GPR.D.surfaceH{ii} = ones(1,n).*tmpDatumH;
end

end
