function [GPR] = HVA(GPR)
for ii = 1:GPR.MD.nFiles
    % Allocation
    offsets = (GPR.Geometry.offset{ii});
    groundTWT = GPR.D.groundTWT{ii};
    nChan = GPR.Geometry.nChan{ii};    
    Rad = GPR.D.Radar{ii}{1};
    [nsamp, ntrcs] = size(Rad);
    % Kill Channels
    kill = [];
    offsets(kill) = [];
    groundTWT(kill,:) = [];
    nChan = nChan - numel(kill);
% Estimate zero-offset travel-time and velocity of free space
    isL1 = 0; isL2 = 1; % Select Norm
%     if isL1
        p = 1;
        tolr = 0.0005;
        tolx = 0.0005;
        maxiter = 50;
%     end
    % Monte Carlo Bootstrapping 
    isBootstrap = 1; 
    % Optionally Several CMP gather can be analyzed within a radius r
    isManyGathers = 1; 
    isDistill = 0; % Monte Carlo with each offset represented only once;
    if isBootstrap
        nMC = 250;
        % Extract Locations
        X = GPR.Geolocation.X{ii};
        Y = GPR.Geolocation.Y{ii};
        r = 50;%2.5; % Bin Radius (m)
    else
        G = [ones(nChan,1),offsets(:)].^2;
    end
    m = zeros(2,ntrcs);mstd = m;MCout =zeros(2,ntrcs,nMC);
    % Least Squares Estimate of Travel-Time and Velocity
    parfor kk = 1:length(Rad(1,:))
        if isBootstrap
        if isManyGathers
        dist = sqrt((X-X(kk)).^2+(Y-Y(kk)).^2);
        binIx = find(dist<r);
        end
        mMC = zeros(2,nMC);
        for ll = 1:nMC
            if isManyGathers
                if isDistill
                    keepBinIx = datasample(binIx,nChan);
                    keepChanIx = datasample(1:nChan,nChan,'replace',false);
                    % Extract  Picks
                    tmpIx = sub2ind(size(groundTWT),keepChanIx,keepBinIx);
                    d = groundTWT(tmpIx).^2;d = d(:);
                    offsetArray = offsets(keepChanIx);offsetArray =offsetArray(:);
                    G = [ones(length(offsetArray),1),offsetArray.^2];
                    % Now Randomly Remove one Channel
                    rmvIx = datasample(1:nChan,1);
                    d(rmvIx) = []; G(rmvIx,:) = [];
                else
                % Many Gathers Meh :/
                nbins = round(.25.*numel(binIx));
                offsetArray = repmat(offsets(:),1,nbins);
                rmvIx = datasample(1:nChan,nbins);
                rmvIx = sub2ind(size(offsetArray),rmvIx,1:nbins);
                MCbin = datasample(binIx,nbins,'replace',false);
                MCix = [1:numel(offsetArray)]';MCix(rmvIx) = [];
                
                % Extract Surface Picks
                d = groundTWT(:,MCbin).^2;d = d(MCix);
                % Build G matrix
                offsetArray = offsetArray(MCix);
                G = [ones(length(offsetArray),1),offsetArray.^2];
                end
                if isL2
                mMC(:,ll) = G\d;
                end
                if isL1
                    mMC(:,ll) = irls(G,d,tolr,tolx,p,maxiter);
                end
            else
                % Single Gather
                rmvIx = datasample(1:nChan,1);
                offsetArray = offsets;
                % Extract Surface Picks
                d = groundTWT(:,kk).^2;d = d(:); d(rmvIx) = [];
                % Build G matrix
                offsetArray = offsetArray(:); offsetArray(rmvIx) = [];
                G = [ones(length(offsetArray),1),offsetArray.^2];
                if isL2
                    mMC(:,ll) = G\d;
                end
                if isL1
                    mMC(:,ll) = irls(G,d,tolr,tolx,p,maxiter);
                end
            end
        end
        m(:,kk) = mean(mMC,2);
        MCout(:,kk,:) = mMC;
%         mstd(:,kk) = std(mMC,[],2);
        else
        % Without Bootstrapping
        d = groundTWT(:,kk).^2;
        if isL2
            m(:,kk) = G\d;
        end
        if isL1
            m(:,kk) = irls(G,d,tolr,tolx,p,maxiter);
        end
        end
    end
    % The least-squares estimate can be thought of as a filter.   
    % Velocity Calibration
    t0 = sqrt(m(1,:));
    vs = sqrt(1./(m(2,:)));
    MCt0 = sqrt(MCout(1,:,:));
    MCvs = sqrt(1./(MCout(2,:,:)));
%     t0std = (sqrt(m(1,:)+mstd(1,:))) - sqrt(m(1,:));
%     vsstd =  (sqrt(1./(m(2,:))+(1./mstd(2,:))))-sqrt(1./(m(2,:)));
%     t0err = t0cal - sqrt(m(1,:));
% 
%     % Height of Antennas Above Ground
%     hest = sqrt(m(1,:)).*sqrt(1./m(2,:))./2;
    h = (t0.*vs)./2;
    MCh = (MCt0.*MCvs)./2;
%     herr = t0err.*va./2;
    % Store Output Vars
    GPR.D.groundT0{ii} = t0;
    GPR.D.groundH{ii} = h;
    GPR.D.groundVRMS{ii} = vs;
%     GPR.D.airVelocity{ii} = sqrt(1./m(2,:));
% Monte Carlo Output
GPR.D.MCt0{ii} = MCt0;
GPR.D.MCh{ii} = MCh;
GPR.D.MCvrms{ii} = MCvs;
end
end