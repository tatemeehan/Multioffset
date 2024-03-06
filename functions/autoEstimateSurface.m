function [GPR] = autoEstimateSurface(GPR)
for ii = 1:GPR.MD.nFiles
% Extract Variables
surfaceTWT = GPR.D.surfaceTWT{ii};
offsets = GPR.Geometry.offset{ii};
nChan = GPR.Geometry.nChan{ii};
Rad = GPR.D.Radar{ii}{1};
[nsamp, ntrcs] = size(Rad);
va = 0.3;

% Estimate zero-offset travel-time and velocity of free space
isSolveVa = 0;
% Monte Carlo Bootstrapping
isBootstrap = 1;
% Optionally Several CMP gather can be analyzed within a radius r
isManyGathers = 0;
%     if isBootstrap
nMC = 50;
% Extract Locations
X = GPR.Geolocation.X{ii};
Y = GPR.Geolocation.Y{ii};
Distance = GPR.Geolocation.Distance{ii};
r = 0.5; % Bin Radius (m)
%     else
%     end
if isSolveVa
    m = zeros(2,ntrcs);
else
    m = zeros(1,ntrcs);
end
% Least Squares Estimate of Travel-Time and Velocity
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary')
%     parfor (kk = 1:length(Rad(1,:)),GPR.MD.nWorkers)

for kk = 1:length(Rad(1,:))
    if isBootstrap
        if isManyGathers
            dist = sqrt((X-X(kk)).^2+(Y-Y(kk)).^2);
            binIx = find(dist<r);
        end
        if isSolveVa
            mMC = zeros(2,nMC);
        else
            mMC = zeros(1,nMC);
        end
        for ll = 1:nMC
            if isManyGathers
                % Many Gathers Meh :/
                nbins = round(.25.*numel(binIx));
                offsetArray = repmat(offsets(:),1,nbins);
                rmvIx = datasample(1:nChan,nbins);
                rmvIx = sub2ind(size(offsetArray),rmvIx,1:nbins);
                MCbin = datasample(binIx,nbins,'replace',false);
                MCix = [1:numel(offsetArray)]';MCix(rmvIx) = [];
                % Extract Surface Picks
                d = surfaceTWT(:,MCbin).^2;d = d(MCix);
                offsetArray = offsetArray(MCix);
                if isSolveVa
                    % Build G matrix
                    G = [ones(length(offsetArray),1),offsetArray.^2];
                    mMC(:,ll) = G\d;
                else
                    % Solve for Single Parameter Intercept
                    G = [ones(length(offsetArray),1)];
                    d = [d-(1./va).^2.*offsetArray(:).^2];
                    mMC(ll) = G\d;
                end
            else
                % Single Gather
                rmvIx = datasample(1:nChan,1);
                offsetArray = offsets;
                % Extract Surface Picks
                d = surfaceTWT(:,kk).^2;d = d(:); d(rmvIx) = [];
                offsetArray = offsetArray(:); offsetArray(rmvIx) = [];
                if isSolveVa
                    % Build G matrix
                    G = [ones(length(offsetArray),1),offsetArray.^2];
                    mMC(:,ll) = G\d;
                else
                    % Solve for Single Parameter Intercept
                    G = [ones(length(offsetArray),1)];
                    d = [d-(1./va).^2.*offsetArray(:).^2];
                    mMC(ll) = G\d;
                end
            end
        end
        m(:,kk) = mean(mMC,2);
    else
        % Without Bootstrapping
        if isSolveVa
            G = [ones(nChan,1),offsets(:)].^2;
            d = surfaceTWT(:,kk).^2;
            m(:,kk) = G\d;
        else
            % Solve for Single Parameter Intercept
            G = [ones(nChan,1)];
            d = [surfaceTWT(:,kk).^2-(1./va).^2.*offsets(:).^2];
            m(kk) = G\d;
        end
    end
end
% The least-squares estimate can be thought of as a filter.
% Smooth Estimates
smoothr = 0.5;%2;
if isSolveVa
    [ m(1,:) ] = nonParametricSmooth( Distance,m(1,:),Distance,smoothr );
    [ m(2,:) ] = nonParametricSmooth( Distance,m(2,:),Distance,smoothr );
else
    [ m(1,:) ] = nonParametricSmooth( Distance,m(1,:),Distance,smoothr );
end
% Velocity Calibration
%     t0cal = sqrt(m(1,:));
% Just use one channel 0.5 m and see the diff
t0cal = sqrt(surfaceTWT(2,:).^2-(0.5.^2./0.3.^2));
%     t0cal = sqrt(m(1,:)).*sqrt(1./m(2,:))./va;
t0err = t0cal - sqrt(m(1,:));
figure();histogram(t0err,'FaceColor','k','normalization','probability');title('travelime to surface');ylabel('probability');xlabel('error (ns)')
%
%     % Height
%     hest = sqrt(m(1,:)).*sqrt(1./m(2,:))./2;
hcal = (t0cal.*va)./2;
herr = t0err.*va./2;
%     figure();histogram(herr,'FaceColor','k','normalization','probability');title('height above surface');ylabel('probability');xlabel('error (m)')
% Store Output Vars
GPR.D.surfaceT0{ii} = t0cal;
GPR.D.surfaceH{ii} = hcal;
%     GPR.D.airVelocity{ii} = sqrt(1./m(2,:));
end
end