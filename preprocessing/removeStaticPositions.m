function [removeTraces] = removeStaticPositions(trhd,v,velocity,nChan)
% removeStaticPositions.m uses GPS velocity stored in trhd and a velocity
% threshold v to determine the indicies of stationary positions. 
% trhd vars used
% trhd(23,:) Channel number
% trhd(18,:) Velocity

if nargin < 3
    velocity = trhd(18,:);
    nChan = length(unique(trhd(23,:)));
end
if nargin < 4
    nChan = length(unique(trhd(23,:)));
end
multiplexNtrcs = length(trhd);
if nChan > 1
    % multichannel
    binIx = find(trhd(23,:) == 1);
else
    % single channel
    binIx = 1:multiplexNtrcs;
end
rmIx = find(velocity(binIx) < v);
staticPlexTraces = binIx(rmIx);
% Infill Multiplex Trace Indicies For Removal
% removeBin = zeros(length(staticPlexTraces),nChan);
% removeBin = [staticPlexTraces(:),staticPlexTraces(:)+1,...
%     staticPlexTraces(:)+2,staticPlexTraces(:)+3];
removeBin = [];
for kk = 1:nChan
removeBin = [removeBin;staticPlexTraces(:)+(kk-1)];
end

% for ii = 1:nChan
%     if ii == nChan
%         removeBin(:,ii) = staticPlexTraces;
%     else
%     removeBin(:,ii) = staticPlexTraces - (nChan - ii);
%     end
% end
    % Indicies of Static Traces
    removeTraces = removeBin(removeBin(:) <= multiplexNtrcs);

end

