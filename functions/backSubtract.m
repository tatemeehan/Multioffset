function [rad] = backSubtract(rad,x,R)
% submed.m is a moving window mean subtaction filter for use as a
% background subtraction method for GPR.
tic
rad2 = rad;
[~, ntrc] = size(rad);
dx = mean(diff(x));
nR = round(R.*dx);
modTrc = 500; % Compute Mean every 25 traces
filterMod = round(nR./modTrc);
for kk = 1:ntrc
    dist = sqrt((x-x(kk)).^2);
    ix = dist<R;
    tmp1 = rad2(:,kk);
    tmp2 = rad2(:,ix);
    % Improved Computation Efficiency
    if kk == 1
        tmp3 = mean(tmp2,2);
        rad(:,kk) = tmp1 - tmp3;
    elseif (mod(kk,filterMod) == 0)
        tmp3 = mean(tmp2,2);
        rad(:,kk) = tmp1 - tmp3;
    else
         rad(:,kk) = tmp1 - tmp3;
    end
end
toc
end

