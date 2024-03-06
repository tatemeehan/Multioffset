function [rad] = backSubtract(rad,x,R)
% submed.m is a moving window mean subtaction filter for use as a
% background subtraction method for GPR.
rad2 = rad;
[~, ntrc] = size(rad);
for kk = 1:ntrc
    dist = sqrt((x-x(kk)).^2);
    ix = dist<R;
    tmp1 = rad2(:,kk);
    tmp2 = rad2(:,ix);
    rad(:,kk) = tmp1 - mean(tmp2,2);
end
end

