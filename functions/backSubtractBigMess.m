function [rad] = backSubtractBigMess(rad, x, R)
rad2 = rad;
[ns, ntrc] = size(rad);
% Precompute distances and means
dist = sqrt((x - x').^2);  % Calculate distances between all x points
isDist = (dist < R);
mean_tmp2 = zeros(size(rad2));
for kk = 1:10:ntrc
    % Precompute mean for all traces where dist < R
    mean_tmp2(:,kk) = mean(rad2(:, isDist(:,kk)), 2);  
end
interpIx = find(~mean_tmp2(1,:));
goodIx = find(mean_tmp2(1,:));
for kk = 1:ns
    mean_tmp2(kk,interpIx) = interp1(x(goodIx),mean_tmp2(kk,goodIx),x(interpIx),"linear","extrap");
end

for kk = 1:ntrc
    tmp1 = rad2(:, kk);
    % Use precomputed mean
    tmp3 = mean_tmp2(:, kk);  % mean for current trace kk
    rad(:, kk) = tmp1 - tmp3;
end
end