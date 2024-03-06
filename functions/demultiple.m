function [out] = demultiple(dat,NF,L,mu)
% demultiple loops throguh the traces and applies the predictive
% deconvolution algorithm (Sacchi 2008).

[m,n] = size(dat);
out = zeros(size(dat));
for kk = 1:n
    w = dat(:,kk);
    [f,o] = predictive(w,NF,L,mu);
    out(:,kk) = o;
end
end

