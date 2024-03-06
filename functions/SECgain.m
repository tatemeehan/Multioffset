function [dataG] = SECgain(data, a, b)
% SECgain - Spreading and Exponential Compensation (SEC) Gain
% Applies the gain function g = exp(at)*t^b to the traces
% Inputs: data - the GPR traces
%            a - the exponential gain coefficient
%            b - the power gain coefficient
% Ouput: dataG - the gained data
% Created by Tate Meehan - 2022
[m,n] = size(data);
dataG = zeros(m,n);
% Defaults
if nargin == 1
    a = 0.001; % Shallow exponential
    b = 2; % t-squared Gain
elseif nargin == 2
    b = 2; % t-squared Gain
end
% SEC Gain Function
g = exp(a.*(1:m)').*(1:m)'.^b;
for kk = 1:n
    dataG(:,kk) = g.*data(:,kk);
end
end

