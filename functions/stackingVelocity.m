function [vStack] = stackingVelocity(v,z)
% stackingVelocity.m computes the stacking velocity from instantaneous v(z)
% Input: v - instantaneous (interval) velocity
%        z - depth of instantaneous velocities
%
% Output: vStack - Stacking Velocity Profile
% Created by Tate Meehan - 2018
%% Linear System of Equations for Dix Inversion

interval = diff(z); interval = [interval(1);interval];
% Compute Weights for Cumulative Average Density
W = tril(ones(length(interval),1)*interval');
W = W./sum(W,2);

% Cumulative Average as Stacking Velocity
vStack = W*v;
end

