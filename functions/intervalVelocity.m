function [vInterval] = intervalVelocity(v,z)
% stackingVelocity.m computes the stacking velocity from instantaneous v(z)
% Input: v - stacking velocity
%        z - depth of stacking velocities
%
% Output: vStack - Instantaneous (Interval) Velocity Profile
%
% Created by Tate Meehan - 2018
%% Linear System of Equations for Dix Inversion

interval = diff(z); interval = [interval(1);interval];
% Compute Weights for Cumulative Average Density
W = tril(ones(length(interval),1)*interval');
W = W./sum(W,2);

% Invert Stacking Velocity with Interval Weights for Instantaneous V
vInterval = W\v;
end