function [S,tau,v] = lmoVelocityCoherence(d,dt,h,vmin,vmax,nv,R,L,stretchMax)
% VelocityCoherence: A program to compute velocity spectra
%                    applying linear move out equation.
%
%  [S,tau,v] = lmoVelocityCoherence(d,dt,h,vmin,vmav,nv,R,L);
%
%  IN   data:      data
%       dt:        sampling interval in secs
%       h:         vector of offsets in meters
%       vmin:      min velocity to scan
%       vmax:      max velocity to scan
%       nv:        number of velocities
%       R:         integer (2,3,4) indicating that the coherence
%                  is computed every R time samples
%       L:         the length of the gate of analysis is 2L+1
%       stretchMax: 0 - 1 stretch threshold
%
%  OUT  S:         Measure of energy - in this case unnormalized
%                  cross-correlation in the gate
%       v,tau:     axis vectors to plot the semblance using, for instance,
%                  imagesc(v,tau,S)
%
%  Reference: Yilmaz, 1987, Seismic Data Processing, SEG Publication
%
%  Example: see va_demo.m
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%
% Modified by Tate Meehan - Boise State University - Geophysics
%   Provided Iterative Stretch Muting - GreenTrACS 2016

[nt, nh] = size(d); %nt no. time samples, nh no. horizontal offsets

% [ntog, nh] = size(d); %nt no. time samples, nh no. horizontal offsets
% nt = 250; % Limit Linear Computation to 50 samples (Direct Wave)
% tauog = (0:R:ntog-1)*dt;        % Stacking Time
% ntauog = length(tauog);         % no. of Time Samples
v = linspace(vmin,vmax,nv); % Allocation for Velocity Estimates
nv = length(v);             % no. of Velocity Guesses
tau = (0:R:nt-1)*dt;        % Stacking Time
ntau = length(tau);         % no. of Time Samples
taper = hamming(2*L+1);     % Create Hamming Weighting Filter
H = taper*ones(1,nh);       % Apply Smoothing Along all Traces; nt = 2L + 1
S = zeros(nt,nv);         % Allocate Coherence Matrix

    for iv = 1:nv
        tmpData = d;
        for it = 1:ntau

            % Apply Linear MoveOut Forward Model
            time = (h/v(iv));%tau(it) - (h/v(iv));

            s = zeros(2*L+1,nh);    % Allocate Velocity Analysis Window
            
            for ig = -L:L
                ts = time + (ig-1)*dt; % Computed Time (+-) Filter Sample Time
                
                % Calculate Percent Stretch
                stretch = (ts-tau(it))./(tau(it)+1e-10);
                
                for ih = 1:nh
                    
                    is = ts(ih)/dt+1; % Sample Time Index
                    i1 = floor(is);   % Round to lower index
                    i2 = i1 + 1;      % Adjacent Index
                    
                    if i1>=1 && i2<=nt % Condition Handle - Prevents Subscript Exceed
                        % Apply Stretch Mute
                        if stretch(ih)>stretchMax
%                             tmpData([i1,i2],ih) = 0;
                        end
                        a = is-i1;      % Weight for Sum
                        s(ig+L+1,ih) = (1.-a)*tmpData(i1,ih) + a*tmpData(i2,ih);
                        %Grab sample with weighted sum interpolation
                    end
                    
                end
            end
            
            s = s.*H;                   % Weighted Amplitude Window
            s1  = sum( (sum(s,2)).^2);  % Stacked Energy by NMO at Estimated Vrms - Weighted by Hamming Window
            s2  = sum( sum(s.^2));      % Unstacked Energy - Weighted by Hamming Window
 
            % Half Difference Between Out Energy and In Energy                          
            S(it,iv) = 0.5.*(s1-s2);           % -inf to inf Coherence
            if S(it,iv)<0
             S(it,iv) = 0;                     % 0 to inf Coherence
            end
        end
    end
    
    S = S/max(max(S));                % Normalization Occurs Externally

end