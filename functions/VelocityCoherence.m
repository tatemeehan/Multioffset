function [GPR] = VelocityCoherence(GPR,vmin,vmax,nv,R,L,stretch)
% VelocityCoherence.m Computes the Radar Velocity Spetra from CMP Gathers.
% Inputs:   GPR  - Radar Data Structure
%           vmin - Minimum Velocity in Analysis (Default 0.17)
%           vmax - Maximum Velocity in Analysis (Default 0.25)
%           nv   - Number of Velocities to Compute (Default 100 = 0.0008 m/ns)
%           R    - Interval of Time Samples to Compute (Default 5)
%           L    - Length of Smoothing Operator (Default 5)
%        stretch - Percentage Threshold of Stretch Mute (Default 1)
%
% Outputs:  GPR  -

% Establish Default Values
if nargin<2
    vmin = 0.15;
    vmax = 0.25;
    nv = 100;
    R = 5;
    L = 5;
    stretch = 1;
elseif nargin < 3
    vmax = 0.25;
    nv = 100;
    R = 5;
    L = 5;
    stretch = 1;
elseif nargin < 4
    nv = 100;
    R = 5;
    L = 5;
    stretch = 1;
elseif nargin < 5
    R = 5;
    L = 5;
    stretch = 1;
elseif nargin < 6
    L = 5;
    stretch = 1;
elseif nargin < 7
    stretch = 1;
end

for ii = 1 %: GPR.MD.nFiles
    dt = GPR.D.dt{ii};
    nt = size(GPR.D.CMP{ii}{1},1);
    offset = GPR.Geometry.offset{ii};
    % Coherence Gain Window
    gt = 50; % 50 nanosecond Gain Window
    gix = floor((gt./dt)./R)+1;
    g = (4.*tukeywin(gix)+1);
    % Store Velocity and Time Axes
    GPR.D.velocityAxis{ii} = linspace(vmin,vmax,nv);   % Velocity Estimates
    GPR.D.stackingTimeAxis{ii} = (0:R:nt-1)*dt;        % Stacking Time
    % Allocation
    ncmp = size(GPR.D.CMP{ii},2);
    C = cell(1,ncmp);%S = C;
    CMP = GPR.D.CMP{ii};
    % Velocity Analysis
    parfor (jj = 1:ncmp,GPR.MD.nWorkers)
        % for jj = 1:length(CMP)
        tmpCMP = single(CMP{jj});
%         tmpCMP = single(AGCgain(tmpCMP,150,2));
        tmpCMP = rmsAmplitude(tmpCMP,50);
%         [S{jj},~,~] = lmoVelocityCoherence(tmpCMP,dt,offset,vmin,vmax,nv,R,L,1);
        [C{jj},~,~] = nmoVelocityCoherence(tmpCMP,dt,offset,vmin,vmax,nv,R,L,stretch);
        % Gain Coherence of Early Reflections with Reduced Data Fold
        C{jj}(1:gix,:) = C{jj}(1:gix,:).*g;
        % Normalize
        C{jj} = C{jj}./max(C{jj}(:));
    end
    % Store Velocity Coherence Spectra
    GPR.D.velocityCoherence{ii} = C;
end
% End of Function
end