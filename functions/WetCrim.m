function [ pF,PhiF ] = WetCrim( Vf,W )
%   WetCrim         Approximates Wet Firn/Sno Density and Porosity 
%                   Using a 3-Phase Mixing Model
%
%   Inputs:         EM Velocity of Firn/Sno,       Vf
%                   Percent Water Content          W
%
%   Outputs:        Density of Firn,               pF
%                   Porosity of Firn,              PhiF
%
%   Defined:
%                   EM Velovity of Ice,            Vi
%                   EM Velocity of Air,            Va
%                   Density of Ice,                pI
%   Units are MKS
%
%   Tate Meehan - Boise State University Geophysics - March 7, 2016
%       Adapted From A.P. Annan, S.W. Cosway, and T. Sigurdsson (1994).
%       GPR For Snow Pack Water Content. Proceedings of the Fifth
%       International Conference on Ground Penetrating Radar, June 1994.

Va = 2.998E8;   % Velocity of Free Space [m/s]
Ki = 3.15;      % Real Dielectric Constant of Ice at Microwave Frequency 
                % Ulaby et al. (1986)
Vi = Va/sqrt(Ki);   % Velocity of Ice [m/s]
pi = 917;        % Pure Ice Density, Ulaby et al.(1986)
Kw = 80.0723;       % Abs Complex Dielectric of Water
Vw = Va/sqrt(Kw);   % Velocity of Water [m/s]

% Condition Units of Velocity [m/ns] to [m/s]
if any(Vf(:) < 1E8)
    Vf = Vf.*(10^9); 
end

Kf = (Va./Vf).^2; % Approximate Dielectric Constant of Firn eqn. E. 82 Ulaby et al. (1986)

% PhiF = (sqrt(Kf)-sqrt(Ki)+W.*sqrt(Ki)-W.*sqrt(Kw))./(1-sqrt(Ki)); 
                                      %   Wet Firn/Sno Porosity (Annan et
                                      %   al., 1994)
% pF = pI.*((Vi.*Va)/(Va-Vi)).*((1/Vf)-(W/Vw)-(1/Va)+(W/Va)); %  Wet Density
% pF = pi.*(1-PhiF);                    %   Firn/Sno Density
% pF = pi.*(1-PhiF-W)+W;                    %   Firn/Sno Density
% pF = pF/1000;                         %   Firn/Sno Percent Density

% 1-23-17
fi = (sqrt(Kf) - sqrt(Kw).*W - (1-W))./(sqrt(Ki)-1);

rhoi = pi.*fi;
rhoW = 1000.*W;
pF = (rhoi+rhoW)./1000;
PhiF = 1-fi-W;

end

