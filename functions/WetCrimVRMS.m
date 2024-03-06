function [ VRMS ] = WetCrimVRMS( rho, W )
%WetCrimVRMS Returns the RMS Velocity Caused by a SnowPack of Density rho
%  
%
Va = 2.998E8;   % Velocity of Free Space [m/s]
rhoi = 917;       % Pure Ice Density, Ulaby et al.(1986)
Ki = 3.15;      % Real Dielectric Constant of Ice at Microwave Frequency 
                % Ulaby et al. (1986)
Kw = 80.0723;       % Abs Complex Dielectric of Water
% Kw = 86;% Bradford & Harper (2005)
Vw = Va/sqrt(Kw);   % Velocity of Water [m/s]
                                
% Handle Units of Density
if rho<=1
    rho = rho.*1000;
end
fi = rho./rhoi;

Kf = ((1 - fi - W) + sqrt(Ki).*fi + sqrt(Kw).*W).^2;

VRMS = (Va./sqrt(Kf)).*10^-9; %[m/ns]
                



end