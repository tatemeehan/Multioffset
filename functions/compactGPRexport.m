function [GPRlite] = compactGPRexport(GPR)
% GPRlite.m Extracts the Most inportant Data from GPR.mat to create a more
% compact/transefferable data package.
% Meta Data
GPRlite.MD = GPR.MD;
% Data
GPRlite.D.TimeAxis = GPR.D.TimeAxis;
GPRlite.D.DepthAxis = GPR.D.DepthAxis;
GPRlite.D.stackingVelocity = GPR.D.stackingVelocity;
GPRlite.D.intervalVelocity = GPR.D.intervalVelocity;
GPRlite.D.Depth = GPR.D.Depth;
GPRlite.D.Density = GPR.D.Density;
GPRlite.D.wDensity = GPR.D.wDensity;
GPRlite.D.LWC = GPR.D.LWC;
GPRlite.D.RadarStack = GPR.D.RadarStack;
GPRlite.D.RadarDepth = GPR.D.RadarDepth;
% Geometry
GPRlite.Geometery = GPR.Geometry;
% Geolocation
GPRlite.Geolocation = GPR.Geolocation;

end