function [x, y, z] = mysph2cart(theta, phi, gain)

%SPH2CART Transform spherical to Cartesian coordinates.
%   [X,Y,Z] = SPH2CART(TH,PHI,R) transforms corresponding elements of
%   data stored in spherical coordinates (azimuth TH, elevation PHI,
%   radius R) to Cartesian coordinates X,Y,Z.  The arrays TH, PHI, and
%   R must be the same size (or any of them can be scalar).  TH and
%   PHI must be in radians.
%
%   TH is the counterclockwise angle in the xy plane measured from the
%   positive x axis.  PHI is the elevation angle from the xy plane.
%
%   Class support for inputs TH,PHI,R:
%      float: double, single
%
%   See also CART2SPH, CART2POL, POL2CART.

%   Copyright 1984-2005 The MathWorks, Inc. 


rcoselev = gain .* cos(phi);
x = rcoselev .* cos(theta);
y = rcoselev .* sin(theta);
z = gain .* sin(phi);%.*sign(y).*sign(x);
