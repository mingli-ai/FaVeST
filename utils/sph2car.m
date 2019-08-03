function x = sph2car(theta,phi)
% sph2car(theta,phi) convert spherical coordinates to Cartesian
% coordinates in R^3.
% 
% INPUTS
%  theta  - lattitude radian from the North pole, range in [0,\pi],
%            size(theta) = [1 N]
%  phi    - longitude radian anticlock-wise from X-axis, range in [0,2\pi)
%            size(phi) = [1 N]
% OUTPUTS
%  x      - 3-by-N coordinates corresponding to (theta,phi), where N =
%            numel(theta)


x = [sin(theta).*cos(phi); sin(theta).*sin(phi); cos(theta)];