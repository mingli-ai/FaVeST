function [w,x_sph,x_ca] = QpS2(L,QN)
% [w,x_sph,x_ca] = QpS(L,QN) computes the quadrature rule that is exact 
% for polynomials of degree L on S^2, and returns the weights, quadrature 
% nodes in both Cartesian and spherical coordinates.
%
% Input
%  L      - Degree for which the QN is exact
%
%  QN     - name of the quadrature rule4
%
% Output
%  w      - weights of QN, size(w) = [1 N]
%
%  x_sph  - spherical coordinates for quadrature nodes, size(x_sph) = [N 3]
%
%  x_ca   - Cartesian coordinates for quadratrue nodes, size(x_ca) = [2 N],
%            first-row phi are in [0,2pi) and the second-row theta in [0,pi].

%L = 100;
%QN = 'SD';

switch QN
    case 'SD'
        qN = @SD;
        [w,x_ca] = qN(2*L);
        w = w';
        x_sph = zeros(size(x_ca,1),2);
         % turn x_ca to its spherical coordinates
        [x_sph(:,2),x_sph(:,1)] = car2sph(x_ca);
        x_sph = x_sph';
        x_ca = x_ca';
    case 'GL'
        qN = @gl;
        [x_sph,w] = qN(L);
        x_ca = sph2car(x_sph(2,:),x_sph(1,:));
end
       
    
end

function [w,y] = SD(L)
% [w,y] = SD(L,d)
% computes the symmetric spherical design of Rob
%
% Inputs:
% L -- degree for which the SD quadrature is exact
% d -- dimension of sphere; by default d = 2
%
% Outputs:
% w -- weights of quadrature rule, size(w) = [Num. weights, 1]
% y -- nodes of quadrature rule, size(y)=[Num. points, 3]

if nargin < 2
    d = 2;
end

if mod(L,2)==0
    L = L+1;
%     disp('L should be odd.');
%     return;
end

loadfpath = 'SD/';
if mod(L,2)~=0
    if L<10
        Ltxt = ['00' num2str(L)];
    elseif L<100
        Ltxt = ['0' num2str(L)];
    else
        Ltxt = num2str(L);
    end
    ld = [loadfpath 'ss' Ltxt '.mat'];
    load(ld,'y');
end
w = (4*pi)/size(y,1)*ones(size(y,1),1);
end


%GL Gauss-Legendre interpolatory quadrature rule
%   X = GL(N) generates the (2N+2)*(N+1) Gauss-Legendre nodes and returns a
%   2x[(2N+2)*(N+1)] matrix X containing their spherical coordinates. The first
%   row contains the longitudes in [0,2pi] and the second row the colatitudes in
%   [0,pi].
%
%   [X,W] = GL(N) in addition generates the quadrature weights W. The resulting
%   quadrature rule is exact up to polynomial degree 2*N.
%
%   Example:
%   [X,W] = GL(1) gives
%
%   X =
%        0    1.5708    3.1416    4.7124         0    1.5708    3.1416    4.7124
%   0.9553    0.9553    0.9553    0.9553    2.1863    2.1863    2.1863    2.1863
%
%   W =
%   1.5708    1.5708    1.5708    1.5708    1.5708    1.5708    1.5708    1.5708
%
%   See also cc, healpix, equidist
%
%   References
%   TODO Add references.
%
%   Copyright (c) 2002, 2016 Jens Keiner, Stefan Kunis, Daniel Potts

% Copyright (c) 2002, 2016 Jens Keiner, Stefan Kunis, Daniel Potts
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51
% Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
function [x,w] = gl(n)

% correctness conditions
if (nargin ~= 1)
  error('gl: Exactly one input argument required.');
end
if (~isscalar(n) || n < 0)
  error('gl: Input argument must be a non-negative number.');
end
if (nargout > 2)
  error('gl: No more than two output arguments allowed.');
end

% Gauss-Legendre nodes and weights (maybe) on [-1,1]
if (nargout == 2)
  [theta,w] = lgwt(n+1,-1,1);
else
  theta = lgwt(n+1,-1,1);
end
% coordinate transformation to [0,pi] for colatitude
theta = acos(theta);
% equispaced nodes for longitude
phi = 2*pi*(0:2*n+1)/(2*n+2);
% tensor product
[x,y] = meshgrid(theta,phi);
% linearisation to 2x((2*n+2)(n+1)) matrix
x = [y(:),x(:)]';
if (nargout == 2)
  % weights
  w = repmat((pi/(n+1))*w,1,2*n+2)';
  w = w(:)';
end
end

% The following function is based on code by Greg von Winckel, 02/25/2004
% See: TODO Add references.
function [x,w] = lgwt(n,a,b)

n = n-1; n1= n + 1; n2 = n + 2;
% n1 uniform nodes in [-1,1]
xu = linspace(-1,1,n1)';
% initial guess for nodes
y=cos((2*(0:n)'+1)*pi/(2*n+2)) + (0.27/n1)*sin(pi*xu*n/n2);
% Gauss-Legendre Vandermonde matrix
L=zeros(n1,n2);
% derivative of that
Lp=zeros(n1,n2);
% We compute the zeros of the n1th Legendre polynomial using the recursion
% relation and the Newton-Raphson method.
y0=2;
% Iterate until new points are uniformly within epsilon of old points.
while (max(abs(y - y0)) > eps)
  L(:,1) = 1; Lp(:,1) = 0; L(:,2) = y; Lp(:,2) = 1;
  for k = 2:n1
    L(:,k+1) = ((2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1))/k;
  end
  Lp = n2 * (L(:,n1) - y.*L(:,n2))./(1-y.^2); y0 = y; y = y0 - L(:,n2)./Lp;
end
% linear map from[-1,1] to [a,b]
x = (a*(1-y)+b*(1+y))/2;
if (nargout == 2)
  % weights
  w = (b-a)./((1-y.^2).*Lp.^2)*(n2/n1)^2;
end
end

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
end


function [theta,phi] = car2sph(xyz)
% [theta,phi]=Cartesian2Polar(xyz)
% converts xyz to polar coordinates, where size(xyz)=[N 3];
theta=acos(xyz(:,3));
x1=xyz(:,1);
x2=xyz(:,2);
x3=xyz(:,3);
phi=zeros(size(x1));
t=phi;
% x3==1
logic_x3=x3==1|x3==-1;
phi(logic_x3)=0;
% x2>=0 & -1<x3<1
logic_x2_1=x2>=0&x3<1&x3>-1;
t(logic_x2_1)=x1(logic_x2_1)./sqrt(1-(x3(logic_x2_1)).^2);
t(t>1)=1;
t(t<-1)=-1;
phi(logic_x2_1)=acos(t(logic_x2_1));
% x2<0 & x3<1
logic_x2_2=x2<0&x3<1&x3>-1;
t(logic_x2_2)=x1(logic_x2_2)./sqrt(1-(x3(logic_x2_2)).^2);
t(t>1)=1;
t(t<-1)=-1;
phi(logic_x2_2)=2*pi-acos(t(logic_x2_2));
end
