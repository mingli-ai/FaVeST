function [F1,F2] = FaVeST_fwd(T,L,X,w)
% Forward FFTs for vector spherical harmonics
% INPUTS:
%  T      - scalar field or component of vector field
%  L,m    - degree and order for vector spherical harmonic
%  X,w    - quadrature rule for evaluating FFT
% OUTPUTS:
%  F1,F2  - Fourier coefficients for divergence-free and curl-free vector spherical harmonics of degree L,m


% scalar Fourier coeff
That1 = fftS2_fwd(T(1,:),X,L+1,w);                 % of size (2L+3)*(L+2)
That2 = fftS2_fwd(T(2,:),X,L+1,w);                 % of size (2L+3)*(L+2)
That3 = fftS2_fwd(T(3,:),X,L+1,w);                 % of size (2L+3)*(L+2)
L1 = L + 2; % number of columns of F1,F2
L0 = L1; % mid line of scalar fft for That1, That2, That3

F1 = zeros(2*L+1,L); % lower triangle stands for m>0
F2 = zeros(2*L+1,L); % lower triangle stands for m>0
Lmid = L+1;          % the mid line for FaVeST: F1, F2
 for l=1:L
     for m=0:l
         if m==0
             % for m=0
            F1(Lmid+m,l) = 1/sqrt(2)*(xi1(l-1,m-1)*(-1*That1(L0+(m-1),(l-1)+1)+1i*That2(L0+(m-1),(l-1)+1)) ...
                                       + xi2(l+1,m-1)*(-1*That1(L0+(m-1),(l+1)+1)+1i*That2(L0+(m-1),(l+1)+1)) ...
                                       + xi3(l-1,m+1)*(-1)^(m+1)*(That1(L0+(m+1),(l-1)+1)+1i*That2(L0+(m+1),(l-1)+1)) ...
                                       + xi4(l+1,m+1)*(-1)^(m+1)*(That1(L0+(m+1),(l+1)+1)+1i*That2(L0+(m+1),(l+1)+1))) ...
                                       + xi5(l-1,m)*That3(L0+m,(l-1)+1) ...
                                       + xi6(l+1,m)*That3(L0+m,(l+1)+1);  
            F2(Lmid+m,l) = -1/sqrt(2)*1i*(mu1(l,m-1)*(-1*That1(L0+(m-1),l+1)+1i*That2(L0+(m-1),l+1)) ...
                                         + mu3(l,m+1)*(-1)^(m+1)*(That1(L0+(m+1),l+1)+1i*That2(L0+(m+1),l+1))) ...
                                         - 1i*mu2(l,m)*That3(L0+m,l+1);

         else
       % for m>0
            F1(Lmid+m,l) = 1/sqrt(2)*(xi1(l-1,m-1)*(-1)^(m-1)*(-1*That1(L0+(m-1),(l-1)+1)+1i*That2(L0+(m-1),(l-1)+1)) ...
                                   + xi2(l+1,m-1)*(-1)^(m-1)*(-1*That1(L0+(m-1),(l+1)+1)+1i*That2(L0+(m-1),(l+1)+1)) ...
                                   + xi3(l-1,m+1)*(-1)^(m+1)*(That1(L0+(m+1),(l-1)+1)+1i*That2(L0+(m+1),(l-1)+1)) ...
                                   + xi4(l+1,m+1)*(-1)^(m+1)*(That1(L0+(m+1),(l+1)+1)+1i*That2(L0+(m+1),(l+1)+1))) ...
                                   + xi5(l-1,m)*(-1)^m*That3(L0+m,(l-1)+1) ...
                                   + xi6(l+1,m)*(-1)^m*That3(L0+m,(l+1)+1);

                          
            F2(Lmid+m,l) = -1/sqrt(2)*1i*(mu1(l,m-1)*(-1)^(m-1)*(-1*That1(L0+(m-1),l+1)+1i*That2(L0+(m-1),l+1)) ...
                                         + mu3(l,m+1)*(-1)^(m+1)*(That1(L0+(m+1),l+1)+1i*That2(L0+(m+1),l+1))) ...
                                         - 1i*mu2(l,m)*(-1)^m*That3(L0+m,l+1);
         end
     end
     for m = -l:-1
         % upper triangle of F1/F2 corresponds to m<0
         F1(Lmid+m,l) = (-1)^(abs(m))*conj(F1(Lmid+(-m),l));
         F2(Lmid+m,l) = (-1)^(abs(m))*conj(F2(Lmid+(-m),l));
     end
     F1(:,l) = F1(:,l)/sqrt(l*(l+1));
     F2(:,l) = F2(:,l)/sqrt(l*(l+1));
 end
end
%% subfunctions
function f_hat = fftS2_fwd(f,X,N,w)
%% evaluate scalar spherical harmonic coefficient by nfsft
% 
%        f_hat = fftS2(f,X,N,w)
%
% INPUT:
%   X     - X  = (x1,...,x_M), points on the 2-sphere, where each
%           x_k=(\theta,\phi) with \theta the longitudes in [0,2pi] and 
%           \phi the colatitudes in [0,pi].
%
%   f     - vector (f(x_1),...,f(x_M))
%
%   N     - bandwidth, see f_hat below
%
%   w     - weight on f(x_j), default = 1
%
% OUTPT:  
%   f_hat - matrix size is (2N+1) x (N+1), w.r.t.  
%           f_hat(l,m) = \sum_{k=1}^M f(x_k) bar{Y}_{l,m}(x_k);
%           for l = 0,..,N and m=-l,..,l.
%%
if length(f) ~= size(X,2)
    error('Input signal f must have the same length as X!');
end

if nargin  == 3
    w = ones(size(f));
end

M = size(X,2);

nfsft_precompute(N,1000);                             % precomputation

plan = nfsft_init_advanced(N,M,NFSFT_NORMALIZED);     % Create plan.

nfsft_set_x(plan,X);                                  % Set nodes.

nfsft_precompute_x(plan);                             % node-dependent precomputation

nfsft_set_f(plan,f.*w);                               % set f

nfsft_adjoint(plan);                                  % adjoint transform          

f_hat = nfsft_get_f_hat(plan);                        % get f_hat
 
nfsft_finalize(plan);                                 % release memory
end
function [xi_1] = xi1(l,m)
    if m>=-l&&m<=l
        xi_1 = c_l(l+1)*sqrt((l+1+m)*(l+m+2)/((2*l+2)*(2*l+1)));  %C_{l,m,1,1}^{l+1,m+1}
    else
        xi_1=0;
    end
end

function [xi_2] = xi2(l,m)
     if m>=-l&&m<=l
        xi_2 = d_l(l-1)*sqrt((l-m-1)*(l-m)/((2*l+1)*(2*l))); %C_{l,m,1,1}^{l-1,m+1}
     else
        xi_2 = 0;
     end
end

function [xi_3] = xi3(l,m)
    if m>=-l&&m<=l
      xi_3 = c_l(l+1)*sqrt((l-m+2)*(l-m+1)/((2*l+2)*(2*l+1))); %C_{l,m,1,1}^{l+1,m-1}
    else
      xi_3 = 0;
    end
end

function [xi_4] = xi4(l,m)
    if m>=-l&&m<=l
        xi_4 =d_l(l-1)*sqrt((l+m)*(l+m-1)/((2*l+1)*(2*l))); %C_{l,m,1,1}^{l-1,m-1}
    else
        xi_4 = 0;
    end
end

function [xi_5] = xi5(l,m)
    if m>=-l&&m<=l
        xi_5 = c_l(l+1)*sqrt((l+m+1)*(l-m+1)/((l+1)*(2*l+1))); %C_{l,m,1,0}^{l+1,m}
    else
        xi_5 = 0;
    end
end

function [xi_6] = xi6(l,m)
    if m>=-l&&m<=l
        xi_6 = -d_l(l-1)*sqrt((l+m)*(l-m)/((2*l+1)*(l))); %C_{l,m,1,0}^{l-1,m}
    else
        xi_6 = 0;
    end
end

function [mu_1] = mu1(l,m)
    if m>=-l&&m<=l
        lambda_l=l*(l+1);
        mu_1 = - sqrt(lambda_l)*sqrt((l+m+1)*(l-m)/(l*(2*l+2))); %C_{l,m,1,1}^{l,m+1}
    else
        mu_1 = 0;
    end
end

function [mu_2] = mu2(l,m)
    if m>=-l&&m<=l
        lambda_l=l*(l+1);
        mu_2=sqrt(lambda_l)*m/sqrt(l*(l+1)); %C_{l,m,1,0}^{l,m}
    else
        mu_2 = 0;
    end
end

function [mu_3] = mu3(l,m)
    if m>=-l&&m<=l
        lambda_l=l*(l+1);
        mu_3=sqrt(lambda_l)*sqrt((l+m)*(l-m+1)/(l*(2*l+2))); %C_{l,m,1,-1}^{l,m-1}
    else
        mu_3 = 0;
    end
end

function y = c_l(l)
    y = (l+1)*sqrt(l/(2*l+1));
end

function y = d_l(l)
    y = l*sqrt((l+1)/(2*l+1));
end

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