function [f] = FaVeST_adj(alm,blm,X)
% Adjoint FFTs for vector spherical harmonics
% INPUTS:
%  alm    - Fourier coefficients for divergence-free part
%  blm    - Fourier coefficients for curl-free part
%  X      - quadrature rule points for evaluating FFT
% OUTPUT:
%  f      - Fourier partial sum
for l=1:size(alm,2)
    alm(:,l) = alm(:,l)/sqrt(l*(l+1));
    blm(:,l) = blm(:,l)/sqrt(l*(l+1));
end
L = size(alm,2);
L1=L-1;
nu_1 = zeros(2*L1+1,L1+1);
for l = 0:L1
    for m = -l:l
       nu_1(m+L1+1,l+1) = alm(m+L+2,l+1)*xi1(l,m);
    end
end
%% 
L2 = L+1;
nu_2 = zeros(2*L2+1,L2+1);
for l = 2:L2
    for m=-l:(l-2)
       nu_2(m+L2+1,l+1) = alm(m+L+2,l-1)*xi2(l,m);
    end
end
%%
L3 = L-1;
nu_3 = zeros(2*L3+1,L3+1);
for l = 0:L3
    for m = -l:l
        nu_3(m+L3+1,l+1) = alm(m+L,l+1)*xi3(l,m);
    end
end
%%
L4 = L+1;
nu_4 = zeros(2*L4+1,L4+1);
for l = 2:L4
    for m = (-l+2):l
         nu_4(m+L4+1,l+1) = alm(m+L,l-1)*xi4(l,m);
    end
end
%%
L5 = L-1;
nu_5 = zeros(2*L5+1,L5+1);
for l = 0:L5
    for m=-l:l
        nu_5(m+L5+1,l+1) = alm(m+L+1,l+1)*xi5(l,m);
    end
end
%%
L6 = L+1;
nu_6 = zeros(2*L6+1,L6+1);
for l = 2:L6
    for m = (-l+1):(l-1)
        nu_6(m+L6+1,l+1) = alm(m+L+1,l-1)*xi6(l,m);
    end
end
%%
eta_1 = zeros(2*L+1,L+1);
for l = 1:L
    for m = -l:(l-1)
        eta_1(m+L+1,l+1) = 1i*blm(m+L+2,l)*mu1(l,m);
    end
end
%%
eta_2 = zeros(2*L+1,L+1);
for l = 1:L
    for m = (-l+1):l
        eta_2(m+L+1,l+1) = 1i*blm(m+L,l)*mu3(l,m);
    end
end
%%
eta_3 = zeros(2*L+1,L+1);
for l = 1:L
    for m = -l:l
        eta_3(m+L+1,l+1) = 1i*blm(m+L+1,l)*mu2(l,m);
    end
end
%%
nu1mnu3 = nu_1-nu_3;
nu2mnu4 = nu_2-nu_4;
eta1meta2 = eta_1-eta_2; 
%-------post process for nu1m3/nu2m4/eta1m2 
nu1mnu3=[nu1mnu3(1:L,:);flipud(conj(nu1mnu3(1:L-1,:)))];
nu2mnu4=[nu2mnu4(1:L+2,:);flipud(conj(nu2mnu4(1:L+1,:)))];
eta1meta2 =[eta1meta2(1:L+1,:);flipud(conj(eta1meta2(1:L,:)));];
%-----------------------------------------------------------------------
nu1pnu3 = 1i*(nu_1+nu_3);
nu2pnu4 = 1i*(nu_2+nu_4);
eta1peta2 = 1i*(eta_1+eta_2); 
%post process for nu1p3/nu2p4/eta1p2
nu1pnu3=[nu1pnu3(1:L,:);flipud(conj(nu1pnu3(1:L-1,:)))];
nu2pnu4=[nu2pnu4(1:L+2,:);flipud(conj(nu2pnu4(1:L+1,:)))];
eta1peta2 =[eta1peta2(1:L+1,:);flipud(conj(eta1peta2(1:L,:)))];
%-----------------------------------------------------------------------
%post process for nu_5/nu_6/eta_3
nu_5=[nu_5(1:L,:);flipud(conj(nu_5(1:L-1,:)))]; 
nu_6=[nu_6(1:L+2,:);flipud(conj(nu_6(1:L+1,:)));];
eta_3 =[eta_3(1:L+1,:);flipud(conj(eta_3(1:L,:)));];
%%
f1 = -1/sqrt(2)*(fftS2_adj(nu1mnu3,X) +...
               fftS2_adj(nu2mnu4,X) +...
               fftS2_adj(eta1meta2,X));
f2 = -1/sqrt(2)*(fftS2_adj(nu1pnu3,X) +...
               fftS2_adj(nu2pnu4,X) +...
               fftS2_adj(eta1peta2,X));
f3 = fftS2_adj(nu_5,X) +...
     fftS2_adj(nu_6,X) +...
     fftS2_adj(eta_3,X);
 
f = [f1,f2,f3];

%% subfunctions
function f = fftS2_adj(f_hat,X)
%% Adjoint FFT of nfsft for scalar case
% 
%  f = fftS2_adj(f_hat,X)
%
% INPUT:
%   f_hat - matrix size is (2N+1) x (N+1), w.r.t.  
%           f(x_k) = \sum_{l=0}^N\sum_{m=-l}^l f_hat(l,m)Y_{l,m}(x_k);          
%
%   X     - X  = (x1,...,x_M), points on the 2-sphere, where each
%           x_k=(\theta,\phi) with \theta the longitudes in [0,2pi] and 
%           \phi the colatitudes in [0,pi].
%
% OUTPT:
%   f     - f(x_k) = \sum_{l=0}^N\sum_{m=-l}^l f_hat(l,m)Y_{l,m}(x_k);

%%
N    = size(f_hat,2)-1;   

M    = size(X,2);

nfsft_precompute(N,1000);                             % precomputation

plan = nfsft_init_advanced(N,M,NFSFT_NORMALIZED);     % Create plan.

nfsft_set_x(plan,X);                                  % Set nodes.

nfsft_precompute_x(plan);                             % node-dependent precomputation

nfsft_set_f_hat(plan,f_hat);                          % set f_hat

nfsft_trafo(plan);                                    % forward transform          

f = nfsft_get_f(plan);                                % get f
 
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



end