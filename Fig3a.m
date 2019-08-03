%% Compute and plot Figure 3(a)
clear,clc
close all
p1 = pwd;
addpath(genpath(p1))
fprintf('****** Reconstruction of Tangent Field by FaVeST ******\n')
% % Generate an simulated tangent field
ivf = 1;
switch ivf
    case 1
        vf = @flow1;
        vf_txt = 'Tangent Field A';
    case 2
        vf = @flow2;
        vf_txt = 'Tangent Field B';
    case 3
        vf = @flow3;
        vf_txt = 'Tangent Field C';
end
L = 30;
QH = 'SD';
fprintf('%s, Quadrature: %s, L: %d\n',vf_txt,QH,L)
[w_gl,x_gl,X_gl] = QpS2(L,QH);
[lam,th,r] = cart2sph(X_gl(1,:)',X_gl(2,:)',X_gl(3,:)');
Y_tar = vf(X_gl);
% % Running FaVeST_fwd and FaVeST_adj
% Fast evaluate Fourier coefficients for divergent-free and curl-free parts
[F1,F2] = FaVeST_fwd(Y_tar',L,x_gl,w_gl);
% Fast compute Fourier summation with the given Fourier coefficients
Y_rec = FaVeST_adj(F1,F2,x_gl);
Y_rec = real(Y_rec);
% Calculate the approximation error
err_pntwise = Y_rec-Y_tar;
err_abs = norm(err_pntwise);
err_rela = err_abs/norm(Y_tar);
fprintf('== Absolute Error: %.4e,  Relative Error: %.4e\n',err_abs,err_rela);
% Plotting
Visualization(lam,th,Y_tar,Y_rec,err_pntwise)