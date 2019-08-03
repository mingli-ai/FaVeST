%% Program for Table 2 about the CPU time vs Degree L

clear,clc
close all

p1 = pwd;
addpath(genpath(p1))

vf = @flow1;  % Choice of vf includes flow1, flow2, flow3  
Lv = 250:250:2250;
num_trial = 5;
time_fwd = zeros(numel(Lv),num_trial);
time_adj = zeros(numel(Lv),num_trial);
fprintf('****** CPU Time for FaVeST ******\n')
for i= 1:size(Lv,2)
    [x_gl,w_gl] = gl(Lv(i));
    X_gl = sph2car(x_gl(2,:),x_gl(1,:));
    Y_tar = vf(X_gl);
    Y_tar = Y_tar';
    fprintf('*** Degree L up to %d ***\n',Lv(i))
    
    fprintf('* Forward FaVeST\n')
    for j = 1:num_trial
        t0_fwd = tic;
        [F1,F2] = FaVeST_fwd(Y_tar,Lv(i),x_gl,w_gl);
        time_fwd(i,j) = toc(t0_fwd);
        fprintf(' - %d-th trial, time: %.4e\n',j,time_fwd(i,j))
    end
    fprintf('* Adjoint FaVeST ***\n')
    for j = 1:num_trial
        t0_adj = tic;
        FaVeST_adj(F1,F2,x_gl);
        time_adj(i,j) = toc(t0_adj);
        fprintf(' - %d-th trial, time: %.4e\n',j,time_adj(i,j))
    end
end

%% take mean
time_fwd_1 = sum(time_fwd(:,2:num_trial),2)/(num_trial-1);
Time_adj_1 = sum(time_adj(:,2:num_trial),2)/(num_trial-1);
Time_total = time_fwd_1+Time_adj_1;
ratio_fwd = time_fwd_1(2:end)./time_fwd_1(1:end-1);
ratio_adj = Time_adj_1(2:end)./Time_adj_1(1:end-1);
ratio_total = Time_total(2:end)./Time_total(1:end-1);

% save data
clear x_gl X_gl w_gl Y_tar
sv_dat = ['favest_cpu_time_' num2str(Lv(end)) '.mat'];
save(sv_dat)

for i= 1:size(Lv,2)
    [x_gl,w_gl] = gl(Lv(i));
    Nv(i)=size(x_gl,2);
end
for i= 1:size(Lv,2)
    Mv(i)=Lv(i)*(Lv(i)+2);
end
%% fit and plot
colors=[0,255,255;  %cyan
        255,0,0;  %red
        153,0,0; % OU Crimson Red
        77,255,166;  % medium aquamarine
        0,179,89; % Green (pigment)
        102,204,25;  %Maya blue
        0,77,230; % Blue (RYB)
        106,0,128; %Patriarch (purple)
        0,213,255]./255; % Capri (blue)

% fwd
[p_fwd,pstr_fwd] = fitpow(Nv(2:end),time_fwd_1(2:end),'N');
lgstr_fwd = ['Forward'];
Nv1 = linspace(0.9*Nv(1),1.1*Nv(end),100);
time_fwd_fit = p_fwd(1)*(Nv1).^(p_fwd(2));
fig_fwd = figure;
loglog(Nv,time_fwd_1,'p',Nv1,time_fwd_fit,'-','MarkerSize',3,'LineWidth',2)
lg_fwd = legend(lgstr_fwd,pstr_fwd);
set(lg_fwd,'interpreter','latex','location','North','Fontsize',14)
grid on
xlim([6*10^4 2*10^7])
ylim([10^(-1) 3*10^2])
xlabel('Number $N$ of quadrature nodes','interpreter','latex','Fontsize',12)
ylabel('CPU time (s)','interpreter','latex','Fontsize',12)
set(gca,'Fontsize',12)
sv_fig_fwd = ['time_favest_fwd_to_2250' '.eps'];
print(sv_fig_fwd,'-depsc2',fig_fwd)

% adj
[p_adj,pstr_adj] = fitpow(Mv(2:end),Time_adj_1(2:end),'M');
lgstr_adj = ['Adjoint'];
Mv1 = linspace(0.9*Mv(1),1.1*Mv(end),100);
time_adj_fit = p_adj(1)*(Mv1).^(p_adj(2));
fig_adj = figure;
loglog(Mv,Time_adj_1,'p',Mv1,time_adj_fit,'-','MarkerSize',3,'LineWidth',2)
lg_adj = legend(lgstr_adj,pstr_adj);
set(lg_adj,'interpreter','latex','location','North','Fontsize',14)
grid on
xlim([10^4 10^7])
ylim([10^(-1) 3*10^2])
xlabel('Number $M$ of coefficients','interpreter','latex','Fontsize',12)
ylabel('CPU time (s)','interpreter','latex','Fontsize',12)
set(gca,'Fontsize',12)
sv_fig_adj = ['time_favest_adj_to_2250' '.eps'];
print(sv_fig_adj,'-depsc2',fig_adj)
    
