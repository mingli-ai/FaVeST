function Visualization(lam,th,Y_tar,Y_rec,error_Y)
%% Plotting the target field, approximated field, error field, respectively;
% Determine the max-norm of the field
% Copywright@Grady Wright, 2018, Boise State University
% Reference:
% Edward J Fuselier and Grady B Wright. 2009. Stability and error estimates for vector eld interpolation and decomposition on the sphere with RBFs.
% SIAM J. Numer. Anal. 47, 5 (2009), 3213?3239.
error_Y_refined = error_Y;
error_Y_refined(error_Y_refined<mean(mean(abs(error_Y)))) = 0;  
% for a better resolution in the following visualization
lon = lam*180/pi;
lat = th*180/pi;

[u_tar,v_tar] = convertTo2D(lam,th,Y_tar);
[u_rec,v_rec] = convertTo2D(lam,th,Y_rec);
[u_error,v_error] = convertTo2D(lam,th,error_Y);
[u_error_refined,v_error_refined] = convertTo2D(lam,th,error_Y_refined);

nm_tar = max(sqrt(sum(u_tar.^2+v_tar.^2,2)));
enm_tar = floor(log10(nm_tar));
fig_flow = figure;
m_proj('ortho','lat',0,'long',0);
m_grid('ytick',[-90 -60 -30 0 30 60 90],'xaxis','middle','xticklabels',[],'yticklabels',[],'linestyle','-','LineWidth',1,'Color',[0.2 0.2 0.2]);
hold on;
m_quiver(lon,lat,u_tar,v_tar,2,'k-');
ti = sprintf('$\\|T\\|_{\\infty} = %1.1f\\times 10^{%d}$',nm_tar/10^(enm_tar),enm_tar);
title(ti,'Interpreter','Latex');
set(gcf,'PaperPosition',[3 3 3 3]);

%------------------------------------------------------------------------------------------------------------------------------------------------
nm_rec = max(sqrt(sum(u_rec.^2+v_rec.^2,2)));
enm_rec = floor(log10(nm_rec));
fig_flow_rec = figure;
m_proj('ortho','lat',0,'long',0);
m_grid('ytick',[-90 -60 -30 0 30 60 90],'xaxis','middle','xticklabels',[],'yticklabels',[],'linestyle','-','LineWidth',1,'Color',[0.2 0.2 0.2]);
hold on;
m_quiver(lon,lat,u_rec,v_rec,2,'b-');
ti = sprintf('$\\|T^{\\rm rec}\\|_{\\infty} = %1.1f\\times 10^{%d}$',nm_rec/10^(enm_rec),enm_rec);
title(ti,'Interpreter','Latex');
set(gcf,'PaperPosition',[3 3 3 3]);

%------------------------------------------------------------------------------------------------------------------------------------------------
nm_error = max(sqrt(sum(u_error.^2+v_error.^2,2)));
enm_error = floor(log10(nm_error));
fig_flow_error = figure;
m_proj('ortho','lat',0,'long',0);
m_grid('ytick',[-90 -60 -30 0 30 60 90],'xaxis','middle','xticklabels',[],'yticklabels',[],'linestyle','-','LineWidth',1,'Color',[0.2 0.2 0.2]);
hold on;
m_quiver(lon,lat,u_error_refined,v_error_refined,2,'r-');
ti = sprintf('$\\|T-T^{\\rm rec}\\|_{\\infty}= %1.1f\\times 10^{%d}$',nm_error/10^(enm_error),enm_error);
title(ti,'Interpreter','Latex');
set(gcf,'PaperPosition',[3 3 3 3]);
end

function [f_u,f_v] = convertTo2D(lam,th,f)


% vectors for translating the field in Cartesian coordinates to a field
% in spherical coordinates.
c2s_u = [-sin(lam) -sin(th).*cos(lam)];
c2s_v = [cos(lam) -sin(th).*sin(lam)];
c2s_w = [zeros(size(lam)) cos(th)];



%project to 2D representation with u and v parts
f_u = c2s_u(:,1).*f(:,1) + c2s_v(:,1).*f(:,2) + c2s_w(:,1).*f(:,3);
f_v = c2s_u(:,2).*f(:,1) + c2s_v(:,2).*f(:,2) + c2s_w(:,2).*f(:,3);

end



