function [f,f_u,f_v] = flow1(X)

x = X(1,:)';
y = X(2,:)';
z = X(3,:)';
[la,th,~] = cart2sph(x,y,z);
% vectors for translating the field in Cartesian coordinates to a field
% in spherical coordinates.
c2s_u = [-sin(la) -sin(th).*cos(la)];
c2s_v = [cos(la) -sin(th).*sin(la)];
c2s_w = [zeros(size(la)) cos(th)];

% Orders of the spherical harmonics for the curl-free field.
mu = [6 4];
nu = [-3 0];
alp = [1 1]/25;

Yx = zeros(size(x)); Yy = Yx; Yz = Yx;
for j=1:length(mu)
   [~,Yxt,Yyt,Yzt] = dsph(mu(j),x,y,z);
   Yx = Yx + alp(j)*Yxt(:,mu(j)+nu(j)+1);
   Yy = Yy + alp(j)*Yyt(:,mu(j)+nu(j)+1);
   Yz = Yz + alp(j)*Yzt(:,mu(j)+nu(j)+1);
end

% Project the vectors onto the sphere.  This will give the curl-free field;
curl_free_dimx = ( (1 - x.^2).*Yx - x.*y.*Yy - x.*z.*Yz);
curl_free_dimy = (-x.*y.*Yx + (1 - y.^2).*Yy - y.*z.*Yz);
curl_free_dimz = (-x.*z.*Yx - y.*z.*Yy + (1 - z.^2).*Yz);

curl_free = [curl_free_dimx,curl_free_dimy,curl_free_dimz];

%project to 2D representation with u and v parts
curl_free_u = c2s_u(:,1).*curl_free_dimx + c2s_v(:,1).*curl_free_dimy + c2s_w(:,1).*curl_free_dimz;
curl_free_v = c2s_u(:,2).*curl_free_dimx + c2s_v(:,2).*curl_free_dimy + c2s_w(:,2).*curl_free_dimz;

% Orders of the spherical harmonics for the div-free field.
%

mu = [1 5];
nu = [0 4];
alp = 1./[-sqrt(3) 3/8*sqrt(385/2)];

Yx = zeros(size(x)); Yy = Yx; Yz = Yx;
for j=1:length(mu)
   [~,Yxt,Yyt,Yzt] = dsph(mu(j),x,y,z);
   Yx = Yx + alp(j)*Yxt(:,mu(j)+nu(j)+1);
   Yy = Yy + alp(j)*Yyt(:,mu(j)+nu(j)+1);
   Yz = Yz + alp(j)*Yzt(:,mu(j)+nu(j)+1);
end

% Project the vectors onto the sphere.  This will give the div-free field;
diver_free_dimx = -z.*Yy + y.*Yz;
diver_free_dimy =  z.*Yx - x.*Yz;
diver_free_dimz = -y.*Yx + x.*Yy;
diver_free = [diver_free_dimx,diver_free_dimy,diver_free_dimz];

diver_free_u = c2s_u(:,1).*diver_free_dimx + c2s_v(:,1).*diver_free_dimy + c2s_w(:,1).*diver_free_dimz;
diver_free_v = c2s_u(:,2).*diver_free_dimx + c2s_v(:,2).*diver_free_dimy + c2s_w(:,2).*diver_free_dimz;

f = curl_free+diver_free;

f_u = curl_free_u+diver_free_u;
f_v = curl_free_v+diver_free_v;

end

