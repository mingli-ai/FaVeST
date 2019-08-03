function [f,f_u,f_v] = flow2(X)

x = X(1,:)';
y = X(2,:)';
z = X(3,:)';
[la,th,~] = cart2sph(x,y,z);
% Curl-free part is generated from a b-spline velocity potential
%
sig = [3 3.5 3 3.5];
psi0 = [-1/7 1/8 -1/8 1/9];

la0 = [-pi/7 0 pi/3 pi/2]; th0 = [pi/5 pi/6 -pi/5 -pi/6];
[x0,y0,z0] = sph2cart(la0,th0,ones(size(la0)));
% vectors for translating the field in Cartesian coordinates to a field
% in spherical coordinates.
c2s_u = [-sin(la) -sin(th).*cos(la)];
c2s_v = [cos(la) -sin(th).*sin(la)];
c2s_w = [zeros(size(la)) cos(th)];

% Cubic b-spline
ucfc = zeros(size(th)); 
vcfc = ucfc;
wcfc = ucfc;

for k=1:length(la0)
   rd = sqrt(2*(1-(x0(k)*x + y0(k)*y + z0(k)*z)));
   t1 = rd;
   id = abs(t1) < 2/sig(k);
   t1(~id) = 0;

   temp = zeros(size(th)); cx = temp; cy = temp; cz = temp;
   temp(id) = rd(id);
   temp(abs(th-th0(k))<5*eps & abs(la-la0(k))<5*eps ) = 1;
   cx(id) = (x(id)-x0(k))./temp(id);
   cy(id) = (y(id)-y0(k))./temp(id);
   cz(id) = (z(id)-z0(k))./temp(id);

   ucfct = zeros(size(th)); 
   vcfct = ucfct;
   wcfct = ucfct;
   temp =  zeros(size(th));
   temp1 = zeros(size(th,1),5);
   coeff = [1 -4 6 -4 1];
   for j=1:5
      temp1(id,j) = coeff(j)*(t1(id)-(3-j)/sig(k)).^2.*sign(t1(id)-(3-j)/sig(k));
   end
   temp(id) = sum(sort(temp1(id,:),2),2);

   ucfct(id) = psi0(k)*sig(k)^3/4*cx(id).*temp(id);
   vcfct(id) = psi0(k)*sig(k)^3/4*cy(id).*temp(id);
   wcfct(id) = psi0(k)*sig(k)^3/4*cz(id).*temp(id);

   ucfc = ucfc + ucfct;
   vcfc = vcfc + vcfct;
   wcfc = wcfc + wcfct;
end

% Project the vectors onto the sphere.  This will give the curl-free field;
curl_free_dimx = ( (1 - x.^2).*ucfc - x.*y.*vcfc - x.*z.*wcfc);
curl_free_dimy = (-x.*y.*ucfc + (1 - y.^2).*vcfc - y.*z.*wcfc);
curl_free_dimz = (-x.*z.*ucfc - y.*z.*vcfc + (1 - z.^2).*wcfc);
curl_free = [curl_free_dimx,curl_free_dimy,curl_free_dimz];

curl_free_u = c2s_u(:,1).*curl_free_dimx + c2s_v(:,1).*curl_free_dimy + c2s_w(:,1).*curl_free_dimz;
curl_free_v = c2s_u(:,2).*curl_free_dimx + c2s_v(:,2).*curl_free_dimy + c2s_w(:,2).*curl_free_dimz;
%
% Div free part is just a Rosby-Haurwitz wave.
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

% Project the vectors onto the sphere.  This will give the curl-free field;
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