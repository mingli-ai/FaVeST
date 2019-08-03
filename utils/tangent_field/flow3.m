function [f,f_u,f_v] = flow3(X)
x = X(1,:)';
y = X(2,:)';
z = X(3,:)';
[la,th,~] = cart2sph(x,y,z);

% Center of the "bell"
lac = -pi/12;
thc = pi/4;
alp = -3;

udf = zeros(size(la)); vdf = udf; ucf = udf; vcf = udf; 

% Do the div-free first
for j=1:size(lac,1)
   t0 = thc(j);
   l0 = lac(j);
   for i=1:max(size(la))
      lai = la(i);
      thi = th(i);
      q = (cos(t0).*cos(lai-l0).*cos(thi) + sin(thi).*sin(t0));
      dql = -cos(t0).*sin(lai-l0);
      dqt = -cos(t0).*cos(lai-l0).*sin(thi) + cos(thi).*sin(t0);
      if q < 1-8*eps && q > -(1-8*eps)
         temp = 1/(1+q)*(-1-3*q+sqrt(2-2*q)*(2+3*q)+log((1-q)/(1+sqrt(2-2*q)-q))*(2-q-3*q^2));
      elseif q <= -(1-8*eps)
         temp = 7/2-log(32);
      else
         temp = 0;
      end
      udf(i) = udf(i) + alp(j)*(-temp.*dqt);
      vdf(i) = vdf(i) + alp(j)*(temp.*dql);
   end
end

% Add jet to the div-free flow.
udf = udf + sin(2*th).^14;

lac = [pi/9 0 pi/10];
thc = [pi/6 pi/4 5*pi/16];
alp = [-7/4 5/2 -3/2];

for j=1:max(size(lac))
   t0 = thc(j);
   l0 = lac(j);
   for i=1:max(size(la))
      lai = la(i);
      thi = th(i);
      q = cos(t0).*cos(lai-l0).*cos(thi) + sin(thi).*sin(t0);
      dql = -cos(t0).*sin(lai-l0);
      dqt = -cos(t0).*cos(lai-l0).*sin(thi) + cos(thi).*sin(t0);
      if q < 1-8*eps && q > -(1-8*eps)
         temp = 1/(1+q)*(-1-3*q+sqrt(2-2*q)*(2+3*q)+log((1-q)/(1+sqrt(2-2*q)-q))*(2-q-3*q^2));
      elseif q <= -(1-8*eps)
         temp = 7/2-log(32);
      else
         temp = 0;
      end
      ucf(i) = ucf(i) + alp(j)*(temp.*dql);
      vcf(i) = vcf(i) + alp(j)*(temp.*dqt);
   end
end

u = udf + ucf;
v = vdf + vcf;
f_u = u;
f_v = v;
f1 = -sin(la).*u-cos(la).*sin(th).*v;
f2 = cos(la).*u-sin(la).*sin(th).*v;
f3 = cos(th).*v;

f=[f1 f2 f3];

end
