% Computes the spherical harmonics of order mu, as well as the Cartesian 
% derivatives D_x(Y), D_y(Y), and D_z(Y) of these functions.  Uses the 
% trivariate polynomial form of the spherical harmonics.
function [Y,Dx,Dy,Dz] = dsph(mu,x,y,z)

nd = length(x);

Dx = zeros(nd,2*mu+1); 
Dy = Dx; Dz = Dx; T = Dx; Tx = Dx; Ty = Dx; Tz = Tx;

c=sqrt((2*mu+1)/pi)/2;

for nu=0:mu
   n = floor((mu-nu)/2);
   k = 0:n;
   p = 0:nu;
   xm = repmat(x,[1 nu+1]);
   ym = repmat(y,[1 nu+1]);
   zm = repmat(z,[1 n+1]);
   pm = repmat(p,[nd 1]);
   km = repmat(k,[nd 1]);
   %C = (-1).^k.*binom(mu,k).*binom(2*(mu-k),mu).*factorial(mu-2*k)./factorial(mu-2*k-nu);
   temp = ones(size(k));
   for j=1:length(k)
      temp(j) = prod((mu-2*k(j)):-1:(mu-2*k(j)-nu+1));
   end
   C = (-1).^k.*binom(mu,k).*binom(2*(mu-k),mu).*temp;
   %D = sqrt(factorial(mu-nu)./factorial(mu+nu))./2.^mu;
   D = sqrt(prod(1./((mu+nu):-1:(mu-nu+1))))./2.^mu;
   %E = factorial(nu)./(factorial(p).*factorial(nu-p));
   E = binom(nu,p);
   A = repmat(E,[nd 1]).*(xm.^pm).*(ym.^(nu-pm));
   B = sum(A.*repmat(sin((nu-p)*pi/2),[nd 1]),2);
   A = sum(A.*repmat(cos((nu-p)*pi/2),[nd 1]),2);

   idp1 = 2:nu+1;
   if nu==0
      Ax = zeros(nd,1);
      Bx = Ax;
   else
      Ax = repmat(E(idp1).*p(idp1),[nd 1]).*(xm(:,idp1).^(pm(:,idp1)-1)).*(ym(:,idp1).^(nu-pm(:,idp1)));
      Bx = sum(Ax.*repmat(sin((nu-p(idp1))*pi/2),[nd 1]),2);
      Ax = sum(Ax.*repmat(cos((nu-p(idp1))*pi/2),[nd 1]),2);
   end

   idm1 = 1:nu;
   if nu==0
      Ay = zeros(nd,1);
      By = Ay;
   else
      Ay = repmat(E(idm1).*(nu-p(idm1)),[nd 1]).*(xm(:,idm1).^(pm(:,idm1))).*(ym(:,idm1).^(nu-pm(:,idm1)-1));
      By = sum(Ay.*repmat(sin((nu-p(1:end-1))*pi/2),[nd 1]),2);
      Ay = sum(Ay.*repmat(cos((nu-p(1:end-1))*pi/2),[nd 1]),2);
   end

%    idp1 = 2:n+1;
   P = D*(repmat(C,[nd 1]).*zm.^(mu-nu-2*km));
%    Px = sum(P(:,idp1).*repmat(k(idp1),[nd 1]),2);
%    Py = 2*y.*Px;
%    Pz = 2*z.*Px;
%    Px = 2*x.*Px;
   idm1 = 1:n;
   id = 1:n+1;
   if mu-nu-1 == 0
      %Pz = Pz + D*sum(repmat(C(1).*(mu-nu),[nd 1]),2);
      Pz = D*sum(repmat(C(1).*(mu-nu),[nd 1]),2);
   else
      if mod(mu-nu,2) == 1
         %Pz = Pz + D*sum(repmat(C(id).*(mu-nu-2*k(id)),[nd 1]).*zm(:,id).^(mu-nu-2*km(:,id)-1),2);
         Pz = D*sum(repmat(C(id).*(mu-nu-2*k(id)),[nd 1]).*zm(:,id).^(mu-nu-2*km(:,id)-1),2);
      else
         %Pz = Pz + D*sum(repmat(C(idm1).*(mu-nu-2*k(idm1)),[nd 1]).*zm(:,idm1).^(mu-nu-2*km(:,idm1)-1),2);
         Pz = D*sum(repmat(C(idm1).*(mu-nu-2*k(idm1)),[nd 1]).*zm(:,idm1).^(mu-nu-2*km(:,idm1)-1),2);
      end
   end
   P = sum(P,2);

   T(:,nu+1) = c*(-1)^nu*P.*(A + 1i*B);
   %Tx(:,nu+1) = c*(-1)^nu*(Px.*(A + i*B) + P.*(Ax + i*Bx));
   %Ty(:,nu+1) = c*(-1)^nu*(Py.*(A + i*B) + P.*(Ay + i*By));
   %Tz(:,nu+1) = c*(-1)^nu*Pz.*(A + i*B);

   Tx(:,nu+1) = c*(-1)^nu*(P.*(Ax + 1i*Bx));
   Ty(:,nu+1) = c*(-1)^nu*(P.*(Ay + 1i*By));
   Tz(:,nu+1) = c*(-1)^nu*Pz.*(A + 1i*B);
end
Dx(:,1:nu) = imag(Tx(:,nu+1:-1:2));
Dx(:,nu+1) = real(Tx(:,1));
Dx(:,nu+2:(2*nu+1)) = real(Tx(:,2:nu+1));

Dy(:,1:nu) = imag(Ty(:,nu+1:-1:2));
Dy(:,nu+1) = real(Ty(:,1));
Dy(:,nu+2:(2*nu+1)) = real(Ty(:,2:nu+1));

Dz(:,1:nu) = imag(Tz(:,nu+1:-1:2));
Dz(:,nu+1) = real(Tz(:,1));
Dz(:,nu+2:(2*nu+1)) = real(Tz(:,2:nu+1));

Y(:,1:nu) = imag(T(:,nu+1:-1:2));
Y(:,nu+1) = real(T(:,1));
Y(:,nu+2:(2*nu+1)) = real(T(:,2:nu+1));
end
