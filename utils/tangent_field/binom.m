% Computes the binomial coefficient m choose k. 
function b = binom(m,k)

nm = length(m);
nk = length(k);

if nm > 1 & nk == 1
   for j=1:nm
      idn = m(j):-1:(m(j)-k+1);
      idd = ones(size(idn));
      idd(1:k) = k:-1:1;
      b(j) = prod(idn./idd);
   end
elseif nk > 1 & nm == 1
   for j=1:nk
      idn = m:-1:(m-k(j)+1);
      idd = ones(size(idn));
      idd(1:k(j)) = k(j):-1:1;
      b(j) = prod(idn./idd);
   end
elseif nk == 1 & nm == 1
   for j=1:nk
      idn = m:-1:(m-k+1);
      idd = ones(size(idn));
      idd(1:k) = k:-1:1;
      b(j) = prod(idn./idd);
   end
else
   b = factorial(m)./(factorial(m-k).*factorial(k));
end
