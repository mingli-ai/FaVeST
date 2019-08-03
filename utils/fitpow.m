function [p, pstr] = fitpow(x, y, ipr)
% [p, pstr] = fitpow(x, y, ipr)
% Find the least squares fit of the form y = p(1)*x^p(2)
% If ipr > 0 print one line summary (defualt ipr = 0)
% If two output arguments, pstr is a string of fit

if nargin < 3
    ipr = '';
end;

x =x (:);
y = y(:);

lx = length(x);
ly = length(y);
if lx ~= ly
    fprintf('Error in fitpow: length(x) = %d, length(y) = % must be equal', lx, ly);
    return;
end;

A = [ones(size(x)) log(x)];
p = A \ log(y);
p(1) = exp(p(1));

if nargout > 1
   if isempty(ipr)
       pstr = sprintf('%5.1e L^{%.1f}', p); % scientific notation
   elseif ipr == '\sigma' 
       pstr = sprintf('%.1f %s^{%.1f}', p(1),ipr,p(2));
   else
       pstr = sprintf('%1.1e$%s^{%.1f}$', p(1),ipr,p(2));
   end
    % compute for 3 significant figures
%     num_sig = 3; % number of significant figures %%%%%%
%     p1_log10 = floor(log10(p(1))); % order of magnitude
%     if p1_log10 < num_sig - 1
%         pstr = sprintf('%3.2f L^{%.1f}', p);
%     elseif p1_log10 == num_sig - 1
%         pstr = sprintf('%3.0f L^{%.1f}', p);
%     else                
%         num_div = p1_log10 - num_sig + 1;
%         p1 = round(p(1)/10^(num_div))*10^num_div;
%         pstr = sprintf('%5.0f L^{%.1f}', p1,p(2));
%     end
else
    pstr = [];
end

% if ipr > 0
%     fprintf('fitpow: y = p(1)*x^p(2): p = %.6e, %.6e\n', p);
% end;
