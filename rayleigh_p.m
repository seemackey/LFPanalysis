function [Rval, pval] = rayleigh_p(n, pval)

% n = number of epochs
% pval = p values for which R values should be calculated


Rbar=0:0.0001:1;
p = [];
for i=1:size(Rbar,2)
  Z = n*Rbar(i)^2;
  p(i) = exp(-Z) * (1 + (2*Z - Z^2) / (4*n) - (24*Z - 132*Z^2 + 76*Z^3 - 9*Z^4) / (288*n^2));
end

for i=1:size(pval,2)
    z=find(abs(pval(i)-p)==min(abs(pval(i)-p)));
    Rval(i)=Rbar(z);
end