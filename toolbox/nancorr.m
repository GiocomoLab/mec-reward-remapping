function [rho,pval] = nancorr(temp1,temp2)

temp1(isnan(temp2)) = []; temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = []; temp1(isnan(temp1)) = [];
if isempty(temp1)
    rho = NaN; pval = NaN;
else
    [rho,pval] = corr(temp1,temp2);
end

return