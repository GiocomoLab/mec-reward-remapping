function b = nanpolyfit(x,y,p)

x = reshape(x,numel(x),1);
y = reshape(y,numel(y),1);

x(isnan(y)) = []; y(isnan(y)) = [];
y(isnan(x)) = []; x(isnan(x)) = [];

b = polyfit(x,y,p);

return