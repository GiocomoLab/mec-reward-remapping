function [temp1,temp2] = remove_nan(temp1,temp2)

temp1(isnan(temp2)) = []; temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = []; temp1(isnan(temp1)) = [];


return
