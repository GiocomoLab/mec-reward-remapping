function [cross_corr,numpoints] = correlation_kh(map1,map2)

% assume that map1 and map2 are the same size

[row,col] = size(map1);
cross_corr = nan(size(2*row-1,2*col-1));
numpoints = nan(size(cross_corr));
for i = 1:2*row-1
    if i <= row
        row_ind1 = row-i+1:row;
        row_ind2 = 1:i;
    else
        row_ind1 = 1:2*row-i;
        row_ind2 = i-row+ 1:row;
    end
    for j = 1:2*col-1
        if j <= col
            col_ind1 = col-j+1:col;
            col_ind2 = 1:j;
        else
            col_ind1 = 1:2*col-j;
            col_ind2 = j-col+ 1:col;
        end
        map1_ij = map1(row_ind1,col_ind1);
        map2_ij = map2(row_ind2,col_ind2);
        cross_corr(i,j) = nancorr(map1_ij(:),map2_ij(:));
        numpoints(i,j) = numel(map1_ij);
    end
end

return