function [tuning_curve] = fill_in_nan_2d(tuning_curve)
% this will first interpolate to fill in nan's, and then will use another
% backup method to fill in anything remaining

%% first, interpolate to fill in NaN's

[nanrow,nancol] = find(isnan(tuning_curve));
nanind = sub2ind(size(tuning_curve),nanrow,nancol);

interp_tuning_curve = interp2(tuning_curve,nanrow,nancol);
tuning_curve(nanind) = interp_tuning_curve;


%% second, if the above method leaves NaNs, fill them in with mean of neighbors
[~,numbins] = size(tuning_curve);

% fill in the NaNs with neigboring values
nan_ind = find(isnan(tuning_curve));
[j,i] = ind2sub(size(tuning_curve),nan_ind);
nan_num= numel(nan_ind);

% fill in the NaNs with neigboring values
for n = 1:nan_num
    ind_i = i(n); ind_j = j(n);
    
    right = tuning_curve(ind_j,min(ind_i+1,numbins));
    left = tuning_curve(ind_j,max(ind_i-1,1));
    down = tuning_curve(min(ind_j+1,numbins),ind_i);
    up = tuning_curve(max(ind_j-1,1),ind_i);
    
    ru = tuning_curve(max(ind_j-1,1),min(ind_i+1,numbins));
    lu = tuning_curve(max(ind_j-1,1),max(ind_i-1,1));
    ld = tuning_curve(min(ind_j+1,numbins),max(ind_i-1,1));
    rd = tuning_curve(max(ind_j-1,1),min(ind_i+1,numbins));
    
    tuning_curve(ind_j,ind_i) = nanmean([left right up down lu ru rd ld]);
    
end

% fill in the NaNs with neigboring values
nan_ind = flip(find(isnan(tuning_curve)));
[j,i] = ind2sub(size(tuning_curve),nan_ind);
nan_num= numel(nan_ind);
for n = 1:nan_num
    ind_i = i(n); ind_j = j(n);
    
    right = tuning_curve(ind_j,min(ind_i+1,numbins));
    left = tuning_curve(ind_j,max(ind_i-1,1));
    down = tuning_curve(min(ind_j+1,numbins),ind_i);
    up = tuning_curve(max(ind_j-1,1),ind_i);
    
    ru = tuning_curve(max(ind_j-1,1),min(ind_i+1,numbins));
    lu = tuning_curve(max(ind_j-1,1),max(ind_i-1,1));
    ld = tuning_curve(min(ind_j+1,numbins),max(ind_i-1,1));
    rd = tuning_curve(max(ind_j-1,1),min(ind_i+1,numbins));
    
    tuning_curve(ind_j,ind_i) = nanmean([left right up down lu ru rd ld]);
    
end