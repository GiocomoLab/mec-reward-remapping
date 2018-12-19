function [X,cpts_all] = spline_1d(y,cpts,s)


%% This is adapted from Uri Eden's spline lab. 

% y = the values of the variable over time
% cpts = control points
% s = tension parameters

% S matrix
S = [-s 2-s s-2 s;2*s s-3 3-2*s -s;-s 0 s 0;0 1 0 0];

% Construct spline regressors
bin1 = cpts(2) - cpts(1);
bin2 = cpts(end) - cpts(end-1);
cpts_all = [cpts(1)-bin1 cpts cpts(end)+bin2];
X = zeros(length(y),length(cpts_all));
num_c_pts = length(cpts_all);  %number of control points in total


%for each timepoint, calculate the corresponding row of the glm input matrix
for i=1:length(y)  
    
    % find the nearest, and next, control point
    nearest_c_pt_index = max(find(cpts_all < y(i)));
    nearest_c_pt_time = cpts_all(nearest_c_pt_index);
    next_c_pt_time = cpts_all(nearest_c_pt_index+1);
    
    % compute the alpha (u here)
    u = (y(i)-nearest_c_pt_time)/(next_c_pt_time-nearest_c_pt_time);
    p=[u^3 u^2 u 1]*S;
    
    % fill in the X matrix, with the right # of zeros on either side 
    X(i,:) = [zeros(1,nearest_c_pt_index-2) p zeros(1,num_c_pts-4-(nearest_c_pt_index-2))];
    
end