function [X] = spline_1d_circ(y,cpts,s)


%% This is adapted from Uri Eden's spline lab.

% y = the values of the variable over time
% cpts = control points
% s = tension parameter

% This has been adapted to include a circular variable that goes from 0 to
% 2*pi

% S matrix
S = [-s 2-s s-2 s;2*s s-3 3-2*s -s;-s 0 s 0;0 1 0 0];

% Construct spline regressors
X = zeros(length(y),length(cpts));
%for each timepoint, calculate the corresponding row of the glm input matrix
for i=1:length(y)
    
    % find the nearest, and next, control point
    if y(i) > 0
        nearest_c_pt_index = max(find(cpts<y(i)));
        nearest_c_pt_time = cpts(nearest_c_pt_index);
    else
        nearest_c_pt_index = numel(cpts);
        nearest_c_pt_time = cpts(nearest_c_pt_index);
    end
    
    if y(i) > cpts(end) || y(i) == 0
        next_c_pt_time = cpts(1);
    else
        next_c_pt_time = cpts(nearest_c_pt_index+1);
    end
    
    
    % compute the alpha (u here)
    u = mod(y(i)-nearest_c_pt_time,pi)/mod(next_c_pt_time-nearest_c_pt_time,pi);
    p=[u^3 u^2 u 1]*S;
    
    % fill in the X matrix, with the right zeros on either side
    temp_x = [p zeros(1,numel(cpts)-4)]; % [p1 p2 p3 p4 0 ... 0] <-- shift this appropriately
    X(i,:) = circshift(temp_x,[0,nearest_c_pt_index-2]);
    
end