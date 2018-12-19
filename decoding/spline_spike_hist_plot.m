function [y,x] = spline_spike_hist_plot(param,cpts,s)


%% This is adapted from Uri Eden's spline lab.

% spiketrain = the # of spikes over time
% cpts = control points
% s = tension parameters

param = reshape(param,numel(param),1);

% S matrix
S = [-s 2-s s-2 s;2*s s-3 3-2*s -s;-s 0 s 0;0 1 0 0];

% Q - max time in the past
Q = cpts(end); % this is in units of time bins

% Construct spline regressors
bin1 = cpts(2) - cpts(1); bin2 = cpts(end)-cpts(end-1);
cpts_all = [cpts(1)-bin1 cpts cpts(end)+bin2];

num_c_pts = length(cpts_all);  %number of control points in total

p_mat = nan(Q,num_c_pts);
for tau = 1:Q
    nearest_c_pt_index = max(find(cpts_all<tau));
    nearest_c_pt_time = cpts_all(nearest_c_pt_index);
    next_c_pt_time = cpts_all(nearest_c_pt_index+1);
    u = (tau-nearest_c_pt_time)/(next_c_pt_time-nearest_c_pt_time);
    p =[u^3 u^2 u 1]*S;
    p_mat(tau,:) = [zeros(1,nearest_c_pt_index-2) p zeros(1,num_c_pts-4-(nearest_c_pt_index-2))];
end

y = p_mat*param;
x = 1:Q;

return




