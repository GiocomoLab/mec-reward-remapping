var_name = {'position','hd','speed'};
numVar = numel(var_name);

variables = sort(best_models);
b0 = final_param(1);
param = final_param(2:end);

total_ind = 0;
% position
numvar = 1;
if ismember(numvar,variables)
    param_ind = size(A{numvar},2);
    total_ind = total_ind + param_ind;
    param1 = param(total_ind - param_ind + 1:total_ind);
    [pos_y,x_spline,y_spline] = spline_2d_plot(param1,ctl_pts_all{numvar},s);
    [~,pos_occupancy] = compute_2d_tuning_curve(x_spline,y_spline,ones(length(x_spline),1),length(x_spline),[x_spline(1) y_spline(1)],[x_spline(end) y_spline(end)]);
    pos_occupancy = pos_occupancy(:); pos_occupancy = pos_occupancy./sum(pos_occupancy);
    scale(numvar) = exp(sum(pos_occupancy.*pos_y(:)));
end


% head direction
numvar = 2;
if ismember(numvar,variables)
    param_ind = size(A{numvar},2);
    total_ind = total_ind + param_ind;
    param1 = param(total_ind - param_ind + 1:total_ind);
    [hd_y,hd_x] = spline_1d_circ_plot(param1,ctl_pts_all{numvar},s);
    [~,hd_occupancy] = compute_1d_tuning_curve(hd,ones(length(hd),1),numel(hd_x),hd_x(1),hd_x(end));
    hd_occupancy = hd_occupancy./sum(hd_occupancy);
    scale(numvar) = exp(sum(hd_occupancy.*hd_y));
end


% speed
numvar = 3;
if ismember(numvar,variables)
    param_ind = size(A{numvar},2);
    total_ind = total_ind + param_ind;
    param1 = param(total_ind - param_ind + 1:total_ind);
    [speed_y,speed_x] = spline_1d_plot(param1,ctl_pts_all{numvar},s);
    [~,speed_occupancy] = compute_1d_tuning_curve(speed,ones(length(speed),1),numel(speed_x),speed_x(1),speed_x(end));
    speed_occupancy = speed_occupancy./sum(speed_occupancy);
    scale(numvar) = exp(sum(speed_occupancy.*speed_y));
end

tuning_curves = {};
var_k = 1;
if ismember(var_k,variables)
    scale_factor_ind = setdiff(variables,var_k); scale_factor = scale(scale_factor_ind);
    tuning_curves{var_k} = exp(pos_y(:))*exp(b0)*prod(scale_factor);
    
    if plotfig
        figure(1)
        subplot(3,3,4)
        imagesc([0 150],[0 150],exp(pos_y)*exp(b0)*prod(scale_factor))
        axis off
        axis tight
        temp = exp(pos_y)*exp(b0)*prod(scale_factor);
        formatMin = sprintf('%0.2f',min(temp(:)));
        formatMax = sprintf('%0.2f',max(temp(:)));
        t = ['min = ' formatMin ' max = ' formatMax];
        title(t)
    end
end

var_k = 2;
if ismember(var_k,variables)
    
    scale_factor_ind = setdiff(variables,var_k); scale_factor = scale(scale_factor_ind);
    tuning_curves{var_k} = exp(hd_y')*exp(b0)*prod(scale_factor);
    
    if plotfig
        figure(1)
        subplot(3,3,5)
        plot(hd_x,exp(hd_y)*exp(b0)*prod(scale_factor),'k','linewidth',2);
        box off
        xlabel('HD angle')
        ylabel('spikes/s')
        axis tight
        vector_y = tuning_curves{var_k};
        title('LNP hd tuning')
    end
end

var_k = 3;
if ismember(var_k,variables)
    
    scale_factor_ind = setdiff(variables,var_k); scale_factor = scale(scale_factor_ind);
    tuning_curves{var_k} = exp(speed_y')*exp(b0)*prod(scale_factor);
    
    if plotfig
        subplot(3,3,6)
        positive_val = speed_x > 0;
        plot(speed_x(positive_val),exp(speed_y(positive_val))*prod(scale_factor)*exp(b0),'k','linewidth',2);
        box off
        xlabel('speed')
        ylabel('spikes/s')
        axis tight
        vector_x = speed_x(positive_val); vector_y = tuning_curves{var_k}(positive_val);
        title('LNP speed tuning')
    end
end
