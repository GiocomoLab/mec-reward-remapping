function [tuning_curves] = plot_tuning(A,variables,parameters,ctl_pts_all,s,plotfig,dt)

%% Description
% Given the variables, A (and thetagrid and shgrid), and the parameters,
% this will return the tuning curves for the cell
% if plotfig = 1, this will also plot the tuning curves

% NOTE: I just use A to compute the correct indexes
tuning_curves = {};
variables = sort(variables);
b0 = parameters(1);
param = parameters(2:end);

% position
total_ind = 0;
scale = [];

if ismember(1,variables)
    param_ind = size(A{1},2);
    total_ind = total_ind + param_ind;
    param1 = param(total_ind - param_ind + 1:total_ind);
    [pos_y,~,~] = spline_2d_plot(param1,ctl_pts_all{1},s);
    scale = [scale; mean(mean(pos_y))];
end

% spike history
if ismember(2,variables)
    param_ind = size(A{2},2);
    total_ind = total_ind + param_ind;
    param1 = param(total_ind - param_ind + 1:total_ind);
    [sh_y,sh_x] = spline_spike_hist_plot(param1,ctl_pts_all{2},s);
    scale = [scale; mean(sh_y)];
end

tuning_curves = {};


if ismember(1,variables)
    scale_factor = scale(ismember(variables,setdiff(variables,1)));
    temp = exp(pos_y)*exp(b0)*prod(exp(scale_factor));
    tuning_curves{end+1} = temp(:);
    if plotfig
        close all
        figure(1)
        subplot(1,2,1)
        imagesc(exp(pos_y)*exp(b0)*prod(exp(scale_factor))/dt);
        colorbar; axis off
    end
end

if ismember(2,variables)
    tuning_curves{end+1} = exp(sh_y);
    if plotfig
        figure(1)
        subplot(1,2,2)
        plot(sh_x*dt*1e3,exp(sh_y),'k','linewidth',2);
        hold on
        plot(sh_x*dt*1e3,ones(size(sh_x)),'--r')
        hold off
        box off
        xlabel('ms in past')
        ylabel('gain')
        axis tight
    end
end

return