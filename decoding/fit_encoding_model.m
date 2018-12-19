function [param_mat,tuning_curves] = fit_encoding_model(posgrid,shgrid,pos_ind,sh_ind,S,train_ind,test_ind,ctl_pts_all,s,dt)

num_cells = size(S,1);
best_models = cell(num_cells,1);
param_mat = zeros(num_cells,1+size(posgrid,2)+size(shgrid,2));
A{1} = posgrid;
final_pval = nan(num_cells,1);
tuning_curves = cell(num_cells,1);

for n = 1:num_cells
    
    fprintf('Fitting cell  %d \n', n);
    
    A{2} = shgrid{n};
    
    spiketrain = S(n,:);
    
    [~, ~, best_models{n}, ~, parameters, ~, final_pval(n)] = forward_search_kfold(A,spiketrain,train_ind,test_ind);
    
    if numel(best_models{n})>0
        
        param_mean = mean(parameters);
        param_mat(n,1) = param_mean(1);
        
        [tuning_curves{n}] = plot_tuning(A,best_models{n},param_mean,ctl_pts_all,s,0,dt);
        
        if ismember(1,best_models{n})
            
            param_mat(n,pos_ind) = param_mean(pos_ind);
            
            if ismember(2,best_models{n})
                param_mat(n,sh_ind) = param_mean(sh_ind);
            end
            
        elseif ismember(2,best_models{n})
            param_mat(n,sh_ind) = param_mean(2:numel(sh_ind)+1);
        end
        
    end
    
    
end

return