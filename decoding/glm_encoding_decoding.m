function [free_decode_info,task_decode_info,encoding_cells_free,encoding_cells_task] = glm_encoding_decoding(free_path,free_cells,task_path,task_cells)

%% random things to set initially

boxSize = 150;
dt = 0.02;
s = 0.5;
bin_p = 9;
posVec = linspace(0,boxSize,bin_p); posVec(1) = -0.01;
spikeHistoryVec = [0:5 7 10]; %in bins
ctl_pts_all{1} = posVec;
ctl_pts_all{2} = spikeHistoryVec;

total_decode_time = 5*60/dt; % in time bins
num_int = 0.5/dt;
num_decode_chunks = round(total_decode_time/num_int);

num_folds = 10;
chunks = 30;


%% load the behavioral data for each session
% and compute the position grid for each session

%%%%%%%%%%%%%% FREE %%%%%%%%%%%%%%%%%%%%
% load free data
[x_free, y_free, ~, ts_free, ~] = load_session_data_and_preprocess(free_path);

% fill in nans and slightly upsample data
[x_us_free,y_us_free,ts_us_free] = upsample_pos_data(x_free,y_free,ts_free,dt);
T_free = numel(ts_us_free);

% compute the position grid (free)
[posgrid_free,~] = spline_2d(x_us_free,y_us_free,posVec,s);

%%%%%%%%%%%%%% TASK %%%%%%%%%%%%%%%%%%%%
% load task data
[x_task, y_task, ~, ts_task, ~] = load_session_data_and_preprocess(task_path);

% slightly upsample data
[x_us_task,y_us_task,ts_us_task] = upsample_pos_data(x_task,y_task,ts_task,dt);
T_task = numel(ts_us_task);

% compute the position grid (free)
[posgrid_task,~] = spline_2d(x_us_task,y_us_task,posVec,s);


%% load the cell data for each session
% also get the spike history splines while we are at it

num_cells = numel(free_cells);
S_free = nan(num_cells,T_free);
S_task = nan(num_cells,T_task);
shgrid_free = cell(num_cells,1);
shgrid_task = cell(num_cells,1);
for k = 1:num_cells
    
    % free
    [cellTS_free] = load_spike_data(free_cells{k});
    S_free(k,:) = hist(cellTS_free,ts_us_free);
    [shgrid,~] = spline_spike_hist(S_free(k,:),spikeHistoryVec,s);
    shgrid_free{k} = shgrid;
    
    % task
    [cellTS_task] = load_spike_data(task_cells{k});
    S_task(k,:) = hist(cellTS_task,ts_us_task);
    [shgrid,~] = spline_spike_hist(S_task(k,:),spikeHistoryVec,s);
    shgrid_task{k} = shgrid;
    
end

num_pos_param = size(posgrid_free,2);
num_sh_param = size(shgrid,2);
pos_ind = 2:2+num_pos_param-1;
sh_ind = 2+num_pos_param:2+num_pos_param+num_sh_param-1;

%% take out decoding chunks for free and task (set aside test data for encoding)

[decoding_ind_free,decoding_start_free,decoding_stop_free] = get_decoding_ind(ts_us_free,num_decode_chunks,num_int);

[decoding_ind_task,decoding_start_task,decoding_stop_task] = get_decoding_ind(ts_us_task,num_decode_chunks,num_int);

%% match the occupancy for the training data (get test/train data for encoding model)

wo_decoding_free = setdiff(1:T_free,decoding_ind_free);
wo_decoding_task = setdiff(1:T_task,decoding_ind_task);

[ds_ind1, ds_ind2] = downsample_match_position(x_us_free(wo_decoding_free),y_us_free(wo_decoding_free),x_us_task(wo_decoding_task),y_us_task(wo_decoding_task));
free_ind = wo_decoding_free(ds_ind1);
task_ind = wo_decoding_task(ds_ind2);

%% split encoding data into testing and training datasets

[train_ind_free,test_ind_free] = get_train_test_ind(num_folds,chunks,free_ind,dt);

[train_ind_task,test_ind_task] = get_train_test_ind(num_folds,chunks,task_ind,dt);

%% fit the encoding model for each session


[param_mat_free,tuning_curves_free] = fit_encoding_model(posgrid_free,shgrid_free,pos_ind,sh_ind,S_free,train_ind_free,test_ind_free,ctl_pts_all,s,dt);
[param_mat_task,tuning_curves_task] = fit_encoding_model(posgrid_task,shgrid_task,pos_ind,sh_ind,S_task,train_ind_task,test_ind_task,ctl_pts_all,s,dt);


encoding_cells_free = find(param_mat_free(:,2) ~= 0);
encoding_cells_task = find(param_mat_task(:,2) ~= 0);

% for each encoding cell, plot the position tuning curve
%{
for k = 1:numel(encoding_cells_free)
    figure(1)
    subplot(2,4,k)
    imagesc(reshape(tuning_curves_free{encoding_cells_free(k)}{1},51,51))
    axis off
    
    figure(2)
    subplot(2,4,k)
    [tuning_curve] = compute_2d_tuning_curve(x_us_free,y_us_free,S_free(encoding_cells_free(k),:),30,[0 0],[150 150],1);
    imagesc(tuning_curve)
    axis off
end
%}

%% make sure that encoding model is correct
% X = [ones(size(S_free,2),1) posgrid_free shgrid_free{1}];

%% decode position for each session

%%%%%% free %%%%

S_free_encoding = S_free(encoding_cells_free,:);
param_free_encoding = param_mat_free(encoding_cells_free,:);
[x_real_free,x_decode_free,y_real_free,y_decode_free] = do_decoding(x_us_free,y_us_free,S_free_encoding,shgrid_free,decoding_start_free,decoding_stop_free,param_free_encoding,pos_ind,sh_ind,posVec,s);

free_decode_info = [x_real_free x_decode_free y_real_free y_decode_free];
free_error = mean(sqrt((x_real_free - x_decode_free).^2 + (y_real_free - y_decode_free).^2))

%%%%%% task %%%%
S_task_encoding = S_task(encoding_cells_task,:);
param_task_encoding = param_mat_task(encoding_cells_task,:);
[x_real_task,x_decode_task,y_real_task,y_decode_task] = do_decoding(x_us_task,y_us_task,S_task_encoding,shgrid_task,decoding_start_task,decoding_stop_task,param_task_encoding,pos_ind,sh_ind,posVec,s);

task_decode_info = [x_real_task x_decode_task y_real_task y_decode_task];
task_error = mean(sqrt((x_real_task - x_decode_task).^2 + (y_real_task - y_decode_task).^2))

return


