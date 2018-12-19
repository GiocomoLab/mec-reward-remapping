function [theta_power,thetapow_speed_info_free,thetapow_speed_info_task,thetafreq_speed_info_free,thetafreq_speed_info_task] ...
    = compare_theta_lfp(free_path,task_path)

% things this piece of code will answer:
% theta power in each session
% theta power,frequence vs running speed


%% load the free and task x, y, hd, speed data
[x_free, y_free, hd_free, ts_free, speed_free] = load_session_data_and_preprocess(free_path);
[x_task, y_task, hd_task, ts_task, speed_task] = load_session_data_and_preprocess(task_path);

%% load lfp, compute theta power

% find the tetrode with the highest theta power

% FREE
% get the list of all csc's for that session
split = strsplit(free_path, '_');
folderpath_free = char(strcat('Z:\Users\WButler\D254 Neuralynx\',split(1),'\',split(2),'_',split(3),'\'));
mat_files_free = dir([folderpath_free,'/*.ncs']);
csc_files_free = {mat_files_free.name}';


split = strsplit(task_path, '_');
folderpath_task = char(strcat('Z:\Users\WButler\D254 Neuralynx\',split(1),'\',split(2),'_',split(3),'\'));
mat_files_task = dir([folderpath_task,'/*.ncs']);
csc_files_task = {mat_files_task.name}';


thetafreq_band = [6 10];
widefreq_band = [1 50];
theta_power_task = nan(numel(csc_files_free),1);
theta_power_free = nan(numel(csc_files_task),1);
for k = 1:numel(csc_files_free)
    
    % free
    csc_num = strsplit(csc_files_free{k},'.ncs');
    csc_num = strsplit(csc_num{1},'CSC'); csc_num = csc_num{2};
    
    
    [lfp_free_k,lfp_ts_free] = load_csc_data_and_preprocess(free_path,csc_num);
    lfp_fs_free = 1/((lfp_ts_free(2) - lfp_ts_free(1))/1e6);
    [theta_power_free(k)] = compute_theta_power(lfp_free_k,lfp_fs_free,thetafreq_band,widefreq_band);
    
     % task
    [lfp_task_k,lfp_ts_task] = load_csc_data_and_preprocess(task_path,csc_num);
    lfp_fs_task = 1/((lfp_ts_task(2) - lfp_ts_task(1))/1e6);
    [theta_power_task(k)] = compute_theta_power(lfp_task_k,lfp_fs_task,thetafreq_band,widefreq_band);
    
end

% compare tetrode with highest theta
[highest_theta_free,highest_theta_free_ind] = max(theta_power_free);
[highest_theta_task,highest_theta_task_ind] = max(theta_power_task);

theta_power = [highest_theta_free highest_theta_task];

csc_num = strsplit(csc_files_free{highest_theta_free_ind},'.ncs');
csc_num = strsplit(csc_num{1},'CSC'); csc_num = csc_num{2};
[lfp_free,lfp_ts_free] = load_csc_data_and_preprocess(free_path,csc_num);
lfp_fs_free = 1/((lfp_ts_free(2) - lfp_ts_free(1))/1e6);

csc_num = strsplit(csc_files_task{highest_theta_task_ind},'.ncs');
csc_num = strsplit(csc_num{1},'CSC'); csc_num = csc_num{2};
[lfp_task,lfp_ts_task] = load_csc_data_and_preprocess(task_path,csc_num);
lfp_fs_task = 1/((lfp_ts_task(2) - lfp_ts_task(1))/1e6);

%% for csc with highest theta, compute the theta power/freq with running speed

% compute theta power - speed relationship
time_window = 2; % seconds
[avg_speed_free,theta_power_free,theta_freq_free] = compute_theta_power_speed(lfp_free,lfp_ts_free,thetafreq_band,widefreq_band,speed_free,ts_free,time_window);
[avg_speed_task,theta_power_task,theta_freq_task] = compute_theta_power_speed(lfp_task,lfp_ts_task,thetafreq_band,widefreq_band,speed_task,ts_task,time_window);

% compute correlation, slope, intercept for power
[corr_free,~] = nancorr(avg_speed_free,theta_power_free);
[corr_task,~] = nancorr(avg_speed_task,theta_power_task);

slope_free = polyfit(avg_speed_free,theta_power_free,1); intercept_free = slope_free(2); slope_free = slope_free(1);
slope_task = polyfit(avg_speed_task,theta_power_task,1); intercept_task = slope_task(2); slope_task = slope_task(1);

thetapow_speed_info_free = [corr_free slope_free intercept_free];
thetapow_speed_info_task = [corr_task slope_task intercept_task];

% compute correlation, slope, intercept for freq
[freq_corr_free,~] = nancorr(avg_speed_free,theta_freq_free);
[freq_corr_task,~] = nancorr(avg_speed_task,theta_freq_task);

freq_slope_free = polyfit(avg_speed_free,theta_freq_free,1); freq_intercept_free = freq_slope_free(2); freq_slope_free = freq_slope_free(1);
freq_slope_task = polyfit(avg_speed_task,theta_freq_task,1); freq_intercept_task = freq_slope_task(2); freq_slope_task = freq_slope_task(1);

thetafreq_speed_info_free = [freq_corr_free freq_slope_free freq_intercept_free];
thetafreq_speed_info_task = [freq_corr_task freq_slope_task freq_intercept_task];

return