function [mvl_theta,pref_angle,theta_index,theta_skipping] = compare_theta_single_cell(free_cell_path,task_cell_path)

% things this piece of code will answer:
% theta locking with lfp (mvl and pref angle)
% theta modulation (theta index)
% theta skipping (comparing peaks)

%% load the free and task data
[x_free, y_free, hd_free, ts_free, speed_free, spikes_free, fr_free, cellTS_free] = load_data_and_preprocess(free_cell_path);
[x_task, y_task, hd_task, ts_task, speed_task, spikes_task, fr_task, cellTS_task] = load_data_and_preprocess(task_cell_path);
%% compute theta locking

% find the right csc and load it
free_tetrode = strsplit(free_cell_path,'TT');
free_tetrode = strsplit(free_tetrode{end},'c'); free_tetrode = free_tetrode{1};
[lfp_free,lfp_ts_free] = load_csc_data_and_preprocess(free_cell_path,free_tetrode);
lfp_fs_free = 1/(median(diff(lfp_ts_free)))/1e6;

task_tetrode = strsplit(task_cell_path,'TT');
task_tetrode = strsplit(task_tetrode{end},'c'); task_tetrode = task_tetrode{1};
[lfp_task,lfp_ts_task] = load_csc_data_and_preprocess(task_cell_path,task_tetrode);
lfp_fs_task = 1/(median(diff(lfp_ts_task)))/1e6;

% filter for theta
theta_window = [6 10];
[lfp_theta_free, theta_phase_free, lfp_fs_free] = theta_filter(lfp_free,lfp_ts_free,theta_window);
[lfp_theta_task, theta_phase_task,lfp_fs_task] = theta_filter(lfp_task,lfp_ts_task,theta_window);

% for every spike, compute the theta phase
[spiketrain_free] = histcounts(cellTS_free,lfp_ts_free);
spike_phase_free = theta_phase_free(spiketrain_free>0);

[spiketrain_task] = histcounts(cellTS_task,lfp_ts_task);
spike_phase_task = theta_phase_task(spiketrain_task>0);

% compute mvl for theta locking
[rvl_free,pref_angle_free] = compute_mvl(spike_phase_free);
[rvl_task,pref_angle_task] = compute_mvl(spike_phase_task);

mvl_theta = [rvl_free rvl_task];
pref_angle_free(pref_angle_free < 0) = pref_angle_free(pref_angle_free < 0) + 2*pi;
pref_angle_task(pref_angle_task < 0) = pref_angle_task(pref_angle_task < 0) + 2*pi;
pref_angle = [pref_angle_free pref_angle_task];

%% compute theta modulation and theta skipping

% compute isi autocorr

[N_free,m_free] = compute_isi_lag(cellTS_free);
[N_task,m_task] = compute_isi_lag(cellTS_task);

% compute theta index
thetafreq_band = [6 10];
[theta_index_free] = compute_theta_index(N_free,m_free,thetafreq_band);
[theta_index_task] = compute_theta_index(N_task,m_task,thetafreq_band);

theta_index = [theta_index_free theta_index_task];

% jitter spikes and then recompute index
%{
numshuff = 100;
theta_index_free_shuff = nan(numshuff,1);
theta_index_task_shuff = nan(numshuff,1);
for k = 1:numshuff
    cellTS_free_k = mod(cellTS_free + (-10 + 20*rand(size(cellTS_free)))*1e6,max(ts_free));
    cellTS_task_k = mod(cellTS_task + (-10 + 20*rand(size(cellTS_task)))*1e6,max(ts_task));
    
    [N_free_k,m_free_k] = compute_isi_lag(cellTS_free_k);
    [N_task_k,m_task_k] = compute_isi_lag(cellTS_task_k);
    
    theta_index_free_shuff(k) = compute_theta_index(N_free_k,m_free_k,thetafreq_band);
    theta_index_task_shuff(k) = compute_theta_index(N_task_k,m_task_k,thetafreq_band);
    k
end
%}

% look at theta skipping
N_free_smooth = conv(N_free,gausswin(5)/sum(gausswin(5)),'same');
N_task_smooth = conv(N_task,gausswin(5)/sum(gausswin(5)),'same');

% find first peak
first_theta_free = find(abs(m_free - 130) < 30);
[~,max_theta_freq_ind] = max(N_free_smooth(first_theta_free));
first_theta_peak_free = N_free_smooth(first_theta_free(max_theta_freq_ind));

first_theta_task = find(abs(m_task - 130) < 30);
[~,max_theta_freq_ind] = max(N_task_smooth(first_theta_task));
first_theta_peak_task = N_task_smooth(first_theta_task(max_theta_freq_ind));

% find second peak
second_theta_free = find(abs(m_free - 260) < 60);
[~,max_theta_freq_ind] = max(N_free_smooth(second_theta_free));
second_theta_peak_free = N_free_smooth(second_theta_free(max_theta_freq_ind));

second_theta_task = find(abs(m_task - 260) < 60);
[~,max_theta_freq_ind] = max(N_task_smooth(second_theta_task));
second_theta_peak_task = N_task_smooth(second_theta_task(max_theta_freq_ind));

theta_skipping_free = (second_theta_peak_free - first_theta_peak_free)/max(second_theta_peak_free,first_theta_peak_free);
theta_skipping_task = (second_theta_peak_task - first_theta_peak_task)/max(second_theta_peak_task,first_theta_peak_task);

theta_skipping = [theta_skipping_free theta_skipping_task];

return