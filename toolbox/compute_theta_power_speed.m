function [avg_speed,theta_power,theta_freq] = compute_theta_power_speed(lfp,lfp_ts,lfp_fs,thetafreq_band,widefreq_band,speed,ts,time_window)

start = max(lfp_ts(1),ts(1)); stop = min(lfp_ts(end),ts(end));
time_vec = start:time_window*1e6:stop;
avg_speed = nan(numel(time_vec)-1,1);
theta_power = nan(numel(time_vec)-1,1);
theta_freq = nan(numel(time_vec)-1,1);

for k = 1:numel(time_vec)-1
    
    speed_ind = ts >= time_vec(k) & ts < time_vec(k+1);
    lfp_ind = lfp_ts >= time_vec(k) & lfp_ts < time_vec(k+1);
    
    
    
    if sum(speed_ind) > 50 && sum(lfp_ind) > lfp_fs*(time_window-1)
    
    avg_speed(k) = nanmean(speed(ts >= time_vec(k) & ts < time_vec(k+1)));
    
    lfp_window = lfp(lfp_ts >= time_vec(k) & lfp_ts < time_vec(k+1));
    
    %{
    L = 2^nextpow2(length(lfp_window));
    lfp_fft = fft(lfp_window,L);
    P2 = abs(lfp_fft/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = lfp_fs*(0:(L/2))/L;
    %}
    
    window = round(numel(lfp_window));
    [P1,f] = pwelch(lfp_window,window,round(window/2),2^nextpow2(window),lfp_fs);
    
    theta_band = f >= thetafreq_band(1) & f <= thetafreq_band(2); theta_f = f(theta_band);
    large_band = f >= widefreq_band(1) & f <= widefreq_band(2);
    theta_power(k) = sum(P1(theta_band))/sum(P1(large_band));
    [~,max_theta_freq_ind] = max(P1(theta_band));
    theta_freq(k) = theta_f(max_theta_freq_ind);
    end
end

nan_ind = [find(isnan(avg_speed)); find(isnan(theta_power)); find(isnan(theta_freq))];
avg_speed(nan_ind) = [];
theta_power(nan_ind) = [];
theta_freq(nan_ind) = [];

return