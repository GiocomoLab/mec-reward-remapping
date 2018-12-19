function [power_position,freq_position] = compute_theta_power_position(lfp,lfp_ts,thetafreq_band,widefreq_band,x,y,ts,time_window)

lfp_fs = 1/((lfp_ts(2) - lfp_ts(1))/1e6);
dt = (ts(2) - ts(1))/1e6;

start = max(lfp_ts(1),ts(1)); stop = min(lfp_ts(end),ts(end));
time_vec = start:time_window*1e6:stop;
avg_x = nan(numel(time_vec)-1,1);
avg_y = nan(numel(time_vec)-1,1);
theta_power = nan(numel(time_vec)-1,1);
theta_freq = nan(numel(time_vec)-1,1);
for k = 1:numel(time_vec)-1
    
    time_ind = ts >= time_vec(k) & ts < time_vec(k+1);
    lfp_ind = lfp_ts >= time_vec(k) & lfp_ts < time_vec(k+1);
    
    if sum(time_ind) > time_window/dt/2 && sum(lfp_ind) > lfp_fs*(time_window)/2
        
        avg_x(k) = nanmean(x(time_ind));
        avg_y(k) = nanmean(y(time_ind));
        
        lfp_window = lfp(lfp_ts >= time_vec(k) & lfp_ts < time_vec(k+1));
        
        window = round(numel(lfp_window));
        [P1,f] = pwelch(lfp_window,window,round(window/2),2^nextpow2(window),lfp_fs);
        
        theta_band = f >= thetafreq_band(1) & f <= thetafreq_band(2); theta_f = f(theta_band);
        large_band = f >= widefreq_band(1) & f <= widefreq_band(2);
        theta_power(k) = sum(P1(theta_band))/sum(P1(large_band));
        [~,max_theta_freq_ind] = max(P1(theta_band));
        theta_freq(k) = theta_f(max_theta_freq_ind);
    end
end

nan_ind = [find(isnan(avg_x)); find(isnan(avg_y)); find(isnan(theta_power)); find(isnan(theta_freq))];
avg_x(nan_ind) = [];
avg_y(nan_ind) = [];
theta_power(nan_ind) = [];
theta_freq(nan_ind) = [];

% compute the position-power and position-frequency tuning curve
numbins = 25;
posvec = linspace(0,150,numbins+1);
power_position = nan(numbins,numbins);
freq_position = nan(numbins,numbins);
for i = 1:numbins
    if i == numbins
        x_ind = find(avg_x >= posvec(i) & avg_x <= posvec(i+1));
    else
        x_ind = find(avg_x >= posvec(i) & avg_x < posvec(i+1));
    end
    for j = 1:numbins
        if j == numbins
            y_ind = find(avg_y >= posvec(j) & avg_y <= posvec(j+1));
        else
            y_ind = find(avg_y >= posvec(j) & avg_y < posvec(j+1));
        end
        
        index = intersect(x_ind,y_ind);
        power_position(numbins-j+1,i) = nanmean(theta_power(index));
        freq_position(numbins-j+1,i) = nanmean(theta_freq(index));
    end
end

[power_position] = fill_in_nan_2d(power_position);
[freq_position] = fill_in_nan_2d(freq_position);
H = fspecial('gaussian',[3 3],0.5); % using default values - size=[3 3] and sigma=0.5
power_position = imfilter(power_position,H);
freq_position = imfilter(freq_position,H);

return