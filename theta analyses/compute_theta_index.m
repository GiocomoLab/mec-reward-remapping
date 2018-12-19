function [theta_index] = compute_theta_index(N,m,thetafreq_band)

% compute sampling rate
fs = 1/((m(2) - m(1))/1e3);

% compute fft
window = round(numel(N));
[P1,f] = pwelch(N,window,round(window/2),2^nextpow2(window),fs);

% smooth the power spectrum
theta_f = find(f >= thetafreq_band(1) & f <= thetafreq_band(2));
[~,max_theta_freq_ind] = max(P1(theta_f));
max_theta = f(theta_f(max_theta_freq_ind));

theta_power = mean(P1(f >= max_theta - 1 & f <= max_theta + 1));
other_power = mean(P1(f >= 1 & f <= 50));
theta_index = theta_power/other_power;

return