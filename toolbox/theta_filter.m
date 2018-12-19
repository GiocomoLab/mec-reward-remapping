function [lfp_theta, theta_phase, lfp_fs] = theta_filter(lfp,lfp_ts,theta_window)

lfp_fs = 1/((lfp_ts(2) - lfp_ts(1))/1e6);
[b,a] = butter(3,theta_window./(lfp_fs/2)); %bandpass between 6 and 10 Hz
lfp_theta = filtfilt(b,a,lfp);
hilb_lfp = hilbert(lfp_theta); % compute hilbert transform
theta_phase = atan2(imag(hilb_lfp),real(hilb_lfp)); %inverse tangent (-pi to pi)
ind = theta_phase <0; theta_phase(ind) = theta_phase(ind)+2*pi; % from 0 to 2*pi

return