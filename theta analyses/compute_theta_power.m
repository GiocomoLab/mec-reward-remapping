function [theta_power] = compute_theta_power(lfp,lfp_fs,thetafreq_band,widefreq_band)

window = round(10*lfp_fs);
[P1,f] = pwelch(lfp,window,round(window/2),2^nextpow2(window),lfp_fs);

theta_band = f >= thetafreq_band(1) & f <= thetafreq_band(2);
large_band = f >= widefreq_band(1) & f <= widefreq_band(2);
theta_power = sum(P1(theta_band))/sum(P1(large_band));

return