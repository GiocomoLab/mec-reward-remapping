function [mvl,direction] = compute_mvl_angle(tuning_curve)

numHDBin = numel(tuning_curve);

% compute the head direction score
exp_term  = exp(-1i*2*pi*(1:numHDBin)/numHDBin);
vector = pi/(numHDBin*sin(pi/numHDBin))*(exp_term*tuning_curve)/sum(tuning_curve);
rvl = sqrt(real(vector)^2 + imag(vector)^2);

mvl = round(rvl*100)/100;
direction = angle(vector);

return