function [tuning_curve,x,sem,occupancy] = compute_1d_tuning_curve(variable,fr,numBin,minVal,maxVal)

%bin it
var_vec = linspace(minVal,maxVal,numBin+1);
tuning_curve = nan(numBin,1);
sem = nan(numBin,1);
occupancy = nan(numBin,1);

% compute mean fr for each bin
for n = 1:numBin
    if n == numBin
        index = variable >= var_vec(n) & variable <= var_vec(n+1);
    else
        index = variable >= var_vec(n) & variable < var_vec(n+1);
    end
    tuning_curve(n) = nanmean(fr(index));
    sem(n) = nanstd(fr(index))./sqrt(sum(index));
    occupancy(n) = sum(index);
    
end

x = var_vec(1:end-1); x = x + (var_vec(2)-var_vec(1))/2;

return