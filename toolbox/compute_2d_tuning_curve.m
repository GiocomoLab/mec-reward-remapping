function [tuning_curve,occupancy_curve] = compute_2d_tuning_curve(variable_x,variable_y,fr,numBin,minVal,maxVal,smoothing)

% this assumes that the 2d environment is a square box, and that the
% variable is recorded along the x- and y-axes

%% define the axes and initialize variables

xAxis = linspace(minVal(1),maxVal(1),numBin+1);
yAxis = linspace(minVal(2),maxVal(2),numBin+1);

% initialize 
tuning_curve = zeros(numBin,numBin);
rawtuning_curve = nan(numBin,numBin);
occupancy_curve = zeros(numBin,numBin);
%% fill out the tuning curve

% find the mean firing rate in each position bin
for i  = 1:numBin
    start_x = xAxis(i); stop_x = xAxis(i+1);
    % find the times the animal was in the bin
    if i == numBin
        x_ind = find(variable_x >= start_x & variable_x <= stop_x);
    else
        x_ind = find(variable_x >= start_x & variable_x < stop_x);
    end
    
    for j = 1:numBin
        
        start_y = yAxis(j); stop_y = yAxis(j+1);
        
        if j == numBin
            y_ind = find(variable_y >= start_y & variable_y <= stop_y);
        else
            y_ind = find(variable_y >= start_y & variable_y < stop_y);
        end
        
        ind = intersect(x_ind,y_ind);
        
        % fill in rate map
        tuning_curve(numBin+1 - j,i) = sum(fr(ind));
        rawtuning_curve(numBin+1 - j,i) = mean(fr(ind));
        occupancy_curve(numBin+1 - j,i) = numel(ind);
    end
end

%%
if smoothing  
    %new smoothing (very similar to Stensola et al. 2015 and other Moser papers)
    %they use 5 x 5 filter for smoothing 3cm bins; Leutgebs use 5 x 5 filter for
    %5 cm bins
    H = fspecial('gaussian',[5 5],1);

    %smooth spikes and occupancy maps separately, then calculate rate map
    smspikes = imfilter(tuning_curve,H);
    smoccu = imfilter(occupancy_curve, H);
    tuning_curve = smspikes./smoccu;
else 
    tuning_curve = tuning_curve./occupancy_curve;
end

