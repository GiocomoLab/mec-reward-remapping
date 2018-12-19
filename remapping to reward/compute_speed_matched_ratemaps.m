function [ratemap1,ratemap2,occupancy1] = compute_speed_matched_ratemaps(filename1,filename2,numiter)

%% load the data

%load video data for first session; must have run 'analyze_session_data' first
[x1, y1, ~, ~, speed1, ~, fr1] = load_data_and_preprocess(filename1);

%load video data for second session
[x2, y2, ~, ~, speed2, ~, fr2] = load_data_and_preprocess(filename2);

%% Downsample both sessions to match on position and speed samplings
numPosBin = 30; boxSize = 150;
ratemap1_all = nan(numiter,numPosBin,numPosBin);
ratemap2_all = nan(numiter,numPosBin,numPosBin);

for i = 1:numiter
    %get indices for downsampling behavioral data
    [ind1, ind2] = downsample_match_speed(x1,y1,speed1,x2,y2,speed2);

    %downsampled rate map for first session
    [ratemap1_all(i,:,:),occupancy1] = compute_2d_tuning_curve(x1(ind1),y1(ind1),fr1(ind1),numPosBin,[0 0],[boxSize boxSize],1);
    [ratemap2_all(i,:,:),occupancy2] = compute_2d_tuning_curve(x2(ind2),y2(ind2),fr2(ind2),numPosBin,[0 0],[boxSize boxSize],1);

end

% compute the average downsampled ratemap
ratemap1 = squeeze(nanmean(ratemap1_all,1));
ratemap2 = squeeze(nanmean(ratemap2_all,1));

return

