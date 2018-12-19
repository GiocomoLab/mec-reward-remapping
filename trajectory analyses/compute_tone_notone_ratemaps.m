function [tone_ratemap,notone_ratemap,tone_speedmap,notone_speedmap, zone_center, mean_fr, mean_speed] = compute_tone_notone_ratemaps(cell_file)

tone_ratemap = {};
notone_ratemap = {};
tone_speedmap = {};
notone_speedmap = {};
zone_center = nan(1,2);
mean_fr = nan(1,2);
mean_speed = nan(1,2);

%load video data for first session; must have run 'analyze_session_data' first
[x, y, ~, ts, speed, spiketrain, fr, cellTS] = load_data_and_preprocess(cell_file);
T = numel(ts); dt = median(diff(ts))/1e6;
boxSize = 150;

%% compute the zone, tone start, and zone entries

task_file = strsplit(cell_file,'_');
task_file = strjoin(task_file(1:end-1),'_');
[zone, ~, tone_start_ts, zone_entry_ts,zone_exit_ts,other_zone_entry_ts,other_zone_exit_ts] = compute_zone(task_file);
zone_center = [(zone(3)-zone(1))/2+zone(1)  (zone(4)-zone(2))/2+zone(2)];

if isempty(tone_start_ts)
    return
end

minx = zone(1); maxx = zone(3);
miny = zone(2); maxy = zone(4);

points = 100;
top_x = linspace(minx,maxx,points); top_y = maxy*ones(1,points);
bottom_x = linspace(minx,maxx,points); bottom_y = miny*ones(1,points);
left_x = minx*ones(1,points); left_y = linspace(miny,maxy,points);
right_x = maxx*ones(1,points); right_y = linspace(miny,maxy,points);
x_boundary = [top_x bottom_x left_x right_x]'; y_boundary = [top_y bottom_y left_y right_y]';

% find all the times inside the zone
x_inside = find(x > minx & x < maxx); y_inside = find(y > miny & y < maxy);
zone_inside = intersect(x_inside,y_inside);
zone_binary = zeros(numel(ts),1);
zone_binary(zone_inside) = 1;
zone_change = diff(zone_binary);
zone_enter = find(zone_change == 1); all_zone_entries_ts = ts(zone_enter);
zone_exit = find(zone_change == -1);all_zone_exits_ts = ts(zone_exit);

%% find all the tone-on times
tone_ind = [];
tone_start_ind = nan(numel(tone_start_ts),1);
zone_entry_ind = nan(numel(tone_start_ts),1);
for k = 1:numel(tone_start_ts)
    [~,tone_start_ind(k)] = min(abs(tone_start_ts(k)-ts));
    [~,zone_entry_ind(k)] = min(abs(zone_entry_ts(k)-ts));
    tone_ind = [tone_ind; (tone_start_ind(k):zone_entry_ind(k))'];
end
tone_ind(tone_ind > T) = [];

% find all non-tone times
no_tone_ind = setdiff(1:T,tone_ind); T_notone = numel(no_tone_ind);

%% make rate map for tone and no-tone times - NO DOWNSAMPLING

% for all the tone times, compute the rate map, and the speeds in each bin.
% for all non tone times, find similar speed points (in each bin)
num_p_bins = 20; num_s_bins = 30; maxSpeed = 150; numiter = 50;

posVec = linspace(0,boxSize,num_p_bins+1);
tone_ratemap = nan(num_p_bins,num_p_bins,numiter);
notone_ratemap = nan(num_p_bins,num_p_bins,numiter);
tone_speedmap = nan(num_p_bins,num_p_bins,numiter);
notone_speedmap = nan(num_p_bins,num_p_bins,numiter);
mean_fr_all = nan(numiter,2);
mean_speed_all = nan(numiter,2);

for k = 1:numiter
    
    tone_ind1 = [];
    no_tone_ind1 = [];

    for i = 1:num_p_bins
        
        indexi = find(x >= posVec(i)  & x < posVec(i+1));
        
        if i == num_p_bins
            indexi = find(x >= posVec(i)  & x <= posVec(i+1));
        end
        
        for j = 1:num_p_bins
            
            indexj = find(y >= posVec(j)  & y < posVec(j+1));
            
            if j == num_p_bins
                indexj = find(y >= posVec(j)  & y <= posVec(j+1));
            end
            
            index = intersect(indexi,indexj);
            
            index_tone = intersect(index,tone_ind);
            index_notone1 = intersect(index,no_tone_ind);
            
            if numel(index) > 0
                
                speed_notone1_ij = speed(index_notone1);
                speed_tone_ij = speed(index_tone);
                [spd_tone_ind1,spd_notone_ind1] = match_1d_dist(speed_tone_ij,speed_notone1_ij,linspace(0,maxSpeed,num_s_bins));
                
                tone_ind1_ij = index_tone(spd_tone_ind1);
                no_tone_ind1_ij = index_notone1(spd_notone_ind1);
                
                tone_ind1 = [tone_ind1; tone_ind1_ij];
                no_tone_ind1 = [no_tone_ind1; no_tone_ind1_ij'];

                
                tone_ratemap(num_p_bins - j+1,i,k) = nanmean(spiketrain(tone_ind1_ij))/dt;
                notone_ratemap(num_p_bins - j+1,i,k) = nanmean(spiketrain(no_tone_ind1_ij))/dt;
                
                tone_speedmap(num_p_bins - j+1,i,k) = nanmean(speed(tone_ind1_ij));
                notone_speedmap(num_p_bins - j+1,i,k) = nanmean(speed(no_tone_ind1_ij));
            end
        end
    end
    
    mean_fr_all(k,:) = [nanmean(fr(tone_ind1)) nanmean(fr(no_tone_ind1))];
    mean_speed_all(k,:) = [nanmean(speed(tone_ind1)) nanmean(speed(no_tone_ind1))];
    
end

tone_ratemap = squeeze(nanmean(tone_ratemap,3));
notone_ratemap = squeeze(nanmean(notone_ratemap,3));

mean_fr = nanmean(mean_fr_all);
mean_speed = nanmean(mean_speed_all);

return

