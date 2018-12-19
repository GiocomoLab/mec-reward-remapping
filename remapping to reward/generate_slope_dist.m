function [slope_all,activity_diff] = generate_slope_dist(ratemaps,zone_centers,numint)


numratemaps = numel(ratemaps);
numbins = size(ratemaps{1},1);
%% generate zone centers
boxSize = 150;

% compute the average distance from nearest all for all zone centers
dist_to_wall = nan(length(zone_centers),1);
for k  = 1:length(zone_centers)
    minx = min(abs(zone_centers(k,1) - 150),abs(zone_centers(k,1) - 0));
    miny = min(abs(zone_centers(k,2) - 150),abs(zone_centers(k,2) - 0));
    dist_to_wall(k) = min(minx,miny);
end
dist_to_wall = unique(dist_to_wall);
min_dist_to_wall = min(dist_to_wall);
max_dist_to_wall = max(dist_to_wall);

% generate grid of potential centers
posvec = linspace(0,boxSize,boxSize+1); binw = posvec(2)-posvec(1);
posvec_x = binw/2:binw:boxSize-binw/2;
posvec_y = flip(posvec_x);
[rows,cols] = ind2sub([length(posvec_x),length(posvec_x)],1:(length(posvec_x))^2);
posvec_allx = posvec_x(cols);
posvec_ally = posvec_y(rows);

% find the distance of each possible center to the nearest wall
dist_to_wall_all = nan(numel(posvec_allx),1);
for k = 1:numel(posvec_allx)
    minx = min(abs(posvec_allx(k) - 150),abs(posvec_allx(k) - 0));
    miny = min(abs(posvec_ally(k) - 150),abs(posvec_ally(k) - 0));
    dist_to_wall_all(k) = min(minx,miny);
end

% reject the ones that are outside the range we observe in the data
null_center_x = posvec_allx;  null_center_y = posvec_ally;
null_center_x(dist_to_wall_all > max_dist_to_wall | dist_to_wall_all < min_dist_to_wall) = [];
null_center_y(dist_to_wall_all > max_dist_to_wall | dist_to_wall_all < min_dist_to_wall) = [];
numcenters = numel(null_center_x);

%% pick a ratemap, pick a zone center, and compute the fr distance tuning curve

% re-compute the position of every bin
boxSize = 150;
posvec = linspace(0,boxSize,numbins+1); binw = posvec(2)-posvec(1);
posvec_x = binw/2:binw:boxSize-binw/2;
posvec_y = flip(posvec_x);
[rows,cols] = ind2sub([numbins,numbins],1:numbins^2);
posvec_allx = posvec_x(cols);
posvec_ally = posvec_y(rows);
numdistbins = 10; minval = 0; maxval = 100;
distvec1 = linspace(minval,maxval,numdistbins+1); bw = distvec1(2) - distvec1(1);
distvec = bw/2+minval:bw:maxval-bw/2;

slope_all = nan(numint,1);
activity_diff = nan(numint,1);
for k = 1:numint
    
    % first, grab and normalize ratemap
    randind = randi(numratemaps);
    ratemap1 = ratemaps{randind};
    ratemap1_norm = (ratemap1(:) - nanmin(ratemap1(:)))./(nanmax(ratemap1(:)) - nanmin(ratemap1(:)));
    
    % second, grab zone and compute distance to center
    dist_to_orig_center = sqrt((null_center_x - zone_centers(randind,1)).^2 + (null_center_y - zone_centers(randind,2)).^2);
    good_centers = find(dist_to_orig_center > 30);
    index = randi(numel(good_centers));
    zone_center1 = [null_center_x(good_centers(index)) null_center_y(good_centers(index))];
    dist_to_center = sqrt((zone_center1(1) - posvec_allx').^2 + (zone_center1(2) - posvec_ally').^2);
    
    % third, compute a fr vs distance tuning curve
    [fr_by_distance,~,~] = compute_1d_tuning_curve(dist_to_center,ratemap1_norm,numdistbins,0,100);
    bestline = polyfit(distvec,fr_by_distance',1); slope_all(k) = bestline(1);
    activity_diff(k) = nanmean(ratemap1_norm(dist_to_center < 30)) - nanmean(ratemap1_norm(dist_to_center > 80));
end


return