function [field_distance_null] = generate_field_distance_null(thresh_maps,zone_centers,numint)


numratemaps = numel(thresh_maps);
numbins = size(thresh_maps{1},1);
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

field_distance_null = nan(numint,1);
nummaps = numel(thresh_maps);
for k = 1:numint
    
    % first, compute the field properties
    index_field = randi(nummaps);
    [~,~,~,field_com] = compute_field_properties(thresh_maps{index_field});
    
    % next, compute the distance of the nearest field to zone
    index_zone = randi(numcenters);
    zone_center1 = [null_center_x(index_zone) null_center_y(index_zone)];
    [min_dist_free_k,~] = min(pdist2(field_com,zone_center1));
    if ~isempty(min_dist_free_k)
        field_distance_null(k) = min_dist_free_k;
    end
end


return