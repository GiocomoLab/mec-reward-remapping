function [fr_by_distance,distance_vector,fr_by_distance_sem] = compute_fr_by_distance_1map(ratemap,zone_center,posvec_allx,posvec_ally,numdistbins,minval,maxval,normalize)

% first, normalize the ratemaps
if normalize
    ratemap_norm = (ratemap(:) - nanmin(ratemap(:)))./(nanmax(ratemap(:)) - nanmin(ratemap(:)));
else
    ratemap_norm = ratemap(:);
end

% second, compute the distance of every bin to the reward zone
dist_to_center = sqrt((zone_center(1) - posvec_allx').^2 + (zone_center(2) - posvec_ally').^2);

% third, compute a fr vs distance tuning curve
[fr_by_distance,distance_vector,fr_by_distance_sem] = compute_1d_tuning_curve(dist_to_center,ratemap_norm,numdistbins,minval,maxval);

return