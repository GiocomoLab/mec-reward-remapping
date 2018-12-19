function [all_dist, all_fr_t, all_fr_f] = fr_by_distance(inds, plot_or_not)

% compute fr by distance from reward zone for both free and task, for a
% subset of cells specified by the index array 'inds'; returns long arrays
% of the distances, fr during task, and fr during free for every ratemap
% bin

load('Z:\Users\WButler\D254 Neuralynx\data for manuscript\all_cells_both_boxes_best_iso_downsampled_speed_pos_only');
%^^ loads 'cell_downsampled_trained2' structure with all 778 cells'
%downsampled ratemaps

load('Z:\Users\WButler\D254 Neuralynx\data for manuscript\all_zones_and_zone_centers');
%^^ loads 'all_zones' and 'all_zone_centers' cell arrays

% compute average normalized ratemaps from set of downsampled ratemaps
for i = inds
    ratemap1 = zeros(30,30); ratemap2 = zeros(30,30);
    for j = 1:50
        ratemap1 = ratemap1 + cell_downsampled_trained2(i).ratemaps_free{j};
        ratemap2 = ratemap2 + cell_downsampled_trained2(i).ratemaps_task{j};
    end
    ratemap1_norm = (ratemap1-nanmin(ratemap1(:)))./(nanmax(ratemap1(:))-nanmin(ratemap1(:)));
    ratemap2_norm = (ratemap2-nanmin(ratemap2(:)))./(nanmax(ratemap2(:))-nanmin(ratemap2(:)));
    cell_downsampled_trained2(i).mean_ratemap1norm = ratemap1_norm;
    cell_downsampled_trained2(i).mean_ratemap2norm = ratemap2_norm;
end

all_dist = []; all_fr_t = []; all_fr_f = [];
% pull out firing rate and distance for every bin in every normalized
% average ratemap
for k = inds
    mapf = cell_downsampled_trained2(k).mean_ratemap1norm;
    mapt = cell_downsampled_trained2(k).mean_ratemap2norm;
    zonecen = all_zone_centers{k};
    for i = 1:30
        for j = 1:30
            d = pdist2([j i],zonecen); % swap x and y for pdist
            all_dist = [all_dist d];
            all_fr_t = [all_fr_t mapt(i,j)];
            all_fr_f = [all_fr_f mapf(i,j)];
        end
    end
end

% remove nans
nan_list = isnan(all_fr_f);
all_dist(nan_list) = []; all_fr_f(nan_list) = []; all_fr_t(nan_list) = [];

if plot_or_not % can plot as ribbons with s.e.m.
    dist_bin_edges = [0:12]*2; % in ratemap bins, multiply by 5 to get cm
    ribbon_mean_dist = zeros(12,7);
    % array of values: mean free, sem free, mean task, sem task, mean diff, sem
    % diff, number of bins
    for i = 1:12
        inds = intersect(find(all_dist>=dist_bin_edges(i)),find(all_dist<dist_bin_edges(i+1)));
        ribbon_mean_dist(i,1) = nanmean(all_fr_f(inds));
        ribbon_mean_dist(i,2) = nanstd(all_fr_f(inds))/sqrt(numel(inds));
        ribbon_mean_dist(i,3) = nanmean(all_fr_t(inds));
        ribbon_mean_dist(i,4) = nanstd(all_fr_t(inds))/sqrt(numel(inds));
        ribbon_mean_dist(i,5) = nanmean(all_fr_t(inds)-all_fr_f(inds));
        ribbon_mean_dist(i,6) = nanstd(all_fr_t(inds)-all_fr_f(inds))/sqrt(numel(inds));
        ribbon_mean_dist(i,7) = numel(inds);
    end
    
    % plot results
    figure()
    subplot(1,2,1)
    plot(dist_bin_edges(1:11)*5,ribbon_mean_dist(1:11,1),'b')
    hold on
    plot(dist_bin_edges(1:11)*5,ribbon_mean_dist(1:11,1)+ribbon_mean_dist(1:11,2),'b')
    plot(dist_bin_edges(1:11)*5,ribbon_mean_dist(1:11,1)-ribbon_mean_dist(1:11,2),'b')
    plot(dist_bin_edges(1:11)*5,ribbon_mean_dist(1:11,3),'r')
    plot(dist_bin_edges(1:11)*5,ribbon_mean_dist(1:11,3)+ribbon_mean_dist(1:11,4),'r')
    plot(dist_bin_edges(1:11)*5,ribbon_mean_dist(1:11,3)-ribbon_mean_dist(1:11,4),'r')
    
    % difference between task and free foraging
    subplot(1,2,2)
    plot(dist_bin_edges(1:11)*5,ribbon_mean_dist(1:11,5),'k')
    hold on
    plot(dist_bin_edges(1:11)*5,ribbon_mean_dist(1:11,5)+ribbon_mean_dist(1:11,6),'k')
    plot(dist_bin_edges(1:11)*5,ribbon_mean_dist(1:11,5)-ribbon_mean_dist(1:11,6),'k')
end