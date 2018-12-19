%% Code to compare the average firing rate as a function of the distance to reward zone center
% This script will compute the average speed-position matched rate maps for
% each cell (or load a saved structure).
% Second, this script will use the ratemap, in conjunction with the reward
% zone center, to compute the avg FR vs distance from reward zone center

%% clear workspace

clear all; close all; clc

%% add the paths
windows = 1;
if windows
    addpath('C:\Users\khardcas\Dropbox\Butler Hardcastle rat project\Code for paper\toolbox')
    addpath('C:\Users\khardcas\Dropbox\Butler Hardcastle rat project\Code for paper\data structures')
    addpath('C:\Users\khardcas\Dropbox\Butler Hardcastle rat project\Code for paper\neuralynx to matlab')
else
    addpath('/Users/kiah/Dropbox/Butler Hardcastle rat project/Code for paper/toolbox')
    addpath('/Users/kiah/Dropbox/Butler Hardcastle rat project/Code for paper/data structures')
    addpath('/Users/kiah/Dropbox/Butler Hardcastle rat project/Code for paper/neuralynx to matlab')
end


%% load the data structure that has all of the file names

load all_cells_both_boxes_best_iso_with_scores_model
numCells = numel(cell_info_trained);

%% get the list of animal names
animal_name = cell(numCells,1);
for k = 1:numCells
    animal_name_k = strsplit(cell_info_trained(k).uniqueCellId,'-');
    animal_name{k} = animal_name_k{1};
end

unique_animal_names = unique(animal_name);

%% compute the ratemaps for every cell, in each environment
ratemap1_all = cell(numCells,1);
ratemap2_all = cell(numCells,1);

for k = 1:numCells
    
    % load data
    [x_free, y_free, ~, ~, ~, ~, fr_free, ~] = load_data_and_preprocess(cell_info_trained(k).best_free);
    [x_task, y_task, ~, ~, ~, ~, fr_task, ~] = load_data_and_preprocess(cell_info_trained(k).best_task);
    
    % compute the tuning curves
    [ratemap1_all{k},~] = compute_2d_tuning_curve(x_free,y_free,fr_free,30,[0 0],[ 150 150],1);
    [ratemap2_all{k},~] = compute_2d_tuning_curve(x_task,y_task,fr_task,30,[0 0],[ 150 150],1);
    k
end

%% next, compute the zone and zone center for every cell, unless this has already been computed. then load it
computezones = 0; % if = 0, don't re-compute all the zones
[zone_all,zone_center] = do_compute_all_zones(cell_info_trained,computezones);

%% get the indices
[grid_cell_ind, other_pos_ind, border_cell_ind, ~, ~, ~, ~] = find_grid_border_cells;

%% compute the spacing for each grid cell

load grid_spacing_orientation
%}

%% get the best grid scores

min_grid_score = nan(size(grid_cell_ind));
max_grid_score  = nan(size(grid_cell_ind));
task_session = cell(size(grid_cell_ind));
for k = 1:numel(grid_cell_ind)
    min_grid_score(k) = min(cell_info_trained(grid_cell_ind(k)).grid_score_free,cell_info_trained(grid_cell_ind(k)).grid_score_task);
    max_grid_score(k) = max(cell_info_trained(grid_cell_ind(k)).grid_score_free,cell_info_trained(grid_cell_ind(k)).grid_score_task);
    task_session_k = strsplit(cell_info_trained(grid_cell_ind(k)).best_task,'_');
    task_session{k} = strjoin(task_session_k(1:2),'_');
end
min_grid_score_thresh = 0.35;
best_grid_ind = find(min_grid_score > min_grid_score_thresh);
best_grids = grid_cell_ind(best_grid_ind);
task_session = task_session(best_grid_ind);

min_spacing = min(spacing_free(best_grid_ind),spacing_task(best_grid_ind));
max_spacing = max(spacing_free(best_grid_ind),spacing_task(best_grid_ind));

%% find the translation for each grid cell

binsize = 150/sqrt(numel(ratemap1_all{1}));
maxshift = round(min_spacing*0.8/binsize);
shift_val = nan(numel(best_grids),2);
shift_hyp= nan(numel(best_grids),1);
center_vals = nan(numel(best_grids),1);
max_vals = nan(numel(best_grids),1);
for k = 1:numel(best_grids)
    
    % first, change orientation of task grid cell to match what was
    % observed in free
    orientation_diff = orientation_task(best_grid_ind(k)) - orientation_free(best_grid_ind(k));
    reorient_free_ratemap = imrotate(ratemap1_all{best_grids(k)},orientation_diff);
    task_ratemap = ratemap2_all{best_grids(k)};
    [row_free,col_free] = size(task_ratemap);
    [row_reorient,col_reorient] = size(reorient_free_ratemap);
    middle_row_reorient_free_ratemap = round(row_reorient/2);
    middle_col_reorient_free_ratemap = round(col_reorient/2);
    free_inside_row = middle_row_reorient_free_ratemap - round(row_free/2):middle_row_reorient_free_ratemap + round(row_free/2)-1;
    free_inside_col = middle_col_reorient_free_ratemap - round(col_free/2):middle_col_reorient_free_ratemap + round(col_free/2)-1;
    middle_reorient_free_ratemap = reorient_free_ratemap(free_inside_row,free_inside_col);
   
    [cross_corr,numpoints] = correlation_kh(middle_reorient_free_ratemap,task_ratemap);
    cross_corr(isnan(cross_corr)) = 0;
    
    % find the nearest local maximum to the center
    % all possible peaks have to be within 75% of the max of center
    inside = round(size(cross_corr,1)/2)-maxshift(k):round(size(cross_corr,1)/2)+maxshift(k);
    cross_corr_inside = cross_corr(inside,inside);
    max_inside = max(max(cross_corr_inside));
    
    cross_corr_thresh = cross_corr; cross_corr_thresh(cross_corr < 0.5*max_inside) = 0;
    peaks_in_corr = imregionalmax(cross_corr_thresh);
    [row,col] = find(peaks_in_corr);
    middle = round(size(cross_corr,1)/2);
    dist_to_center = sqrt((row - middle).^2 + ((col - middle).^2));
    [shift_hyp(k),ind] = min(dist_to_center);
    ns = middle - row(ind); ew = col(ind)-middle;
    shift_val(k,:) = [ns ew]*binsize;
    
    %{
    figure(1)
    subplot(1,3,1)
    imagesc(middle_reorient_free_ratemap)
    axis off
    title('ENV1')
    subplot(1,3,2)
    imagesc(task_ratemap)
    axis off
    title('ENV2')
    subplot(1,3,3)
    imagesc(cross_corr)
    title(['cross-corr, shift = ',num2str(shift_val(k,:))])
    hold on
    plot(middle,middle,'*k','markersize',10)
    plot(col(ind),row(ind),'*r','markersize',10)
    hold off
    axis off
    keyboard
    %}
    
end

hyp = sqrt(shift_val(:,1).^2 + shift_val(:,2).^2);

[n,x] = hist(hyp./min_spacing,20);
figure(1)
bar(x,n,'k')
ylabel('# of cells')
xlabel('translation (cm) / spacing (cm)')
set(gca,'fontsize',15)



%% look at the shift per animal

zone_center_grid = zone_center(best_grids,:);
count = 0;
angle_diff = [];
iter = 500;
for k = 1:numel(unique_animal_names)
    
    ind = intersect(find(ismember(animal_name,unique_animal_names{k})),best_grids);
    grid_ind = find(ismember(best_grids,ind));
    if numel(ind) > 0
        figure(2)
        count = count + 1;
        subplot(1,7,count)
        axis_ind = nan(numel(grid_ind,1));
        angle1_all= nan(numel(grid_ind,1));
        for j = 1:numel(ind)
            plot([0 shift_val(grid_ind(j),2) ],[0 shift_val(grid_ind(j),1) ],'k','linewidth',2)
            hold on
            
            % find angle
            hyp1 = sqrt(shift_val(grid_ind(j),1).^2 + shift_val(grid_ind(j),2).^2);
            if hyp1 > 0
            angle1 = acosd(shift_val(grid_ind(j),2)/hyp1)*sign(shift_val(grid_ind(j),1));
            angle1(angle1 < 0) = 360+angle1(angle1<0);
            angle1_all(j) = angle1;
            else
                angle1_all(j) = NaN;
            end
        end

        zone_center_k = [zone_center_grid(grid_ind(1),1)-75 zone_center_grid(grid_ind(1),2)-75];
        
        hyp_zc = sqrt(zone_center_k(1).^2 + zone_center_k(2).^2);
        angle_zc = acosd(zone_center_k(1)/hyp_zc)*sign(zone_center_k(2));
        angle_zc(angle_zc < 0) = 360+angle_zc(angle_zc<0);
        angle1_all1 = angle1_all; angle1_all1(angle1_all1 == 0) = NaN;
        
        
        plot((zone_center_grid(grid_ind(1),1)-75 - 10)*ones(1,2),[zone_center_grid(grid_ind(1),2)-85 zone_center_grid(grid_ind(1),2)-65],'r')
        hold on
        plot((zone_center_grid(grid_ind(1),1)-75 + 10)*ones(1,2),[zone_center_grid(grid_ind(1),2)-85 zone_center_grid(grid_ind(1),2)-65],'r')
        plot([zone_center_grid(grid_ind(1),1)-85 zone_center_grid(grid_ind(1),1)-65],(zone_center_grid(grid_ind(1),2)-75-10)*ones(1,2),'r')
        plot([zone_center_grid(grid_ind(1),1)-85 zone_center_grid(grid_ind(1),1)-65],(zone_center_grid(grid_ind(1),2)-75+10)*ones(1,2),'r')
        hold off
        axis([-100 100 -100 100])
    end
end


angle_diff(angle_diff > 180) = angle_diff(angle_diff > 180)-360;
figure()
[n,x] = hist(angle_diff,30);
bar(x,n,'k')
box off
hold on
plot([-90 -90],[0 max(n)],'--r','linewidth',2)
plot([90 90],[0 max(n)],'--r','linewidth',2)
hold off
title(['fraction moving towards: ',num2str(numel(find(abs(angle_diff < 90)))),'/',num2str(sum(isfinite(angle_diff)))])


%% compute which cells exhibit a significant shift
iter = 100;
null_shift_val_free = nan(numel(best_grids),iter);
null_shift_val_task = nan(numel(best_grids),iter);
minbin =33;
for k = 1:numel(best_grids)
    
    % FREE
    % get the session information
    [x_free, y_free, ~, ~, ~, ~, fr_free, ~] = load_data_and_preprocess(cell_info_trained(best_grids(k)).best_free);
    
    % split into segments
    T = numel(x_free);
    index = 1:floor(T/minbin)*minbin;
    indexmat = reshape(index,minbin,numel(index)/minbin);
    [~,numpoints] = size(indexmat);
    
    tic
    for j = 1:iter
        
        first_col = datasample(1:numpoints,round(numpoints/2),'Replace',false);
        second_col = setdiff(1:numpoints,first_col);
        
        first = indexmat(:,first_col); first = first(:);
        second = indexmat(:,second_col); second = second(:);
        
        [ratemap1_j_first,~] = compute_2d_tuning_curve(x_free(first),y_free(first),fr_free(first),30,[0 0],[ 150 150],1);
        [ratemap1_j_second,~] = compute_2d_tuning_curve(x_free(second),y_free(second),fr_free(second),30,[0 0],[ 150 150],1);
        
        [cross_corr,~] = correlation_kh(ratemap1_j_first,ratemap1_j_second);
        cross_corr(isnan(cross_corr)) = 0;
        
        % find the nearest local maximum to the center
        % all possible peaks have to be within 75% of the max of center
        inside = round(size(cross_corr,1)/2)-maxshift(k):round(size(cross_corr,1)/2)+maxshift(k);
        cross_corr_inside = cross_corr(inside,inside);
        max_inside = max(max(cross_corr_inside));
        
        cross_corr_thresh = cross_corr; cross_corr_thresh(cross_corr < 0.75*max_inside) = 0;
        peaks_in_corr = imregionalmax(cross_corr_thresh);
        [row,col] = find(peaks_in_corr);
        middle = round(size(cross_corr,1)/2);
        dist_to_center = sqrt((row - middle).^2 + ((col - middle).^2));
        [shift_hyp(k),ind] = min(dist_to_center);
        ns = middle - row(ind); ew = col(ind)-middle;
        shift_val = [ns ew]*binsize;
        null_shift_val_free(k,j) = sqrt(shift_val(1).^2 + shift_val(2).^2);
    end
    toc
    
    
   % FREE
    % get the session information
    [x_task, y_task, ~, ~, ~, ~, fr_task, ~] = load_data_and_preprocess(cell_info_trained(best_grids(k)).best_task);
    
    % split into segments
    T = numel(x_task);
    index = 1:floor(T/minbin)*minbin;
    indexmat = reshape(index,minbin,numel(index)/minbin);
    [~,numpoints] = size(indexmat);
    
    for j = 1:iter
        
        first_col = datasample(1:numpoints,round(numpoints/2),'Replace',false);
        second_col = setdiff(1:numpoints,first_col);
        
        first = indexmat(:,first_col); first = first(:);
        second = indexmat(:,second_col); second = second(:);
        
        [ratemap1_j_first,~] = compute_2d_tuning_curve(x_task(first),y_task(first),fr_task(first),30,[0 0],[ 150 150],1);
        [ratemap1_j_second,~] = compute_2d_tuning_curve(x_task(second),y_task(second),fr_task(second),30,[0 0],[ 150 150],1);
        
        [cross_corr,~] = correlation_kh(ratemap1_j_first,ratemap1_j_second);
        cross_corr(isnan(cross_corr)) = 0;
        
        % find the nearest local maximum to the center
        % all possible peaks have to be within 75% of the max of center
        inside = round(size(cross_corr,1)/2)-maxshift(k):round(size(cross_corr,1)/2)+maxshift(k);
        cross_corr_inside = cross_corr(inside,inside);
        max_inside = max(max(cross_corr_inside));
        
        cross_corr_thresh = cross_corr; cross_corr_thresh(cross_corr < 0.75*max_inside) = 0;
        peaks_in_corr = imregionalmax(cross_corr_thresh);
        [row,col] = find(peaks_in_corr);
        middle = round(size(cross_corr,1)/2);
        dist_to_center = sqrt((row - middle).^2 + ((col - middle).^2));
        [shift_hyp(k),ind] = min(dist_to_center);
        ns = middle - row(ind); ew = col(ind)-middle;
        shift_val = [ns ew]*binsize;
        null_shift_val_task(k,j) = sqrt(shift_val(1).^2 + shift_val(2).^2);
        j
    end
    
    k
    
end


%% compute the number of cells that shift
sig_translate = nan(size(best_grids));
for k = 1:numel(best_grids)
    
    hyp_k = hyp(k);
    
    null_dist = [null_shift_val_free(k,:) null_shift_val_task(k,:)];
    
    sig_translate(k) = sum(null_dist >= hyp_k)/numel(null_dist);
    
end





