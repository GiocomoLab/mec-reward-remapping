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

%% compute the avg speed-position matched downsampled rate maps for every cell

computemaps = 0; % if = 0, don't re-compute all the maps
[ratemap1_all,ratemap2_all] = do_compute_speed_matched_ratemaps(cell_info_trained,computemaps);

%% next, compute the zone and zone center for every cell, unless this has already been computed. then load it

computezones = 0; % if = 0, don't re-compute all the zones
[zone_all,zone_center] = do_compute_all_zones(cell_info_trained,computezones);

%% compute the firing rate by distance for every cell

% go through every ratemap and zone center combo, and compute the average
% norm firing rate (avg spatially) vs the reward zone
[fr_by_distance_free,fr_by_distance_task,distvec] = compute_fr_by_distance(ratemap1_all,ratemap2_all,zone_center);
fr_diff = fr_by_distance_task - fr_by_distance_free;

%% compute the slopes
slope_free = nan(numCells,1); slope_task = nan(numCells,1); slope_diff = nan(numCells,1);
for k = 1:numCells
    slope_free1 = polyfit(distvec,fr_by_distance_free(k,:),1); slope_free(k) = slope_free1(1);
    slope_task1 = polyfit(distvec,fr_by_distance_task(k,:),1); slope_task(k) = slope_task1(1);
    slope_diff1 = polyfit(distvec,fr_diff(k,:),1); slope_diff(k) = slope_diff1(1);
end

%% get the indices
[grid_cell_ind, other_pos_ind, border_cell_ind, ~, ~, ~, ~] = find_grid_border_cells;

%% do a bunch of plotting 

% grid cells
figure(1)
subplot(3,3,1)
top = nanmean(fr_by_distance_free(grid_cell_ind,:))+ nanstd(fr_by_distance_free(grid_cell_ind,:))/sqrt(numel(grid_cell_ind));
bottom = nanmean(fr_by_distance_free(grid_cell_ind,:)) - nanstd(fr_by_distance_free(grid_cell_ind,:))/sqrt(numel(grid_cell_ind));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[1 0 0]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_by_distance_free(grid_cell_ind,:)),'r','linewidth',2)
alpha(h,0.3)
top = nanmean(fr_by_distance_task(grid_cell_ind,:))+ nanstd(fr_by_distance_task(grid_cell_ind,:))/sqrt(numel(grid_cell_ind));
bottom = nanmean(fr_by_distance_task(grid_cell_ind,:)) - nanstd(fr_by_distance_task(grid_cell_ind,:))/sqrt(numel(grid_cell_ind));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[0 0 1]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_by_distance_task(grid_cell_ind,:)),'b','linewidth',2)
alpha(h,0.3)
hold off
box off
xlabel('distance from reward zone center')
ylabel('normalized firing rate difference')
legend('free','task')


subplot(3,3,4)
top = nanmean(fr_diff(grid_cell_ind,:))+ nanstd(fr_diff(grid_cell_ind,:))/sqrt(numel(grid_cell_ind));
bottom = nanmean(fr_diff(grid_cell_ind,:)) - nanstd(fr_diff(grid_cell_ind,:))/sqrt(numel(grid_cell_ind));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[0 0 0]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_diff(grid_cell_ind,:)),'k','linewidth',2)
alpha(h,0.3)
plot(distvec,zeros(size(distvec)),'--r','linewidth',2)
hold off
box off
xlabel('distance from reward zone center')
ylabel('normalized firing rate difference')
legend('difference','baseline firing')

subplot(3,3,7)
[n,x] = hist(slope_free(grid_cell_ind),20);
[nn,xx] = hist(slope_task(grid_cell_ind),20);
[nnn,xxx] = hist(slope_diff(grid_cell_ind),20);
%s1 = bar(x,n,'r');
%s2 = bar(xx,nn,'b');
s3 = bar(xxx,nnn,'k');
hold on
alpha(0.3)
plot([0 0],[0 max([n nn nnn])],'--r','linewidth',2)
hold off
legend('free','task','difference')
box off
xlabel('slope value')
ylabel('# of cells')
p1 = signrank(slope_free(grid_cell_ind));
p2 = signrank(slope_task(grid_cell_ind));
p3 = signrank(slope_diff(grid_cell_ind));
title(['p_f = ',num2str(p1),' p_s = ',num2str(p2),' p_d = ',num2str(p3)])

% border cells
subplot(3,3,2)
top = nanmean(fr_by_distance_free(border_cell_ind,:))+ nanstd(fr_by_distance_free(border_cell_ind,:))/sqrt(numel(border_cell_ind));
bottom = nanmean(fr_by_distance_free(border_cell_ind,:)) - nanstd(fr_by_distance_free(border_cell_ind,:))/sqrt(numel(border_cell_ind));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[1 0 0]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_by_distance_free(border_cell_ind,:)),'r','linewidth',2)
alpha(h,0.3)
top = nanmean(fr_by_distance_task(border_cell_ind,:))+ nanstd(fr_by_distance_task(border_cell_ind,:))/sqrt(numel(border_cell_ind));
bottom = nanmean(fr_by_distance_task(border_cell_ind,:)) - nanstd(fr_by_distance_task(border_cell_ind,:))/sqrt(numel(border_cell_ind));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[0 0 1]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_by_distance_task(border_cell_ind,:)),'b','linewidth',2)
alpha(h,0.3)
box off
xlabel('distance from reward zone center')
ylabel('normalized firing rate difference')
legend('free','task')

subplot(3,3,5)
top = nanmean(fr_diff(border_cell_ind,:))+ nanstd(fr_diff(border_cell_ind,:))/sqrt(numel(border_cell_ind));
bottom = nanmean(fr_diff(border_cell_ind,:)) - nanstd(fr_diff(border_cell_ind,:))/sqrt(numel(border_cell_ind));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[0 0 0]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_diff(border_cell_ind,:)),'k','linewidth',2)
alpha(h,0.3)
plot(distvec,zeros(size(distvec)),'--r','linewidth',2)
hold off
box off
xlabel('distance from reward zone center')
ylabel('normalized firing rate difference')
legend('difference','baseline firing')

subplot(3,3,8)
[n,x] = hist(slope_free(border_cell_ind),20);
[nn,xx] = hist(slope_task(border_cell_ind),20);
[nnn,xxx] = hist(slope_diff(border_cell_ind),20);
s3 = bar(xxx,nnn,'k');
hold on
alpha(0.3)
plot([0 0],[0 max([n nn nnn])],'--r','linewidth',2)
hold off
legend('free','task','difference')
box off
xlabel('slope value')
ylabel('# of cells')
p1 = signrank(slope_free(border_cell_ind));
p2 = signrank(slope_task(border_cell_ind));
p3 = signrank(slope_diff(border_cell_ind));
title(['p_f = ',num2str(p1),' p_s = ',num2str(p2),' p_d = ',num2str(p3)])

% non-grid non-border cells
subplot(3,3,3)
top = nanmean(fr_by_distance_free(other_pos_ind,:))+ nanstd(fr_by_distance_free(other_pos_ind,:))/sqrt(numel(other_pos_ind));
bottom = nanmean(fr_by_distance_free(other_pos_ind,:)) - nanstd(fr_by_distance_free(other_pos_ind,:))/sqrt(numel(other_pos_ind));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[1 0 0]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_by_distance_free(other_pos_ind,:)),'r','linewidth',2)
alpha(h,0.3)
top = nanmean(fr_by_distance_task(other_pos_ind,:))+ nanstd(fr_by_distance_task(other_pos_ind,:))/sqrt(numel(other_pos_ind));
bottom = nanmean(fr_by_distance_task(other_pos_ind,:)) - nanstd(fr_by_distance_task(other_pos_ind,:))/sqrt(numel(other_pos_ind));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[0 0 1]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_by_distance_task(other_pos_ind,:)),'b','linewidth',2)
alpha(h,0.3)
box off
xlabel('distance from reward zone center')
ylabel('normalized firing rate difference')
legend('free','task')

subplot(3,3,6)
top = nanmean(fr_diff(other_pos_ind,:))+ nanstd(fr_diff(other_pos_ind,:))/sqrt(numel(other_pos_ind));
bottom = nanmean(fr_diff(other_pos_ind,:)) - nanstd(fr_diff(other_pos_ind,:))/sqrt(numel(other_pos_ind));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[0 0 0]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_diff(other_pos_ind,:)),'k','linewidth',2)
alpha(h,0.3)
plot(distvec,zeros(size(distvec)),'--r','linewidth',2)
hold off
box off
xlabel('distance from reward zone center')
ylabel('normalized firing rate difference')
legend('difference','baseline firing')

subplot(3,3,9)
[n,x] = hist(slope_free(other_pos_ind),20);
[nn,xx] = hist(slope_task(other_pos_ind),20);
[nnn,xxx] = hist(slope_diff(other_pos_ind),20);
s3 = bar(xxx,nnn,'k');
hold on
alpha(0.3)
plot([0 0],[0 max([n nn nnn])],'--r','linewidth',2)
hold off
legend('free','task','difference')
box off
xlabel('slope value')
ylabel('# of cells')
p1 = signrank(slope_free(other_pos_ind));
p2 = signrank(slope_task(other_pos_ind));
p3 = signrank(slope_diff(other_pos_ind));
title(['p_f = ',num2str(p1),' p_s = ',num2str(p2),' p_d = ',num2str(p3)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate estimate of slopes in env1
iter = 500;
med_slope_grid_env1 = nan(iter,1);
med_activity_diff_grid_env1= nan(iter,1);
for k = 1:iter
    [slope_all,act_diff] = generate_slope_dist(ratemap1_all(grid_cell_ind),zone_center(grid_cell_ind,:),numel(grid_cell_ind));
    med_slope_grid_env1(k) = median(slope_all);
    med_activity_diff_grid_env1(k) = median(act_diff);
    k
end

figure(1)
[n,x] = hist(med_slope_grid_env1,20);
bar(x,n,'k')
hold on
plot(median(slope_free(grid_cell_ind))*ones(1,2),[0 max(n)],'--r','linewidth',2)
hold off
box off
ylabel('# of iterations')
xlabel('slope in ENV1')

figure(2)
[n,x] = hist(med_activity_diff_grid_env1,20);
bar(x,n,'k')
hold on
plot(median(act_diff_free(grid_cell_ind))*ones(1,2),[0 max(n)],'--r','linewidth',2)
hold off
box off
ylabel('# of iterations')
xlabel('slope in ENV1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% randomly sample 102 nonspatial cells to get estimate of grid cell like dist - is effect bigger?

iter = 5000;
null_dist = nan(iter,1);
for k = 1:iter
    null_dist(k) = nanmedian(datasample(slope_diff(other_pos_ind),numel(grid_cell_ind)));
end

figure(1)
[n,x] = hist(null_dist,20);
bar(x,n,'k')
hold on
plot(median(slope_diff(grid_cell_ind))*ones(1,2),[0 max(n)],'--r','linewidth',2)
plot(median(slope_diff(other_pos_ind))*ones(1,2),[0 max(n)],'--b','linewidth',2)
hold off
box off
ylabel('# of iterations')
xlabel('median slope')
set(gca,'fontsize',15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% look at the same thing for head direction cells

mvl_grid_free = nan(numel(grid_cell_ind),1);
mvl_grid_task = nan(numel(grid_cell_ind),1);

for k = 1:numel(grid_cell_ind)
    mvl_grid_free(k) = cell_info_trained(grid_cell_ind(k)).hd_score_free;
    mvl_grid_task(k) = cell_info_trained(grid_cell_ind(k)).hd_score_task;
end

mean_mvl = (mvl_grid_free + mvl_grid_task)/2;
[n,x] = hist(mean_mvl,20);

figure(3)
subplot(1,2,1)
plot(mvl_grid_free,mvl_grid_task,'ok','markersize',7)
hold on
plot([0 1],0.2*ones(1,2),'--r','linewidth',2)
plot(0.2*ones(1,2),[0 1],'--r','linewidth',2)
hold off
box off
ylabel('MVL (ENV1)')
xlabel('MVL (ENV2)')

subplot(1,2,2)
bar(x,n,'k')
ylabel('# of cells')
xlabel('MVL')
box off
set(gca,'fontsize',15)

high_hd_thresh = 0.2;
low_hd_thresh = 0.2;
grid_hd = grid_cell_ind(mean_mvl > high_hd_thresh);
grid_nonhd = grid_cell_ind(mean_mvl <= low_hd_thresh);


figure(1)
subplot(3,2,1)
top = nanmean(fr_by_distance_free(grid_hd,:))+ nanstd(fr_by_distance_free(grid_hd,:))/sqrt(numel(grid_hd));
bottom = nanmean(fr_by_distance_free(grid_hd,:)) - nanstd(fr_by_distance_free(grid_hd,:))/sqrt(numel(grid_hd));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[0 0 1]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_by_distance_free(grid_hd,:)),'b','linewidth',2)
alpha(h,0.3)
top = nanmean(fr_by_distance_task(grid_hd,:))+ nanstd(fr_by_distance_task(grid_hd,:))/sqrt(numel(grid_hd));
bottom = nanmean(fr_by_distance_task(grid_hd,:)) - nanstd(fr_by_distance_task(grid_hd,:))/sqrt(numel(grid_hd));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[1 0 0]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_by_distance_task(grid_hd,:)),'r','linewidth',2)
alpha(h,0.3)
hold off
box off
xlabel('distance from reward zone center')
ylabel('normalized firing rate difference')
legend('free','task')


subplot(3,2,3)
top = nanmean(fr_diff(grid_hd,:))+ nanstd(fr_diff(grid_hd,:))/sqrt(numel(grid_hd));
bottom = nanmean(fr_diff(grid_hd,:)) - nanstd(fr_diff(grid_hd,:))/sqrt(numel(grid_hd));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[0 0 0]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_diff(grid_hd,:)),'k','linewidth',2)
alpha(h,0.3)
plot(distvec,zeros(size(distvec)),'--r','linewidth',2)
hold off
box off
xlabel('distance from reward zone center')
ylabel('normalized firing rate difference')
legend('difference','baseline firing')

subplot(3,2,5)
[n,x] = hist(slope_free(grid_hd),20);
[nn,xx] = hist(slope_task(grid_hd),20);
[nnn,xxx] = hist(slope_diff(grid_hd),20);
s3 = bar(xxx,nnn,'k');
hold on
alpha(0.3)
plot([0 0],[0 max([n nn nnn])],'--r','linewidth',2)
hold off
legend('free','task','difference')
box off
xlabel('slope value')
ylabel('# of cells')
p1 = signrank(slope_free(grid_hd));
p2 = signrank(slope_task(grid_hd));
p3 = signrank(slope_diff(grid_hd));
title(['p_f = ',num2str(p1),' p_s = ',num2str(p2),' p_d = ',num2str(p3)])

subplot(3,2,2)
top = nanmean(fr_by_distance_free(grid_nonhd,:))+ nanstd(fr_by_distance_free(grid_nonhd,:))/sqrt(numel(grid_nonhd));
bottom = nanmean(fr_by_distance_free(grid_nonhd,:)) - nanstd(fr_by_distance_free(grid_nonhd,:))/sqrt(numel(grid_nonhd));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[0 0 1]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_by_distance_free(grid_nonhd,:)),'b','linewidth',2)
alpha(h,0.3)
top = nanmean(fr_by_distance_task(grid_nonhd,:))+ nanstd(fr_by_distance_task(grid_nonhd,:))/sqrt(numel(grid_nonhd));
bottom = nanmean(fr_by_distance_task(grid_nonhd,:)) - nanstd(fr_by_distance_task(grid_nonhd,:))/sqrt(numel(grid_nonhd));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[1 0 0]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_by_distance_task(grid_nonhd,:)),'r','linewidth',2)
alpha(h,0.3)
hold off
box off
xlabel('distance from reward zone center')
ylabel('normalized firing rate difference')
legend('free','task')


subplot(3,2,4)
top = nanmean(fr_diff(grid_nonhd,:))+ nanstd(fr_diff(grid_nonhd,:))/sqrt(numel(grid_nonhd));
bottom = nanmean(fr_diff(grid_nonhd,:)) - nanstd(fr_diff(grid_nonhd,:))/sqrt(numel(grid_nonhd));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[0 0 0]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_diff(grid_nonhd,:)),'k','linewidth',2)
alpha(h,0.3)
plot(distvec,zeros(size(distvec)),'--r','linewidth',2)
hold off
box off
xlabel('distance from reward zone center')
ylabel('normalized firing rate difference')
legend('difference','baseline firing')

subplot(3,2,6)
[n,x] = hist(slope_free(grid_nonhd),20);
[nn,xx] = hist(slope_task(grid_nonhd),20);
[nnn,xxx] = hist(slope_diff(grid_nonhd),20);
s3 = bar(xxx,nnn,'k');
hold on
alpha(0.3)
plot([0 0],[0 max([n nn nnn])],'--r','linewidth',2)
hold off
legend('free','task','difference')
box off
xlabel('slope value')
ylabel('# of cells')
p1 = signrank(slope_free(grid_nonhd));
p2 = signrank(slope_task(grid_nonhd));
p3 = signrank(slope_diff(grid_nonhd));
title(['p_f = ',num2str(p1),' p_s = ',num2str(p2),' p_d = ',num2str(p3)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the above, but for each animal

num_animals = numel(unique_animal_names);
slope_pos = nan(num_animals,1);
slope_other_animal = nan(num_animals,2);
slope_grid_animal = nan(num_animals,2);
for k = 1:num_animals
    
    animal_ind = find(strcmp(unique_animal_names{k},animal_name));
    
    figure(1)
    subplot(1,3,1)
    errorbar(distvec,nanmean(fr_diff(intersect(grid_cell_ind,animal_ind),:)),nanstd(fr_diff(intersect(grid_cell_ind,animal_ind),:))/sqrt(numel(intersect(grid_cell_ind,animal_ind))),'linewidth',2)
    xlabel('distance from reward zone center')
    ylabel('normalized firing rate difference')
    hold on
    box off
    
    subplot(1,3,2)
    errorbar(distvec,nanmean(fr_diff(intersect(border_cell_ind,animal_ind),:)),nanstd(fr_diff(intersect(border_cell_ind,animal_ind),:))/sqrt(numel(intersect(border_cell_ind,animal_ind))),'linewidth',2)
    xlabel('distance from reward zone center')
    ylabel('normalized firing rate difference')
    hold on
    box off
    
    subplot(1,3,3)
    errorbar(distvec,nanmean(fr_diff(intersect(other_pos_ind,animal_ind),:)),nanstd(fr_diff(intersect(other_pos_ind,animal_ind),:))/sqrt(numel(intersect(other_pos_ind,animal_ind))),'linewidth',2)
    xlabel('distance from reward zone center')
    ylabel('normalized firing rate difference')
    hold on
    box off
    
    slope = polyfit(distvec,nanmean(fr_diff(intersect(other_pos_ind,animal_ind),:)),1); slope_pos(k) = slope(1);
    
    cells_k = intersect(grid_cell_ind,animal_ind);
    figure(2)
    subplot(2,4,k)
    top = nanmean(fr_by_distance_free(cells_k,:))+ nanstd(fr_by_distance_free(cells_k,:))/sqrt(numel(cells_k));
    bottom = nanmean(fr_by_distance_free(cells_k,:)) - nanstd(fr_by_distance_free(cells_k,:))/sqrt(numel(cells_k));
    h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[0 0 1]);
    set(h,'EdgeColor','None')
    hold on
    plot(distvec,nanmean(fr_by_distance_free(cells_k,:)),'b','linewidth',2)
    alpha(h,0.3)
    top = nanmean(fr_by_distance_task(cells_k,:))+ nanstd(fr_by_distance_task(cells_k,:))/sqrt(numel(cells_k));
    bottom = nanmean(fr_by_distance_task(cells_k,:)) - nanstd(fr_by_distance_task(cells_k,:))/sqrt(numel(cells_k));
    h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[1 0 0]);
    set(h,'EdgeColor','None')
    hold on
    plot(distvec,nanmean(fr_by_distance_task(cells_k,:)),'r','linewidth',2)
    alpha(h,0.3)
    hold off
    box off
    xlabel('distance from reward zone center')
    ylabel('normalized firing rate difference')
    title(['# of cells = ',num2str(numel(cells_k))])
    
    slope_grid_animal_k = polyfit(distvec,nanmean(fr_by_distance_task(cells_k,:)),1);
    slope_grid_animal(k,2) = slope_grid_animal_k(1);
    slope_grid_animal_k = polyfit(distvec,nanmean(fr_by_distance_free(cells_k,:)),1);
    slope_grid_animal(k,1) = slope_grid_animal_k(1);
    
    cells_k = intersect(other_pos_ind,animal_ind);
    figure(3)
    subplot(2,4,k)
    top = nanmean(fr_by_distance_free(cells_k,:))+ nanstd(fr_by_distance_free(cells_k,:))/sqrt(numel(cells_k));
    bottom = nanmean(fr_by_distance_free(cells_k,:)) - nanstd(fr_by_distance_free(cells_k,:))/sqrt(numel(cells_k));
    h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[0 0 1]);
    set(h,'EdgeColor','None')
    hold on
    plot(distvec,nanmean(fr_by_distance_free(cells_k,:)),'b','linewidth',2)
    alpha(h,0.3)
    top = nanmean(fr_by_distance_task(cells_k,:))+ nanstd(fr_by_distance_task(cells_k,:))/sqrt(numel(cells_k));
    bottom = nanmean(fr_by_distance_task(cells_k,:)) - nanstd(fr_by_distance_task(cells_k,:))/sqrt(numel(cells_k));
    h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[1 0 0]);
    set(h,'EdgeColor','None')
    hold on
    plot(distvec,nanmean(fr_by_distance_task(cells_k,:)),'r','linewidth',2)
    alpha(h,0.3)
    hold off
    box off
    xlabel('distance from reward zone center')
    ylabel('normalized firing rate difference')
    title(['# of cells = ',num2str(numel(cells_k))])
    
    slope_other_animal_k = polyfit(distvec,nanmean(fr_by_distance_task(cells_k,:)),1);
    slope_other_animal(k,2) = slope_other_animal_k(1);
    slope_other_animal_k = polyfit(distvec,nanmean(fr_by_distance_free(cells_k,:)),1);
    slope_other_animal(k,1) = slope_other_animal_k(1);
    
    
end

figure(2)
subplot(2,4,8)
plot(slope_grid_animal(:,1),ones(1,7),'*b')
hold on
plot(slope_grid_animal(:,2),2*ones(1,7),'*r')
plot([0 0],[0 3],'--k')
hold off
ylabel('free, task')
xlabel('slope per animal')

figure(3)
subplot(2,4,8)
plot(slope_other_animal(:,1),ones(1,7),'*b')
hold on
plot(slope_other_animal(:,2),2*ones(1,7),'*r')
plot([0 0],[0 3],'--k')
hold off
ylabel('free, task')
xlabel('slope per animal')

