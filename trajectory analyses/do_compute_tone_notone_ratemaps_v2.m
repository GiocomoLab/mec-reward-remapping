%% This function will run compute_tone_notone_ratemaps over several cells

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

%% get the indices
[grid_cell_ind, other_pos_ind, border_cell_ind, ~, ~, ~, ~] = find_grid_border_cells;
%cells_to_test = [grid_cell_ind; other_pos_ind];
cells_to_test = 1:numCells;

%% for each spatial session, compute the tone times, speed and position matched non-tone times, and compare
tone_ratemap = cell(numel(cells_to_test),1);
notone_ratemap = cell(numel(cells_to_test),1);
tone_speedmap = cell(numel(cells_to_test),1);
notone_speedmap = cell(numel(cells_to_test),1);
mean_fr = nan(numel(cells_to_test),2);
mean_speed = nan(numel(cells_to_test),2);
ratemap = cell(numel(cells_to_test),1);
zone_center1 = nan(numel(cells_to_test),2);

for k = 1:numel(cells_to_test)
    tic
    % get the path name for the spatial session
    cell_name = cell_info_trained(cells_to_test(k)).best_task;
    [tone_ratemap{k},notone_ratemap{k},tone_speedmap{k},notone_speedmap{k},...
        zone_center1(k,:),mean_fr(k,:),mean_speed(k,:)] = compute_tone_notone_ratemaps(cell_name);
    toc
    k
end

save('do_compute_tone_notone_ratemaps_v2_output.mat');
keyboard

%% load the zone center
load all_zones

%% compute avg firing per distance tuning curve

boxSize = 150; numbins = 20;
posvec = linspace(0,boxSize,numbins+1); binw = posvec(2)-posvec(1);
posvec_x = binw/2:binw:boxSize-binw/2;
posvec_y = flip(posvec_x);
[rows,cols] = ind2sub([numbins,numbins],1:numbins^2);
posvec_allx = posvec_x(cols);
posvec_ally = posvec_y(rows);
numdistbins = 10; minval = 20; maxval = 100;
distvec1 = linspace(minval,maxval,numdistbins+1); bw = distvec1(2) - distvec1(1);
distvec = bw/2+minval:bw:maxval-bw/2;

fr_by_distance_tone = nan(numel(cells_to_test),numdistbins);
fr_by_distance_notone = nan(numel(cells_to_test),numdistbins);
fr_by_distance_diff = nan(numel(cells_to_test),numdistbins);

slope_distance_tone = nan(numel(cells_to_test),1);
slope_distance_notone = nan(numel(cells_to_test),1);
slope_distance_diff = nan(numel(cells_to_test),1);
normmap = 0;
for k = 1:numel(cells_to_test)
    
    if isfinite(zone_center1(k,1))
        [fr_by_distance_tone_k,~,] = compute_fr_by_distance_1map(tone_ratemap{k},zone_center(k,:),posvec_allx,posvec_ally,numdistbins,minval,maxval,normmap);
        [fr_by_distance_notone_k] = compute_fr_by_distance_1map(notone_ratemap{k},zone_center(k,:),posvec_allx,posvec_ally,numdistbins,minval,maxval,normmap);
        
        fr_by_distance_tone(k,:) = (fr_by_distance_tone_k - min(fr_by_distance_tone_k))./range(fr_by_distance_tone_k);
        fr_by_distance_notone(k,:) = (fr_by_distance_notone_k - min(fr_by_distance_notone_k))./range(fr_by_distance_notone_k);
        fr_by_distance_diff(k,:) = fr_by_distance_tone(k,:) - fr_by_distance_notone(k,:);
        
        slope_distance_tone_k = polyfit(distvec,fr_by_distance_tone(k,:),1); slope_distance_tone(k) = slope_distance_tone_k(1);
        slope_distance_notone_k = polyfit(distvec,fr_by_distance_notone(k,:),1); slope_distance_notone(k) = slope_distance_notone_k(1);
        slope_distance_diff_k = polyfit(distvec,fr_by_distance_diff(k,:),1); slope_distance_diff(k) = slope_distance_notone_k(1);
    end
    k
end

%% make a plot of each cell 

normmap = 0;
for k = 1:numel(cells_to_test)
    
    if isfinite(zone_center1(k,1))
        
        [fr_by_distance_tone_k] = compute_fr_by_distance_1map(tone_ratemap{k},zone_center(k,:),posvec_allx,posvec_ally,numdistbins,minval,maxval,normmap);
        [fr_by_distance_notone_k] = compute_fr_by_distance_1map(notone_ratemap{k},zone_center(k,:),posvec_allx,posvec_ally,numdistbins,minval,maxval,normmap);
        
        fr_by_distance_tone_norm = (fr_by_distance_tone_k - min(fr_by_distance_tone_k))./range(fr_by_distance_tone_k);
        fr_by_distance_notone_norm = (fr_by_distance_notone_k - min(fr_by_distance_notone_k))./range(fr_by_distance_notone_k);
        
        if slope_distance_tone(k) < 0 && slope_distance_notone(k) < 0
            
            figure(1)
            
            subplot(1,4,1)
            imagesc([0 150],[0 150],ratemapt_all{k})
            axis off
            hold on
            plot([zone_all(k,1) zone_all(k,3)],150-zone_all(k,2)*ones(1,2),'r','linewidth',1.5)
            plot([zone_all(k,1) zone_all(k,3)],150-zone_all(k,4)*ones(1,2),'r','linewidth',1.5)
            plot(zone_all(k,1)*ones(1,2),150-[zone_all(k,2) zone_all(k,4)],'r','linewidth',1.5)
            plot(zone_all(k,3)*ones(1,2),150 - [zone_all(k,2) zone_all(k,4)],'r','linewidth',1.5)
            hold off
            
            
            subplot(1,4,2)
            imagesc([0 150],[0 150],tone_ratemap{k})
            axis off
            hold on
            plot([zone_all(k,1) zone_all(k,3)],150-zone_all(k,2)*ones(1,2),'r','linewidth',1.5)
            plot([zone_all(k,1) zone_all(k,3)],150-zone_all(k,4)*ones(1,2),'r','linewidth',1.5)
            plot(zone_all(k,1)*ones(1,2),150-[zone_all(k,2) zone_all(k,4)],'r','linewidth',1.5)
            plot(zone_all(k,3)*ones(1,2),150 - [zone_all(k,2) zone_all(k,4)],'r','linewidth',1.5)
            hold off
            minmap = max(min(tone_ratemap{k}(:)),min(notone_ratemap{k}(:)));
            maxmap = min(max(tone_ratemap{k}(:)),max(notone_ratemap{k}(:)));
            caxis([minmap maxmap])
            
            
            subplot(1,4,3)
            imagesc([0 150],[0 150],notone_ratemap{k})
            axis off
            hold on
            plot([zone_all(k,1) zone_all(k,3)],150-zone_all(k,2)*ones(1,2),'r','linewidth',1.5)
            plot([zone_all(k,1) zone_all(k,3)],150-zone_all(k,4)*ones(1,2),'r','linewidth',1.5)
            plot(zone_all(k,1)*ones(1,2),150-[zone_all(k,2) zone_all(k,4)],'r','linewidth',1.5)
            plot(zone_all(k,3)*ones(1,2),150 - [zone_all(k,2) zone_all(k,4)],'r','linewidth',1.5)
            hold off
            caxis([minmap maxmap])
            
            subplot(1,4,4)
            plot(distvec,fr_by_distance_notone_norm,'b','linewidth',2)
            hold on
            plot(distvec,fr_by_distance_tone_norm,'r','linewidth',2)
            hold off
            box off
            xlabel('distance from zone center')
            ylabel('norm firing rate')
            
            keyboard
            
        end
        
    end
    k
end

%% plot
grids = grid_cell_ind;
nongrids = other_pos_ind;

figure(1)
subplot(2,3,1)
top = nanmean(fr_by_distance_tone(grids,:)) + nanstd(fr_by_distance_tone(grids,:))./sqrt(numel(grids));
bottom = nanmean(fr_by_distance_tone(grids,:)) - nanstd(fr_by_distance_tone(grids,:))./sqrt(numel(grids));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[1 0 0]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_by_distance_tone(grids,:)),'r','linewidth',2)
box off
alpha(h,0.3)
top = nanmean(fr_by_distance_notone(grids,:)) + nanstd(fr_by_distance_notone(grids,:))./sqrt(numel(grids));
bottom = nanmean(fr_by_distance_notone(grids,:)) - nanstd(fr_by_distance_notone(grids,:))./sqrt(numel(grids));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[0 0 1]);
set(h,'EdgeColor','None')
plot(distvec,nanmean(fr_by_distance_notone(grids,:)),'b','linewidth',2)
box off
alpha(h,0.3)
hold off

subplot(2,3,2)
top = nanmean(fr_by_distance_tone(grids,:)) - nanmean(fr_by_distance_notone(grids,:)) + nanstd(fr_by_distance_tone(grids,:) - fr_by_distance_notone(grids,:))./sqrt(numel(grids));
bottom = nanmean(fr_by_distance_tone(grids,:)) - nanmean(fr_by_distance_notone(grids,:)) - nanstd(fr_by_distance_tone(grids,:) - fr_by_distance_notone(grids,:))./sqrt(numel(grids));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[0 0 0]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_by_distance_tone(grids,:)) - nanmean(fr_by_distance_notone(grids,:)),'k','linewidth',2)
box off
alpha(h,0.3)
plot(distvec,zeros(size(distvec)),'--r','linewidth',2)
hold off

subplot(2,3,3)
[n,x] = hist(slope_distance_tone(grids) - slope_distance_notone(grids),15);
bar(x,n,'k')
box off
hold on
plot([0 0],[0 max(n)],'--r','linewidth',2)
hold off
ylabel('# of cells')
xlabel('slope difference')
title(['p = ',num2str(signrank(slope_distance_tone(grids),slope_distance_notone(grids)))])

subplot(2,3,4)
top = nanmean(fr_by_distance_tone(nongrids,:)) + nanstd(fr_by_distance_tone(nongrids,:))./sqrt(numel(nongrids));
bottom = nanmean(fr_by_distance_tone(nongrids,:)) - nanstd(fr_by_distance_tone(nongrids,:))./sqrt(numel(nongrids));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[1 0 0]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_by_distance_tone(nongrids,:)),'r','linewidth',2)
box off
alpha(h,0.3)
top = nanmean(fr_by_distance_notone(nongrids,:)) + nanstd(fr_by_distance_notone(nongrids,:))./sqrt(numel(nongrids));
bottom = nanmean(fr_by_distance_notone(nongrids,:)) - nanstd(fr_by_distance_notone(nongrids,:))./sqrt(numel(nongrids));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[0 0 1]);
set(h,'EdgeColor','None')
plot(distvec,nanmean(fr_by_distance_notone(nongrids,:)),'b','linewidth',2)
box off
alpha(h,0.3)
hold off

subplot(2,3,5)
top = nanmean(fr_by_distance_tone(nongrids,:)) - nanmean(fr_by_distance_notone(nongrids,:)) + nanstd(fr_by_distance_tone(nongrids,:) - fr_by_distance_notone(nongrids,:))./sqrt(numel(nongrids));
bottom = nanmean(fr_by_distance_tone(nongrids,:)) - nanmean(fr_by_distance_notone(nongrids,:)) - nanstd(fr_by_distance_tone(nongrids,:) - fr_by_distance_notone(nongrids,:))./sqrt(numel(nongrids));
h = fill([distvec fliplr(distvec)],[top fliplr(bottom)],[0 0 0]);
set(h,'EdgeColor','None')
hold on
plot(distvec,nanmean(fr_by_distance_tone(nongrids,:)) - nanmean(fr_by_distance_notone(nongrids,:)),'k','linewidth',2)
box off
alpha(h,0.3)
plot(distvec,zeros(size(distvec)),'--r','linewidth',2)
hold off

subplot(2,3,6)
[n,x] = hist(slope_distance_diff(nongrids),15);
bar(x,n,'k')
box off
hold on
plot([0 0],[0 max(n)],'--r','linewidth',2)
hold off
ylabel('# of cells')
xlabel('slope difference')
title(['p = ',num2str(signrank(slope_distance_tone(nongrids),slope_distance_notone(nongrids)))])

figure(2)
subplot(2,2,1)
[n,x] = hist(slope_distance_tone(grids),20);
bar(x,n,'r')
hold on
plot([0 0],[0 max(n)],'--r','linewidth',2)
hold off
box off
title(['p = ',num2str(signrank(slope_distance_tone(grids)))])
axis([-0.02 0.015 0 inf])

subplot(2,2,2)
[n,x] = hist(slope_distance_notone(grids),20);
bar(x,n,'b')
hold on
plot([0 0],[0 max(n)],'--r','linewidth',2)
hold off
box off
title(['p = ',num2str(signrank(slope_distance_notone(grids)))])
axis([-0.02 0.015 0 inf])

subplot(2,2,3)
[n,x] = hist(slope_distance_tone(nongrids),20);
bar(x,n,'r')
hold on
plot([0 0],[0 max(n)],'--r','linewidth',2)
hold off
box off
title(['p = ',num2str(signrank(slope_distance_tone(nongrids)))])
axis([-0.02 0.015 0 inf])

subplot(2,2,4)
[n,x] = hist(slope_distance_notone(nongrids),20);
bar(x,n,'b')
hold on
plot([0 0],[0 max(n)],'--r','linewidth',2)
hold off
box off
title(['p = ',num2str(signrank(slope_distance_notone(nongrids)))])
axis([-0.02 0.015 0 inf])