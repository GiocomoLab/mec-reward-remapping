function [fr_by_distance_free,fr_by_distance_task,distvec,fr_by_distance_free_sem,fr_by_distance_task_sem] ...
    = compute_fr_by_distance(ratemap1_all,ratemap2_all,zone_center)

% ratemap1_all is free
% ratemap2_all is task

[numbins,~] = size(ratemap1_all{1}); numCells = length(ratemap1_all);

% compute the location of every position bin
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
fr_by_distance_free = nan(numCells,numdistbins);
fr_by_distance_task = nan(numCells,numdistbins);
fr_by_distance_free_sem = nan(numCells,numdistbins);
fr_by_distance_task_sem = nan(numCells,numdistbins);

normalize = 1;

for k = 1:numCells
    
    [fr_by_distance_free(k,:),~,fr_by_distance_free_sem(k,:)] = compute_fr_by_distance_1map(ratemap1_all{k},zone_center(k,:),posvec_allx,posvec_ally,numdistbins,minval,maxval,normalize);
    [fr_by_distance_task(k,:),~,fr_by_distance_task_sem(k,:)] = compute_fr_by_distance_1map(ratemap2_all{k},zone_center(k,:),posvec_allx,posvec_ally,numdistbins,minval,maxval,normalize);
    
end

return