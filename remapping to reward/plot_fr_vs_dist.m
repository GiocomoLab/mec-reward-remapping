figure(1)
subplot(2,3,1)
fr_temp = fr_by_distance_free(grid_cell_ind,:); fr_temp = fr_temp(:);
dist_temp = distvec_free(grid_cell_ind,:); dist_temp = dist_temp(:);
[y,x,sem,~] = compute_1d_tuning_curve(dist_temp,fr_temp,numdistbins,nanmin(dist_temp),nanmax(dist_temp));
top = y+sem; bottom = y - sem;
h = fill([x fliplr(x)],[top' flip(bottom')],[0 0 0]);
set(h,'EdgeColor','None')
hold on
plot(x,y,'k','linewidth',2)
alpha(h,0.3)
hold off
box off
xlabel('distance from most occupied spot')
ylabel('normalized firing rate')
box off

subplot(2,3,4)
[n,x] = hist(slope_free(grid_cell_ind),20);
plot(x,n,'k','linewidth',2)
hold on
plot([0 0],[0 max(n)],'--r','linewidth',2)
hold off
box off
xlabel('slope value')
ylabel('# of cells')
p1 = signrank(slope_free(grid_cell_ind));
title(['p_f = ',num2str(p1)])

% border cells
subplot(2,3,2)
fr_temp = fr_by_distance_free(border_cell_ind,:); fr_temp = fr_temp(:);
dist_temp = distvec_free(border_cell_ind,:); dist_temp = dist_temp(:);
[y,x,sem,~] = compute_1d_tuning_curve(dist_temp,fr_temp,numdistbins,nanmin(dist_temp),nanmax(dist_temp));
top = y+sem; bottom = y - sem;
h = fill([x fliplr(x)],[top' flip(bottom')],[0 0 0]);
set(h,'EdgeColor','None')
hold on
plot(x,y,'k','linewidth',2)
alpha(h,0.3)
hold off
box off
xlabel('distance from most occupied spot')
ylabel('normalized firing rate')
box off

subplot(2,3,5)
[n,x] = hist(slope_free(border_cell_ind),20);
plot(x,n,'k','linewidth',2)
hold on
plot([0 0],[0 max(n)],'--r','linewidth',2)
hold off
box off
xlabel('slope value')
ylabel('# of cells')
p1 = signrank(slope_free(border_cell_ind));
title(['p_f = ',num2str(p1)])

% non-grid non-border cells
subplot(2,3,3)
fr_temp = fr_by_distance_free(other_pos_ind,:); fr_temp = fr_temp(:);
dist_temp = distvec_free(other_pos_ind,:); dist_temp = dist_temp(:);
[y,x,sem,~] = compute_1d_tuning_curve(dist_temp,fr_temp,numdistbins,nanmin(dist_temp),nanmax(dist_temp));
top = y+sem; bottom = y - sem;
h = fill([x fliplr(x)],[top' flip(bottom')],[0 0 0]);
set(h,'EdgeColor','None')
hold on
plot(x,y,'k','linewidth',2)
alpha(h,0.3)
hold off
box off
xlabel('distance from most occupied spot')
ylabel('normalized firing rate')
box off

subplot(2,3,6)
[n,x] = hist(slope_free(other_pos_ind),20);
plot(x,n,'k','linewidth',2)
hold on
plot([0 0],[0 max(n)],'--r','linewidth',2)
hold off
box off
xlabel('slope value')
ylabel('# of cells')
p1 = signrank(slope_free(other_pos_ind));
title(['p_f = ',num2str(p1)])