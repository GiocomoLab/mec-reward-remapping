clear all; clc; close all

load single_cell_theta

figure(1)
subplot(1,2,1)
plot(mvl_theta(:,1),mvl_theta(:,2),'ok')
hold on
plot([0 0.8],[0 0.8],'--r','linewidth',2)
hold off
box off
set(gca,'fontsize',20)
ylabel('theta locking (task)')
xlabel('theta locking (free)')
title(['p = ',num2str(signrank(mvl_theta(:,1),mvl_theta(:,2)))])

subplot(1,2,2)
either = union(find(mvl_theta(:,1) > 0.3),find(mvl_theta(:,2) > 0.3));
both = intersect(find(mvl_theta(:,1) > 0.3),find(mvl_theta(:,2) > 0.3));
plot(mvl_theta(either,1),mvl_theta(either,2),'ok')
hold on
plot([0 0.8],[0 0.8],'--r','linewidth',2)
hold off
box off
set(gca,'fontsize',20)
ylabel('theta locking (task)')
xlabel('theta locking (free)')
title(['p = ',num2str(signrank(mvl_theta(either,1),mvl_theta(either,2)))])


figure(2)
angle_diff = pref_angle(both,1) - pref_angle(both,2);
angle_diff(angle_diff > pi) = 2*pi - angle_diff(angle_diff > pi);
angle_diff(angle_diff < -pi) = -(2*pi + angle_diff(angle_diff <- pi));
[n,x] = hist(angle_diff,30);
bar(x,n,'k')
hold on
plot([0 0],[0 max(n)],'--r','linewidth',2)
hold off
box off
set(gca,'fontsize',20)
ylabel('# of unique cells')
xlabel('preferred angle')
title(['p = ',num2str(signrank(angle_diff))])
axis([-pi pi 0 max(n)])


figure(3)
subplot(1,2,1)
either = union(find(theta_index(:,1) > 5),find(theta_index(:,2) > 5));
plot(theta_index(:,1),theta_index(:,2),'ok')
hold on
plot([0 20],[0 20],'--r','linewidth',2)
hold off
box off
set(gca,'fontsize',20)
ylabel('theta index (task)')
xlabel('theta locking (free)')
title(['p = ',num2str(signrank(theta_index(:,1),theta_index(:,2)))])

subplot(1,2,2)
plot(theta_index(either,1),theta_index(either,2),'ok')
hold on
plot([0 20],[0 20],'--r','linewidth',2)
hold off
box off
set(gca,'fontsize',20)
ylabel('theta locking (task)')
xlabel('theta locking (free)')
title(['p = ',num2str(signrank(theta_index(either,1),theta_index(either,2)))])


figure(4)
subplot(1,2,1)
plot(theta_skipping(:,1),theta_skipping(:,2),'ok')
hold on
plot([-1 1],[-1 1],'--r','linewidth',2)
hold off
box off
set(gca,'fontsize',20)
ylabel('theta skipping (task)')
xlabel('theta skipping (free)')
title(['p = ',num2str(signrank(theta_skipping(:,1),theta_skipping(:,2)))])

subplot(1,2,2)
plot(theta_skipping(either,1),theta_skipping(either,2),'ok')
hold on
plot([0 1],[0 1],'--r','linewidth',2)
hold off
box off
set(gca,'fontsize',20)
ylabel('theta skipping (task)')
xlabel('theta skipping (free)')
title(['p = ',num2str(signrank(theta_skipping(either,1),theta_skipping(either,2)))])