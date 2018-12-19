%% clear workspace and load data
close all
clc
clear all
load do_compare_theta_lfp_output_v2

%% make the figures
figure(1)
plot(theta_power(:,1),theta_power(:,2),'ok')
hold on
plot([0 0.7],[0 0.7],'--r','linewidth',2)
hold off
box off
title(['p = ',num2str(signrank(theta_power(:,1),theta_power(:,2)))])
set(gca,'fontsize',15)
xlabel('theta power in free')
ylabel('theta power in task')

figure(2)
subplot(1,3,1)
plot(thetapow_speed_info_free(:,1),thetapow_speed_info_task(:,1),'ok')
hold on
plot([0 0.7],[0 0.7],'--r','linewidth',2)
hold off
box off
title(['p = ',num2str(signrank(thetapow_speed_info_free(:,1),thetapow_speed_info_task(:,1)))])
set(gca,'fontsize',15)
xlabel('theta-speed corr in free')
ylabel('theta-speed corr in task')

subplot(1,3,2)
plot(thetapow_speed_info_free(:,2),thetapow_speed_info_task(:,2),'ok')
hold on
plot([0 0.01],[0 0.01],'--r','linewidth',2)
hold off
box off
title(['p = ',num2str(signrank(thetapow_speed_info_free(:,2),thetapow_speed_info_task(:,2)))])
set(gca,'fontsize',15)
xlabel('theta-speed slope in free')
ylabel('theta-speed slope in task')

subplot(1,3,3)
plot(thetapow_speed_info_free(:,3),thetapow_speed_info_task(:,3),'ok')
hold on
plot([0 0.6],[0 0.6],'--r','linewidth',2)
hold off
box off
title(['p = ',num2str(signrank(thetapow_speed_info_free(:,3),thetapow_speed_info_task(:,3)))])
set(gca,'fontsize',15)
xlabel('theta-speed intercept in free')
ylabel('theta-speed intercept in task')

figure(3)
subplot(1,3,1)
plot(thetafreq_speed_info_free(:,1),thetafreq_speed_info_task(:,1),'ok')
hold on
plot([0 0.7],[0 0.7],'--r','linewidth',2)
hold off
box off
title(['p = ',num2str(signrank(thetafreq_speed_info_free(:,1),thetafreq_speed_info_task(:,1)))])
set(gca,'fontsize',15)
xlabel('theta freq-speed corr in free')
ylabel('theta freq-speed corr in task')

subplot(1,3,2)
plot(thetafreq_speed_info_free(:,2),thetafreq_speed_info_task(:,2),'ok')
hold on
plot([0 0.06],[0 0.06],'--r','linewidth',2)
hold off
box off
title(['p = ',num2str(signrank(thetafreq_speed_info_free(:,2),thetafreq_speed_info_task(:,2)))])
set(gca,'fontsize',15)
xlabel('theta freq-speed slope in free')
ylabel('theta freq-speed slope in task')

subplot(1,3,3)
plot(thetafreq_speed_info_free(:,3),thetafreq_speed_info_task(:,3),'ok')
hold on
plot([7 8],[7 8],'--r','linewidth',2)
hold off
box off
title(['p = ',num2str(signrank(thetafreq_speed_info_free(:,3),thetafreq_speed_info_task(:,3)))])
set(gca,'fontsize',15)
xlabel('theta freq-speed intercept in free')
ylabel('theta freq-speed intercept in task')

