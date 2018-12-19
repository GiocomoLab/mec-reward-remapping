function [best_models,tuning_curves,normal_tuning_curves] = compute_glm(cell_path)

%% load the data
[x, y, hd, ts, speed, spiketrain,~,~] = load_data_and_preprocess(cell_path);
dt = [diff(ts) median(diff(ts))]/1e6; % add a last point so the vectors work out
boxsize = 150;
%% compute zone vector for every time point

% re-compute the head direction for every time point so it's in normal
% units
hd = mod(pi/2 - hd,2*pi);
%% build the input features for position, hd, speed, zone_hd
isnan_x = find(isnan(x));
isnan_y = find(isnan(y));
isnan_h = find(isnan(hd));
isnan_s = find(isnan(speed));
too_fast = find(speed > quantile(speed,0.99));
isnan_all = unique([isnan_x isnan_y isnan_h isnan_s too_fast]);

spiketrain(isnan_all) = [];
x(isnan_all) = [];
y(isnan_all) = [];
speed(isnan_all) = [];
hd(isnan_all) = [];
dt(isnan_all) = [];

%% compute the input matrices for head direction and gaze
fprintf('Making input features\n');
%%%%%%%%% POSITION %%%%%%%%%
s = 0.5;
bin_x = 12;
pos_vec = linspace(0,boxsize,bin_x); pos_vec(1) = -0.01;
[posgrid,pos_x_pts] = spline_2d(x,y,pos_vec,s);
A{1} = posgrid;
ctl_pts_all{1} = pos_vec;

%%%%%%%%% HEAD DIRECTION %%%%%%%%%
bin_h = 10;
hd_vec = linspace(0,2*pi,bin_h+1); hd_vec = hd_vec(1:end-1);
[hdgrid] = spline_1d_circ(hd,hd_vec,s);
A{2} = hdgrid;
ctl_pts_all{2} = hd_vec;

%%%%%%%%% SPEED %%%%%%%%%
spdVec = linspace(0,max(speed),8); spdVec(1) = -0.1;
[speedgrid,~] = spline_1d(speed,spdVec,s);
A{3} = speedgrid;
ctl_pts_all{3} = spdVec;

%% compute the normal tuning curves for comparison
[position_tuning_curve,~] = compute_2d_tuning_curve_dt(x,y,spiketrain,dt,30,[0 0],[150 150]);
[hd_tuning_curve,hd_x,hd_sem] = compute_1d_tuning_curve(hd,spiketrain./dt,18,0,2*pi);
[speed_tuning_curve,spd_x,spd_sem] = compute_1d_tuning_curve(speed,spiketrain./dt,10,0,max(speed));
normal_tuning_curves{1} = position_tuning_curve;
normal_tuning_curves{2} = [hd_tuning_curve  hd_x'];
normal_tuning_curves{3} = [speed_tuning_curve spd_x'];

plotfig = 0;
if plotfig
    figure(1)
    subplot(3,3,1)
    imagesc([0 150],[0 150],position_tuning_curve)
    axis off
    
    subplot(3,3,2)
    errorbar(hd_x,hd_tuning_curve,hd_sem,'k','linewidth',2);
    axis([0 2*pi -inf inf])
    box off
    xlabel('hd')
    ylabel('fr')
    
    subplot(3,3,3)
    errorbar(spd_x,speed_tuning_curve,spd_sem,'k','linewidth',2);
    axis([0 max(speed) -inf inf])
    box off
    xlabel('speed')
    ylabel('fr')
    
end

%% run the forward search

%%%%%%% COMPUTE TEST AND TRAIN INDICES %%%%%
numFolds = 10; % will select model using 10 folds
T = numel(spiketrain); numPts = round(20*1/(median(diff(ts))/1e6)); % 20ish seconds
numSeg = ceil(T/(numFolds*numPts));
oneSeg = ones(numPts,1)*(1:numFolds);
new_ind = repmat(oneSeg(:),numSeg,1);
new_ind = new_ind(1:T);

test_ind = cell(numFolds,1);
train_ind = cell(numFolds,1);
for k = 1:numFolds
    test_ind{k} = find(new_ind == k);
    train_ind{k} = setdiff(1:T,test_ind{k});
end

%%%%%%%% FORWARD SEARCH PROCEDURE %%%%%%%%%
[allModelFits, allModelTrainFits, best_models, bestModelFits, parameters, pvals, final_pval] = forward_search_kfold(A,spiketrain,dt',train_ind,test_ind);

%%%% compute the tuning curves %%%%%
if final_pval < 0.05
    final_param = parameters{end};
    plot_tuning
else
    best_models = [];
    tuning_curves = [];
end

if plotfig
    subplot(3,3,7:9)
    errorbar(1:3,mean(allModelFits{1}),std(allModelFits{1})./sqrt(numFolds),'.k','linewidth',2)
    hold on
    plot([1 3],[0 0],'--r','linewidth',2)
    hold off
    box off
    ylabel('model fit')
    xlabel('models')
    
end

if plotfig
    if windows
        temp = strsplit(cell_path,'\'); cell_name = temp{end};
    else
        temp = strsplit(cell_path,'/'); cell_name = temp{end};
    end
    filename = cell_name;
    fig1 = figure(1);
    set(fig1,'units','pix','pos',[0,0,800,600])
    epsfig = hgexport('factorystyle');
    epsfig.Format = 'jpg';
    if windows
        hgexport(fig1,['C:\Users\khardcas\Dropbox\Rat data\randomfigures\',filename],hgexport('factorystyle'),'Format','jpeg')
    else
        hgexport(fig1,['/Users/kiah/Desktop/rat-task-analysis/model output figures/',filename],hgexport('factorystyle'),'Format','jpeg')
    end
end

close all

return

