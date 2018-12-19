%% this script will load the free foraging and spatial task sessions and compare the activity


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

%% load the structure

load all_cells_both_boxes_best_iso_with_scores_model
numCells = numel(cell_info_trained);

mvl_theta = nan(numCells,2);
pref_angle = nan(numCells,2);
theta_index = nan(numCells,2);
theta_skipping = nan(numCells,2);

for k = 1:numCells
    tic
    [mvl_theta(k,:),pref_angle(k,:),theta_index(k,:),theta_skipping(k,:)] = compare_theta_single_cell(cell_info_trained(k).best_free,cell_info_trained(k).best_task);
    toc
    k
end
save('single_cell_theta.mat')





