%% this script will load the free foraging and spatial task sessions and compare the activity


%% clear workspace

clear all; close all; clc

%% load the structure that tells me which cells to run
addpath('C:\Users\khardcas\Dropbox\Butler Hardcastle rat project\Code for paper\data structures')
load all_control_cells_both_boxes_best_iso_with_scores_model

%% for each cell, load the data if it doesn't already exist on my dropbox

numCells = length(control_info);
best_models_free = cell(numCells,1);
best_models_spatial = cell(numCells,1);
tuning_curves_free = cell(numCells,1);
tuning_curves_spatial = cell(numCells,1);
normal_tuning_free = cell(numCells,1);
normal_tuning_spatial = cell(numCells,1);
for k = 1:numCells
    
    [best_models_free{k},tuning_curves_free{k},normal_tuning_free{k}] = compute_glm(control_info(k).best_free);
    [best_models_spatial{k},tuning_curves_spatial{k},normal_tuning_spatial{k}] = compute_glm(control_info(k).best_task);
    
    k
end

save('glm_output_control.mat')

keyboard

%% look at variables encoded
reward_zone_cells_free = nan(numCells,1);
pos_cells_free = nan(numCells,1);
hd_cells_free = nan(numCells,1);
spd_cells_free = nan(numCells,1);
reward_zone_cells_spatial = nan(numCells,1);
pos_cells_spatial = nan(numCells,1);
hd_cells_spatial = nan(numCells,1);
spd_cells_spatial = nan(numCells,1);
encoding_cells_free = zeros(numCells,1);
encoding_cells_spatial = zeros(numCells,1);
for k = 1:numCells
    pos_cells_free(k) = ismember(1,best_models_free{k});
    pos_cells_spatial(k) = ismember(1,best_models_spatial{k});
    
    hd_cells_free(k) = ismember(2,best_models_free{k});
    hd_cells_spatial(k) = ismember(2,best_models_spatial{k});
    
    spd_cells_free(k) = ismember(3,best_models_free{k});
    spd_cells_spatial(k) = ismember(3,best_models_spatial{k});
    
    reward_zone_cells_free(k) = ismember(4,best_models_free{k});
    reward_zone_cells_spatial(k) = ismember(4,best_models_spatial{k});
    
    encoding_cells_free(k) = ~isempty(best_models_free{k});
    encoding_cells_spatial(k) = ~isempty(best_models_spatial{k});
end



%% load this info to the structure

load glm_output_control
for i = 1:numCells
control_info(i).model_free = best_models_free{i};
control_info(i).model_task = best_models_spatial{i};
end









