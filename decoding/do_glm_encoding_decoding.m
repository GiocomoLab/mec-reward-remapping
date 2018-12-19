%% clear workspace
clear all; close all; clc

%% add paths

addpath('C:\Users\khardcas\Dropbox\Butler Hardcastle rat project\Code for paper\toolbox')
addpath('C:\Users\khardcas\Dropbox\Butler Hardcastle rat project\Code for paper\data structures')
addpath('C:\Users\khardcas\Dropbox\Butler Hardcastle rat project\Code for paper\neuralynx to matlab')


%% load the data structure
load All_cells_in_experimental_animals
num_cells = length(Allexperimentalcells);

%% find pairs of good sessions for each animal
all_free_trained_sessions = {};
all_free_dates = {};
all_free_cells = {};
all_free_cell_id = {};

all_task_trained_sessions = {};
all_task_dates = {};
all_task_cells = {};
all_task_cell_id = {};

for k = 1:num_cells
    
    cell_name = Allexperimentalcells(k).cellName;
    words = strsplit(cell_name,'_');
    session_name = strjoin(words(1:end-1),'_');
    date_k = strjoin(words(1:2),'_');
    trained_status = is_cell_from_trained(cell_name);
    
    if trained_status
        
        if strcmp('free',Allexperimentalcells(k).sessionType)
            all_free_cells{end+1} = cell_name;
            all_free_trained_sessions{end+1} = session_name;
            all_free_cell_id{end+1} = Allexperimentalcells(k).uniqueCellId;
            all_free_dates{end+1} = date_k;
            
        elseif strcmp('spatial',Allexperimentalcells(k).sessionType)
            all_task_cells{end+1} = cell_name;
            all_task_trained_sessions{end+1} = session_name;
            all_task_cell_id{end+1} = Allexperimentalcells(k).uniqueCellId;
            all_task_dates{end+1} = date_k;
        end
    end
end

[unique_sessions_free, unique_sessions_free_ind] = unique(all_free_trained_sessions);
[unique_sessions_task, unique_sessions_task_ind] = unique(all_task_trained_sessions);

unique_dates_free = all_free_dates(unique_sessions_free_ind);
unique_dates_task = all_task_dates(unique_sessions_task_ind);

[overlapping_dates,date_free_ind,date_task_ind] = intersect(all_free_dates,all_task_dates);
free_cells_per_session = cell(numel(overlapping_dates),1);
task_cells_per_session = cell(numel(overlapping_dates),1);
cell_id_per_session = cell(numel(overlapping_dates),1);
num_cells_per_session = nan(numel(overlapping_dates),1);

% find set of paired sessions
paired_session_free = all_free_trained_sessions(date_free_ind);
paired_session_task = all_task_trained_sessions(date_task_ind);

% for each task/free session on the same day, find the set of
% overlapping cells.

for k = 1:numel(overlapping_dates)
    
    free_ind = find(ismember(all_free_dates,overlapping_dates{k}));
    task_ind = find(ismember(all_task_dates,overlapping_dates{k}));
    
    cells_free_k = all_free_cell_id(free_ind);
    cells_task_k = all_task_cell_id(task_ind);
    
    [overlap_cells,free_ind_overlap,task_ind_overlap] = intersect(cells_free_k,cells_task_k);
    
    free_cells_per_session{k} = all_free_cells(free_ind(free_ind_overlap));
    task_cells_per_session{k} = all_task_cells(task_ind(task_ind_overlap));
    
    cell_id_per_session{k} = overlap_cells;
    
    num_cells_per_session(k) = numel(overlap_cells);
    
end

% find sessions with no overlapping cells and remove them
no_cells = find(num_cells_per_session==0);
num_cells_per_session(no_cells) = [];
overlapping_dates(no_cells) = [];
cell_id_per_session(no_cells) = [];
paired_session_free(no_cells) = [];
paired_session_task(no_cells) = [];
free_cells_per_session(no_cells) = [];
task_cells_per_session(no_cells) = [];

% for each pair of sessions, compute the degree of overlap
num_cell_overlap = nan(numel(overlapping_dates),numel(overlapping_dates));

for k = 1:numel(overlapping_dates)
    
    cells_k = cell_id_per_session{k};
    
    for j = 1:numel(overlapping_dates)
        
        cells_j = cell_id_per_session{j};
        
        num_cell_overlap(k,j) = numel(intersect(cells_k,cells_j));
        
    end
end

% pick the best set of sessions: find the one with the most cells, reject
% sessions that have 50% or more of overlap, find the next best session,
% repeat

final_task_session_list = {};
final_free_session_list = {};

final_task_cell_list = {};
final_free_cell_list = {};

temp_free_session_list = paired_session_free;
temp_task_session_list = paired_session_task;
temp_num_cell_overlap = num_cell_overlap;
temp_num_cells_per_session = num_cells_per_session;

temp_free_cell_list = free_cells_per_session;
temp_task_cell_list = task_cells_per_session;

final_num_cells = [];

while ~isempty(temp_task_session_list) && ~isempty(temp_free_session_list)
    
    % of all the possibilities, find the session with the highest number of
    % cells
    [~,highest_session_ind] = max(temp_num_cells_per_session);
    final_free_session_list{end+1} = temp_free_session_list{highest_session_ind};
    final_task_session_list{end+1} = temp_task_session_list{highest_session_ind};
    
    final_free_cell_list{end+1} = temp_free_cell_list{highest_session_ind};
    final_task_cell_list{end+1} = temp_task_cell_list{highest_session_ind};
    final_num_cells = [final_num_cells; numel(temp_free_cell_list{highest_session_ind})];
    
    % find the sessions that overlap too much and delete them
    overlap = temp_num_cell_overlap(highest_session_ind,:);
    overlap_proportion = overlap'./temp_num_cells_per_session;
    too_much_overlap = find(overlap_proportion > 0.5);
    
    temp_num_cells_per_session(too_much_overlap) = [];
    temp_num_cell_overlap(too_much_overlap,:) = [];
    temp_num_cell_overlap(:,too_much_overlap) = [];
    temp_free_session_list(too_much_overlap) = [];
    temp_task_session_list(too_much_overlap) = [];
    temp_free_cell_list(too_much_overlap) = [];
    temp_task_cell_list(too_much_overlap) = [];
end

%% for each pair of sessions, decode!
mincellnum = 5;
final_free_session_list = final_free_session_list(final_num_cells >= mincellnum);
final_free_cell_list = final_free_cell_list(final_num_cells >= mincellnum);
final_task_session_list = final_task_session_list(final_num_cells >= mincellnum);
final_task_cell_list = final_task_cell_list(final_num_cells >= mincellnum);
final_num_cells = final_num_cells(final_num_cells >=mincellnum);

num_sessions = numel(final_free_session_list);

free_decode_info = cell(num_sessions,1);
task_decode_info = cell(num_sessions,1);

encoding_cells_free = cell(num_sessions,1);
encoding_cells_task = cell(num_sessions,1);

for k = 1:num_sessions
    
    [free_decode_info{k},task_decode_info{k},encoding_cells_free{k},encoding_cells_task{k}] = glm_encoding_decoding(final_free_session_list{k},final_free_cell_list{k},final_task_session_list{k},final_task_cell_list{k});
    save('decoding_info_v2.mat')
end
