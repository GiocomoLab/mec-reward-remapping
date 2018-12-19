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

%% get the list of all possible spatial and free foraging sessions
load All_cells_in_experimental_animals
numCells = numel(Allexperimentalcells);

session_names_free = {};
session_names_task = {};
 
session_dates_free = {};
session_dates_task = {};

for k = 1:numCells
    cellname = Allexperimentalcells(k).cellName;
    cellname_split = strsplit(cellname,'_');
    filedate = strjoin(cellname_split(1:end-2),'_');
    if is_cell_from_trained(filedate)
        sessiontype = Allexperimentalcells(k).sessionType;
        if strcmp(sessiontype,'free')
            session_names_free{end+1} = strjoin(cellname_split(1:end-1),'_');
            session_dates_free{end+1} = filedate;
        elseif strcmp(sessiontype,'spatial')
            session_names_task{end+1} = strjoin(cellname_split(1:end-1),'_');
            session_dates_task{end+1} = filedate;
        end
    end
end
[unique_session_names_task,unique_task_ind] = unique(session_names_task);
[unique_session_names_free,unique_free_ind]  = unique(session_names_free);
unique_free_date = session_dates_free(unique_free_ind);


%% load the data structure
load('behavioral_sessions')

% get the phase II trained dates
animal_list_unique = {'Halo 1','Halo 2','Halo 3','Halo 5','Halo 6','Halo 7','Harlan 18'};
trained_dates = {'2017-10-18','2017-10-18','2017-10-24','2018-01-30',...
    '2018-01-22','2018-02-02','2017-07-04'};

%% get list of sessions before, after training
task_sessions = {};
free_sessions = {};
task_dates = {};
trained = [];
animal_list = {};
for i = 1:length(behavioral_sessions)
    if behavioral_sessions(i).past_phase1 % only use phase 2 and beyond
        task_sessions_k =  strcat(behavioral_sessions(i).animal,'_',behavioral_sessions(i).session_id);
        words = strsplit(task_sessions_k,'_');
        if ~strcmp(words{2},'19')
            task_dates{end+1} = words{3};
            if strcmp(words{1},'halo')
                task_sessions{end+1} = char(strcat('Halo',{' '},words{2},'_',strjoin(words(3:end),'_')));
            elseif strcmp(words{1},'harlan')
                task_sessions{end+1} = char(strcat('Harlan',{' '},words{2},'_',strjoin(words(3:end),'_')));
            end
            animal_k = strsplit(task_sessions{end},'_');
            animal_list{end+1} = animal_k{1};
            
            % find the matching free session
            task_date = strcat(animal_list{end},'_',task_dates{end});
            free_sessions{end+1} = char(unique_session_names_free(strcmp(task_date,unique_free_date)));
            if behavioral_sessions(i).trained
                trained = [trained; 1];
            else
                trained = [trained; 0];
            end
        end
    end
end

unique_session_names_task = {};
unique_session_names_free = {};
for k = 1:numel(task_sessions)
    if trained(k) == 1 && ~isempty(free_sessions{k})
        unique_session_names_task{end+1} = task_sessions{k};
        unique_session_names_free{end+1} = free_sessions{k};
    end
end

numSessions = numel(unique_session_names_free);


%% run the code
theta_power = nan(numSessions,2);
thetapow_speed_info_free = nan(numSessions,3);
thetapow_speed_info_task = nan(numSessions,3);
thetafreq_speed_info_free = nan(numSessions,3);
thetafreq_speed_info_task = nan(numSessions,3);

for k = 1:numSessions
    
    % get the path names
    
    [theta_power(k,:),thetapow_speed_info_free(k,:),thetapow_speed_info_task(k,:),thetafreq_speed_info_free(k,:),thetafreq_speed_info_task(k,:)] ...
        = compare_theta_lfp(unique_session_names_free{k},unique_session_names_task{k});
    
    k
end

save('do_compare_theta_lfp_output_v2.mat')



