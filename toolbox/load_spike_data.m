function [cellTS] = load_spike_data(filename)

%load video data for session; must have run 'analyze_session_data' first
split = strsplit(filename, '_');
folderpath = strcat('Z:\Users\WButler\D254 Neuralynx\',split(1),'\',split(2),'_',split(3),'\');

%Load spikes 
cell1 = load(strcat(folderpath{1},filename));
cellTS = cell1.cellTS;


return