function [lfp,lfp_ts] = load_csc_data_and_preprocess(filename,tetrode)

%% load the csc data for that session
csc_file_name = ['CSC',tetrode];

split = strsplit(filename, '_');
folderpath = strcat('Z:\Users\WButler\D254 Neuralynx\',split(1),'\',split(2),'_',split(3),'\');
matlabpath = strcat(folderpath,csc_file_name);

[cscTimestamps, ~, ~,~, Samples, ~]  = Nlx2MatCSC([char(matlabpath),'.ncs'],[1 1 1 1 1], 1, 1, [] );

%% do some processing.. reshape it and stuff
% reshape CSC
lfp = Samples(:);

% redo the timestamps
lfp_ts = nan(size(lfp));
for j = 1:numel(cscTimestamps)
    if j < numel(cscTimestamps)
        csc_dt = (cscTimestamps(j+1) - cscTimestamps(j))/512;
    else
        csc_dt = (cscTimestamps(j) - cscTimestamps(j-1))/512;
    end
    lfp_ts(512*(j-1)+1:512*j) = cscTimestamps(j):csc_dt:cscTimestamps(j)+csc_dt*511;
end


return