function [ratemap1_all,ratemap2_all,occupancy] = do_compute_speed_matched_ratemaps(cell_info_trained,computemaps)

numCells = numel(cell_info_trained);

if computemaps
    ratemap1_all = cell(numCells,1);
    ratemap2_all = cell(numCells,1);
    numiter = 50;
    for k = 1:numCells
        [ratemap1_all{k},ratemap2_all{k},occupancy] = compute_speed_matched_ratemaps(cell_info_trained(k).best_free,cell_info_trained(k).best_task,numiter);
    end
    save('C:\Users\khardcas\Dropbox\Butler Hardcastle rat project\Code for paper\data structures\pos_speed_downsampled_ratemaps.mat','ratemap1_all','ratemap2_all','occupancy');
else
    load_struct = 1;
    if load_struct
        % if the ratemaps already exist, then just load them and reformat to
        % match what's above
        load('all_cells_both_boxes_best_iso_downsampled_speed_pos_only')
        ratemap1_all = cell(numCells,1);
        ratemap2_all = cell(numCells,1);
        for k = 1:numCells
            temp1 = cell_downsampled_trained2(k).ratemaps_free;
            temp2 = cell_downsampled_trained2(k).ratemaps_task;
            ratemap1_all{k} = mean(cat(3,temp1{:}),3);
            ratemap2_all{k} = mean(cat(3,temp2{:}),3);
        end
    else
        load pos_speed_downsampled_ratemaps
    end
end

return