function [zone_all,zone_center] = do_compute_all_zones(cell_info_trained,computezones)

numCells = numel(cell_info_trained);

if computezones
    zone_all = nan(numCells,4);
    zone_center = nan(numCells,2);
    animalnum_all = nan(numCells,1);
    for k = 1:numCells
        animalnum = strsplit(cell_info_trained(k).best_task,'_');
        animalnum = strsplit(animalnum{1}); animalnum = str2double(animalnum{2});
        animalnum_all(k) = animalnum;
        [zone_all(k,:),~,~,~] = compute_zone(cell_info_trained(k).best_task);
        zone_center(k,:) = [(zone_all(k,3)-zone_all(k,1))/2+zone_all(k,1) (zone_all(k,4)-zone_all(k,2))/2+zone_all(k,2)];
        k
    end
    unique_animal_name = unique(animalnum_all);
    for k = 1:numel(unique_animal_name)
        animal_ind = find(animalnum_all == unique_animal_name(k));
        zone_all(animal_ind,:) = repmat(nanmean(zone_all(animal_ind,:)),numel(animal_ind),1);
        zone_center(animal_ind,:) = repmat(nanmean(zone_center(animal_ind,:)),numel(animal_ind),1);
    end
    save('C:\Users\khardcas\Dropbox\Butler Hardcastle rat project\Code for paper\data structures\all_zones.mat','zone_center','zone_all');
else
    load all_zones
end


return