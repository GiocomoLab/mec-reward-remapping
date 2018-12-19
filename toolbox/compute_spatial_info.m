function [spatial_info] = compute_spatial_info(tuning,occupancy)

mean_rate = sum(sum(tuning.*occupancy));

% make sure that the input is a column vector
tuning_vec = reshape(tuning,numel(tuning),1);
occupancy_vec = reshape(occupancy,numel(occupancy),1);

spatial_info = 0;
for k = 1:numel(tuning_vec)
    if occupancy_vec(k) > 0 && tuning_vec(k) > 0
        spatial_info = spatial_info + (occupancy_vec(k) * (tuning_vec(k)/mean_rate) * log2( tuning_vec(k) / mean_rate )); 
        %spatial_info = spatial_info + (pos_occupancy_vec(k) * (pos_tuning_vec(k)) * log2( pos_tuning_vec(k) / mean_rate )); 
    end
end

return