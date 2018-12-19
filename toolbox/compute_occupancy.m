function [occupancy] = compute_occupancy(filename)

[x, y, ~, ~, ~, ~, fr] = load_data_and_preprocess(filename);

[~,occupancy_curve] = compute_2d_tuning_curve(x,y,fr,30,[min(x) min(y)],[max(x) max(y)],1);

occupancy = sum(occupancy_curve(:)>0)/900;