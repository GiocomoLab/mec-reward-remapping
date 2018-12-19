function [hd_tuning, mvl, mean_direction,x] = compute_hd_tuning_curve(filename)

[~, ~, angle, ~, ~, ~, fr] = load_data_and_preprocess(filename);

[hd_tuning,~,x] = compute_1d_tuning_curve_periodic(angle,fr,60);

[mvl,mean_direction] = compute_mvl_angle(hd_tuning);