function [autocorrelogram, tuning_curve] = compute_autocorrelogram(filename)

[x, y, ~, ~, ~, ~, fr] = load_data_and_preprocess(filename);

[tuning_curve,~] = compute_2d_tuning_curve(x,y,fr,30,[min(x) min(y)],[max(x) max(y)],1);

autocorrelogram = correlation(tuning_curve,tuning_curve);

return