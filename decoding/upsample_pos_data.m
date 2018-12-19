function [x_us,y_us,ts_us] = upsample_pos_data(x,y,ts,dt)


% fill in nans
x = fill_in_all_nan(x,ts);
y = fill_in_all_nan(y,ts);

% upsample data
dt_micro = dt*1e6;
ts_us = ts(1):dt_micro:ts(numel(ts)); % upsampled timestamps
ts_isi = diff(ts); median_isi = median(ts_isi);
break_isi = find(ts_isi > median_isi*2);
for k = 1:numel(break_isi)
    start_break_ts = ts(break_isi(k)); end_break_ts = ts(break_isi(k)+1);
    [~,start_break_ts_us_ind] = min(abs(start_break_ts - ts_us));
    [~,end_break_ts_us_ind] = min(abs(end_break_ts - ts_us));
    ts_us(start_break_ts_us_ind:end_break_ts_us_ind) = [];
end
x_us = interp1(ts,x,ts_us);
y_us = interp1(ts,y,ts_us);

return