function [x, y, angle, ts, speed] = load_session_data_and_preprocess(filename)

%load video data for session; must have run 'analyze_session_data' first
split = strsplit(filename, '_');
folderpath = strcat('Z:\Users\WButler\D254 Neuralynx\',split(1),'\',split(2),'_',split(3),'\');
matlabpath = strcat(folderpath,split(1),'_',split(2),'_',split(3));
file1 = load(matlabpath{1});
x = file1.x; y = file1.y; angle = file1.angle; ts = file1.Timestamps;
pos_nan = isnan(x);

% Compute the running speed
speed = nan(size(x));
for i = 1:numel(x)
    if i == 1
        speed(i) = sqrt((x(i) - x(i+1)).^2 + (y(i) - y(i+1)).^2)/(ts(i+1) - ts(i))*1e6;
    elseif i == numel(x)
        speed(i) = sqrt((x(i-1) - x(i)).^2 + (y(i-1) - y(i)).^2)/(ts(i) - ts(i-1))*1e6;
    else
        speed(i) = sqrt((x(i-1) - x(i+1)).^2 + (y(i-1) - y(i+1)).^2)/(ts(i+1) - ts(i-1))*1e6;
    end
end
too_fast = 200; too_fast_ind = speed > too_fast;
x(too_fast_ind) = interp1(ts(~too_fast_ind),x(~too_fast_ind),ts(too_fast_ind),'linear');
y(too_fast_ind) = interp1(ts(~too_fast_ind),y(~too_fast_ind),ts(too_fast_ind),'linear');
speed(too_fast_ind) = NaN;
[speed] = fill_in_all_nan(speed,ts);
smooth_window = 5;
speed = conv(speed,gausswin(smooth_window)/sum(gausswin(smooth_window)),'same');
speed(pos_nan) = NaN;

% take away stationary times
too_slow = find(speed < 2);
x(too_slow) = NaN; y(too_slow) = NaN; angle(too_slow) = NaN; speed(too_slow) = NaN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the length of time in each bin
Fs = 1/((ts(2)-ts(1))/1e6);
diffTS = diff(ts);    TimeStamps_mid = ts(1:end-1)+diffTS/2;
diff_TimeStamps_mid = diff(TimeStamps_mid);
diff_TimeStamps_mid = [diffTS(1) diff_TimeStamps_mid diffTS(end)];

return