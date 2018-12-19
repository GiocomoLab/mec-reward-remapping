function [tuning_curve,hdOccupancyMap, hdVec] = compute_1d_tuning_curve_periodic(hd,fr,numbin)

hdVec = linspace(0,2*pi,numbin+1);
hdRateMap = nan(numbin,1);
hdOccupancyMap = nan(numbin,1);

% compute mean fr for each direction bin
for n = 1:numbin
    start = hdVec(n); stop = hdVec(n+1);
    hdRateMap(n) = nanmean(fr(hd >= start & hd < stop));
    hdOccupancyMap(n) = numel(fr(hd >= start & hd <= stop));
    if n == numbin
        hdRateMap(n) = nanmean(fr(hd >= start & hd <= stop));
        hdOccupancyMap(n) = numel(fr(hd >= start & hd <= stop));
    end
end

% expanded version
pad = 10;
hdRateMap_big = [hdRateMap(end-pad+1:end); hdRateMap; hdRateMap(1:pad)];

% smooth with a Gaussian filter
filter = gaussmf(-2:2,[2 0]); filter = filter/sum(filter);
hdRateMap_big_smooth = conv(hdRateMap_big,filter,'same');
tuning_curve = hdRateMap_big_smooth(pad+1:pad+numbin);

return