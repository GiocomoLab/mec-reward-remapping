function [downsampleindices1, downsampleindices2] = downsample_match_position(x1,y1,x2,y2)
% matches two sessions based on speed x position
numposBin = 15; %uses the same number of bins as Hardcastle et al (2017)

xAxis = linspace(0,150,numposBin+1);
yAxis = linspace(0,150,numposBin+1);

%initialize matrices
map1 = zeros(numposBin, numposBin);
map2 = zeros(numposBin, numposBin);

downsampleindices1 = [];
downsampleindices2 = [];
%% Determine number of observations in each bin in session1
for i = 1:numposBin
    for j = 1:numposBin
        ind1 = find(x1 >= xAxis(i) & x1 < xAxis(i+1) & y1 >= yAxis(j) & y1 < yAxis(j+1));
        map1(i,j)= numel(ind1);
        
        %same thing for map 2
        ind2 = find(x2 >= xAxis(i) & x2 < xAxis(i+1) & y2 >= yAxis(j) & y2 < yAxis(j+1));
        map2(i,j) = numel(ind2);
        
        % Number of observations to keep in each map
        numToKeep = min(numel(ind1),numel(ind2));
        
        downsampleindices1 = [downsampleindices1 datasample(ind1, numToKeep, 'Replace', false)];
        downsampleindices2 = [downsampleindices2 datasample(ind2, numToKeep, 'Replace', false)];
    end
end


downsampleindices1 = sort(downsampleindices1);
downsampleindices2 = sort(downsampleindices2);
return