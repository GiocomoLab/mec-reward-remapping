function [interx] = fill_in_consec_nan(x,Timestamps,maxconsec)

non_detect = find(isnan(x)); non_detect = reshape(non_detect,numel(non_detect),1);
interx = x;

if numel(non_detect) > 0
    k = [true;diff(non_detect)~=1];
    s=cumsum(k);
    xx = histc(s,1:s(end));
    xx(:,2)=1:length(xx);
    goodseqs = xx(xx(:,1)<=maxconsec,2);
    
    %lists indexes of "good" non-detects (that are in sets of <= maxconsec)
    goodnons = non_detect(ismember(s, goodseqs));
    list = 1:length(x);
    idx = ismember(list, goodnons);
    
    %interpolates x non_detects linearly
    interx(idx) = interp1(Timestamps(~idx),interx(~idx),Timestamps(idx),'linear');
    
end

return