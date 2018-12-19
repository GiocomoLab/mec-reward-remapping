function [new_x] = fill_in_all_nan(x,ts)

non_detect = isnan(x); non_detect = reshape(non_detect,numel(non_detect),1);
new_x = x;
if sum(non_detect) > 0
    new_x(non_detect) = interp1(ts(~non_detect),x(~non_detect),ts(non_detect),'linear');
    if sum(isnan(new_x)) > 0
        last_non_nan = new_x(find(isfinite(new_x),1,'last'));
        new_x(isnan(new_x)) = last_non_nan;
    end
end

return