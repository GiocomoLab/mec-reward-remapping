function [slopes, intercepts] = do_fr_by_distance(indices, task_or_free) 
slopes = nan(numel(indices),1);
intercepts = nan(numel(indices),1);

for i = 1:numel(indices)
    [all_dist, all_fr_t, all_fr_f] = fr_by_distance(indices(i), 0);
    if strcmp(task_or_free,'task')
        m = fitlm(all_dist,all_fr_t);
        slope = m.Coefficients.Estimate(2);
        slopes(i) = slope;
        intercept = m.Coefficients.Estimate(1);
        intercepts(i) = intercept;
    elseif strcmp(task_or_free,'free')
        m = fitlm(all_dist,all_fr_f);
        slope = m.Coefficients.Estimate(2);
        slopes(i) = slope;
        intercept = m.Coefficients.Estimate(1);
        intercepts(i) = intercept;
    elseif strcmp(task_or_free,'task minus free')
        m = fitlm(all_dist,[all_fr_t-all_fr_f]);
        slope = m.Coefficients.Estimate(2);
        slopes(i) = slope;
        intercept = m.Coefficients.Estimate(1);
        intercepts(i) = intercept;
    end
end