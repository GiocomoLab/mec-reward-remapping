function [x_ind1,y_ind1] = match_1d_dist(x,y,vector)

x_ind1 = [];
y_ind1 = [];

x_vector = zeros(numel(vector)-1,1);
y_vector = zeros(numel(vector)-1,1);
for k = 1:numel(vector)-1
    
    x_ind = find(x >= vector(k) & x < vector(k+1));
    y_ind = find(y >= vector(k) & y < vector(k+1));
    
    if k == numel(vector)-1
        x_ind = find(x >= vector(k) & x <= vector(k+1));
        y_ind = find(y >= vector(k) & y <= vector(k+1));
    end
    
    if numel(x_ind) > 0 && numel(y_ind) > 0
        if numel(x_ind) < numel(y_ind)
            
            % keep all x indices
            x_ind1 = [x_ind1; x_ind'];
            
            % keep all y indices
            y_ind_samp = randperm(numel(y_ind));
            y_ind_samp = y_ind(y_ind_samp(1:numel(x_ind)));
            y_ind1 = [y_ind1; y_ind_samp'];
            
            x_vector(k) = numel(x_ind);
            y_vector(k) = numel(y_ind_samp);
         
            
        elseif numel(x_ind) > numel(y_ind)
            
            % keep all y indices
            y_ind1 = [y_ind1; y_ind'];
            
            % keep all y indices
            x_ind_samp = randperm(numel(x_ind));
            x_ind_samp = x_ind(x_ind_samp(1:numel(y_ind)));
            x_ind1 = [x_ind1; x_ind_samp'];
            
            x_vector(k) = numel(x_ind_samp);
            y_vector(k) = numel(y_ind);
            
        else
            
            y_ind1 = [y_ind1; y_ind'];
            x_ind1 = [x_ind1; x_ind'];
            
            x_vector(k) = numel(x_ind);
            y_vector(k) = numel(y_ind);
            
        end
    end
    
end

return