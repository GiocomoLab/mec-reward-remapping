function [x_real,x_decode,y_real,y_decode] = do_decoding(x,y,S,shgrid,decoding_start,decoding_stop,param,pos_ind,sh_ind,posVec,s)

% this is going to be memory-less decoding, which will be easier to do

num_cells = size(S,1);

%% compute the position occupancy
posvec = 0:3:150;
pos_occupancy = zeros(numel(posvec),numel(posvec));
for m = 1:numel(posvec)-1
    if m < numel(posvec)-1
        x_index = find(x >= posvec(m) & x < posvec(m+1));
    else
        x_index = find(x >= posvec(m) & x < posvec(m+1));
    end
    for n = 1:numel(posvec)-1
        if n < numel(posvec)-1
            y_index = find(y >= posvec(n) & y < posvec(n+1));
        else
            y_index = find(y >= posvec(n) & y <= posvec(n+1));
        end
        pos_occupancy(n,m) = numel(intersect(x_index,y_index));
    end
end
pos_occupancy_vec = pos_occupancy(:);
pos_occupancy_prob = pos_occupancy_vec./sum(pos_occupancy_vec);

%% make pos_states (vector of all possible positions)
x_all = repmat(posvec,numel(posvec),1); x_all_vec = x_all(:);
y_all = repmat(posvec',1,numel(posvec)); y_all_vec = y_all(:);
states = [x_all_vec'; y_all_vec'];
num_poss_pos = numel(x_all);

%% compute the position rate (estimated rate for every position)
[possib_posgrid,~] = spline_2d(x_all_vec,y_all_vec,posVec,s);
possib_pos_rate = exp(possib_posgrid*param(:,pos_ind)'); % spikes/bin

%% do the decoding

decode_num = numel(decoding_start);
x_decode = nan(decode_num,1);
y_decode = nan(decode_num,1);
x_real = nan(decode_num,1);
y_real = nan(decode_num,1);

for m = 1:decode_num
    
    decode_ind_m = decoding_start(m):decoding_stop(m)-1;
    
    num_int = decoding_stop(m) - decoding_start(m);
    
    spikes_k = sum(S(:,decode_ind_m),2); %spikes/numint
    
    sh_input = nan(num_int,num_cells);
    for j = 1:num_cells
        sh_grid_j = shgrid{j};
        sh_input(:,j) = sh_grid_j(decode_ind_m,:)*(param(j,sh_ind)');
    end
    
    possib_rate_base = exp(param(:,1))';
    possib_rate_sh = sum(exp(sh_input),1); % spikes/numint
    possib_rate_nop = repmat(possib_rate_base.*possib_rate_sh,num_poss_pos,1);
    
    possib_rate = possib_rate_nop.*(possib_pos_rate);
    
    possib_f = sum(-possib_rate+repmat(spikes_k',num_poss_pos,1).*log(possib_rate),2);
    likelihood = exp(possib_f);
    
    % compute the posterior using just the data from previous time
    % bin
    posterior_notnorm = likelihood.*pos_occupancy_prob;
    posterior = posterior_notnorm./sum(posterior_notnorm);
    
    [~,possib_index] = max(posterior);
    x_decode(m) = x_all_vec(possib_index);
    y_decode(m) = y_all_vec(possib_index);
    
    x_real(m) = nanmean(x(decode_ind_m));
    y_real(m) = nanmean(y(decode_ind_m));
    
end

return
