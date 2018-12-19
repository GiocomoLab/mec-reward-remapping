%% estimate the state transition matrices of the behavior
% this needs to be on the right timescale...

numint_dt = round(num_int/dt);

x_numint = x_us(1:numint_dt:length(x_us));
y_numint = y_us(1:numint_dt:length(y_us));

x_k1 = [x_numint(2:end); y_numint(2:end)];
x_k = [x_numint(1:end-1); y_numint(1:end-1)];
F = x_k' \ x_k1';
state_err = x_k1 - F*x_k;
Q = cov(state_err');
inv_Q = inv(Q);
det_Q = det(Q);

%% compute a parameter matrix for all the "encoding cells"

base = 1;
pnum = (bin_p+2)^2;
shnum = numel(spikeHistoryVec)+2;
param_num = base+pnum+shnum;
pos_ind = 1:pnum;
sh_ind = pnum+1:pnum+shnum;

all_param = zeros(numel(encoding_cells),param_num);
for k = 1:numCell
    
    param_k = nanmean(output(k).parameters);
    models = output(k).bestModels;
    
    all_param(k,1) = param_k(1);
    param_k = param_k(2:end);
    indices = 0;
    if ismember(1,models)
        indices = 1:pnum;
        all_param(k,pos_ind+1) = param_k(indices);
    end
    
    if ismember(2,models)
        indices = [indices(end)+1:indices(end)+shnum];
        all_param(k,sh_ind+1) = param_k(indices);
    end
    
end

%% first, compute the position occupancy
posvec = 0:3:150;
pos_occupancy = zeros(numel(posvec),numel(posvec));
for m = 1:numel(posvec)-1
    if m < numel(posvec)-1
        x_index = find(x_us >= posvec(m) & x_us < posvec(m+1));
    else
        x_index = find(x_us >= posvec(m) & x_us < posvec(m+1));
    end
    for n = 1:numel(posvec)-1
        if n < numel(posvec)-1
            y_index = find(y_us >= posvec(n) & y_us < posvec(n+1));
        else
            y_index = find(y_us >= posvec(n) & y_us <= posvec(n+1));
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
posVec = linspace(0,boxSize,bin_p); posVec(1) = -0.01;
[possib_posgrid,~] = spline_2d(x_all_vec,y_all_vec,posVec,s);
possib_pos_rate = exp(possib_posgrid*all_param(encoding_cells,pos_ind+1)');

%% do the decoding!!!
count = 0;
% time, neurons, position
decoded_x = nan(numel(decoding_ind_ds),1);
decoded_y = nan(numel(decoding_ind_ds),1);

for m = 1:num_decode_chunks
    
    decode_chunk_m = (decoding_chunk_start(m):decoding_chunk_start(m)+length_decode_chunk-1)';
    decoding_ind_ds_m = decode_chunk_m(num_int:num_int:numel(decode_chunk_m));
    
    for k = 1:numel(decoding_ind_ds_m)
        
        count = count + 1;
        
        spikes_k = sum(S(encoding_cells,decoding_ind_ds_m(k)-num_int+1:decoding_ind_ds_m(k)),2);
        
        base_input = repmat(all_param(encoding_cells,1)',num_int,1);
        
        count_j = 0;
        sh_input = nan(num_int,numel(encoding_cells));
        for j = encoding_cells'
            count_j = count_j + 1;
            sh_grid_j = shgrid_all{j};
            sh_input(:,count_j) = sh_grid_j(decoding_ind_ds_m(k)-num_int+1:decoding_ind_ds_m(k),:)*(all_param(j,sh_ind+1)');
        end
        
        possib_rate_nop = sum(exp(base_input + sh_input),1);
        possib_rate_nop = repmat(possib_rate_nop,num_poss_pos,1);
        
        possib_rate = possib_rate_nop.*possib_pos_rate;
        possib_rate_sum = sum(possib_rate,2);
        
        possib_f = sum(-possib_rate+repmat(spikes_k',num_poss_pos,1).*log(possib_rate),2);
        likelihood = exp(possib_f);
        
        %if k == 1
        if 1
            % compute the posterior using just the data from previous time
            % bin
            posterior_notnorm = likelihood.*pos_occupancy_prob;
            posterior_new = posterior_notnorm./sum(posterior_notnorm);
        else
            
            prior = new_state_prob*posterior_old;
            posterior_notnorm = likelihood.*prior;
            posterior_new = posterior_notnorm./sum(posterior_notnorm);
        end
        
        [~,possib_index] = max(posterior_new);
        decoded_x(count) = x_all_vec(possib_index);
        decoded_y(count) = y_all_vec(possib_index);
        posterior_old = posterior_new;
        
        %{
        figure(1)
        imagesc([0 150],[0 150],reshape(posterior_new,numel(posvec),numel(posvec))); colorbar
        hold on
        plot(decoded_x(count),decoded_y(count),'*r')
        plot(x_us(decoding_ind_ds(count)),y_us(decoding_ind_ds(count)),'*g')
        hold off
        keyboard
        %}
        % count
    end
end

x_error = abs(x_us(decoding_ind_ds)' - decoded_x);
y_error = abs(y_us(decoding_ind_ds)' - decoded_y);
error = sqrt(x_error.^2 + y_error.^2);

%{

figure(1)
subplot(3,1,1)
plot(decoded_x,'k','linewidth',1.5)
hold on
plot(x_us(decoding_ind_ds),'r','linewidth',2);
hold off
box off
ylabel('x position')
title(['total x-error = ',num2str(nanmean(x_error))])

subplot(3,1,2)
plot(decoded_y,'k','linewidth',1.5)
hold on
plot(y_us(decoding_ind_ds),'b','linewidth',2);
hold off
box off
ylabel('y position')
title(['total y-error = ',num2str(nanmean(y_error))])

subplot(3,1,3)
plot(x_error,'b','linewidth',1.5)
hold on
plot(y_error,'r','linewidth',1.5)
plot(error,'k','linewidth',1.5)
hold off
box off
ylabel('total error (cm)')
xlabel('time points')
title(['total error = ',num2str(nanmean(error))])

filename = 'decoding_fig1_memoryless';
fig1 = figure(1);
set(fig1,'units','pix','pos',[0,0,1500,500])
epsfig = hgexport('factorystyle');
epsfig.Format = 'jpg';
if windows
    hgexport(fig1,[data_path,'\',filename],hgexport('factorystyle'),'Format','jpeg')
else
    hgexport(fig1,[data_path,'/',filename],hgexport('factorystyle'),'Format','jpeg')
end



figure(2)
subplot(1,4,2)
[error_map,~] = compute_2d_tuning_curve(x_us(decoding_ind_ds),y_us(decoding_ind_ds),error,20,0,boxSize);
imagesc(error_map)
axis off
colorbar
caxis([0 100])
title('avg error across position')

subplot(1,4,4)
[spd_error_tc,error_spd_vector] = compute_1d_tuning_curve(speed_us(decoding_ind_ds),error,10,0,40);
plot(error_spd_vector,spd_error_tc,'k','linewidth',2)
hold on
line_fit = polyfit(error_spd_vector',spd_error_tc,1);
plot(error_spd_vector,error_spd_vector*line_fit(1) + line_fit(2),'--r','linewidth',1.5)
[rho,pval] = corr(speed_us(decoding_ind_ds)',error);
hold off
box off
xlabel('speed (cm/s)')
ylabel('avg error')
title(['rho = ',num2str(rho),' p = ',num2str(pval)])

subplot(1,4,3)
imagesc(flipud(reshape(sum(possib_pos_rate,2),numel(posvec),numel(posvec))))
axis off
colorbar
title('summed position rate maps')

subplot(1,4,1)
plot(x_us(decoding_ind),y_us(decoding_ind),'.k')
hold on
scatter(x_us(decoding_ind_ds),y_us(decoding_ind_ds),30,error,'filled')
hold off
axis off
colorbar
caxis([0 100])
axis([0 150 0 150])

filename = 'decoding_fig2_memoryless';
fig1 = figure(2);
set(fig1,'units','pix','pos',[0,0,1500,300])
epsfig = hgexport('factorystyle');
epsfig.Format = 'jpg';
if windows
    hgexport(fig1,[data_path,'\',filename],hgexport('factorystyle'),'Format','jpeg')
else
    hgexport(fig1,[data_path,'/',filename],hgexport('factorystyle'),'Format','jpeg')
end

close all
%}
