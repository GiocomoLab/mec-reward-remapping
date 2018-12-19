function [testFit,trainFit,param_mean,paramMat] = fit_model_kfold(A,spiketrain,dt,test_ind,train_ind)

%% Description
% This code will section the data into 10 different portions. Each portion
% is drawn from across the entire recording session. It will then
% fit the model to 9 sections, and test the model performance on the
% remaining section. This procedure will be repeated 10 times, with all
% possible unique testing sections. The fraction of variance explained, the
% mean-squared error, the log-likelihood increase, and the mean square
% error will be computed for each test data set. In addition, the learned
% parameters will be recorded for each section.

%4/9/18 cm added paramMat as an output for future use in created s.d. on ln
%plots.

%% Initialize matrices and section the data for k-fold cross-validation

numFolds = length(test_ind);

% initialize matrices
testFit = nan(numFolds,1);
trainFit = nan(numFolds,1); 
numCol = size(A,2);
paramMat = nan(numFolds,numCol);

%% perform k-fold cross validation
for k = 1:numFolds
    %fprintf('\t\t- Cross validation fold %d of %d\n', k, numFolds);
    
    test_spikes = spiketrain(test_ind{k}); %test spiking
    test_A = A(test_ind{k},:);
    test_dt = dt(test_ind{k});
    
    % training data
    train_spikes = spiketrain(train_ind{k});
    train_A = A(train_ind{k},:);
    train_dt = dt(train_ind{k});
    
    options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'HessianFcn','objective','Display','off');
    
    data{1} = train_A; data{2} = train_spikes; data{3} = train_dt;
    if k == 1
        init_param = 1e-3*randn(numCol, 1);
    else
        init_param = param;
    end
    [param] = fminunc(@(param) glm_model(param,data),init_param,options);
    
    % save the parameters
    paramMat(k,:) = param;
    
    %%%%%%%%%%%%% TEST DATA %%%%%%%%%%%%%%%%%%%%%%%
    
    % compute llh increase from "mean firing rate model" - NO SMOOTHING
    r = exp(test_A * param).*test_dt; n = reshape(test_spikes,numel(test_spikes),1); meanFR_test = nanmean(test_spikes); 
    
    log_llh_test_model = nansum(r-n.*log(r)+log(factorial(n)))/sum(n);
    log_llh_test_mean = nansum(meanFR_test-n.*log(meanFR_test)+log(factorial(n)))/sum(n);
    log_llh_test = (-log_llh_test_model + log_llh_test_mean);
    log_llh_test = log(2)*log_llh_test;
    
    % fill in all the relevant values for the test fit cases
    testFit(k) = log_llh_test;
    
    %%%%%%%%%%%%% TRAINING DATA %%%%%%%%%%%%%%%%%%%%%%%
    
    % compute log-likelihood
    
    r_train = exp(train_A * param).*train_dt; n_train = reshape(train_spikes,numel(train_spikes),1); meanFR_train = nanmean(train_spikes);   
    log_llh_train_model = nansum(r_train-n_train.*log(r_train)+log(factorial(n_train)))/sum(n_train);
    log_llh_train_mean = nansum(meanFR_train-n_train.*log(meanFR_train)+log(factorial(n_train)))/sum(n_train);
    log_llh_train = (-log_llh_train_model + log_llh_train_mean);
    log_llh_train = log(2)*log_llh_train;
    
    trainFit(k) = log_llh_train;

end

param_mean = nanmean(paramMat);

return
