function [testFit_all, trainFit_all, bestModels, bestModel_testfit, best_parameters, pvals, final_pval] = forward_search_kfold(A,spiketrain,train_ind,test_ind)

% do the forward search

num_var = length(A); % number of variables to search over

% do the forward-search method to identify variables encoded:

% var_vec = list of variables to add in
% variables = the variables currently in the model
% num_var = total number of variables to search over
% allModelFits = the model fit for every variable, on every iteration
% through the forward search procedure
% bestModelFits = the model fit for the best model (for every iteration)
% parametesr = parameters of the best model for each iteration
% pvals = the pvalue of the model comparison for each model

numFolds = length(test_ind);

% initialize values for while loop
var_vec = 1:num_var;
variables = [];
baseModel = -1*ones(numFolds,1);
testFit_all = nan(numFolds,num_var,num_var+1);
trainFit_all = nan(numFolds,num_var,num_var+1);
best_parameters = [];
bestModel_testfit= [];
bestModels = [];
pval = 0;
pvals = [];
loop_number = 0;

while pval < 0.05 && numel(variables) < num_var && numel(var_vec) > 0
    
    testFit = nan(numFolds,num_var);
    trainFit = nan(numFolds,num_var);
    parameters = cell(num_var,1);
    loop_number = loop_number + 1;

    for m = var_vec

        
        %fprintf('Fitting model  %d \n', m);
        
        % create matrix of variables in model currently
        X  = ones(length(spiketrain),1);
        temp_var = [variables m]; temp_var = sort(temp_var);
        for l = temp_var % this is always in order, for simplicity
            X = [X A{l}];
        end
        
        [testFit(:,m),trainFit(:,m),parameters{m}] = fit_model_kfold(X,spiketrain,test_ind,train_ind);
        
    end
    
        % save all of the model fits
        testFit_all(:,:,loop_number) = testFit;
        trainFit_all(:,:,loop_number) = trainFit;
    
        % choose the best model
        [~,topModel_ind] = max(nanmean(testFit));
        topModel = testFit(:,topModel_ind);
        
        % take out NaNs in case there aren't enough spikes or something
        baseModel1 = baseModel; baseModel1(isnan(topModel)) = [];
        topModel1 = topModel; topModel1(isnan(topModel)) = [];
        
        pval = signrank(topModel1,baseModel1,'tail','right');
        
        if pval < 0.05
            bestModel_testfit = [bestModel_testfit topModel];
            bestModels = [bestModels topModel_ind];
            best_parameters = parameters{topModel_ind};
        end
        
        pvals = [pvals pval];
        
        % find the variables in the model so far
        variables = [variables topModel_ind]; 
        variables = sort(variables);
        
        % find the new variables to try adding into the model
        %if loop_number == 1
        %    bad_var = find(mean(testFit) < 0);
        %end
        bad_var = [];
        var_vec = setdiff(1:num_var,union(variables,bad_var));

        baseModel = topModel;
end

% check that the final model is sig better than zero
final_pval = signrank(topModel,zeros(size(topModel)),'tail','right');

if final_pval > 0.05
    bestModels = [];
    best_parameters = {};
end


return




