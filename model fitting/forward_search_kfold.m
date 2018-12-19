function [allModelFits, allModelTrainFits, bestModels, bestModelFits, parameters, pvals, final_pval] = forward_search_kfold(A,spiketrain,dt,train_ind,test_ind)

% do the forward search

num_var = length(A); % number of variables to search over

% do the forward-search method to identify variables encoded:

% var_vec = list of variables to add in
% variables = the variables currently in the model
% num_var = total number of variables to search over
% allModelFits = the model fit for every variable, on every iteration
% through the forward search procedure
% bestModelFits = the model fit for the best model (for every iteration)
% parameters = parameters of the best model for each iteration
% pvals = the pvalue of the model comparison for each model

numFolds = length(test_ind);

% initialize values for while loop
var_vec = 1:num_var;
variables = [];
baseModel = -5*ones(numFolds,1);
allModelFits = {};
allModelTrainFits = {};
bestModelFits= [];
bestModels = [];
parameters = {};
pval = 0;
pvals = [];
loop_number = 0;

while pval < 0.05 && numel(variables) < num_var && numel(var_vec) > 0
    
    testFit = nan(numFolds,num_var);
    trainFit = nan(numFolds,num_var);
    param_mean = cell(num_var,1);
    loop_number = loop_number + 1;
    
    for m = var_vec
        
        
        fprintf('Fitting model  %d \n', m);
        
        % create matrix of variables in model currently
        X  = ones(length(spiketrain),1);
        temp_var = [variables m]; temp_var = sort(temp_var);
        for l = temp_var % this is always in order, for simplicity
            X = [X A{l}];
        end
        
        [testFit(:,m),trainFit(:,m),param_mean{m}] = fit_model_kfold(X,spiketrain,dt,test_ind,train_ind);
        
    end
    
    % save all of the model test fits
    allModelFits{end+1} = testFit;
    allModelTrainFits{end+1} = trainFit;
    
    % choose the best model
    [~,topModel_ind] = max(nanmean(testFit));
    topModel = testFit(:,topModel_ind);
    
    % remove NaNs to do the significance test
    topModel_noNaN = topModel;
    baseModel_noNaN = baseModel;
    topModel_noNaN(isnan(baseModel_noNaN)) = []; baseModel_noNaN(isnan(baseModel_noNaN)) = [];
    baseModel_noNaN(isnan(topModel_noNaN)) = []; topModel_noNaN(isnan(topModel_noNaN)) = [];
    
    pval = signrank(topModel_noNaN,baseModel_noNaN,'tail','right');
    
    if pval < 0.05
        bestModelFits = [bestModelFits topModel ];
        bestModels = [bestModels topModel_ind ];
        parameters{end+1} = param_mean{topModel_ind};
    end
    
    pvals = [pvals pval];
    
    % find the variables in the model so far
    variables = [variables topModel_ind]; variables = sort(variables);
    
    % find the new variables to try adding into the modeld
    var_vec = setdiff(1:num_var,variables);
    
    baseModel = topModel;
    
    
end

% check that the final model is sig better than zero
topModel_noNaN = topModel;
topModel_noNaN(isnan(topModel_noNaN)) = [];
final_pval = signrank(topModel_noNaN,zeros(size(topModel_noNaN)),'tail','right');

return




