function [train_ind,test_ind] = get_train_test_ind(numFolds,chunks,indices,dt)

T = numel(indices);
numPts = chunks*round(1/dt); % cut into 30 second chunks
numSeg = ceil(T/(numFolds*numPts));
oneSeg = ones(numPts,1)*(1:numFolds);
new_ind = repmat(oneSeg(:),numSeg,1);
new_ind = new_ind(1:T);

test_ind = cell(numFolds,1);
train_ind = cell(numFolds,1);
for k = 1:numFolds
    test_ind{k} = indices(new_ind == k);
    train_ind{k} = setdiff(indices,test_ind{k});
end