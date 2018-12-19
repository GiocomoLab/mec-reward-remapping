function [f, df, hessian] = glm_model(param,data)

X = data{1}; % subset of A
n = data{2}; % number of spikes
dt = data{3};
n = reshape(n,numel(n),1);

% compute the firing rate
u = X * param;
rate = exp(u).*dt;

%% compute f, the gradient, and the hessian 

% L2 regularizer weight
alpha = 1;
f_l2 = 0.5*alpha*sum(param.^2);

f = sum(rate-n.*log(rate)) + f_l2;
df = real(X' * (rate - n) + alpha*param);
rX = bsxfun(@times,rate,X);       
hessian = rX'*X + alpha*eye(numel(param));

