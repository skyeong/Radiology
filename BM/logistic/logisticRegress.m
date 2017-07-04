function [prob, theta] = logisticRegress(X,y)
%% ================== Part 1: Compute Cost and Gradient ===================
% Setup the data matrix
[m,n] = size(X);

% Add intercept term to x and X_test
X = [ones(m,1) X];

% Initialize fitting parameters
initial_theta = zeros(n+1,1);


% compute and display initial cost and gradient
[cost,grad] = costFunction(initial_theta, X, y);


%% ================== Part 2: Optimizing using fminunc ===================
options = optimset('GradObj','on','MaxIter',1000,'Display','off');


% Run fminunc to obtain the optimal theta
% This function will return theta and the cost
[theta,cost] = fminunc(@(t)(costFunction(t,X,y)), initial_theta, options);


%% ================== Part 3: Prediction and Accuracies ===================
prob = sigmoid(X*theta);
