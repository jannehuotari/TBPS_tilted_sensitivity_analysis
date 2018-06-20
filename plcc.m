% This function calculates the Pearson's correlation coefficient between
% arratys X and Y
function [rho] = plcc(X,Y)
n_samples = size(X,1); % Amount of samples of variance in parameters
n_inputs = size(X,2); % Amount of parameters

for i = 1:n_inputs
    rho(i) = abs(corr(X(:,i),Y, 'Type','Pearson')); % Pearson's correlation coefficient between -1 and 1. We are not interested in direction of correlation,
    % so absolute value is used...
end
end
