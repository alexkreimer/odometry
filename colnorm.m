function n = colnorm(X)
% L2 over columns of X

n = sqrt(sum(X.*X,1));