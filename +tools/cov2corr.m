
% COV2CORR  Function to convert covariance matrix to correlation matrix.
% Author:  Timothy Sipkens, 2019-10-29
%=========================================================================%

function R = cov2corr(Sigma)

D = 1./sqrt(diag(Sigma));
R = diag(D) * Sigma * diag(D);
    
end