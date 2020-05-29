
% SPARSIFYX  Remove data points where x-gradients are below a threshold.
% Author: Timothy Sipkens, 2020-04-02
% 
% Note: Ix should be vectorized.
%=========================================================================%

function [A] = sparsifyx(A,Ix)

idx_remove = Ix<(max(Ix).*1e-2);

A(idx_remove,:) = [];

end

