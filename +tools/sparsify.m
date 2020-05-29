
% SPARSIFY  Remove data points where x-gradients are below a threshold.
% Author: Timothy Sipkens, 2020-04-02
% 
% Note: Ix should be vectorized.
%=========================================================================%

function [A,Ix,It,idx_remove] = sparsify(A,Ix,It)

idx_remove = abs(Ix)<(max(abs(Ix)).*0.25);

A(:,idx_remove) = [];
It(idx_remove) = [];
Ix(idx_remove) = [];

end

