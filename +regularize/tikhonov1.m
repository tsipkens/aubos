
% TIKHONOV1  Builds appropriate Lpr and inverts the system for 1D (axial slice) case.
%  
%  AUTHOR: Timothy Sipkens, 2021-02-04

function [x, Lpr] = tikhonov1(K, b, Le, lambda)

L_tk = regularize.tikhonov_lpr(2, size(K, 2), size(K, 2));
% L_tk(end,end) = L_tk(end,end) + 1;  % adjust last element boundary condition (to no slope) 

A = [Le * K; lambda.*L_tk];
b = [Le * b; sparse(zeros(size(L_tk,1), 1))];

x = lsqlin(A, b);

end

