
% GEN1  General, one-dimensional optical flow operator (ignore axial deflections).
% Timothy Sipkens, 2020-06-18
%=========================================================================%

function [O] = gen1(dims)

dim1 = dims(1);
dim2 = dims(2);


%-- Generate differential operator ---------%
I1 = speye(dim1, dim1);
E1 = full(sparse(1:dim1-1, 2:dim1, 1, dim1, dim1));
D1 = E1 - I1;

I2 = speye(dim2, dim2);

O = kron(I2,D1);
%-------------------------------------------%


end


