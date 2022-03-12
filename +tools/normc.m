
% NORMC  Normalize the columns of a matrix. 
%   A local version to act as a drop-in replacement of 
%   MATLAB's original function. 
%  
%  AUHTOR: Timothy Sipkens, 2022-03-12

function A = normc(A)

A = A ./ sqrt(sum(A .^ 2));

end

