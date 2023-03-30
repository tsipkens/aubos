
% NORMC  A function to normalize columns of vector. 
%  A simple utility to substitute for Matlab equivalent. 
%  
%  AUTHOR: Timothy Sipkens, 2023-03-30

function v = normc(v)

v = v ./ sqrt(sum(v .^ 2));

end
