
% TOTAL_VAR  Computes the total variation prior.
%  
%  AUTHOR: Timothy Sipkens, 2021-04-27

function L = total_var(x, aso)

if isa(aso, 'Aso2')
    Dx = aso.Dx * x;
    Dr = aso.Dr * x;
    
else  % assume dimensions supplied
    [Dx, Dr] = gradient(reshape(x, aso));
    
end

L = sum(abs(Dx(:)) + abs(Dr(:)));

end

