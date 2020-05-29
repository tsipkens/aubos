
% ASO  A class to handle spatial information for axis-symmetric objects (ASOs).
% Author: Timothy Sipkens, 2020-05-20
%=========================================================================%


classdef Aso
    
    properties
        dr    = [];     % annuli width
        r     = [];     % annuli centers
        re    = [];     % annuli edges
        
        N     = [];     % number of annuli
        
        phi   = [];     % unscaled basis functions
    end
    
    
    
    methods
        function obj = Aso(R,N)
            obj.N = N; % number of annuli
            
            obj.re = linspace(0,R,N+1)'; % linearily space edges from 0 -> R
            obj.r  = (obj.re(2:end) + obj.re(1:(end-1))) ./ 2; % annuli centers
            obj.dr = obj.re(2:end) - obj.re(1:(end-1)); % annuli width
        end
        
        
        
        
    end
end

