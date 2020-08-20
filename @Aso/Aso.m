
% ASO  A class to handle spatial information for a 1D axis-symmetric object (ASO).
% Such an object has no axial dependence, which is useful in demonstrating
% and testing kernels. 
% 
% Author: Timothy Sipkens, 2020-05-20
%=========================================================================%


classdef Aso
    
    properties
        dr    = [];     % annuli width
        r     = [];     % annuli centers
        re    = [];     % annuli edges
        
        R     = [];     % outer radius of object
        Nr    = [];     % number of annuli
        
        grad  = [];     % gradient operator
    end
    
    
    
    methods
        function aso = Aso(R, Nr)
            aso.Nr = Nr; % number of annuli
            aso.R = R; % outer radius
            
            aso.re = linspace(0, R, Nr+1)'; % linearily space edges from 0 -> R
            aso.r  = (aso.re(2:end) + aso.re(1:(end-1))) ./ 2; % annuli centers
            aso.dr = aso.re(2:end) - aso.re(1:(end-1)); % annuli width
            
            % Evaluate the radial gradient, assuming no slope at outer radius.
            D = (eye(aso.Nr+1, aso.Nr+1) - diag(ones(aso.Nr, 1), 1));
            D(end, :) = []; % remove final row
            aso.grad = D ./  aso.dr; % divide by element area
        end
        
        
        
        %== UNIFORM =======================================================%
        %   Evaluates kernel/operator for a uniform basis representation of an ASO.
        %   Not recommended due to noise properties.
        %   Timothy Sipkens, 2020-06-10
        %
        % Inputs:
        %   aso     Axis-symmetric object
        %   m       Set of slopes for the rays
        %   x0      Intersect with line through center of aso
        function K = uniform(aso,m,x0)
            K = kernel.uniform1(aso,m,x0);
        end
        
        
        
        %== LINEAR =======================================================%
        %   Bridge function from Aso to kernel.
        %   Timothy Sipkens, 2020-06-10
        %
        % Inputs:
        %   aso     Axis-symmetric object
        %   m       Set of slopes for the rays
        %   x0      Intersect with line through center of aso
        function K = linear(aso, m, x0)
            K = kernel.linear1(aso,m,x0);
        end
        
        
        
        %== SURF =========================================================%
        %   Plot the axis-symmetric object as a surface. 
        %   Timothy Sipkens, 2020-06-09
        function [h,t,r,z0] = surf(aso, x, f_grid)
            
            if ~exist('f_grid','var'); f_grid = []; end
            if isempty(f_grid); f_grid = 1; end
            
            [t,i] = meshgrid(linspace(0,2*pi,64), 1:(aso.Nr+1));
                % mesh of angles and integers for radial positions
            
            r = aso.re(i); % radii to plot
            x0 = r.*cos(t); % x values for plot
            y0 = r.*sin(t); % y values for plot
            z0 = x(i); % z value for plot, given by input data
            
            h = surf(x0,y0,z0); % generate the surface plot
            h.EdgeColor = 'none'; % remove default edges on plot
            axis image; % such that circles are properly round
            axis off; % remove axis to clean visualization
            
            if f_grid % overlay radial element edges
                hold on;
                plot3(x0', y0', z0', 'k'); % add circles corresponding to element edges
                hold off;
            end
            
            if nargout==0; clear h; end % suppress output if none required
            
        end
        
        
        
        %== PLOT =========================================================%
        %   Plot the axis-symmetric object as a series of annuli.
        %   Timothy Sipkens, 2020-06-09
        %   Note: Works best with monotonically increasing/decreasing z0.
        function h = plot(aso, x)
            [t,i] = meshgrid(linspace(0,2*pi,64), 1:(aso.Nr+1));
            
            x0 = aso.re(i).*cos(t);
            y0 = aso.re(i).*sin(t);
            z0 = x(i);
            
            h = contourf(x0,y0,z0,sort(x),'edgecolor','none');
            axis image;
            
            hold on;
            viscircles(ones(aso.Nr+1, 1) * [0,0], aso.re, ...
                'EnhanceVisibility', 0, 'LineWidth', 0.1, ...
                'Color', 'k', 'LineStyle', '-');
            hold off;
            
            if nargout==0; clear h; end % suppress output if none required
        end
        
        
        %== SRAYS ========================================================%
        %   Plot aso as a surface, with rays overlaid.
        %   Timothy Sipkens, 2020-06-09
        %   Note: Works best with monotonically increasing/decreasing z0.
        function h = srays(aso, x, m, x0)
            
            [h,t0,r0,z0] = aso.surf(x); % generate surface plot
            
            x1 = linspace(-aso.R, aso.R, 150);
            y1 = (m').*x1 + x0';
            r1 = sqrt(x1.^2 + y1.^2);
            t1 = atan2(y1, x1);
            t1(t1<0) = t1(t1<0) + 2*pi;
            
            z1 = interp2(t0, r0, z0, t1, r1, 'linear'); 
            hold on;
            plot3(x1, y1, z1, 'r');
            quiver(-aso.R, -aso.R, aso.R/3, 0, ...
                'MaxHeadSize', 0.8, 'Color', 'r');
            hold off;
            
            if nargout==0; clear h; end % suppress output if none required
        end
    end
end

