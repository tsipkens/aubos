
% ASO  A class to handle spatial information for 1D axisymmetric objects.
%  Such an object has no axial dependence, which is useful in demonstrating
%  and testing kernels. This class provides the basics for discretizing the
%  space and plotting functions on the space. 
%  
%  A = Aso(NR) creates an discrete, axisymmetric object with NR
%  annuli/elements. Uses default radius of R = 1.
%  
%  A = Aso(NR,R) creates an discrete, axisymmetric object with NR
%  annuli/elements and an outer radius of R. 
%  
%  ------------------------------------------------------------------------
% 
%  AUTHOR: Timothy Sipkens, 2020-05-20

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
        %== ASO ==========================================================%
        %   Constructor method.
        function aso = Aso(Nr, R)
            if nargin==0; return; end % return empty object
            
            if ~exist('R', 'var'); R = []; end
            if isempty(R); R = 1; end  % use R = 1 by default
            
            aso.Nr = Nr; % number of annuli
            aso.R = R; % outer radius
            
            aso.re = linspace(0, R, Nr+1)'; % linearily space edges from 0 -> R
            aso.r  = (aso.re(2:end) + aso.re(1:(end-1))) ./ 2; % annuli centers
            aso.dr = aso.re(2:end) - aso.re(1:(end-1)); % annuli width
            
            % Evaluate the radial gradient, assuming no slope at outer radius.
            D = speye(aso.Nr+1, aso.Nr+1) - ...
                spdiags(ones(aso.Nr+1, 1), 1, aso.Nr+1, aso.Nr+1);
            D(end, :) = []; % remove final row
            aso.grad = D ./  aso.dr; % divide by element area
        end
        %=================================================================%
        
        
        %== INTERPC ====================================================%
        %   Interpolate discrete function defined on an ASO to Cartesian
        %   coordinates. Assumes a linear basis.
        %   Timothy Sipkens, 2021-02-12
        function [f] = interpc(aso, yi, zi, bet)
            
            % Get position in cylindrical coordinates
            ri = sqrt(yi.^2 + zi.^2); % radial position
            
            % Interpolate r-gradient and convert to Cartesian coords.
            f = interp1(aso.re, bet, ri, 'linear', 0);
            
        end
        %=================================================================%
        
        
        %== GRADIENTC ====================================================%
        %   Cartesian gradients, interpolated from the ASO object.
        %   Assumes a linear basis.
        %   Timothy Sipkens, 2021-02-11
        function [Dy, Dz] = gradientc(aso, yi, zi, bet)
            
            % Get position in cylindrical coordinates
            ri = sqrt(yi.^2 + zi.^2); % radial position
            q = atan2(yi, -zi); % angle, in coord. system from Sipkens et al.
            
            % Setup grid for interpolation
            Dri = aso.grad * bet;  % reshape radial gradient
            Dri = [Dri; 0]; % append constant slope data for last axial position
            
            % Interpolate r-gradient and convert to Cartesian coords.
            Dro = interp1(aso.re, Dri, ri, 'linear', 0);
            Dy = -sin(q) .* Dro; % get the y-gradient based on the angle
            Dz = cos(q) .* Dro; % get the z-gradient based on the angle
            
        end
        %=================================================================%
        
        
        %== SURF =========================================================%
        %   Plot the axis-symmetric object as a surface. 
        %   Timothy Sipkens, 2020-06-09
        function [h,t,r,z0] = surf(aso, bet, f_grid)
            
            if ~exist('f_grid','var'); f_grid = []; end
            if isempty(f_grid); f_grid = 1; end
            
            [t,i] = meshgrid(linspace(0,2*pi,64), 1:(aso.Nr+1));
                % mesh of angles and integers for radial positions
            
            r = aso.re(i); % radii to plot
            x0 = r.*cos(t); % x values for plot
            y0 = r.*sin(t); % y values for plot
            z0 = bet(i); % z value for plot, given by input data
            
            h = surf(x0,y0,z0); % generate the surface plot
            h.EdgeColor = 'none'; % remove default edges on plot
            axis image; % such that circles are properly round
            axis off; % remove axis to clean visualization
            
            % Add circles corresponding to element edges.
            if f_grid
                hold on;
                plot3(x0', y0', z0', ...
                    'Color', 'k', 'LineWidth', 0.03);
                hold off;
            end
            
            if nargout==0; clear h; end % suppress output if none required
            
        end
        %=================================================================%
        
        
        
        %== PLOT =========================================================%
        %   Plot the axis-symmetric object as a series of annuli.
        %   Timothy Sipkens, 2020-06-09
        %   Note: Works best with monotonically increasing/decreasing z0.
        function h = plot(aso, bet, f_edges, f_img)
            
            if ~exist('f_edges', 'var'); f_edges = []; end
            if isempty(f_edges); f_edges = 1; end
            
            if ~exist('f_img', 'var'); f_img = []; end
            if isempty(f_img); f_img = 0; end
            
            if f_img
                x0 = linspace(-aso.R, aso.R, 201);
                y0 = linspace(-aso.R, aso.R, 205);
                [x1, y1] = meshgrid(x0, y0);
                z0 = aso.interpc(x1, y1, bet);
                h = imagesc(x0, y0, z0);
            else
                [t,i] = meshgrid(linspace(0,2*pi,64), 1:(aso.Nr+1));
                x0 = aso.re(i) .* cos(t);
                y0 = aso.re(i) .* sin(t);
                z0 = bet(i);
                
                h = contourf(x0, y0, z0, sort(bet), ...
                    'edgecolor','none');
            end
            axis image;
            
            if f_edges
                hold on;
                viscircles(ones(aso.Nr+1, 1) * [0,0], aso.re, ...
                    'EnhanceVisibility', 0, 'LineWidth', 0.1, ...
                    'Color', 'k', 'LineStyle', '-');
                hold off;
            else  % plot only outer circle
                hold on;
                viscircles([0,0], aso.R, ...
                    'EnhanceVisibility', 0, 'LineWidth', 0.1, ...
                    'Color', 'k', 'LineStyle', '-');
                hold off;
            end
            
            if nargout==0; clear h; end % suppress output if none required
        end
        %=================================================================%
        
        
        %== SRAYS ========================================================%
        %   Plot aso as a surface, with rays overlaid.
        %   Timothy Sipkens, 2020-06-09
        %   Note: Works best with monotonically increasing/decreasing z0.
        function h = srays(aso, bet, my, y0)
            
            [h,t0,r0,z0] = aso.surf(bet); % generate surface plot
            
            x1 = linspace(-aso.R, aso.R, 150);
            y1 = (my').*x1 + y0';
            r1 = sqrt(x1.^2 + y1.^2);
            t1 = atan2(y1, x1);
            t1(t1<0) = t1(t1<0) + 2*pi;
            
            z1 = interp2(t0, r0, z0, t1, r1, 'linear'); 
            hold on;
            plot3(x1, y1, z1 + 1e-3, ...
                'Color', [1, 0.3, 0.3], ...
                'LineWidth', 0.8);
            quiver(-aso.R, -aso.R, aso.R/3, 0, ...
                'MaxHeadSize', 0.8, ...
                'Color', 'k', ...
                'LineWidth', 0.8);
            hold off;
            
            if nargout==0; clear h; end % suppress output if none required
        end
        
        
        %== PRAYS ========================================================%
        %   Plot aso as a surface, with rays overlaid.
        %   Timothy Sipkens, 2020-06-09
        %   NOTE: Same as 'srays' but operated in 2D.
        function h = prays(aso, bet, my, y0, f_edges)
            
            if ~exist('f_edges', 'var'); f_edges = []; end
            if isempty(f_edges); f_edges = 1; end
            
            g = aso.plot(bet, f_edges); % generate surface plot
            
            x1 = linspace(-aso.R, aso.R, 200);
            
            % Used to restrict points to those in the ASO. 
            xmax = real(-(my') .* y0' + ...
                sqrt((my').^2 .* (y0').^2 - ...
                (1 + (my').^2) .* ((y0').^2 - aso.R.^2))) ./ ...
                ((1 + (my').^2));
            xmin = real(-(my') .* y0' - ...
                sqrt((my').^2 .* (y0').^2 - ...
                (1 + (my').^2) .* ((y0').^2 - aso.R.^2))) ./ ...
                ((1 + (my').^2));
            
            y1 = (my') .* max(min(x1, xmax), xmin) + y0';
             
            hold on;
            plot(max(min(x1, xmax), xmin)', y1', ...
                'Color', [1, 0.3, 0.3], ...
                'LineWidth', 0.8);
            quiver(-aso.R, -aso.R, aso.R/3, 0, ...
                'MaxHeadSize', 0.8, ...
                'Color', 'k', ...
                'LineWidth', 0.8);
            hold off;
            
            xlim(aso.R .* [-1.2, 1.2]);
            ylim(aso.R .* [-1.2, 1.2]);
            
            if nargout==0; clear h; end % suppress output if none required
        end
        %=================================================================%
    end
end

