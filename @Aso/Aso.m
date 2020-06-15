
% ASO  A class to handle spatial information for an axis-symmetric object (aso).
% Author: Timothy Sipkens, 2020-05-20
%=========================================================================%


classdef Aso
    
    properties
        dr    = [];     % annuli width
        r     = [];     % annuli centers
        re    = [];     % annuli edges
        
        R     = [];     % outer radius of object
        N     = [];     % number of annuli
    end
    
    
    
    methods
        function aso = Aso(R,N)
            aso.N = N; % number of annuli
            aso.R = R; % outer radius
            
            aso.re = linspace(0, R, N+1)'; % linearily space edges from 0 -> R
            aso.r  = (aso.re(2:end) + aso.re(1:(end-1))) ./ 2; % annuli centers
            aso.dr = aso.re(2:end) - aso.re(1:(end-1)); % annuli width
        end
        
        
        
        %== GRAD =========================================================%
        %   Radial gradient operator, assuming no slope at outer radius.
        %   Timothy Sipkens, 2020-06-10
        function D = grad(aso)
            D = (eye(aso.N+1, aso.N+1) - diag(ones(aso.N, 1), 1));
            D(end, :) = []; % remove final row
            D = D ./  aso.dr; % divide by element area
        end
        
        
        
        %== UNIFORM =======================================================%
        %   Evaluates kernel/operator for a uniform basis representation of an ASO.
        %   Not recommended due to noise properties.
        %   Timothy Sipkens, 2020-06-10
        %
        % Inputs:
        %   aso     Axis-symmetric object
        %   m       Set of slopes for the rays
        %   u0      Intersect with line through center of aso
        function K = uniform(aso,m,u0)
            
            rj = aso.re(1:(end-1)); % r_j
            rju  = aso.re(2:end); % r_{j+1}
            
            %-{
            Ka = @(m,u0,r) sqrt(r.^2 - u0.^2 ./ (1+m.^2));
            Kb = @(m,u0,r) log(r + Ka(m,u0,r)); % function for indefinite integral
            
            K = real(2 .* u0 .* ( ... % real(.) removes values outside integral bounds
                Kb(m,u0,rju) - ...
                Kb(m,u0,rj)))';
                % uniform basis kernel function at specified m and u0
            
            K = -K*aso.grad; % gradient is implemented as seperate operator (better noise characteristics)
            %}
            
            %{
            Ka = @(m,u0,r) 1 ./ sqrt(r.^2 - u0.^2 ./ (1+m.^2));
            
            K = real(2 .* u0 .* ( ... % real(.) removes values outside integral bounds
                Ka(m,u0,rj) - Ka(m,u0,rju)))';
                % uniform basis kernel function at specified m and u0
            %}
            
            K(abs(K)<300*eps) = 0; % remove numerical noise
            
        end
        
        
        
        %== LINEAR =======================================================%
        %   Evaluates kernel/operator for a linear basis representation of an ASO.
        %   Timothy Sipkens, 2020-06-10
        %
        % Inputs:
        %   aso     Axis-symmetric object
        %   m       Set of slopes for the rays
        %   u0      Intersect with line through center of aso
        function K = linear(aso,m,u0)
            
            if aso.N<3; error('Aso does not have enough annuli for linear basis.'); end
            
            rjd = aso.re(1:(end-2)); % r_{j-1}
            rj  = aso.re(2:(end-1)); % r_j
            rju = aso.re(3:end);     % r_{j+1}
            
            Kb = @(m,u0,r) log(r + sqrt(r.^2 - u0.^2 ./ (1+m.^2)));
            Kc = @(m,u0,r1,r2,r3) 1 ./ (r2 - r1) .* (Kb(m,u0,r3));
                % function for indefinite integral
            
            K = real(2 .* u0 .* ( ... % real(.) removes values outside integral bounds
                [ ...
                 zeros(1,length(m)); ...
                 Kc(m,u0,rjd,rj,rj) - Kc(m,u0,rjd,rj,rjd); ... % integral over rise
                 Kc(m,u0,rj(end),rju(end),rju(end)) - Kc(m,u0,rj(end),rju(end),rj(end)) ... % incline, last element
                ] + [
                 Kc(m,u0,rj(1),rjd(1),rj(1)) - Kc(m,u0,rj(1),rjd(1),rjd(1)); ... % decline, last element
                 Kc(m,u0,rju,rj,rju) - Kc(m,u0,rju,rj,rj); ... % integral over decline
                 zeros(1,length(m))
                ]))';
            
            K(abs(K)<10*eps) = 0; % remove numerical noise
            
        end
        
        
        
        %== SURF =========================================================%
        %   Plot the axis-symmetric object as a surface. 
        %   Timothy Sipkens, 2020-06-09
        function [h,t,r,z0] = surf(aso,x)
            [t,i] = meshgrid(linspace(0,2*pi,64), 1:(aso.N+1));
            
            r = aso.re(i); % plot radii
            x0 = r.*cos(t);
            y0 = r.*sin(t);
            z0 = x(i);
            
            h = surf(x0,y0,z0);
            h.EdgeColor = 'none';
            axis image;
            
            hold on;
            plot3(x0', y0', z0', 'k');
            hold off;
            
            if nargout==0; clear h; end % suppress output if none required
            
        end
        
        
        
        %== PLOT =========================================================%
        %   Plot the axis-symmetric object as a series of annuli.
        %   Timothy Sipkens, 2020-06-09
        %   Note: Works best with monotonically increasing/decreasing z0.
        function h = plot(aso,x)
            [t,i] = meshgrid(linspace(0,2*pi,64), 1:(aso.N+1));
            
            x0 = aso.re(i).*cos(t);
            y0 = aso.re(i).*sin(t);
            z0 = x(i);
            
            h = contourf(x0,y0,z0,sort(x),'edgecolor','none');
            axis image;
            
            hold on;
            viscircles(ones(aso.N+1, 1) * [0,0], aso.re, ...
                'EnhanceVisibility', 0, 'LineWidth', 0.1, ...
                'Color', 'k', 'LineStyle', '-');
            hold off;
            
            if nargout==0; clear h; end % suppress output if none required
        end
        
        
        %== SRAYS ========================================================%
        %   Plot aso as a surface, with rays overlaid.
        %   Timothy Sipkens, 2020-06-09
        %   Note: Works best with monotonically increasing/decreasing z0.
        function h = srays(aso,x,m,u0)
            
            [h,t0,r0,z0] = aso.surf(x); % generate surface plot
            
            x1 = linspace(-aso.R, aso.R, 150);
            y1 = (m').*x1 + u0';
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

