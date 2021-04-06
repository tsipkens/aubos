
% ASO2  A class to handle spatial information for 2D axisymmetric objects.
% Author: Timothy Sipkens, 2020-05-20
%=========================================================================%

classdef Aso2
    
    properties
        %-- Axial properties ---------------------------------------------%        
        dx    = [];     % axial element width
        x     = [];     % axial element centers
        xe    = [];     % axial element edges
        xe2   = [];     % axial edges (for all elements in 2D grid),
                        % to be used with 2D phantoms defined on 
                        % a linear basis for radius and uniform basis for
                        % axial position
        
        X     = [];     % axial extent / overall width
        Nx    = [];     % number of edges in the axial direction
        %-----------------------------------------------------------------%
        
        
        %-- Radial properties (same as Aso) ------------------------------%
        dr    = [];     % annuli width
        r     = [];     % annuli centers
        re    = [];     % annuli edges
        re2   = [];     % annuli edges (for all elements in 2D grid),
                        % to be used with 2D phantoms defined on 
                        % a linear basis for radius and uniform basis for
                        % axial position
        
        R     = [];     % outer radius of object
        Nr    = [];     % number of annuli
        %-----------------------------------------------------------------%
        
        
        N     = [];     % total number of elements (Nr*Nv)
        edges = [];     % edges of all the elements (N x 4 array)
        
        Dr    = [];     % r-gradient operator
        Dx    = [];     % y-gradient operator
    end
    
    
    
    methods
        function aso = Aso2(Nr, R, Nx, X)
            if nargin==0; return; end % return empty object
            
            %-- Radial positions and discretization ----------------------%
            aso.Nr = Nr; % number of annuli
            aso.R = R; % outer radius
            
            aso.re = linspace(0, R, Nr+1)'; % linearily space edges from 0 -> R
            aso.r  = (aso.re(2:end) + aso.re(1:(end-1))) ./ 2; % annuli centers
            aso.dr = aso.re(2:end) - aso.re(1:(end-1)); % annuli width
            %-------------------------------------------------------------%
            
            %-- Axial positions and discretization -----------------------%
            aso.Nx = Nx; % number of annuli
            aso.X = X; % outer radius
            
            aso.xe = linspace(0, X, Nx+1)'; % linearily space edges from 0 -> R
            aso.x  = (aso.xe(2:end) + aso.xe(1:(end-1))) ./ 2; % annuli centers
            aso.dx = aso.xe(2:end) - aso.xe(1:(end-1)); % annuli width
            %-------------------------------------------------------------%
            
            %-- Consolidated element information -------------------------%
            aso.N = aso.Nx * aso.Nr; % total number of elements
            
            % vector list of element edges (for 2D linear basis)
            [aso.xe2, aso.re2] = meshgrid(aso.xe(1:(end-1)), aso.re);
            aso.xe2 = aso.xe2(:); aso.re2 = aso.re2(:);
            
            % vectorized element edges for overall grid
            aso.edges = [repmat(aso.re(1:(end-1)),[Nx,1]), ...
                repmat(aso.re(2:end),[Nx,1]), ...
                reshape(repmat(aso.xe(1:(end-1)),[1,Nr]), [Nx*Nr,1]), ...
                reshape(repmat(aso.xe(2:end),[1,Nr]), [Nx*Nr,1])];
            %-------------------------------------------------------------%
            
            %-- Compute differential operators ---------------------------%
            I1 = speye(aso.Nr+1, aso.Nr+1);
            E1 = sparse(1:aso.Nr+1-1, 2:aso.Nr+1, 1, aso.Nr+1, aso.Nr+1);
            D1 = (E1 - I1) ./ [aso.dr; aso.dr(end)];
            D1(end,end) = 0; % applies no slope at final radial position
            
            I2 = speye(aso.Nx,aso.Nx);
            E2 = sparse(1:aso.Nx-1, 2:aso.Nx, 1, aso.Nx, aso.Nx);
            D2 = (E2 - I2) ./ aso.dx;
            D2(end,end) = 1; D2(end,end-1) = -1; % applies same slope as previous y-position
            
            aso.Dr = kron(I2, D1);
            aso.Dx = kron(D2, I1);
            %-------------------------------------------------------------%
            
            
        end
        
        
        
        %== RESHAPE ======================================================%
        %   Reshape, assuming x is evaluated at aso.re.
        %   Timothy Sipkens, 2020-06-11
        function x = reshape(aso,x)
            x = reshape(x, [aso.Nr + 1, aso.Nx]);
        end
        %=================================================================%
        
        
        
        %== GRADIENTC ====================================================%
        %   Cartesian gradients, interpolated from the ASO object.
        %   Assumes a linear basis.
        %   Timothy Sipkens, 2020-06-11
        function [Dx, Dy, Dz] = gradientc(aso, xi, yi, zi, bet)
            
            % Get position in cylindrical coordinates
            ri = sqrt(yi.^2 + zi.^2); % radial position
            q = atan2(yi, -zi); % angle, in coord. system from Sipkens et al.
            
            % Setup grid for interpolation
            [x0, r0] = meshgrid(aso.xe, aso.re); % grid for input to interpolation
            Dr0 = aso.reshape(aso.Dr * bet); % reshape radial gradient
            Dr0 = [Dr0, Dr0(:,end)]; % append constant slope data for last axial position
            Dx0 = aso.reshape(aso.Dx * bet); % reshape axial gradient
            Dx0 = [Dx0, Dx0(:,end)]; % append constant slope data for last axial position
            
            % Interpolate r-gradient and convert to Cartesian coords.
            Dri = interp2(x0, r0, Dr0, xi, ri, 'linear', 0);
            Dy = sin(q) .* Dri; % get the y-gradient based on the angle
            Dz = -cos(q) .* Dri; % get the z-gradient based on the angle
            
            % Interpolate x-gradient
            Dx = interp2(x0, r0, Dx0, xi, ri, 'linear', 0);
            
        end
        %=================================================================%
        
        
        
        %== INTERPC ====================================================%
        %   Interpolate discrete function defined on an ASO to Cartesian
        %   coordinates. Assumes a linear basis and uniform axial basis.
        %   Timothy Sipkens, 2020-06-11
        function [f] = interpc(aso, xi, yi, zi, bet)
            
            % Get position in cylindrical coordinates
            ri = sqrt(yi.^2 + zi.^2); % radial position
            
            % Setup grid for interpolation
            [x0, r0] = meshgrid(aso.xe(1:end-1), aso.re); % grid for input to interpolation
            bet = reshape(bet, [aso.Nr+1, aso.Nx]);
            
            % Interpolate r-gradient and convert to Cartesian coords.
            f = interp2(x0, r0, bet, xi, ri, 'linear' ,0);
            
        end
        %=================================================================%
        
        
        
        
        %== SURF =========================================================%
        %   Plot the axis-symmetric object as a surface. 
        %   Timothy Sipkens, 2020-06-09
        function [h,x0,y0,z0] = surf(aso, bet, f_grid)
            
            if ~exist('f_grid','var'); f_grid = []; end
            if isempty(f_grid)
                if bet>1e3; f_grid = 1; else; f_grid = 0; end
            end
            
            [iv,ir] = meshgrid(1:aso.Nx, 1:(aso.Nr+1));
            
            y1 = aso.xe(iv);
            x1 = aso.re(ir);
            
            x0 = [flipud(y1);  y1(2:end,:)]; % plot axial positions
            y0 = [-flipud(x1); x1(2:end,:)]; % plot radii
            
            z1 = reshape(bet, [aso.Nr+1, aso.Nx]);
            z0 = [flipud(z1); z1(2:end,:)];
            
            h = surf(x0,y0,z0);
            h.EdgeColor = 'none';
            axis image;
            set(gca,'YDir','normal');
            
            if f_grid % overlay element grid
                hold on;
                plot3(x0', y0', z0', ...
                    'Color', 'k', 'LineWidth', 0.03);
                plot3(x0, y0, z0, ...
                    'Color', 'k', 'LineWidth', 0.03);
                hold off;
            end
            
            xlabel('Axial, x');
            ylabel('Radial, r');
            
            if nargout==0; clear h; end % suppress output if none required
            
        end
        
        
        
        %== PLOT =========================================================%
        %   Plot the axis-symmetric object as a surface. 
        %   Timothy Sipkens, 2020-06-09
        function [h,x0,y0,z0] = plot(aso, bet, f_grid)
            
            if ~exist('f_grid','var'); f_grid = []; end
            if isempty(f_grid)
                if bet>1e3; f_grid = 1; else; f_grid = 0; end
            end
            
            [iv,ir] = meshgrid(1:aso.Nx, 1:(aso.Nr+1));
            
            y1 = aso.xe(iv);
            x1 = aso.re(ir);
            
            x0 = [flipud(y1);  y1(2:end,:)]; % plot axial positions
            y0 = [-flipud(x1); x1(2:end,:)]; % plot radii
            
            z1 = reshape(bet, [aso.Nr+1, aso.Nx]);
            z0 = [flipud(z1); z1(2:end,:)];
            
            h = imagesc(x0(1,:),y0(:,1),z0);
            axis image;
            set(gca,'YDir','normal');
            
            if f_grid % overlay element grid
                hold on;
                plot(x0, y0, 'k');
                plot(x0, y0, 'k');
                hold off;
            end
            
            xlabel('Axial, x');
            ylabel('Radial, r');
            
            if nargout==0; clear h; end % suppress output if none required
            
        end
        
        
        
        %== SRAYS ========================================================%
        %   Plot ASO as a surface, with rays overlaid.
        %   Timothy Sipkens, 2020-06-09
        %   Note: Works best with monotonically increasing/decreasing z0.
        function h = srays(aso, bet, mx, x0, f_grid)
            
            if ~exist('f_grid','var'); f_grid = []; end
            
            [h,x1,y1,z1] = aso.surf(bet,f_grid); % generate surface plot
            
            y2 = linspace(-aso.R, aso.R, 150);
            x2 = (mx(1:(5*aso.Nr):end)').*y2 + x0(1:(5*aso.Nr):end)';
            
            y2 = y2.*ones(size(x2)); % if necessary, will out y1 to have correct dimension
            
            z2 = interp2(x1, y1, z1, x2, y2, 'linear'); 
            hold on;
            plot3(x2', y2', z2', 'r');
            hold off;
            
            % Format z-axis characteristics
            ax = gca;
            ax.ZAxis.Visible = 'off';
            ax.GridLineStyle = 'none';
            ax.Color = 'none';
            
            if nargout==0; clear h; end % suppress output if none required
        end
        
        
        
        %== SRAYS ========================================================%
        %   Plot ASO as a surface, with rays overlaid.
        %   Timothy Sipkens, 2020-06-09
        %   Note: Same as 'srays' but operates in 2D instead of 3D.
        function h = prays(aso, bet, mx, x0, f_grid)
            
            if ~exist('f_grid','var'); f_grid = []; end
            
            h = aso.plot(bet,f_grid); % generate surface plot
            
            y2 = linspace(-aso.R, aso.R, 150);
            x2 = (mx(1:(5*aso.Nr):end)').*y2 + x0(1:(5*aso.Nr):end)';
            
            y2 = y2.*ones(size(x2)); % if necessary, will out y1 to have correct dimension
            
            hold on;
            plot(x2', y2', 'r');
            hold off;
            
            % Format z-axis characteristics
            ax = gca;
            ax.GridLineStyle = 'none';
            ax.Color = 'none';
            
            axis image;
            xlim([min(aso.xe), max(aso.xe)]);
            
            if nargout==0; clear h; end % suppress output if none required
        end
        
        
        
        %== PLOT_SLICE ===================================================%
        %   Plot a slice through the ASO.
        %   Plots slide closest to x0, rather than interpolating.
        %   Timothy Sipkens, 2020-10-19
        function [h,idx] = plot_slice(aso, bet, x0)
            
            %-- Parse inputs, assign x0 index ----------------------------%
            if ~exist('x0','var'); x0 = []; end
            if isempty(x0) % if empty, use midpoint
                x0 = aso.xe(floor((aso.Nx + 1) / 2));
            end
            [~,idx] = min(abs(x0 - aso.xe)); % find index closest to x0
            %-------------------------------------------------------------%
            
            
            bet = aso.reshape(bet); % reshape beta to rectangle
            
            % if x0 is not exact match, display actual slice location
            if aso.xe(idx)~=x0
                disp(['Showing closest slice: x = ', ...
                    num2str(aso.xe(idx))]);
                disp(' ');
            end
            
            % output plot
            h = plot(aso.re, bet(:, idx));
            
            if nargout==0; clear h; end
            
        end
    end
end

