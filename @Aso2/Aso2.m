
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
        function aso = Aso2(R, Nr, X, Nx)
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
            E1 = sparse(1:aso.Nr+1-1,2:aso.Nr+1,1,aso.Nr+1,aso.Nr+1);
            D1 = E1-I1;
            D1(end,end) = 0; % applies no slope at final radial position
            
            I2 = speye(aso.Nx,aso.Nx);
            E2 = sparse(1:aso.Nx-1,2:aso.Nx,1,aso.Nx,aso.Nx);
            D2 = E2-I2;
            D2(end,end) = 1; D2(end,end-1) = -1; % applies same slope as previous y-position
            
            aso.Dr = kron(I2,D1);
            aso.Dx = kron(D2,I1);
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
        %   Assumes a linear radial basis and uniform axial basis.
        %   Timothy Sipkens, 2020-06-11
        function [Dx, Dy, Dz] = gradientc(aso, xi, yi, zi, f)
            
            % Get position in cylindrical coordinates
            ri = sqrt(yi.^2 + zi.^2); % radial position
            q = atan2(yi, -zi); % angle, in coord. system from Sipkens et al.
            
            % Setup grid for interpolation
            [x0, r0] = meshgrid(aso.xe, aso.re); % grid for input to interpolation
            Dri = aso.reshape(aso.Dr * f); % reshape radial gradient
            Dri = [Dri, Dri(:,end)]; % append constant slope data for last axial position
            Dxi = aso.reshape(aso.Dx * f); % reshape axial gradient
            Dxi = [Dxi, Dxi(:,end)]; % append constant slope data for last axial position
            
            % Interpolate r-gradient and convert to Cartesian coords.
            Dro = interp2(x0, r0, Dri, xi, ri, 'linear' ,0);
            Dy = sin(q) .* Dro; % get the y-gradient based on the angle
            Dz = -cos(q) .* Dro; % get the z-gradient based on the angle
            
            % Interpolate x-gradient
            Dx = interp2(x0, r0, Dxi, xi, ri, 'linear', 0);
            
        end
        %=================================================================%
        
        
        
        %== LINEAR =======================================================%
        %   Evaluates kernel/operator for a linear basis representation of an ASO.
        %   Linear basis is applied in radial direction. 
        %   Uniform basis is applied in axial direction.
        %   Timothy Sipkens, 2020-06-10
        %
        % Inputs:
        %   aso     Axis-symmetric object
        %   mu      Set of slopes for the rays
        %   y0      Intersect with line through center of aso
        %
        % Outputs:
        %   K       Overall radial deflection kernel
        %   Ky      Overall axial deflection kernel
        %           (often noisy due to uniform elements in axial direction)
        function [K,Kx] = linear(aso, my, y0, mx, x0)
            
            if aso.N<3; error('Aso does not have enough annuli for linear basis.'); end
            
            mx(abs(mx) < 1e-12) = 1e-12; % avoid division by zero in rv
            
            % edges of annuli and axial elements
            rjd0 = aso.re(1:(end-2));
            rj0  = aso.re(2:(end-1));
            rju0 = aso.re(3:end);
            xj  = aso.xe(1:(end-1));
            xju = aso.xe(2:end);
            
            % functions for indefinite integral
            A = @(mx,x0,ry,ryu,r3) 1 ./ (ryu - ry) .* ...
                log(abs(r3 + sqrt(r3.^2 - x0.^2 ./ (1+mx.^2)))); % main part of integrand
            Ad = @(mx,ryd,ry,r1,r2) (r2 - r1) ./ (ry - ryd) ...
                .* mx .* sqrt(1 + mx.^2); % second part of integrand
            
            
            N_beams = max([length(my), length(y0), ...
                length(mx), length(x0)]); % any of the values could be scalar, so take max
            K = spalloc(N_beams, aso.Nx * (aso.Nr + 1), ...
                round(1e-4 * N_beams * aso.Nx * (aso.Nr + 1)));
                % initialize K, assume 0.5% full
            
            Kx = spalloc(N_beams, aso.Nx * (aso.Nr + 1), ...
                round(1e-5 * N_beams * aso.Nx * (aso.Nr + 1)));
            
            
            disp('Looping through axial slices...');
            tools.textbar(0);
            for ii=1:(length(aso.xe)-1) % loop through and append axial slices
                
                % Find z intersection with axial elements.
                % This is used to determine whether to involve front of back
                % portion of integral.
                zu  = (xj(ii) - x0) ./ mx;
                zuu = (xju(ii) - x0) ./ mx;
                
                % flag if ray intersects elements at all
                f0 = or(or(or(and(zu > -aso.R, zu < aso.R), ... % if zv is in the bounds
                    and(zuu > -aso.R, zuu < aso.R)), ... % if zvu is in the bounds
                    and(zuu > 0, zu < 0)), ...
                    and(zuu < 0, zu > 0)); % if there is a change of sign
                idx_a = 1:length(f0);
                idx_a = idx_a(f0); % indices of rays that intersect element
                
                
                
                % Find radial intersection with axial elements.
                % This will be a strictly positive number.
                % This could be a bound for either term of the integrand.
                rx = sqrt(1 ./ (1+my(idx_a).^2) .* ...
                    ((1+my(idx_a).^2) ./ mx(idx_a) .* ...
                    (xj(ii) - x0(idx_a)) + my(idx_a).*y0(idx_a)) .^ 2 + ...
                    y0(idx_a).^2); % lower edge, vector over u0 and v0
                rxu = sqrt(1 ./ (1+my(idx_a).^2) .* ...
                    ((1+my(idx_a).^2) ./ mx(idx_a) .* ...
                    (xju(ii) - x0(idx_a)) + my(idx_a).*y0(idx_a)) .^ 2 + ...
                    y0(idx_a).^2); % upper edge, vector over u0 and v0
                
                
                
                %== FIRST INTEGRAND ======================================%
                % flag where intersect occurs
                fx  = zu(idx_a)  < (-my(idx_a) .* y0(idx_a) ...
                    ./ (1 + my(idx_a).^2)); % flags if intersect is in first integrand
                fxu = zuu(idx_a) < (-my(idx_a) .* y0(idx_a) ...
                    ./ (1 + my(idx_a).^2)); % flags if upper intersect is in first integrand
                
                % reverse if slope is negative (i.e. rvu is encountered first)
                r1 = rx; r1(mx(idx_a)<0) = rxu(mx(idx_a)<0);
                r2 = rxu; r2(mx(idx_a)<0) = rx(mx(idx_a)<0);
                f1 = fx; f1(mx(idx_a)<0) = fxu(mx(idx_a)<0);
                f2 = fxu; f2(mx(idx_a)<0) = fx(mx(idx_a)<0);
                
                % modified element widths
                % (adjusted for intersections with axial elements)
                rjd = rjd0 .* ones(size(rx)); % repeat for relevant rv elements
                rjd(:,f1) = min(r1(f1), rjd(:,f1)); % adjust for lower axial bound
                rjd(:,f2) = max(r2(f2), rjd(:,f2)); % adjust for upper axial bound
                
                rj = rj0 .* ones(size(rx));
                rj(:,f1) = min(r1(f1), rj(:,f1));
                rj(:,f2) = max(r2(f2), rj(:,f2));
                
                rju = rju0 .* ones(size(rx));
                rju(:,f1) = min(r1(f1), rju(:,f1));
                rju(:,f2) = max(r2(f2), rju(:,f2));
                
                % evaluate lower kernel (up to midpoint in ASO)
                K1 = (y0(idx_a) .* ( ...
                    or(fx,fxu) .* ([ ...
                     zeros(1,size(rj,2)); ...
                     A(my(idx_a), y0(idx_a), rjd0, rj0, rj) - ...
                     A(my(idx_a), y0(idx_a), rjd0, rj0, rjd) + ... % integral over rise
                     Ad(my(idx_a), rjd0, rj0, rjd, rj); ... % added component of integrand
                     A(my(idx_a), y0(idx_a), rj0(end,:), rju0(end,:), rju(end,:)) - ...
                     A(my(idx_a), y0(idx_a), rj0(end,:), rju0(end,:), rj(end,:)) + ...
                     Ad(my(idx_a), rj0(end,:), rju0(end,:), rj(end,:), rju(end,:)) ... % incline, last element
                    ] + [
                     A(my(idx_a), y0(idx_a), rj0(1,:), rjd0(1,:), rj(1,:)) - ...
                     A(my(idx_a), y0(idx_a), rj0(1,:), rjd0(1,:), rjd(1,:)) - ...
                     Ad(my(idx_a), rj0(1,:), rjd0(1,:), rj(1,:), rjd(1,:)); ...
                     A(my(idx_a), y0(idx_a), rju0, rj0, rju) - ...
                     A(my(idx_a), y0(idx_a), rju0, rj0, rj) - ...
                     Ad(my(idx_a), rju0, rj0, rju, rj); ... % integral over decline
                     zeros(1,size(rj,2))
                    ])))';
                
                
                
                %== SECOND INTEGRAND =====================================%
                % flag where intersect occurs
                fx  = ~fx; % flags if intersect is in first integrand
                fxu = ~fxu; % flags if upper intersect is in first integrand
                
                % reverse if slope is negative (i.e. rvu is encountered first)
                r1 = rxu; r1(mx(idx_a)<0) = rx(mx(idx_a)<0);
                r2 = rx; r2(mx(idx_a)<0) = rxu(mx(idx_a)<0);
                f1 = fxu; f1(mx(idx_a)<0) = fx(mx(idx_a)<0);
                f2 = fx; f2(mx(idx_a)<0) = fxu(mx(idx_a)<0);
                
                % modified element width
                % (adjusted for intersections with axial elements)
                rjd = rjd0 .* ones(size(rx)); % repeat for relevant rv elements
                rjd(:,f1) = min(r1(f1), rjd(:,f1)); % adjust for lower axial bound
                rjd(:,f2) = max(r2(f2), rjd(:,f2)); % adjust for upper axial bound
                
                rj = rj0 .* ones(size(rx));
                rj(:,f1) = min(r1(f1), rj(:,f1));
                rj(:,f2) = max(r2(f2), rj(:,f2));
                
                rju = rju0 .* ones(size(rx));
                rju(:,f1) = min(r1(f1), rju(:,f1));
                rju(:,f2) = max(r2(f2), rju(:,f2));
                
                % evaluate second integrand (beyond midpoint in ASO)
                K2 = (y0(idx_a) .* ( ...
                    or(fx,fxu) .* ([ ...
                     zeros(1, size(rj,2)); ...
                     A(my(idx_a), y0(idx_a), rjd0, rj0, rj) - ...
                     A(my(idx_a), y0(idx_a), rjd0, rj0, rjd) - ... % integral over rise
                     Ad(my(idx_a), rjd0,rj0, rjd, rj); ... % added component of integrand
                     A(my(idx_a), y0(idx_a), rj0(end,:), rju0(end,:), rju(end,:)) - ...
                     A(my(idx_a), y0(idx_a), rj0(end,:), rju0(end,:), rj(end,:)) - ...
                     Ad(my(idx_a), rj0(end,:), rju0(end,:), rj(end,:), rju(end,:)) ... % incline, last element
                    ] + [
                     A(my(idx_a), y0(idx_a), rj0(1,:), rjd0(1,:), rj(1,:)) - ...
                     A(my(idx_a), y0(idx_a), rj0(1,:), rjd0(1,:), rjd(1,:)) + ...
                     Ad(my(idx_a), rj0(1,:), rjd0(1,:), rj(1,:), rjd(1,:)); ...
                     A(my(idx_a), y0(idx_a), rju0, rj0, rju) - ...
                     A(my(idx_a), y0(idx_a), rju0, rj0, rj) + ...
                     Ad(my(idx_a), rju0, rj0, rju, rj); ... % integral over decline
                     zeros(1,size(rj,2))
                    ])))';
                
                
                K0 = K1 + K2;
                K0(isnan(K0)) = 0; % remove NaN values that result when modified element width is zero
                K0(abs(K0)<100*eps) = 0; % remove numerical noise
                K0 = sparse(K0);
                
                if any(any(K(idx_a, ((ii-1)*(aso.Nr+1)+1):(ii*(aso.Nr+1)))~=0))
                    disp('HELP');
                end
                K(idx_a, ((ii-1)*(aso.Nr+1)+1):(ii*(aso.Nr+1))) = K0;
                
                
                %== Evaluate axial deflections ===========================%
                %-{
                Kx0 = ([zeros(1, size(rj,2)); ...
                     (rx - rjd0) ./ (rj0 - rjd0) .* and(rjd0 < rx, rx < rj0) - ...
                     (rxu - rjd0) ./ (rj0 - rjd0) .* and(rjd0 < rxu, rxu < rj0); ...
                     (rx - rjd0(end,:)) ./ (rj0(end,:) - rjd0(end,:)) .* and(rjd0(end,:) < rx, rx < rj0(end,:)) - ...
                     (rxu - rjd0(end,:)) ./ (rj0(end,:) - rjd0(end,:)) .* and(rjd0(end,:) < rxu, rxu < rj0(end,:)) % if crosses in rise
                    ] + [(rx - rju0(1,:)) ./ (rj0(1,:) - rju0(1,:)) .* and(rju0(1,:) < rx, rx < rj0(1,:)) - ...
                     (rxu - rju0(1,:)) ./ (rj0(1,:) - rju0(1,:)) .* and(rju0(1,:) < rxu, rxu < rj0(1,:)); ...
                     (rx - rju0) ./ (rj0 - rju0) .* and(rju0 < rx, rx < rj0) - ...
                     (rxu - rju0) ./ (rj0 - rju0) .* and(rju0 < rxu, rxu < rj0); ...
                     zeros(1,size(rj,2)) ... % if crosses in decline
                    ])';
                Kx(idx_a, ((ii-1)*(aso.Nr+1)+1):(ii*(aso.Nr+1))) = sparse(Kx0);
                % + remove derivative below
                %}
                
                %{
                Kx0 = ([zeros(1, size(rj,2)); ...
                     rx .* (rx ./ 2 - rjd0) ./ (rj0 - rjd0) .* and(rjd0 < rx, rx < rj0) - ...
                     rxu .* (rxu ./ 2 - rjd0) ./ (rj0 - rjd0) .* and(rjd0 < rxu, rxu < rj0); ...
                     rx .* (rx ./ 2 - rjd0(end,:)) ./ (rj0(end,:) - rjd0(end,:)) .* and(rjd0(end,:) < rx, rx < rj0(end,:)) - ...
                     rxu .* (rxu ./ 2 - rjd0(end,:)) ./ (rj0(end,:) - rjd0(end,:)) .* and(rjd0(end,:) < rxu, rxu < rj0(end,:)) % if crosses in rise
                    ] + [rx .* (rx ./ 2 - rju0(1,:)) ./ (rj0(1,:) - rju0(1,:)) .* and(rju0(1,:) < rx, rx < rj0(1,:)) - ...
                     rxu .* (rxu ./ 2 - rju0(1,:)) ./ (rj0(1,:) - rju0(1,:)) .* and(rju0(1,:) < rxu, rxu < rj0(1,:)); ...
                     rx .* (rx ./ 2 - rju0) ./ (rj0 - rju0) .* and(rju0 < xy, rx < rj0) - ...
                     rxu .* (rxu ./ 2 - rju0) ./ (rj0 - rju0) .* and(rju0 < rxu, rxu < rj0); ...
                     zeros(1,size(rj,2)) ... % if crosses in decline
                    ])';
                Kx(idx_a, ((ii-1)*(aso.Nr+1)+1):(ii*(aso.Nr+1))) = sparse(Ky0);
                % + add derivative below
                %}
                
                tools.textbar(ii/(length(aso.xe)-1));
            end
            
            % Ky = Ky * aso.Dx;
            
            
        end
        
        
        
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
                plot3(x0', y0', z0', 'k');
                plot3(x0, y0, z0, 'k');
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
            
            if nargout==0; clear h; end % suppress output if none required
        end
    end
end

