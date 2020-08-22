
% ASO2  A class to handle spatial information for 2D axisymmetric objects.
% Author: Timothy Sipkens, 2020-05-20
%=========================================================================%

classdef Aso2
    
    properties
        %-- Axial properties ---------------------------------------------%        
        dy    = [];     % axial element width
        y     = [];     % axial element centers
        ye    = [];     % axial element edges
        ye2   = [];     % axial edges (for all elements in 2D grid),
                        % to be used with 2D phantoms defined on 
                        % a linear basis for radius and uniform basis for
                        % axial position
        
        Y     = [];     % axial extent / overall width
        Ny    = [];     % number of edges in the axial direction
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
        Dy    = [];     % y-gradient operator
    end
    
    
    
    methods
        function aso = Aso2(R,Nr,Y,Ny)
            %-- Radial positions and discretization ----------------------%
            aso.Nr = Nr; % number of annuli
            aso.R = R; % outer radius
            
            aso.re = linspace(0, R, Nr+1)'; % linearily space edges from 0 -> R
            aso.r  = (aso.re(2:end) + aso.re(1:(end-1))) ./ 2; % annuli centers
            aso.dr = aso.re(2:end) - aso.re(1:(end-1)); % annuli width
            %-------------------------------------------------------------%
            
            %-- Axial positions and discretization -----------------------%
            aso.Ny = Ny; % number of annuli
            aso.Y = Y; % outer radius
            
            aso.ye = linspace(0, Y, Ny+1)'; % linearily space edges from 0 -> R
            aso.y  = (aso.ye(2:end) + aso.ye(1:(end-1))) ./ 2; % annuli centers
            aso.dy = aso.ye(2:end) - aso.ye(1:(end-1)); % annuli width
            %-------------------------------------------------------------%
            
            %-- Consolidated element information -------------------------%
            aso.N = aso.Ny * aso.Nr; % total number of elements
            
            % vector list of element edges (for 2D linear basis)
            [aso.ye2, aso.re2] = meshgrid(aso.ye(1:(end-1)), aso.re);
            aso.ye2 = aso.ye2(:); aso.re2 = aso.re2(:);
            
            % vectorized element edges for overall grid
            aso.edges = [repmat(aso.re(1:(end-1)),[Ny,1]), ...
                repmat(aso.re(2:end),[Ny,1]), ...
                reshape(repmat(aso.ye(1:(end-1)),[1,Nr]), [Ny*Nr,1]), ...
                reshape(repmat(aso.ye(2:end),[1,Nr]), [Ny*Nr,1])];
            %-------------------------------------------------------------%
            
            %-- Compute differential operators ---------------------------%
            I1 = speye(aso.Nr+1, aso.Nr+1);
            E1 = sparse(1:aso.Nr+1-1,2:aso.Nr+1,1,aso.Nr+1,aso.Nr+1);
            D1 = E1-I1;
            D1(end,end) = 0; % applies no slope at final radial position
            
            I2 = speye(aso.Ny,aso.Ny);
            E2 = sparse(1:aso.Ny-1,2:aso.Ny,1,aso.Ny,aso.Ny);
            D2 = E2-I2;
            D2(end,end) = 1; D2(end,end-1) = -1; % applies same slope as previous y-position
            
            aso.Dr = kron(I2,D1);
            aso.Dy = kron(D2,I1);
            %-------------------------------------------------------------%
            
            
        end
        
        
        
        %== RESHAPE ======================================================%
        %   Reshape, assuming x is evaluated at aso.re.
        %   Timothy Sipkens, 2020-06-11
        function x = reshape(aso,x)
            x = reshape(x, [aso.Nr + 1, aso.Ny]);
        end
        %=================================================================%
        
        
        
        %== GRADIENTC ====================================================%
        %   Cartesian gradients, interpolated from the ASO object.
        %   Assumes a linear radial basis and uniform axial basis.
        %   Timothy Sipkens, 2020-06-11
        function [Dx, Dy, Dz] = gradientc(aso,xi,yi,zy,f)
            
            % Get position in cylindrical coordinates
            ri = sqrt(xi.^2 + zy.^2); % radial position
            q = atan2(xi, -zy); % angle, in coord. system from Sipkens et al.
            
            % Setup grid for interpolation
            [y0, r0] = meshgrid(aso.ye, aso.re); % grid for input to interpolation
            Dri = aso.reshape(aso.Dr * f); % reshape radial gradient
            Dri = [Dri, Dri(:,end)]; % append constant slope data for last axial position
            Dyi = aso.reshape(aso.Dy * f); % reshape axial gradient
            Dyi = [Dyi, Dyi(:,end)]; % append constant slope data for last axial position
            
            % Interpolate r-gradient and convert to Cartesian coords.
            Dro = interp2(y0, r0, Dri, yi, ri, 'linear' ,0);
            Dx = sin(q) .* Dro; % get the x-gradient based on the angle
            Dz = -cos(q) .* Dro; % get the y-gradient based on the angle
            
            % Interpolate y-gradient
            Dy = interp2(y0, r0, Dyi, yi, ri, 'linear', 0);
            
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
        function [K,Ky] = linear(aso,mx,x0,my,y0)
            
            if aso.N<3; error('Aso does not have enough annuli for linear basis.'); end
            
            my(abs(my) < 1e-12) = 1e-12; % avoid division by zero in rv
            
            % edges of annuli and axial elements
            rjd0 = aso.re(1:(end-2));
            rj0  = aso.re(2:(end-1));
            rju0 = aso.re(3:end);
            yj  = aso.ye(1:(end-1));
            yju = aso.ye(2:end);
            
            % functions for indefinite integral
            A = @(mx,x0,ry,ryu,r3) 1 ./ (ryu - ry) .* ...
                log(abs(r3 + sqrt(r3.^2 - x0.^2 ./ (1+mx.^2)))); % main part of integrand
            Ad = @(mx,ryd,ry,r1,r2) (r2 - r1) ./ (ry - ryd) ...
                .* mx .* sqrt(1 + mx.^2); % second part of integrand
            
            
            N_beams = max([length(mx), length(x0), ...
                length(my), length(y0)]); % any of the values could be scalar, so take max
            K = spalloc(N_beams, aso.Ny * (aso.Nr + 1), ...
                round(1e-4 * N_beams * aso.Ny * (aso.Nr + 1)));
                % initialize K, assume 0.5% full
            
            Ky = spalloc(N_beams, aso.Ny * (aso.Nr + 1), ...
                round(1e-5 * N_beams * aso.Ny * (aso.Nr + 1)));
            
            
            disp('Looping through axial slices...');
            tools.textbar(0);
            for ii=1:(length(aso.ye)-1) % loop through and append axial slices
                
                % Find z intersection with axial elements.
                % This is used to determine whether to involve front of back
                % portion of integral.
                zv  = (yj(ii) - y0) ./ my;
                zvu = (yju(ii) - y0) ./ my;
                
                % flag if ray intersects elements at all
                f0 = or(or(or(and(zv > -aso.R, zv < aso.R), ... % if zv is in the bounds
                    and(zvu > -aso.R, zvu < aso.R)), ... % if zvu is in the bounds
                    and(zvu > 0, zv < 0)), ...
                    and(zvu < 0, zv > 0)); % if there is a change of sign
                idx_a = 1:length(f0);
                idx_a = idx_a(f0); % indices of rays that intersect element
                
                
                
                % Find radial intersection with axial elements.
                % This will be a strictly positive number.
                % This could be a bound for either term of the integrand.
                ry = sqrt(1 ./ (1+mx(idx_a).^2) .* ...
                    ((1+mx(idx_a).^2) ./ my(idx_a) .* ...
                    (yj(ii) - y0(idx_a)) + mx(idx_a).*x0(idx_a)) .^ 2 + ...
                    x0(idx_a).^2); % lower edge, vector over u0 and v0
                ryu = sqrt(1 ./ (1+mx(idx_a).^2) .* ...
                    ((1+mx(idx_a).^2) ./ my(idx_a) .* ...
                    (yju(ii) - y0(idx_a)) + mx(idx_a).*x0(idx_a)) .^ 2 + ...
                    x0(idx_a).^2); % upper edge, vector over u0 and v0
                
                
                
                %== FIRST INTEGRAND ======================================%
                % flag where intersect occurs
                fy  = zv(idx_a)  < (-mx(idx_a) .* x0(idx_a) ...
                    ./ (1 + mx(idx_a).^2)); % flags if intersect is in first integrand
                fyu = zvu(idx_a) < (-mx(idx_a) .* x0(idx_a) ...
                    ./ (1 + mx(idx_a).^2)); % flags if upper intersect is in first integrand
                
                % reverse if slope is negative (i.e. rvu is encountered first)
                r1 = ry; r1(my(idx_a)<0) = ryu(my(idx_a)<0);
                r2 = ryu; r2(my(idx_a)<0) = ry(my(idx_a)<0);
                f1 = fy; f1(my(idx_a)<0) = fyu(my(idx_a)<0);
                f2 = fyu; f2(my(idx_a)<0) = fy(my(idx_a)<0);
                
                % modified element widths
                % (adjusted for intersections with axial elements)
                rjd = rjd0 .* ones(size(ry)); % repeat for relevant rv elements
                rjd(:,f1) = min(r1(f1), rjd(:,f1)); % adjust for lower axial bound
                rjd(:,f2) = max(r2(f2), rjd(:,f2)); % adjust for upper axial bound
                
                rj = rj0 .* ones(size(ry));
                rj(:,f1) = min(r1(f1), rj(:,f1));
                rj(:,f2) = max(r2(f2), rj(:,f2));
                
                rju = rju0 .* ones(size(ry));
                rju(:,f1) = min(r1(f1), rju(:,f1));
                rju(:,f2) = max(r2(f2), rju(:,f2));
                
                % evaluate lower kernel (up to midpoint in ASO)
                K1 = (x0(idx_a) .* ( ...
                    or(fy,fyu) .* ([ ...
                     zeros(1,size(rj,2)); ...
                     A(mx(idx_a), x0(idx_a), rjd0, rj0, rj) - ...
                     A(mx(idx_a), x0(idx_a), rjd0, rj0, rjd) + ... % integral over rise
                     Ad(mx(idx_a), rjd0, rj0, rjd, rj); ... % added component of integrand
                     A(mx(idx_a), x0(idx_a), rj0(end,:), rju0(end,:), rju(end,:)) - ...
                     A(mx(idx_a), x0(idx_a), rj0(end,:), rju0(end,:), rj(end,:)) + ...
                     Ad(mx(idx_a), rj0(end,:), rju0(end,:), rj(end,:), rju(end,:)) ... % incline, last element
                    ] + [
                     A(mx(idx_a), x0(idx_a), rj0(1,:), rjd0(1,:), rj(1,:)) - ...
                     A(mx(idx_a), x0(idx_a), rj0(1,:), rjd0(1,:), rjd(1,:)) - ...
                     Ad(mx(idx_a), rj0(1,:), rjd0(1,:), rj(1,:), rjd(1,:)); ...
                     A(mx(idx_a), x0(idx_a), rju0, rj0, rju) - ...
                     A(mx(idx_a), x0(idx_a), rju0, rj0, rj) - ...
                     Ad(mx(idx_a), rju0, rj0, rju, rj); ... % integral over decline
                     zeros(1,size(rj,2))
                    ])))';
                
                
                
                %== SECOND INTEGRAND =====================================%
                % flag where intersect occurs
                fy  = ~fy; % flags if intersect is in first integrand
                fyu = ~fyu; % flags if upper intersect is in first integrand
                
                % reverse if slope is negative (i.e. rvu is encountered first)
                r1 = ryu; r1(my(idx_a)<0) = ry(my(idx_a)<0);
                r2 = ry; r2(my(idx_a)<0) = ryu(my(idx_a)<0);
                f1 = fyu; f1(my(idx_a)<0) = fy(my(idx_a)<0);
                f2 = fy; f2(my(idx_a)<0) = fyu(my(idx_a)<0);
                
                % modified element width
                % (adjusted for intersections with axial elements)
                rjd = rjd0 .* ones(size(ry)); % repeat for relevant rv elements
                rjd(:,f1) = min(r1(f1), rjd(:,f1)); % adjust for lower axial bound
                rjd(:,f2) = max(r2(f2), rjd(:,f2)); % adjust for upper axial bound
                
                rj = rj0 .* ones(size(ry));
                rj(:,f1) = min(r1(f1), rj(:,f1));
                rj(:,f2) = max(r2(f2), rj(:,f2));
                
                rju = rju0 .* ones(size(ry));
                rju(:,f1) = min(r1(f1), rju(:,f1));
                rju(:,f2) = max(r2(f2), rju(:,f2));
                
                % evaluate second integrand (beyond midpoint in ASO)
                K2 = (x0(idx_a) .* ( ...
                    or(fy,fyu) .* ([ ...
                     zeros(1, size(rj,2)); ...
                     A(mx(idx_a), x0(idx_a), rjd0, rj0, rj) - ...
                     A(mx(idx_a), x0(idx_a), rjd0, rj0, rjd) - ... % integral over rise
                     Ad(mx(idx_a), rjd0,rj0, rjd, rj); ... % added component of integrand
                     A(mx(idx_a), x0(idx_a), rj0(end,:), rju0(end,:), rju(end,:)) - ...
                     A(mx(idx_a), x0(idx_a), rj0(end,:), rju0(end,:), rj(end,:)) - ...
                     Ad(mx(idx_a), rj0(end,:), rju0(end,:), rj(end,:), rju(end,:)) ... % incline, last element
                    ] + [
                     A(mx(idx_a), x0(idx_a), rj0(1,:), rjd0(1,:), rj(1,:)) - ...
                     A(mx(idx_a), x0(idx_a), rj0(1,:), rjd0(1,:), rjd(1,:)) + ...
                     Ad(mx(idx_a), rj0(1,:), rjd0(1,:), rj(1,:), rjd(1,:)); ...
                     A(mx(idx_a), x0(idx_a), rju0, rj0, rju) - ...
                     A(mx(idx_a), x0(idx_a), rju0, rj0, rj) + ...
                     Ad(mx(idx_a), rju0, rj0, rju, rj); ... % integral over decline
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
                Ky0 = ([zeros(1, size(rj,2)); ...
                     (ry - rjd0) ./ (rj0 - rjd0) .* and(rjd0 < ry, ry < rj0) - ...
                     (ryu - rjd0) ./ (rj0 - rjd0) .* and(rjd0 < ryu, ryu < rj0); ...
                     (ry - rjd0(end,:)) ./ (rj0(end,:) - rjd0(end,:)) .* and(rjd0(end,:) < ry, ry < rj0(end,:)) - ...
                     (ryu - rjd0(end,:)) ./ (rj0(end,:) - rjd0(end,:)) .* and(rjd0(end,:) < ryu, ryu < rj0(end,:)) % if crosses in rise
                    ] + [(ry - rju0(1,:)) ./ (rj0(1,:) - rju0(1,:)) .* and(rju0(1,:) < ry, ry < rj0(1,:)) - ...
                     (ryu - rju0(1,:)) ./ (rj0(1,:) - rju0(1,:)) .* and(rju0(1,:) < ryu, ryu < rj0(1,:)); ...
                     (ry - rju0) ./ (rj0 - rju0) .* and(rju0 < ry, ry < rj0) - ...
                     (ryu - rju0) ./ (rj0 - rju0) .* and(rju0 < ryu, ryu < rj0); ...
                     zeros(1,size(rj,2)) ... % if crosses in decline
                    ])';
                Ky(idx_a, ((ii-1)*(aso.Nr+1)+1):(ii*(aso.Nr+1))) = sparse(Ky0);
                % + remove derivative below
                %}
                
                %{
                Ky0 = ([zeros(1, size(rj,2)); ...
                     ry .* (ry ./ 2 - rjd0) ./ (rj0 - rjd0) .* and(rjd0 < ry, ry < rj0) - ...
                     ryu .* (ryu ./ 2 - rjd0) ./ (rj0 - rjd0) .* and(rjd0 < ryu, ryu < rj0); ...
                     ry .* (ry ./ 2 - rjd0(end,:)) ./ (rj0(end,:) - rjd0(end,:)) .* and(rjd0(end,:) < ry, ry < rj0(end,:)) - ...
                     ryu .* (ryu ./ 2 - rjd0(end,:)) ./ (rj0(end,:) - rjd0(end,:)) .* and(rjd0(end,:) < ryu, ryu < rj0(end,:)) % if crosses in rise
                    ] + [ry .* (ry ./ 2 - rju0(1,:)) ./ (rj0(1,:) - rju0(1,:)) .* and(rju0(1,:) < ry, ry < rj0(1,:)) - ...
                     ryu .* (ryu ./ 2 - rju0(1,:)) ./ (rj0(1,:) - rju0(1,:)) .* and(rju0(1,:) < ryu, ryu < rj0(1,:)); ...
                     ry .* (ry ./ 2 - rju0) ./ (rj0 - rju0) .* and(rju0 < ry, ry < rj0) - ...
                     ryu .* (ryu ./ 2 - rju0) ./ (rj0 - rju0) .* and(rju0 < ryu, ryu < rj0); ...
                     zeros(1,size(rj,2)) ... % if crosses in decline
                    ])';
                Ky(idx_a, ((ii-1)*(aso.Nr+1)+1):(ii*(aso.Nr+1))) = sparse(Ky0);
                % + add derivative below
                %}
                
                tools.textbar(ii/(length(aso.ye)-1));
            end
            
            % Ky = Ky * aso.Dy;
            
            
        end
        
        
        
        %== SURF =========================================================%
        %   Plot the axis-symmetric object as a surface. 
        %   Timothy Sipkens, 2020-06-09
        function [h,x0,y0,z0] = surf(aso,x,f_grid)
            
            if ~exist('f_grid','var'); f_grid = []; end
            if isempty(f_grid)
                if x>1e3; f_grid = 1; else; f_grid = 0; end
            end
            
            [iv,ir] = meshgrid(1:aso.Ny, 1:(aso.Nr+1));
            
            x1 = aso.ye(iv);
            y1 = aso.re(ir);
            
            x0 = [flipud(x1);  x1(2:end,:)]; % plot axial positions
            y0 = [-flipud(y1); y1(2:end,:)]; % plot radii
            
            z1 = reshape(x, [aso.Nr+1, aso.Ny]);
            z0 = [flipud(z1); z1(2:end,:)];
            
            h = surf(x0,y0,z0);
            h.EdgeColor = 'none';
            
            if f_grid % overlay element grid
                hold on;
                plot3(x0', y0', z0', 'k');
                plot3(x0, y0, z0, 'k');
                hold off;
            end
            
            xlabel('Axial, v');
            ylabel('Radial, r');
            
            if nargout==0; clear h; end % suppress output if none required
            
        end
        
        
        
        %== PLOT =========================================================%
        %   Plot the axis-symmetric object as a surface. 
        %   Timothy Sipkens, 2020-06-09
        function [h,x0,y0,z0] = plot(aso,x,f_grid)
            
            if ~exist('f_grid','var'); f_grid = []; end
            if isempty(f_grid)
                if x>1e3; f_grid = 1; else; f_grid = 0; end
            end
            
            [iv,ir] = meshgrid(1:aso.Ny, 1:(aso.Nr+1));
            
            x1 = aso.ye(iv);
            y1 = aso.re(ir);
            
            x0 = [flipud(x1);  x1(2:end,:)]; % plot axial positions
            y0 = [-flipud(y1); y1(2:end,:)]; % plot radii
            
            z1 = reshape(x, [aso.Nr+1, aso.Ny]);
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
            
            xlabel('Axial, v');
            ylabel('Radial, r');
            
            if nargout==0; clear h; end % suppress output if none required
            
        end
        
        
        
        %== SRAYS ========================================================%
        %   Plot aso as a surface, with rays overlaid.
        %   Timothy Sipkens, 2020-06-09
        %   Note: Works best with monotonically increasing/decreasing z0.
        function h = srays(aso,x,mv,v0,f_grid)
            
            if ~exist('f_grid','var'); f_grid = []; end
            
            [h,x0,y0,z0] = aso.surf(x,f_grid); % generate surface plot
            
            y1 = linspace(-aso.R, aso.R, 150);
            x1 = (mv(1:10:end)').*y1 + v0(1:10:end)';
            
            y1 = y1.*ones(size(x1)); % if necessary, will out y1 to have correct dimension
            
            z1 = interp2(x0, y0, z0, x1, y1, 'linear'); 
            hold on;
            plot3(x1', y1', z1', 'r');
            hold off;
            
            if nargout==0; clear h; end % suppress output if none required
        end
    end
end

