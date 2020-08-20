
% ASO2 A class to handle spatial information for two-dimensional axis-symmetric objects (asos).
% Author: Timothy Sipkens, 2020-05-20
%=========================================================================%


classdef Aso2
    
    properties
        %-- Axial properties ---------------------------------------------%        
        dy    = [];     % axial element width
        y     = [];     % axial element centers
        ye    = [];     % axial element edges
        
        Y     = [];     % axial extent / overall width
        Ny    = [];     % number of edges in the axial direction
        %-----------------------------------------------------------------%
        
        
        %-- Radial properties (same as Aso) ------------------------------%
        dr    = [];     % annuli width
        r     = [];     % annuli centers
        re    = [];     % annuli edges
        
        R     = [];     % outer radius of object
        Nr    = [];     % number of annuli
        %-----------------------------------------------------------------%
        
        N     = [];     % total number of elements (Nr*Nv)
        edges = [];     % edges of all the elements (N x 4 array)
    end
    
    
    
    methods
        function aso = Aso2(R,Nr,Y,Ny)
            aso.Nr = Nr; % number of annuli
            aso.R = R; % outer radius
            
            aso.re = linspace(0, R, Nr+1)'; % linearily space edges from 0 -> R
            aso.r  = (aso.re(2:end) + aso.re(1:(end-1))) ./ 2; % annuli centers
            aso.dr = aso.re(2:end) - aso.re(1:(end-1)); % annuli width
            
            
            aso.Ny = Ny; % number of annuli
            aso.Y = Y; % outer radius
            
            aso.ye = linspace(0, Y, Ny+1)'; % linearily space edges from 0 -> R
            aso.y  = (aso.ye(2:end) + aso.ye(1:(end-1))) ./ 2; % annuli centers
            aso.dy = aso.ye(2:end) - aso.ye(1:(end-1)); % annuli width
            
            
            aso.N = aso.Ny * aso.Nr; % total number of elements
            aso.edges = [repmat(aso.re(1:(end-1)),[Ny,1]), ...
                repmat(aso.re(2:end),[Ny,1]), ...
                reshape(repmat(aso.ye(1:(end-1)),[1,Nr]), [Ny*Nr,1]), ...
                reshape(repmat(aso.ye(2:end),[1,Nr]), [Ny*Nr,1])];
                % vectorized element edges for overall grid
        end
        
        
        
        %== RESHAPE ======================================================%
        %   Timothy Sipkens, 2020-06-11
        function x = reshape(aso,x)
            x = reshape(x,[aso.Nr + 1,aso.Ny]);
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
        %   u0      Intersect with line through center of aso
        %
        % Outputs:
        %   K       Overall radial deflection kernel
        %   Kv      Overall axial deflection kernel
        %           (often noisy due to uniform elements in axial direction)
        function [K,Kv] = linear(aso,mu,u0,mv,v0)
            
            if aso.N<3; error('Aso does not have enough annuli for linear basis.'); end
            
            mv(abs(mv)<1e-12) = 1e-12; % avoid division by zero in rv
            
            % edges of annuli and axial elements
            rjd0 = aso.re(1:(end-2));
            rj0  = aso.re(2:(end-1));
            rju0 = aso.re(3:end);
            vj  = aso.ye(1:(end-1));
            vju = aso.ye(2:end);
            
            % functions for indefinite integral
            A = @(mu,u0,rv,rvu,r3) 1 ./ (rvu - rv) .* ...
                log(abs(r3 + sqrt(r3.^2 - u0.^2 ./ (1+mu.^2)))); % main part of integrand
            Ad = @(mu,rvd,rv,r1,r2) (r2 - r1) ./ (rv - rvd) ...
                .* mu .* sqrt(1 + mu.^2); % second part of integrand
            
            
            N_beams = max([length(mu), length(u0), ...
                length(mv), length(v0)]); % any of the values could be scalar, so take max
            K = spalloc(N_beams, aso.Ny * (aso.Nr + 1), ...
                round(1e-4 * N_beams * aso.Ny * (aso.Nr + 1)));
                % initialize K, assume 0.5% full
            
            Kv = spalloc(N_beams, aso.Ny * (aso.Nr + 1), ...
                round(1e-5 * N_beams * aso.Ny * (aso.Nr + 1)));
            
            
            disp('Looping through axial slices...');
            tools.textbar(0);
            for ii=1:(length(aso.ye)-1) % loop through and append axial slices
                
                % Find z intersection with axial elements.
                % This is used to determine whether to involve front of back
                % portion of integral.
                zv  = (vj(ii) - v0) ./ mv;
                zvu = (vju(ii) - v0) ./ mv;
                
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
                rv = sqrt(1 ./ (1+mu(idx_a).^2) .* ...
                    ((1+mu(idx_a).^2) ./ mv(idx_a) .* ...
                    (vj(ii) - v0(idx_a)) + mu(idx_a).*u0(idx_a)) .^ 2 + ...
                    u0(idx_a).^2); % lower edge, vector over u0 and v0
                rvu = sqrt(1 ./ (1+mu(idx_a).^2) .* ...
                    ((1+mu(idx_a).^2) ./ mv(idx_a) .* ...
                    (vju(ii) - v0(idx_a)) + mu(idx_a).*u0(idx_a)) .^ 2 + ...
                    u0(idx_a).^2); % upper edge, vector over u0 and v0
                
                
                
                %== FIRST INTEGRAND ======================================%
                % flag where intersect occurs
                fv  = zv(idx_a)  < (-mu(idx_a) .* u0(idx_a) ...
                    ./ (1 + mu(idx_a).^2)); % flags if intersect is in first integrand
                fvu = zvu(idx_a) < (-mu(idx_a) .* u0(idx_a) ...
                    ./ (1 + mu(idx_a).^2)); % flags if upper intersect is in first integrand
                
                % reverse if slope is negative (i.e. rvu is encountered first)
                r1 = rv; r1(mv(idx_a)<0) = rvu(mv(idx_a)<0);
                r2 = rvu; r2(mv(idx_a)<0) = rv(mv(idx_a)<0);
                f1 = fv; f1(mv(idx_a)<0) = fvu(mv(idx_a)<0);
                f2 = fvu; f2(mv(idx_a)<0) = fv(mv(idx_a)<0);
                
                % modified element widths
                % (adjusted for intersections with axial elements)
                rjd = rjd0 .* ones(size(rv)); % repeat for relevant rv elements
                rjd(:,f1) = min(r1(f1), rjd(:,f1)); % adjust for lower axial bound
                rjd(:,f2) = max(r2(f2), rjd(:,f2)); % adjust for upper axial bound
                
                rj = rj0 .* ones(size(rv));
                rj(:,f1) = min(r1(f1), rj(:,f1));
                rj(:,f2) = max(r2(f2), rj(:,f2));
                
                rju = rju0 .* ones(size(rv));
                rju(:,f1) = min(r1(f1), rju(:,f1));
                rju(:,f2) = max(r2(f2), rju(:,f2));
                
                % evaluate lower kernel (up to midpoint in ASO)
                K1 = (u0(idx_a) .* ( ...
                    or(fv,fvu) .* ([ ...
                     zeros(1,size(rj,2)); ...
                     A(mu(idx_a), u0(idx_a), rjd0, rj0, rj) - ...
                     A(mu(idx_a), u0(idx_a), rjd0, rj0, rjd) + ... % integral over rise
                     Ad(mu(idx_a), rjd0, rj0, rjd, rj); ... % added component of integrand
                     A(mu(idx_a), u0(idx_a), rj0(end,:), rju0(end,:), rju(end,:)) - ...
                     A(mu(idx_a), u0(idx_a), rj0(end,:), rju0(end,:), rj(end,:)) + ...
                     Ad(mu(idx_a), rj0(end,:), rju0(end,:), rj(end,:), rju(end,:)) ... % incline, last element
                    ] + [
                     A(mu(idx_a), u0(idx_a), rj0(1,:), rjd0(1,:), rj(1,:)) - ...
                     A(mu(idx_a), u0(idx_a), rj0(1,:), rjd0(1,:), rjd(1,:)) - ...
                     Ad(mu(idx_a), rj0(1,:), rjd0(1,:), rj(1,:), rjd(1,:)); ...
                     A(mu(idx_a), u0(idx_a), rju0, rj0, rju) - ...
                     A(mu(idx_a), u0(idx_a), rju0, rj0, rj) - ...
                     Ad(mu(idx_a), rju0, rj0, rju, rj); ... % integral over decline
                     zeros(1,size(rj,2))
                    ])))';
                
                
                
                %== SECOND INTEGRAND =====================================%
                % flag where intersect occurs
                fv  = ~fv; % flags if intersect is in first integrand
                fvu = ~fvu; % flags if upper intersect is in first integrand
                
                % reverse if slope is negative (i.e. rvu is encountered first)
                r1 = rvu; r1(mv(idx_a)<0) = rv(mv(idx_a)<0);
                r2 = rv; r2(mv(idx_a)<0) = rvu(mv(idx_a)<0);
                f1 = fvu; f1(mv(idx_a)<0) = fv(mv(idx_a)<0);
                f2 = fv; f2(mv(idx_a)<0) = fvu(mv(idx_a)<0);
                
                % modified element width
                % (adjusted for intersections with axial elements)
                rjd = rjd0 .* ones(size(rv)); % repeat for relevant rv elements
                rjd(:,f1) = min(r1(f1), rjd(:,f1)); % adjust for lower axial bound
                rjd(:,f2) = max(r2(f2), rjd(:,f2)); % adjust for upper axial bound
                
                rj = rj0 .* ones(size(rv));
                rj(:,f1) = min(r1(f1), rj(:,f1));
                rj(:,f2) = max(r2(f2), rj(:,f2));
                
                rju = rju0 .* ones(size(rv));
                rju(:,f1) = min(r1(f1), rju(:,f1));
                rju(:,f2) = max(r2(f2), rju(:,f2));
                
                % evaluate second integrand (beyond midpoint in ASO)
                K2 = (u0(idx_a) .* ( ...
                    or(fv,fvu) .* ([ ...
                     zeros(1, size(rj,2)); ...
                     A(mu(idx_a), u0(idx_a), rjd0, rj0, rj) - ...
                     A(mu(idx_a), u0(idx_a), rjd0, rj0, rjd) - ... % integral over rise
                     Ad(mu(idx_a), rjd0,rj0, rjd, rj); ... % added component of integrand
                     A(mu(idx_a), u0(idx_a), rj0(end,:), rju0(end,:), rju(end,:)) - ...
                     A(mu(idx_a), u0(idx_a), rj0(end,:), rju0(end,:), rj(end,:)) - ...
                     Ad(mu(idx_a), rj0(end,:), rju0(end,:), rj(end,:), rju(end,:)) ... % incline, last element
                    ] + [
                     A(mu(idx_a), u0(idx_a), rj0(1,:), rjd0(1,:), rj(1,:)) - ...
                     A(mu(idx_a), u0(idx_a), rj0(1,:), rjd0(1,:), rjd(1,:)) + ...
                     Ad(mu(idx_a), rj0(1,:), rjd0(1,:), rj(1,:), rjd(1,:)); ...
                     A(mu(idx_a), u0(idx_a), rju0, rj0, rju) - ...
                     A(mu(idx_a), u0(idx_a), rju0, rj0, rj) + ...
                     Ad(mu(idx_a), rju0, rj0, rju, rj); ... % integral over decline
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
                Kv0 = ([zeros(1, size(rj,2)); ...
                     (rv - rjd0) ./ (rj0 - rjd0) .* and(rjd0 < rv, rv < rj0) - ...
                     (rvu - rjd0) ./ (rj0 - rjd0) .* and(rjd0 < rvu, rvu < rj0); ...
                     (rv - rjd0(end,:)) ./ (rj0(end,:) - rjd0(end,:)) .* and(rjd0(end,:) < rv, rv < rj0(end,:)) - ...
                     (rvu - rjd0(end,:)) ./ (rj0(end,:) - rjd0(end,:)) .* and(rjd0(end,:) < rvu, rvu < rj0(end,:)) % if crosses in rise
                    ] + [(rv - rju0(1,:)) ./ (rj0(1,:) - rju0(1,:)) .* and(rju0(1,:) < rv, rv < rj0(1,:)) - ...
                     (rvu - rju0(1,:)) ./ (rj0(1,:) - rju0(1,:)) .* and(rju0(1,:) < rvu, rvu < rj0(1,:)); ...
                     (rv - rju0) ./ (rj0 - rju0) .* and(rju0 < rv, rv < rj0) - ...
                     (rvu - rju0) ./ (rj0 - rju0) .* and(rju0 < rvu, rvu < rj0); ...
                     zeros(1,size(rj,2)) ... % if crosses in decline
                    ])';
                Kv(idx_a, ((ii-1)*(aso.Nr+1)+1):(ii*(aso.Nr+1))) = sparse(Kv0);
                
                
                tools.textbar(ii/(length(aso.ye)-1));
            end
            
            
            
            
            
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

