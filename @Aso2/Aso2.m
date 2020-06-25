
% ASO2 A class to handle spatial information for two-dimensional axis-symmetric objects (asos).
% Author: Timothy Sipkens, 2020-05-20
%=========================================================================%


classdef Aso2
    
    properties
        %-- Axial properties ---------------------------------------------%        
        dv    = [];     % axial element width
        v     = [];     % axial element centers
        ve    = [];     % axial element edges
        
        V     = [];     % axial extent / overall width
        Nv    = [];     % number of edges in the axial direction
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
        function aso = Aso2(R,Nr,V,Nv)
            aso.Nr = Nr; % number of annuli
            aso.R = R; % outer radius
            
            aso.re = linspace(0, R, Nr+1)'; % linearily space edges from 0 -> R
            aso.r  = (aso.re(2:end) + aso.re(1:(end-1))) ./ 2; % annuli centers
            aso.dr = aso.re(2:end) - aso.re(1:(end-1)); % annuli width
            
            
            aso.Nv = Nv; % number of annuli
            aso.V = V; % outer radius
            
            aso.ve = linspace(0, V, Nv+1)'; % linearily space edges from 0 -> R
            aso.v  = (aso.ve(2:end) + aso.ve(1:(end-1))) ./ 2; % annuli centers
            aso.dv = aso.ve(2:end) - aso.ve(1:(end-1)); % annuli width
            
            
            aso.N = aso.Nv * aso.Nr; % total number of elements
            aso.edges = [repmat(aso.re(1:(end-1)),[Nv,1]), ...
                repmat(aso.re(2:end),[Nv,1]), ...
                reshape(repmat(aso.ve(1:(end-1)),[1,Nr]), [Nv*Nr,1]), ...
                reshape(repmat(aso.ve(2:end),[1,Nr]), [Nv*Nr,1])];
                % vectorized element edges for overall grid
        end
        
        
        
        %== RESHAPE ======================================================%
        %   Timothy Sipkens, 2020-06-11
        function x = reshape(aso,x)
            x = reshape(x,[aso.Nr + 1,aso.Nv]);
        end
        %=================================================================%
        
        
        
        %== LINEAR =======================================================%
        %   Evaluates kernel/operator for a linear basis representation of an ASO.
        %   Timothy Sipkens, 2020-06-10
        %
        % Inputs:
        %   aso     Axis-symmetric object
        %   mu      Set of slopes for the rays
        %   u0      Intersect with line through center of aso
        function K = linear(aso,mu,u0,mv,v0)
            
            if aso.N<3; error('Aso does not have enough annuli for linear basis.'); end
            
            mv(abs(mv)<1e-12) = 1e-12; % avoid division by zero in rv
            
            % edges of annuli and axial elements
            rjd0 = aso.re(1:(end-2));
            rj0  = aso.re(2:(end-1));
            rju0 = aso.re(3:end);
            vj  = aso.ve(1:(end-1));
            vju = aso.ve(2:end);
            
            % functions for indefinite integral
            Kb = @(mu,u0,r) log(r + sqrt(r.^2 - u0.^2 ./ (1+mu.^2)));
            Kc = @(mu,u0,r1,r2,r3) 1 ./ (r2 - r1) .* (Kb(mu,u0,r3));
            K1 = @(mu,r1,r2,r3,r4) r4 .* (r4 - r3) ./ (r2 - r1) .* ...
                 (mu .* sqrt(1+mu.^2));
            
             
            N_beams = max([length(mu), length(u0), ...
                length(mv), length(v0)]); % any of the values could be scalar, so take max
            K = spalloc(N_beams, aso.Nv * (aso.Nr + 1), ...
                round(0.05 * N_beams * aso.Nv * (aso.Nr + 1)));
                % initialize K, assume 5% full
            
            
            disp('Looping through axial slices...');
            tools.textbar(0);
            for ii=1:(length(aso.ve)-1) % loop through and append axial slices
                
                % Find radial intersection with axial elements.
                % This could be a bound for either term of the integrand.
                rv = sqrt(1 ./ (1+mu.^2) .* ...
                    ((1+mu.^2) ./ mv .* (vj(ii) - v0) + mu.*u0) .^ 2 + ...
                    u0.^2); % lower edge
                rvu = sqrt(1 ./ (1+mu.^2) .* ...
                    ((1+mu.^2) ./ mv .* (vju(ii) - v0) + mu.*u0) .^ 2 + ...
                    u0.^2); % upper edge
                
                
                % Find z intersection with axial elements.
                % This is used to determine whether to involve front of back
                % portion of integral.
                zv  = (vj(ii) - v0) ./ abs(mv);
                zvu = (vju(ii) - v0) ./ abs(mv);
                fv  = zv <= (abs(mu) .* u0); % flags if intersect is in first integrand
                fvu = zvu < (abs(mu) .* u0); % flags if intersect (upper) is in first integrand
                
                
                % modified element width
                % (adjusted for intersections with axial elements)
                rjd = rjd0 .* ones(size(rv));
                rjd(:,fv) = min(rv(fv), rjd(:,fv)); % adjust for lower axial bound
                rjd(:,fvu) = max(rvu(fvu), rjd(:,fvu)); % adjust for upper axial bound
                
                rj = rj0 .* ones(size(rv));
                rj(:,fv) = min(rv(fv), rj(:,fv));
                rj(:,fvu) = max(rvu(fvu), rj(:,fvu));
                
                rju = rju0 .* ones(size(rv));
                rju(:,fv) = min(rv(fv), rju(:,fv));
                rju(:,fvu) = max(rvu(fvu), rju(:,fvu));
                
                % evaluate lower kernel
                Kd = real(u0 .* ( ... % real(.) removes values outside integral bounds
                    or(fv,fvu) .* ([ ...
                     zeros(1,size(rj,2)); ...
                     Kc(mu,u0,rjd0,rj0,rj) - ...
                     Kc(mu,u0,rjd0,rj0,rjd) + ...
                     K1(mu,rjd0,rj0,rjd,rj); ... % integral over rise
                     Kc(mu,u0,rj0(end,:),rju0(end,:),rju(end,:)) - ...
                     Kc(mu,u0,rj0(end,:),rju0(end,:),rj(end,:)) + ...
                     K1(mu,rj0(end,:),rju0(end,:),rj(end,:),rju(end,:)) ... % incline, last element
                    ] + [
                     Kc(mu,u0,rj0(1,:),rjd0(1,:),rj(1,:)) - ...
                     Kc(mu,u0,rj0(1,:),rjd0(1,:),rjd(1,:)) + ...
                     K1(mu,rj0(1,:),rjd0(1,:),rj(1,:),rjd(1,:)); ... % decline, first element
                     Kc(mu,u0,rju0,rj0,rju) - ...
                     Kc(mu,u0,rju0,rj0,rj) + ...
                     K1(mu,rju0,rj0,rju,rj); ... % integral over decline
                     zeros(1,size(rj,2))
                    ])))';
                
                
                % modified element width
                % (adjusted for intersections with axial elements)
                rjd = rjd0 .* ones(size(rv));
                rjd(:,~fv) = max(rv(~fv), rjd(:,~fv)); % adjust for lower axial bound
                rjd(:,~fvu) = min(rvu(~fvu), rjd(:,~fvu)); % adjust for upper axial bound
                
                rj = rj0 .* ones(size(rv));
                rj(:,~fv) = max(rv(~fv), rj(:,~fv));
                rj(:,~fvu) = min(rvu(~fvu), rj(:,~fvu));
                
                rju = rju0 .* ones(size(rv));
                rju(:,~fv) = max(rv(~fv), rju(:,~fv));
                rju(:,~fvu) = min(rvu(~fvu), rju(:,~fvu));
                
                % evaluate upper kernel
                Ku = real(u0 .* ( ... % real(.) removes values outside integral bounds
                    or(~fv,~fvu) .* ([ ...
                     zeros(1,size(rj,2)); ...
                     Kc(mu,u0,rjd0,rj0,rj) - ...
                     Kc(mu,u0,rjd0,rj0,rjd) - ...
                     K1(mu,rjd0,rj0,rjd,rj); ... % integral over rise
                     Kc(mu,u0,rj0(end,:),rju0(end,:),rju(end,:)) - ...
                     Kc(mu,u0,rj0(end,:),rju0(end,:),rj(end,:)) - ...
                     K1(mu,rj0(end,:),rju0(end,:),rj(end,:),rju(end,:)) ... % incline, last element
                    ] + [
                     Kc(mu,u0,rj0(1,:),rjd0(1,:),rj(1,:)) - ...
                     Kc(mu,u0,rj0(1,:),rjd0(1,:),rjd(1,:)) - ...
                     K1(mu,rj0(1,:),rjd0(1,:),rj(1,:),rjd(1,:)); ... % decline, first element
                     Kc(mu,u0,rju0,rj0,rju) - ...
                     Kc(mu,u0,rju0,rj0,rj) - ...
                     K1(mu,rju0,rj0,rju,rj); ... % integral over decline
                     zeros(1,size(rj,2))
                    ])))';
                
                
                K0 = Kd + Ku;
                K0(isnan(K0)) = 0; % remove NaN values that result when modified element width is zero
                K0(abs(K0)<100*eps) = 0; % remove numerical noise
                K0 = sparse(K0);
                
                K(:, ((ii-1)*(aso.Nr+1)+1):(ii*(aso.Nr+1))) = K0;
                
                tools.textbar(ii/(length(aso.ve)-1));
            end
            
            
        end
        
        
        
        %== SURF =========================================================%
        %   Plot the axis-symmetric object as a surface. 
        %   Timothy Sipkens, 2020-06-09
        function [h,x0,y0,z0] = surf(aso,x,f_grid)
            
            if ~exist('f_grid','var'); f_grid = []; end
            if isempty(f_grid); f_grid = 1; end
            
            [iv,ir] = meshgrid(1:aso.Nv, 1:(aso.Nr+1));
            
            x1 = aso.ve(iv);
            y1 = aso.re(ir);
            
            x0 = [flipud(x1);  x1(2:end,:)]; % plot axial positions
            y0 = [-flipud(y1); y1(2:end,:)]; % plot radii
            
            z1 = reshape(x, [aso.Nr+1, aso.Nv]);
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
        
        
        
        %== SRAYS ========================================================%
        %   Plot aso as a surface, with rays overlaid.
        %   Timothy Sipkens, 2020-06-09
        %   Note: Works best with monotonically increasing/decreasing z0.
        function h = srays(aso,x,mv,v0,f_grid)
            
            if ~exist('f_grid','var'); f_grid = []; end
            
            [h,x0,y0,z0] = aso.surf(x,f_grid); % generate surface plot
            
            y1 = linspace(-aso.R, aso.R, 150);
            x1 = (mv').*y1 + v0';
            
            y1 = y1.*ones(size(x1)); % if necessary, will out y1 to have correct dimension
            
            z1 = interp2(x0, y0, z0, x1, y1, 'linear'); 
            hold on;
            plot3(x1', y1', z1', 'r');
            hold off;
            
            if nargout==0; clear h; end % suppress output if none required
        end
    end
end

