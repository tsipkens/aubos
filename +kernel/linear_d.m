
% LINEAR_D  Evaluates kernel/operator for a linear basis representation of a 1D ASO.
% Timothy Sipkens, 2020-06-10
% 
% Inputs:
%   aso_re  Axis-symmetric object or edges of the annuli
%   my      Set of slopes for the rays
%   y0      Intersect with line through center of aso
%=========================================================================%

function [K, Kx] = linear_d(aso_re, y0, my, x0, mx)

%-- Parse inputs ---------------------------------------------------------%
if isa(aso_re,'Aso'); re = aso_re.re; % if input is an Aso
else; re = aso_re; end % if an input is edges
if isa(aso_re,'Aso2'); aso2 = aso_re; re = aso_re.re; end

if ~exist('y0', 'var'); y0 = []; end
if isempty(y0); y0 = re'; end

if ~exist('my', 'var'); my = []; end
if isempty(my); my = zeros(size(y0)); end
%-------------------------------------------------------------------------%


Nr = length(re) - 1;


if Nr<3; error(' Not have enough annuli for linear basis.'); end



if nargin<4  % consider 1D case
    % get range of rj
    rjd = re(1:(end-2)); % r_{j-1}
    rj  = re(2:(end-1)); % r_j
    rju = re(3:end);     % r_{j+1}

    % functions for indefinite integral
    Kb = @(m,y0,r) log(abs(r + sqrt(r.^2 - y0.^2 ./ (1+m.^2))));
    Kc = @(m,y0,r1,r2,r3) 1 ./ (r2 - r1) .* Kb(m, y0, r3 + eps);
        % the + eps allows for finite value of kernel when r3 = x0
    
    K = (2 ./ (1 + my .^ 2) .* y0 .* ( ... % real(.) removes values outside integral bounds
        [ ...
         zeros(1,max(length(my),length(y0))); ... % max allows for either m or u0 to be a scalar
         Kc(my,y0,rjd,rj,rj) - ...
         Kc(my,y0,rjd,rj,rjd); ... % integral over rise
         Kc(my,y0,rj(end),rju(end),rju(end)) - ...
         Kc(my,y0,rj(end),rju(end),rj(end)) ... % incline, last element
        ] + [
         Kc(my,y0,rj(1),rjd(1),rj(1)) - ...
         Kc(my,y0,rj(1),rjd(1),rjd(1)); ... % decline, first element
         Kc(my,y0,rju,rj,rju) - ...
         Kc(my,y0,rju,rj,rj); ... % integral over decline
         zeros(1,max(length(my),length(y0)))
        ]))';

    K(abs(K)<1e3*eps) = 0; % remove numerical noise
    K = sparse(K); % convert to a sparse matrix






else  % consider 2D case
    tools.textheader('Building NRAP kernel');
    disp(' (Linear, 2D, Direct)');

    if aso2.N<3; error(' Aso does not have enough annuli for linear basis.'); end

    mx(abs(mx) < 1e-12) = 1e-12; % avoid division by zero in rv

    % edges of annuli and axial elements
    rjd0 = aso2.re(1:(end-2));
    rj0  = aso2.re(2:(end-1));
    rju0 = aso2.re(3:end);
    xj  = aso2.xe(1:(end-1));
    xju = aso2.xe(2:end);
    
    % functions for indefinite integral
    A = @(my,y0,ry,ryu,r3) y0 ./ (ryu - ry) .* ...
        log(abs(r3 + sqrt(r3.^2 - y0.^2 ./ (1+my.^2))));  % main part of integrand
    Ad = @(my,ryd,ry,r1,r2) my .* ...
        (r1 - r2) ./ (ry - ryd);  % second part of integrand


    N_beams = max([length(my), length(y0), ...
        length(mx), length(x0)]); % any of the values could be scalar, so take max
    K = spalloc(N_beams, aso2.Nx * (aso2.Nr + 1), ...
        round(1e-4 * N_beams * aso2.Nx * (aso2.Nr + 1)));
        % initialize K, assume 0.5% full

    Kx = spalloc(N_beams, aso2.Nx * (aso2.Nr + 1), ...
        round(1e-5 * N_beams * aso2.Nx * (aso2.Nr + 1)));
    
    o = [x0 + eps; y0; 0.*x0];
    v = [mx; my; ones(size(my))]; v = normc(v);
    d = (xj' - o(1,:)') ./ (v' * [1 0 0]'); % distance of start of ASO
    % c = [x0; 0.*x0] + bsxfun(@times, v, d'); % adjust ray position to start of ASO
    
    
    % Find radial intersection with axial elements.
    % This will be a strictly positive number.
    % This could be a bound for either term of the integrand.
    rx0 = sqrt(1 ./ (1+my.^2) .* ...
        ((1+my.^2) ./ mx .* ...
        (xj - x0) + my.*y0) .^ 2 + ...
        y0.^2); % lower edge, vector over u0 and v0
    rxu0 = sqrt(1 ./ (1+my.^2) .* ...
        ((1+my.^2) ./ mx .* ...
        (xju - x0) + my.*y0) .^ 2 + ...
        y0.^2); % upper edge, vector over u0 and v0
    
    
    
    
    % Find z intersection with axial elements.
    % This is used to determine whether to involve front of back
    % portion of integral.
    zu0  = (xj - x0) ./ mx;
    zuu0 = (xju - x0) ./ mx;
    
    
    %== COMMON PROPERTIES ======================================%
    % Flag where intersect occurs.
    fx0  = zu0 < (-my .* y0 ...
        ./ (1 + my.^2));  % flags if intersect is in first integrand
    fxu0 = zuu0 < (-my .* y0 ...
        ./ (1 + my.^2));  % flags if upper intersect is in first integrand
    
    % Flag for +ive/-ive slope.
    f_mx0 = repmat(mx<0, [aso2.Nx, 1]);
    
    % FOR: first integrand.
    r1_a0 = rx0; r1_a0(f_mx0) = rxu0(f_mx0);
    r2_a0 = rxu0; r2_a0(f_mx0) = rx0(f_mx0);
    f1_a0 = fx0; f1_a0(f_mx0) = fxu0(f_mx0);
    f2_a0 = fxu0; f2_a0(f_mx0) = fx0(f_mx0);
    
    % FOR: second integrand.
    % Note that fx* flags are reversed from first integrand.
    r1_b0 = rxu0; r1_b0(f_mx0) = rx0(f_mx0);
    r2_b0 = rx0; r2_b0(f_mx0) = rxu0(f_mx0);
    f1_b0 = ~fxu0; f1_b0(f_mx0) = ~fx0(f_mx0);
    f2_b0 = ~fx0; f2_b0(f_mx0) = ~fxu0(f_mx0);
    
    
    
    disp(' Looping through axial slices:');
    tools.textbar([0,aso2.Nx]);
    Kcol = [];
    for ii=1:aso2.Nx % loop through and append axial slices
        
        % Get z intersection with axial elements.
        zu  = zu0(ii,:);
        zuu = zuu0(ii,:);
        
        
        % Flag if ray intersects elements at all.
        f0 = or(or(or(and(zu > -aso2.R, zu < aso2.R), ... % if zv is in the bounds
            and(zuu > -aso2.R, zuu < aso2.R)), ... % if zvu is in the bounds
            and(zuu > 0, zu < 0)), ...
            and(zuu < 0, zu > 0)); % if there is a change of sign
        idx_a = 1:length(f0);
        idx_a = idx_a(f0); % indices of rays that intersect element
        
        
        % Get radial intersection with axial elements.
        rx = rx0(ii, idx_a); % lower edge, vector over u0 and v0
        rxu = rxu0(ii, idx_a); % upper edge, vector over u0 and v0



        %== FIRST INTEGRAND ======================================%
        % Get where intersect occurs.
        fx  = fx0(ii, idx_a); % flags if intersect is in first integrand
        fxu = fxu0(ii, idx_a); % flags if upper intersect is in first integrand
        
        idx_b = idx_a(or(fx,fxu));  % get subset of rays for this integrand
        
        r1 = r1_a0(ii, idx_b);
        r2 = r2_a0(ii, idx_b);
        f1 = f1_a0(ii, idx_b);
        f2 = f2_a0(ii, idx_b);
        
        % Modified element widths.
        % (adjusted for intersections with axial elements)
        rjd = rjd0 .* ones([1,length(idx_b)]); % repeat for relevant rv elements
        rjd(:,f1) = min(r1(f1), rjd(:,f1)); % adjust for lower axial bound
        rjd(:,f2) = max(r2(f2), rjd(:,f2)); % adjust for upper axial bound

        rj = rj0 .* ones([1,length(idx_b)]);
        rj(:,f1) = min(r1(f1), rj(:,f1));
        rj(:,f2) = max(r2(f2), rj(:,f2));

        rju = rju0 .* ones([1,length(idx_b)]);
        rju(:,f1) = min(r1(f1), rju(:,f1));
        rju(:,f2) = max(r2(f2), rju(:,f2));
        
        % evaluate lower kernel (up to midpoint in ASO)
        K1 = (1 ./ (1+my(idx_b).^2) .* ( ...
            ([ ...
             zeros(1,size(rj,2)); ...
             A(my(idx_b), y0(idx_b), rjd0, rj0, rj) - ...
             A(my(idx_b), y0(idx_b), rjd0, rj0, rjd) + ... % integral over rise
             Ad(my(idx_b), rjd0, rj0, rjd, rj); ... % added component of integrand
             A(my(idx_b), y0(idx_b), rj0(end,:), rju0(end,:), rju(end,:)) - ...
             A(my(idx_b), y0(idx_b), rj0(end,:), rju0(end,:), rj(end,:)) + ...
             Ad(my(idx_b), rj0(end,:), rju0(end,:), rj(end,:), rju(end,:)) ... % incline, last element
            ] + [
             A(my(idx_b), y0(idx_b), rj0(1,:), rjd0(1,:), rj(1,:)) - ...
             A(my(idx_b), y0(idx_b), rj0(1,:), rjd0(1,:), rjd(1,:)) - ...
             Ad(my(idx_b), rj0(1,:), rjd0(1,:), rj(1,:), rjd(1,:)); ...
             A(my(idx_b), y0(idx_b), rju0, rj0, rju) - ...
             A(my(idx_b), y0(idx_b), rju0, rj0, rj) - ...
             Ad(my(idx_b), rju0, rj0, rju, rj); ... % integral over decline
             zeros(1,size(rj,2))
            ])))';


        
        %== SECOND INTEGRAND =====================================%
        idx_c = idx_a(or(~fx,~fxu));
        
        r1 = r1_b0(ii, idx_c);
        r2 = r2_b0(ii, idx_c);
        f1 = f1_b0(ii, idx_c);
        f2 = f2_b0(ii, idx_c);
        
        % Modified element width.
        % (adjusted for intersections with axial elements)
        rjd = rjd0 .* ones([1,length(idx_c)]); % repeat for relevant rv elements
        rjd(:,f1) = min(r1(f1), rjd(:,f1)); % adjust for lower axial bound
        rjd(:,f2) = max(r2(f2), rjd(:,f2)); % adjust for upper axial bound

        rj = rj0 .* ones([1,length(idx_c)]);
        rj(:,f1) = min(r1(f1), rj(:,f1));
        rj(:,f2) = max(r2(f2), rj(:,f2));

        rju = rju0 .* ones([1,length(idx_c)]);
        rju(:,f1) = min(r1(f1), rju(:,f1));
        rju(:,f2) = max(r2(f2), rju(:,f2));
        
        % evaluate second integrand (beyond midpoint in ASO)
        K2 = (1 ./ (1+my(idx_c).^2) .* ( ...
            ([ ...
             zeros(1, size(rj,2)); ...
             A(my(idx_c), y0(idx_c), rjd0, rj0, rj) - ...
             A(my(idx_c), y0(idx_c), rjd0, rj0, rjd) - ... % integral over rise
             Ad(my(idx_c), rjd0,rj0, rjd, rj); ... % added component of integrand
             A(my(idx_c), y0(idx_c), rj0(end,:), rju0(end,:), rju(end,:)) - ...
             A(my(idx_c), y0(idx_c), rj0(end,:), rju0(end,:), rj(end,:)) - ...
             Ad(my(idx_c), rj0(end,:), rju0(end,:), rj(end,:), rju(end,:)) ... % incline, last element
            ] + [
             A(my(idx_c), y0(idx_c), rj0(1,:), rjd0(1,:), rj(1,:)) - ...
             A(my(idx_c), y0(idx_c), rj0(1,:), rjd0(1,:), rjd(1,:)) + ...
             Ad(my(idx_c), rj0(1,:), rjd0(1,:), rj(1,:), rjd(1,:)); ...
             A(my(idx_c), y0(idx_c), rju0, rj0, rju) - ...
             A(my(idx_c), y0(idx_c), rju0, rj0, rj) + ...
             Ad(my(idx_c), rju0, rj0, rju, rj); ... % integral over decline
             zeros(1,size(rj,2))
            ])))';
        
        K0 = zeros(length(idx_a), aso2.Nr+1);
        K0(idx_b - idx_a(1) + 1, :) = K1;
        K0(idx_c - idx_a(1) + 1, :) = K0(idx_c - idx_a(1) + 1, :) + K2;
        K0(isnan(K0)) = 0; % remove NaN values that result when modified element width is zero
        K0(abs(K0)<1e3*eps) = 0;  % supress numerical noise
        K0 = sparse(K0);
        
        K(idx_a, ((ii-1)*(aso2.Nr+1)+1):(ii*(aso2.Nr+1))) = K0;
        
        %== Evaluate axial deflections ===========================%
        %-{
        Kx0 = ([zeros(1, length(idx_a)); ...
             (rx - rjd0) ./ (rj0 - rjd0) .* and(rjd0 < rx, rx < rj0) - ...
             (rxu - rjd0) ./ (rj0 - rjd0) .* and(rjd0 < rxu, rxu < rj0); ...
             (rx - rjd0(end,:)) ./ (rj0(end,:) - rjd0(end,:)) .* and(rjd0(end,:) < rx, rx < rj0(end,:)) - ...
             (rxu - rjd0(end,:)) ./ (rj0(end,:) - rjd0(end,:)) .* and(rjd0(end,:) < rxu, rxu < rj0(end,:)) % if crosses in rise
            ] + [(rx - rju0(1,:)) ./ (rj0(1,:) - rju0(1,:)) .* and(rju0(1,:) < rx, rx < rj0(1,:)) - ...
             (rxu - rju0(1,:)) ./ (rj0(1,:) - rju0(1,:)) .* and(rju0(1,:) < rxu, rxu < rj0(1,:)); ...
             (rx - rju0) ./ (rj0 - rju0) .* and(rju0 < rx, rx < rj0) - ...
             (rxu - rju0) ./ (rj0 - rju0) .* and(rju0 < rxu, rxu < rj0); ...
             zeros(1, length(idx_a)) ... % if crosses in decline
            ])';
        Kx0(abs(Kx0)<1e3*eps) = 0;  % supress numerical noise
        Kx(idx_a, ((ii-1)*(aso2.Nr+1)+1):(ii*(aso2.Nr+1))) = sparse(Kx0);
        % + remove derivative below
        %}
        
        %{
        Kx0 = (1 ./ sqrt(1 + mx(idx_a).^2) .* ([zeros(1, length(idx_a)); ...
             rx .* (rx ./ 2 - rjd0) ./ (rj0 - rjd0) .* and(rjd0 < rx, rx < rj0) - ...
             rxu .* (rxu ./ 2 - rjd0) ./ (rj0 - rjd0) .* and(rjd0 < rxu, rxu < rj0); ...
             rx .* (rx ./ 2 - rjd0(end,:)) ./ (rj0(end,:) - rjd0(end,:)) .* and(rjd0(end,:) < rx, rx < rj0(end,:)) - ...
             rxu .* (rxu ./ 2 - rjd0(end,:)) ./ (rj0(end,:) - rjd0(end,:)) .* and(rjd0(end,:) < rxu, rxu < rj0(end,:)) % if crosses in rise
            ] + [rx .* (rx ./ 2 - rju0(1,:)) ./ (rj0(1,:) - rju0(1,:)) .* and(rju0(1,:) < rx, rx < rj0(1,:)) - ...
             rxu .* (rxu ./ 2 - rju0(1,:)) ./ (rj0(1,:) - rju0(1,:)) .* and(rju0(1,:) < rxu, rxu < rj0(1,:)); ...
             rx .* (rx ./ 2 - rju0) ./ (rj0 - rju0) .* and(rju0 < rx, rx < rj0) - ...
             rxu .* (rxu ./ 2 - rju0) ./ (rj0 - rju0) .* and(rju0 < rxu, rxu < rj0); ...
             zeros(1, length(idx_a)) ... % if crosses in decline
            ]))';
        Kx(idx_a, ((ii-1)*(aso2.Nr+1)+1):(ii*(aso2.Nr+1))) = sparse(Kx0);
        % + add derivative below
        %}

        tools.textbar([ii, aso2.Nx]);
    end

    % Kx = Kx * aso2.Dx;
    tools.textheader;
end


end
