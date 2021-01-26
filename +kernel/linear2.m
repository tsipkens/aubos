
% LINEAR2  Evaluates kernel/operator for a linear basis representation of a 2D ASO.
%  Linear basis is applied in radial direction. 
%  Uniform basis is applied in axial direction.
% 
%  INPUTS:
%   aso2    Axis-symmetric object
%   my      Set of slopes in y-direction for the rays
%   y0      Intersect with line through center of aso
%   mx      Set of slopes in x-direction for the rays
%   x0      Intersect with line through center of aso
%
%  OUTPUTS:
%   K       Overall radial deflection kernel
%   Kx      Overall axial deflection kernel
%           (often noisy due to uniform elements in axial direction)
%  
%  AUTHOR: Timothy Sipkens, 2020-06-10

function [K, Kx] = linear2(aso2, my, y0, mx, x0)


tools.textheader('Building 2D linear kernel');

if aso2.N<3; error('Aso does not have enough annuli for linear basis.'); end

mx(abs(mx) < 1e-12) = 1e-12; % avoid division by zero in rv

% edges of annuli and axial elements
rjd0 = aso2.re(1:(end-2));
rj0  = aso2.re(2:(end-1));
rju0 = aso2.re(3:end);
xj  = aso2.xe(1:(end-1));
xju = aso2.xe(2:end);

% functions for indefinite integral
A = @(mx,x0,ry,ryu,r3) 1 ./ (ryu - ry) .* ...
    log(abs(r3 + sqrt(r3.^2 - x0.^2 ./ (1+mx.^2)))); % main part of integrand
Ad = @(mx,ryd,ry,r1,r2) (r2 - r1) ./ (ry - ryd) ...
    .* mx .* sqrt(1 + mx.^2); % second part of integrand


N_beams = max([length(my), length(y0), ...
    length(mx), length(x0)]); % any of the values could be scalar, so take max
K = spalloc(N_beams, aso2.Nx * (aso2.Nr + 1), ...
    round(1e-4 * N_beams * aso2.Nx * (aso2.Nr + 1)));
    % initialize K, assume 0.5% full

Kx = spalloc(N_beams, aso2.Nx * (aso2.Nr + 1), ...
    round(1e-5 * N_beams * aso2.Nx * (aso2.Nr + 1)));


disp('Looping through axial slices...');
tools.textbar(0);
for ii=1:(length(aso2.xe)-1) % loop through and append axial slices
    
    % Find z intersection with axial elements.
    % This is used to determine whether to involve front of back
    % portion of integral.
    zu  = (xj(ii) - x0) ./ mx;
    zuu = (xju(ii) - x0) ./ mx;
    
    % flag if ray intersects elements at all
    f0 = or(or(or(and(zu > -aso2.R, zu < aso2.R), ... % if zv is in the bounds
        and(zuu > -aso2.R, zuu < aso2.R)), ... % if zvu is in the bounds
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
    
    if any(any(K(idx_a, ((ii-1)*(aso2.Nr+1)+1):(ii*(aso2.Nr+1)))~=0))
        disp('HELP');
    end
    K(idx_a, ((ii-1)*(aso2.Nr+1)+1):(ii*(aso2.Nr+1))) = K0;
    
    
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
    Kx(idx_a, ((ii-1)*(aso2.Nr+1)+1):(ii*(aso2.Nr+1))) = sparse(Kx0);
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
    
    tools.textbar(ii/(length(aso2.xe)-1));
end

% Ky = Ky * aso.Dx;


tools.textheader;

end