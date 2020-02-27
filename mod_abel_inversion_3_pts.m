%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementing modified Abel inversion using 3 points method according to https://doi.org/10.1103/PhysRevSTAB.18.081302
% The differences between modified Abel inversion and normal Abel inversion is that the modified one takes into account
%   if the object is moving to longitudinal direction and/or the probe is not crossing the object perpendicularly.
% Inputs:
%   * z, y: one dimensional axis in real unit (or in your own unit); z is longitudinal and y is transverse
%   * F: the 2D half-image to be inverted; axis of the cylinder is at index (1,:) of the image
%   * a: a variable to take into account if the object is moving (opt, def=0)
% Outputs:
%   * zeta, r: one dimensional axis in real unit (or in your own unit) of the output image
%   * f: the 2D half-image resulted from the modified Abel inversion
% Author: Muhammad F. Kasim (University of Oxford, 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zeta, r, f] = mod_abel_inversion_3_pts(z, y, F, a)
    if (nargin < 4)
        a = 0;
    end
    
    cutoff = 0.8; % should be 0.5
    
    % extract the dimension parameters
    dz = z(2) - z(1);
    dr = y(2) - y(1);
    [Nr, Nz] = size(F);
    kz = [-floor(Nz/2):ceil(Nz/2)-1] * 2 * pi / Nz / dz;
    
    % input map in z-frequency domain
    fftF = fftshift(fft(F, [], 2), 2);
    
    % construct the matrix B which can be evaluated without k
    [I, J] = meshgrid([0:Nr-1]);
    I2J2 = ((I.*I - J.*J) .* (I > J));
    Bij0 = zeros(Nr);
    Bij0_a = log((I + 0.5 + sqrt((I+.5).^2 - J.^2)) ./ I);
    Bij0_b = log((I + 0.5 + sqrt((I+.5).^2 - J.^2)) ./ (I - 0.5 + sqrt((I-.5).^2 - J.^2)));
    Bij0((I == J) & (I ~= 0)) = Bij0_a((I == J) & (I ~= 0));
    Bij0(I > J) = Bij0_b(I > J);
    Bij0(1,1) = Bij0(2,2); % handle the singularity (can be commented out to ignore the singularity at all)
    
    Bij1 = zeros(Nr);
    Bij1_a = sqrt((I+.5).^2 - J.^2) - I .* Bij0;
    Bij1_b = sqrt((I+.5).^2 - J.^2) - sqrt((I-.5).^2 - J.^2) - I .* Bij0;
    Bij1((I == J) & (I ~= 0)) = Bij1_a((I == J) & (I ~= 0));
    Bij1(I > J) = Bij1_b(I > J);
    Bij1(1,1) = Bij1(2,2); % handle the singularity (can be commented out to ignore the singularity at all)
    
    % the shifting matrix operators
    Pneg = circshift(eye(Nr), [0, -1]);
    Pneg(1,:) = Pneg(3,:);
    Ppos = circshift(eye(Nr), [0, 1]);
    Ppos(Nr,:) = zeros(1, Nr);
    
    % iterate for each value of frequency
    fftf = zeros(size(fftF));
    for (i = [1:Nz])
        ka = kz(i) * a;
        if (abs(ka * Nr * dr) >= 1/(1-cutoff) * (2*pi/Nz/dz*a*Nr*dr)) ka = 0; end;
        
        Fy = fftF(:,i);
        
        % construct the integral matrices
        Cij = cosh(ka * dr * I2J2.^.5);
        Sij = sinh(ka * dr * I2J2.^.5) ./ (I2J2 - 9999 * (I <= J)) * I * ka * dr;
        Sij = Sij .* (I > J) + (ka * dr)^2 * I .* (I == J);
        CijBij0 = Cij .* Bij0;
        CijBij1 = Cij .* Bij1;
        SijBij1 = Sij .* Bij1;
        W = (-.5 * CijBij0 + CijBij1 - .5 * SijBij1) * Pneg + ...
            ( .5 * CijBij0 + CijBij1 + .5 * SijBij1) * Ppos + ...
            (-2  * CijBij1);
        
        % decimate the interpolated data to fit within the output size
        fftf(:,i) = W * Fy / (-pi * dr);
    end
    
    % get the output
    f = ifft(ifftshift(fftf, 2), [], 2);
    zeta = z;
    r = y;
end
