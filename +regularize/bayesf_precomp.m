
% BAYESF_PRECOMP Evaluate the Bayes factor with pre-comptued GSVD.
% This is useful in the context that Lpr0 does not change, but lambda does.
% Author: Timothy Sipkens, 2019-12-21
%=========================================================================%

function [B,F,C] = bayesf_precomp(A,b,x,Lpr0,lambda,S1,S2,order)

if ~exist('order','var'); order = []; end
if isempty(order); order = 0; end
    % if not Tikhonov and Lpr0 is full rank, use order = 0


r = length(x)-order;
lx = size(A,2); % n in Thompson and Kay
lb = size(A,1); % s in Thompson and Kay

%-- If the gsvd was not pre-computed --%
if ~exist('S1','var'); S1 = []; S2 = []; end
if isempty(S1)
    [~,~,~,S1,S2] = gsvd(full(A),full(Lpr0));
end

F = -1/2*(norm(A*x-b)^2 + norm(lambda.*Lpr0*x)^2); % fit

h = max(S1,[],2); % picks off diagonal elements of each row
the = diag(S2);

ct = 2*sum(log(the((lx-r+1):lb))); % contributes to the constant term
det_po = sum(log(h((lx-r+1):lb).^2 + lambda^2.*the((lx-r+1):lb).^2));

C = (r+lb-lx)*log(lambda) - det_po/2 + ct/2;
    % measurement credence
    % ignores det_b contributions

B = F+C;

end

