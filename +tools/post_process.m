
% POST_PROCESS  Post-process and reformat multiple reconstructions.

function [out] = post_process(f, sz, C0, varargin)

% Update size based on f. 
dim_half = max(sum(reshape(f, sz)));
sz = [dim_half, sum(f) / dim_half];

for ii = 1:length(varargin)
    name0 = inputname(ii + nargin - length(varargin));
    name = parsename(name0);
    
    NAMES{ii} = name;
    BETA(ii, :) = varargin{ii};
    
    DEL(ii, :) = varargin{ii} - varargin{1};
    
    ERROR(ii) = norm(varargin{ii}(f) - varargin{1}(f)) ...
         / sum(f) / mean(varargin{1}(f));
    
    SELF_SIM(ii) = ssim( ...
        reshape(varargin{ii}(f), sz) ./ C0, ...
        reshape(varargin{ 1}(f), sz) ./ C0);
    
end

% Relative error, phrased relative to the 2nd entry.
for ff=1:length(NAMES)
    REL_ERROR(ff) = (ERROR(ff) - ERROR(2)) ./ ERROR(2);
end

% Transpose for table.
ERROR = ERROR';
SELF_SIM = SELF_SIM';
REL_ERROR = REL_ERROR' .* 100;

out = table(ERROR, SELF_SIM, REL_ERROR, ...
    BETA, DEL, 'RowNames', upper(NAMES));

end



function out = parsename(txt)

% Remove text before first underscore.
out = split(txt, '_');
out = [out, repmat({'-'}, [length(out), 1])]';  % add underscores back
out = out(:);  out = out(2:(end-1));
out = [out{2:end}];

if isempty(out)
    out = 'TRUTH';
elseif out(end) == 'a'
    out = out(1:(end-1));
end

end

