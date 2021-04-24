
% POST_PROCESS  Post-process and reformat multiple reconstructions.

function [out] = post_process(f, sz, varargin)

% Update size based on f. 
dim_half = max(sum(reshape(f, sz)));
sz = [dim_half, sum(f) / dim_half];

for ii = 1:length(varargin)
    name = parsename(inputname(ii + 2));
    
    NAMES{ii} = name;
    BETA(ii, :) = varargin{ii};
    
    DEL(ii, :) = varargin{ii} - varargin{1};
    
    ERROR(ii) = norm(varargin{ii}(f) - varargin{1}(f)) ...
         / sum(f) / mean(varargin{1}(f));
    
    SELF_SIM(ii) = ssim( ...
        reshape(varargin{ii}(f), sz), ...
        reshape(varargin{1}(f), sz));
    
end

% Relative error, phrased relative to the last entry.
for ff=1:length(NAMES)
    REL_ERROR(ff) = (ERROR(end) - ERROR(ff)) ./ ERROR(ff);
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

