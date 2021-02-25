
% TEXTBAR  Print out a text-based progress bar.
% 
%  DESCRIPTION:
%   tools.textbar() or textbar(0) initializes a textbar without
%       a trailing fraction: 
%           0% |                    |
% 
%   tools.textbar(pct) displays a textbar with the provided percentage, 
%       removing text corresponding to the previous textbar.  
%       For example, textbar(0.1) yields:
%           10% |██                  |
% 
%   tools.textbar([i(1), i(2)]) displays a textbar with an appended 
%       fraction, again removing text corresponding to the previous textbar. 
%       In this case, i(1) corresponds to the current iterable, and 
%       i(2) to the total. For example, textbar([1, 10]) yields:  
%           10% |██                  | 1/10
% 
%   tools.textbar([i(1), i(2), i(3), i(4), ...]) displays a textbar considering
%       multiple iterables. These are given as pairs of integers, starting
%       from the most recent loop and moving outward. Integer pairs can be
%       added for as many iterables are present. For example, 
%       textbar([1, 3, 2, 5]) yields:
%           27% |█████▌              | 4/15
%    
%   tools.textbar(..., f_back) displays a textbar with control over whether  
%       the previous textbar is removed, where f_back=1 will remove the
%       previous text and f_back = 0 will not. By default, f_back=1, that
%       is, textbar(0.4, 1) is the same as textbar(0.4). By contrast, 
%       the code `for ii=1:3; textbar([ii, 3], 0); end` will yield:
%            33% |███████             | 1/3
%            67% |█████████████▌      | 2/3
%           100% |████████████████████| 3/3
%           
%  NOTE:
%   If this function is in a tools package, append `tools.` to the 
%   function call. For example: `tools.textbar(0)` to initialize bar.
%  
%  AUTHOR: Timothy Sipkens, 2020-11-02
%  INSPIRED BY: Samuel Grauer, 2017-11-16 + tqdm for Python

function textbar(i, f_back)

%--- Initialization ------------------------------------------------------%
%-- Parse inputs ------------------------------------------%
if ~exist('i', 'var'); i = 0; end

if ~exist('f_back', 'var'); f_back = []; end
if isempty(f_back); f_back = 1; end  % by default, remove
%----------------------------------------------------------%


% If initializing textbar.
if (isempty(i) || i(1) < 0); pct = 0;
    
% If percent is given directly.
elseif length(i)==1; pct = i; i(2) = 0;  

% If integer pair is given.
elseif length(i)==2; pct = i(1) / i(2); 

% If a series of integer pairs is given, 
% get global index and total elements.
else
    j = prod(i(2:2:end));
    k = i(1); % initialize global index
    for ii=3:2:length(i)% loop over dimensions
        k = k + (i(ii) - 1) * prod(i(2:2:(ii-1))); % update global index
    end
    i = [max(k, 0), j];
    pct = i(1) / i(2);
end


% If initializing progress bar, do not erase.
if pct==0; f_back = 0; end


%-- Other parameters ----%
n_dot = 20;  % number of elements in progress bar
n_xtra = 10;  % padded elements  for percent / general spacing
n_frac = 2 * length(num2str(i(2))) + 1;  % number of extra elements due to fraction
n_str = n_dot + n_xtra + n_frac;
%-------------------------------------------------------------------------%


%--- Print progress ------------------------------------------------------%
if f_back; str_back = repmat(char(8), [1, n_str]); % if continuing textbar
else; str_back = ''; end  % if initializing textbar


% Format percentage leading bar.
str_p00 = num2str(100 * pct, '%.0f');   % formatted percent
str_p00 = [repmat(' ', [1, 3-length(str_p00)]), str_p00];  % pad with necessary spaces


% Format text for middle of the bar.
nc = ceil(pct * n_dot);  % number of completed elements
str_p01 = repmat('█', [1, nc]);  % completed portion of bar
if ((nc - pct * n_dot) > 0.5); str_p01(end) = '▌'; end  % allow for half blocks
str_p02 = repmat(' ', [1, n_dot-nc]);  % uncompleted portion of bar


% Format fraction trailing bar.
% Note that this rounds the first number, i(1), down.
% This accomodates non-integer i(1), which can be used for partial steps.
if i(2)==0; str_p03 = '';  % if only percentage given, don't show fraction
else  % otherwise, format fraction of end of bar
    str_p03 = [num2str(floor(i(1)), '%.0f'), ...  
        '/', num2str(i(2), '%.0f')];
end
str_p03 = [str_p03, repmat(' ', [1, n_frac-length(str_p03)])];  % pad with necessary spaces


% Compile formatted sting.
str_out = [' ', str_p00, '%%', ' |', str_p01, str_p02, '| ', str_p03];
str_all = [str_out, repmat(' ', [1, n_str-length(str_out)]), '\n'];

if pct<1
    fprintf([str_back, '[', 8,  str_all, ']', 8]);  % give orange format
else
    fprintf([str_back, str_all]);
end
%-------------------------------------------------------------------------%

end


