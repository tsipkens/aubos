
% IDX2TXT  Converts an index to a three character string.
% Author: Timothy Sipkens, 2020-04-01
%=========================================================================%

function [txt] = idx2txt(idx)

ni = floor(log10(idx))+1;
if ni==1; txt = ['00',num2str(idx)];
elseif ni==2; txt = ['0',num2str(idx)];
elseif ni==3; txt = num2str(idx);
end

end

