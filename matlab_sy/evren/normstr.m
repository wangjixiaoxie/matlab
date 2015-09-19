function out_str=normstr(in_str);
% out_str=normstr(in_str);
% replace '_' in in_str with '\_' so that it looks right
% when put on a plot
%
out_str=in_str;
out_str = strrep(in_str,'_','\_');
return;
