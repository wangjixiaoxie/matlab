function [vals,dn]=get_toff(trigs);
% [vals,dn]=get_toff(trigs);
%

vals=[];dn=[];
for ii=1:length(trigs)
	vals=[vals;trigs(ii).toffset];
	dn=[dn;trigs(ii).toffset*0+fn2datenum(trigs(ii).fn)];
end
return;
