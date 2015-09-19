function [avfv]=get_avfv(fv);
%[avfv]=get_avfv(fv);
%
tmp=0*fv(1).fdat;
for ii=1:length(fv)
	tmp=tmp+sqrt(abs(fv(ii).fdat));
end
avfv=tmp./length(fv);
return;
