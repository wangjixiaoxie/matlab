function vec=get_fvday(fv);
% vec=get_fvday(fv);
%
vec=zeros([length(fv),1]);
for ii=1:length(fv)
	[tmp1,tmp2,tmp3,tmp4,dn]=fn2date(fv(ii).fn);
	vec(ii)=dn;
end
return;
