function [fvind,vals]=find_fvind(fv,dn,fileind,noteind);
% fvind=find_fvind(fv,dn,fileind,noteind);
%

fvind=[];
vals=getvals(fv,1,'ind');
vals=[vals(:,1),zeros([length(vals),1]),vals(:,3)];
for ii=1:length(fv)
	fn=fv(ii).fn;
	pp=findstr(fn,'.cbin');
	ppp=findstr(fn,'.');
	tmp=find(ppp<pp);ppp=ppp(tmp(end));
	vals(ii,2)=str2num(fn(ppp+1:pp-1));
end

for ii=1:length(dn)
	pp=find((vals(:,1)==dn(ii))&...
		(vals(:,2)==fileind(ii))&...
		(vals(:,3)==noteind(ii)));
	if (length(pp)==1)
		fvind=[fvind,pp];
	elseif (length(pp)>1)
		disp(['hey ind=',num2str(ii)]);
	end
end
return;
