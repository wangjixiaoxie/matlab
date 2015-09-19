function [vals,cnts]=get_dates(ff,DISP);
%[vals,cnt]=get_dates(ff,DISP);
% looks through ff=dir('*.cbin'); and returns all the possible
% datenums can also return the counts for each one
if (~exist('DISP'))
	DISP=0;
end
vals=zeros([length(ff),1]);
for ii=1:length(ff)
	vals(ii)=fn2datenum(ff(ii).name);
end
vals=sort(fix(vals));
vals2=unique(vals);
cnts=zeros([length(vals2),2]);
for ii=1:length(vals2)
	pp=find(vals==vals2(ii));
	cnts(ii,:)=[vals2(ii),length(pp)];
end

if (DISP==1)
	for ii=1:size(vals2,1)
		dv=datevec(cnts(ii,1));
		disp([num2str(dv(2)),'/',num2str(dv(3)),'/',...
                      num2str(dv(1)),' - ',num2str(cnts(ii,2))]);
	end
end

vals=vals2;
return;
