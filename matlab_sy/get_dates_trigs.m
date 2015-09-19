function outvals=get_dates_trigs(ff,DISP);
%get_dates(ff,DISP);
% looks through ff=dir('*.cbin'); and returns all the possible
% datenums can also return the counts for each one
if (~exist('DISP'))
	DISP=0;
end

vals=zeros([length(ff),2]);
for ii=1:length(ff)
	vals(ii,1)=fn2datenum(ff(ii).name);
	rd=readrecf(ff(ii).name);
	if (isfield(rd,'ttimes'))
		vals(ii,2)=length(rd.ttimes);
	else
		vals(ii,2)=0;
	end
end

outvals=[];
[tmpy,tmpi]=sort(fix(vals(:,1)));
vals=vals(tmpi,:);
vals2=unique(fix(vals(:,1)));
cnts=zeros([length(vals2),3]);
for ii=1:length(vals2)
	pp=find(fix(vals(:,1))==vals2(ii));
	cnts(ii,:)=[vals2(ii),length(pp),sum(vals(pp,2))];
	if (DISP==1)
		%for ii=1:size(vals2,1)
			dv=datevec(cnts(ii,1));
			disp([num2str(dv(2)),'/',num2str(dv(3)),'/',...
                      	num2str(dv(1)),' - ',num2str(cnts(ii,2)),...
				' - ',num2str(cnts(ii,3))]);
            outvals=[outvals;cnts(ii,:)],
		%end
	end
end
return;
