function outvals=plt_ampday(vals,sym,start_t,PLTIT);
% outvals=plt_ampday(vals,sym,start_t,PLTIT);

if (~exist('sym'))
    sym='bs';
else 
    if (length(sym)<1)
        sym='bs-';
    end
end

if (~exist('start_t'))
    start_t=0;
else
    if (length(start_t)==0)
        start_t=0;
    end
end

if (~exist('PLTIT'))
    PLTIT=1;
end

days = vals(:,1);

outvals=[];
if (PLTIT==1)
    hold on;
end

dvals=unique(fix(days));
for jj=1:length(dvals)
	day=dvals(jj);
	pos=find(fix(days)==day);
	if (length(pos)==0)
		continue;
	end

        outvals=[outvals;day-start_t,mean(vals(pos,2)),std(vals(pos,2))];
        if (PLTIT==1)
            errorbar(day-start_t,mean(vals(pos,2)),std(vals(pos,2)),sym);
        end
end

return;
