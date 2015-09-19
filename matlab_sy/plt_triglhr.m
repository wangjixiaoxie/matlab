function plt_triglhr(trigs,vals,sym,dy);
%plt_triglhr(trigs,vals,sym,dy);

if (~exist('sym'))
    sym='bs';
end

days  = zeros([length(trigs),1]);
hours = zeros([length(trigs),1]);
for ii = 1:length(trigs)
    [hr,day]=fn2date(trigs(ii).fn);
    days(ii)  = day;
    hours(ii) = hr;
end

if (exist('dy'))
    pp=find(days==dy);
else
    pp=[1:length(days)];
end
if (length(pp)<2)
    disp(['no data']);
    return;
end

trigs=trigs(pp);
hours=hours(pp);

hold on;
for hh=min(floor(hours)):max(floor(hours))
    pp=find(floor(hours)==hh);
    if length(pp)>1
        tmp=sum(vals(pp,:));
        plot(hh,tmp(1)./tmp(2),sym);
    elseif length(pp)==1
        tmp=vals(pp,:);
        plot(hh,tmp(1)./tmp(2),sym);
    end
end
return;