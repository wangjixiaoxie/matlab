function outvals=plt_trigday(trigs,tvals,sym,start_t,PLTIT,TIMERNG);
%outvals=plt_trigday(trigs,tvals,sym,start_t,PLTIT,TIMERNG);
%

%mndy=[];

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

days  = zeros([length(trigs),1]);
hours = 0*days;
for ii = 1:length(trigs)
    [h,d,m,y,dn]  = fn2date(trigs(ii).fn);
    days(ii)  = dn;
    hours(ii) = h;
end

outvals=[];
if (PLTIT==1)
    hold on;
end

for day=min(fix(days)):max(fix(days))
    pos=find(fix(days)==day);
    if (length(pos)>0)
        if (exist('TIMERNG'))
            pp=find((hours(pos)>=TIMERNG(1))&(hours(pos)<=TIMERNG(2)));
            pos=pos(pp);
        end
    end

    if (length(pos)>0)
        vals=zeros([length(pos),size(tvals,2)]);
        for ii=1:length(pos)
            vals(ii,:)=tvals(pos(ii),:);
        end



        if (size(vals,1)>1)
            tmpsum=sum(vals);
        else
            tmpsum=vals(1,:);
        end
        outvals=[outvals;fix(days(pos(1)))-start_t,tmpsum,length(vals)];
        if (PLTIT==1)
            plot(fix(days(pos(1)))-start_t,tmpsum(1)./tmpsum(2),sym);
        end
    end
end

return;