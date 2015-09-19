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

if (~exist('PLTIT'))
    PLTIT=1;
end

days   = zeros([length(trigs),1]);
months = zeros([length(trigs),1]);
hours  = zeros([length(trigs),1]);
years  = zeros([length(trigs),1]);
for ii = 1:length(trigs)
    [h,d,m,y]  = fn2date(trigs(ii).fn);
    days(ii)   = d;
    months(ii) = m;
    hours(ii)  = h;
    years(ii)  = y;
end

firstval=1;
outvals=[];
if (PLTIT==1)
    hold on;
end
for year = min(years):max(years)
    for month=min(months):max(months)
        for day=min(days):max(days)
            if (exist('TIMERNG'))
                pp=find((days==day)&(months==month)&(years==year)&(hours>=TIMERNG(1))&(hours<=TIMERNG(2)));
            else
                pp=find((days==day)&(months==month)&(years==year));
            end

            if (length(pp)>1)
                vals=zeros([length(pp),size(tvals,2)]);
                for ii=1:length(pp)
                    vals(ii,:)=tvals(pp(ii),:);
                end
                h=trigs(pp(ii)).ftime;
                minu = (h-fix(h))*60;
                dtvec=[year,month,day,fix(h),fix(minu),fix((minu-fix(minu))*60)];

                daten = fix(datenum(dtvec));
                if (exist('start_t'))
                    if (length(start_t)>0)
                        daten=daten-start_t;
                    end
                end
                if (firstval==1)
                    startdate=daten;
                    firstval=0;
                end

                tmpsum=sum(vals);
                outvals=[outvals;daten,tmpsum,length(vals)];

                if (PLTIT==1)
                    plot(fix(daten),tmpsum(1)./tmpsum(2),sym);
                end
            end
        end
    end
end
return;