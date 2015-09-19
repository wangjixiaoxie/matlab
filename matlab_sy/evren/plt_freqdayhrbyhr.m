function outvals=plt_freqdayhrbyhr(fv,tbins,fbins,sym,start_t,TIMERNG,NN,MM);
%outvals=plt_freqdayhrbyhr(fv,tbins,fbins,sym,start_t,TIMERNG,NN,MM);

%mndy=[];

if (~exist('sym'))
    sym='bs';
else
    if (length(sym)<1)
        sym='bs-';
    end
end

if (~exist('NN'))
    NN=3;
end

if (~exist('MM'))
    MM=2;
end


days   = zeros([length(fv),1]);
months = zeros([length(fv),1]);
hours  = zeros([length(fv),1]);
for ii = 1:length(fv)
    days(ii)   = fv(ii).day;
    months(ii) = fv(ii).month;
    hours(ii)  = fv(ii).hour;
end

firstval=1;
outvals=[];
hold on;
for month=min(months):max(months)
    for day=min(days):max(days)
        if (exist('TIMERNG'))  
            pp=find((days==day)&(months==month)&(hours>=TIMERNG(1))&(hours<=TIMERNG(2)));
        else
            pp=find((days==day)&(months==month));
        end

        if (length(pp)>1)
            for hh=min(fix(hours(pp))):max(fix(hours(pp)))
                ppp=find(fix(hours(pp))==hh);
                if (length(ppp)>0)
                    vals=zeros([length(ppp),1]);
                    for ii=1:length(ppp)
                        vals(ii)=fv(pp(ppp(ii))).mxvals(MM,NN+1);
                    end
                    [h,d,m,y]=fn2date(fv(pp(ppp(ii))).fn);
                    minu = (h-fix(h))*60;
                    dtvec=[y,m,d,fix(hh),0,0];
                    stdv=std(vals);

                    daten = datenum(dtvec);
                    if (firstval==1)
                        if (~exist('start_t'))
                            startdate=daten;
                        else
                            startdate=start_t;
                        end
                        firstval=0;
                    end
                    daten=daten-startdate;

                    outvals=[outvals;daten,median(vals),stdv,stdv./sqrt(length(vals)),length(vals)];

                    plot(daten,median(vals),sym);
                    plot(daten+[0,0],median(vals)+[-1,1]*std(vals)/sqrt(size(vals,1)),[sym(1),'-']);
                end
            end
        end
    end
end
return;