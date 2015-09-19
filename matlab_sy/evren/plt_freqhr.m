function outvals=plt_freqhr(fv,tbins,fbins,sym,dy,NN,MM);
%outvals=plt_freqhr(fv,tbins,fbins,sym,dy,NN,MM);

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

days  = zeros([length(fv),1]);
hours = zeros([length(fv),1]);
for ii = 1:length(fv)
    days(ii)  = fv(ii).day;
    hours(ii) = fv(ii).hour;
end

if (exist('dy'))
    pp=find(days==dy);
    fv=fv(pp);
    hours=hours(pp);
end

outvals=[];

hours=fix(hours);

hold on;
for hh=min(hours):max(hours)
    pp=find(hours==hh);
    if (length(pp)>1)
        vals=zeros([length(pp),1]);
        for ii=1:length(pp)
            vals(ii)=fv(pp(ii)).mxvals(MM,NN+1);
        end
        outvals=[outvals;hh,median(vals),std(vals),std(vals)./sqrt(length(vals))];

        plot(hh,median(vals),sym);
        plot(hh+[0,0],median(vals)+[-1,1]*std(vals),[sym(1),'-']);
    end
end
return;