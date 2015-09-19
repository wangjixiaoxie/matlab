%modified from Evren's original code in order to include divider

function outvals=plt_freqdaymodtw(fv,tbins,fbins,BASELINE,sym,start_t,PLTIT,TIMERNG,NN,MM,FREQRATIO,FLIP,PRNG,VALRNG);
%outvals=plt_freqday(fv,tbins,fbins,sym,start_t,PLTIT,TIMERNG,NN,MM,PRNG,VALRNG);
%PRNG = percentile range
%VALRNG = limit data point range
%mndy=[];

%formatting

lw=1
msize=2
errlw=1
mucol=[0.4 0.4 1]
mufillcol=[0.82 0.82 0.82]
acfillcol=[0.92 0.96 0.98]
accol='k'
col{1}='r';
col{2}=accol;
meanlw=2
errht=50
acmark='o'
mumark='o'


if(~exist('FREQRATIO'))
    FREQRATIO=1;
end

if(~exist('FLIP'))
    FLIP=1;
end


if (~exist('PRNG'))    
    PRNG=0;
end
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
    MM=1;
end

if (~exist('PLTIT'))
    PLTIT=1;
end

days   = zeros([length(fv),1]);
hours  = zeros([length(fv),1]);
for ii = 1:length(fv)
    [tmp1,tmp2,tmp3,tmp4,dn]=fn2datev2(fv(ii).fn);
    days(ii)  = dn;
    hours(ii) = (dn-fix(dn))*24;
end

firstval=1;
outvals=[];
if (PLTIT==1)
    hold on;
end

for day=min(fix(days)):max(fix(days))
    pos=find(fix(days)==day);
    if (length(pos)>0)
        if (exist('TIMERNG'))
            pp=find((hours(pos)>=TIMERNG(1))&(hours(pos)<=TIMERNG(2)));
            pos = pos(pp);
        end
    end
    if (length(pos)>1)
        vals=zeros([length(pos),1]);
        for ii=1:length(pos)
            vals(ii)=fv(pos(ii)).mxvals(MM,NN+1);
        end
	if (exist('VALRNG'))
		outpp=find((vals>=VALRNG(1))&(vals<=VALRNG(2)));
        	disp([num2str(length(vals)-length(outpp)),....
                                 ' outliers removed from ',...
                                 num2str(length(vals)),' points!A']);
		vals=vals(outpp);
	end
        %look for outliers
        outpp = find((vals<fbins(NN,1))|(vals>fbins(NN,2)));
	if (length(outpp)>0)
          	disp([num2str(length(outpp)),' outliers removed from ',...
                               num2str(length(vals)),' points!']);
	end
        vals(outpp)=[];
        stdv=std(vals);

        daten=day;
        if (exist('start_t'))
            if (length(start_t)>0)
                daten=daten-start_t;
            end
        end
        rng=prctile(vals,[PRNG,100-PRNG]);
        outvals=[outvals;daten,median(vals),stdv,stdv./sqrt(length(vals)),length(vals),rng];

        if (PLTIT==1)
            plot(daten,((median(vals)./FREQRATIO)-BASELINE)*FLIP,'Marker',acmark,'Color',accol,'Markersize',msize,'MarkerFaceColor',accol);
            plot(daten+[0,0],(((median(vals)+[-1,1]*std(vals))./FREQRATIO)-BASELINE)*FLIP,'k','Linewidth',errlw);
        end
    end
end
return;
