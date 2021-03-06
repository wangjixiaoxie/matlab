function avls=calcrat(avls)
avls.acmean=[];
avls.mumean=[];
for ntind=1:length(avls.NT)
    for muind=1:length(avls.mulist)
        muindvl=avls.mulist(muind);
    if(muindvl)
        
        acvls1=avls.adjvls{ntind}{avls.aclist(muind,1)}(:,2);
        if(avls.aclist(muind,2)<1)
            acvls2=[];
        else
            acvls2=avls.adjvls{ntind}{avls.aclist(muind,2)}(:,2);
        end

        initsd=avls.initsd{ntind}
        initmean=avls.initmean{ntind}
        acvls=[acvls1 ;acvls2]
        acprevls=[acvls1]
        acpstvls=[acvls2]
        muvls=[avls.adjvls{ntind}{muindvl}(:,2)]
        avls.acstdv(ntind,muindvl)=std(acvls)
        avls.acstderr(ntind,muindvl)=std(acvls)/sqrt(length(acvls));
        avls.mustdv(ntind,muindvl)=std(muvls);
        avls.mustderr(ntind,muindvl)=std(muvls)/sqrt(length(muvls));
        acmean=mean(acvls);
        acpremean=mean(acprevls);
        acpstmean=mean(acpstvls);
        acn=length(acvls);
        acpren=length(acprevls);
        acpstn=length(acpstvls);
        mumean=mean(muvls);
        mun=length(muvls);
        avls.acmean(ntind,muindvl)=acmean;
        avls.acpremean(ntind,muindvl)=acpremean;
        avls.acpstmean(ntind,muindvl)=acpstmean;
        avls.mumean(ntind,muindvl)=mumean;
        avls.acn(ntind,muindvl)=acn;
        avls.acpren(ntind,muindvl)=acpren;
        avls.acpstn(ntind,muindvl)=acpstn;
        avls.mun(ntind,muindvl)=mun;
        accv=std(acvls)/mean(acvls);
        mucv=std(muvls)/mean(muvls);
        avls.cvfac(ntind, muindvl)=accv/mucv
        avls.acz(ntind,muindvl)=(acmean-initmean)/initsd;
        avls.acprez(ntind,muindvl)=(acpremean-initmean)/initsd;
        avls.acpstz(ntind,muindvl)=(acpstmean-initmean)/initsd;
        avls.muz(ntind,muindvl)=(mumean-initmean)/initsd;
        avls.acerracz(ntind,muindvl)=avls.acstderr(ntind,muindvl)/initsd;
        avls.mustderrz(ntind,muindvl)=avls.mustderr(ntind,muindvl)/initsd;
        
    end
    end
end