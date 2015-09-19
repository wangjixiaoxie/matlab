function avls=calcrat(avls)
avls.acmean=[];
avls.mumean=[];
for ntind=1:length(avls.NT)
    for muind=1:length(avls.mulist)
        muindvl=avls.mulist(muind);
        if(muindvl==0)
            acindvl=avls.aclist(muind);
        else
            acindvl=muindvl;
        end
        
        
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
        avls.acstdv(ntind,acindvl)=std(acvls)
        avls.acstderr(ntind,acindvl)=std(acvls)/sqrt(length(acvls));
         acmean=mean(acvls);
        acpremean=mean(acprevls);
        acpstmean=mean(acpstvls);
        acn=length(acvls);
        acpren=length(acprevls);
        acpstn=length(acpstvls);
        accv=std(acvls)/mean(acvls);
        avls.acz(ntind,acindvl)=(acmean-initmean)/initsd;
        avls.acn(ntind,acindvl)=acn;
        avls.acprez(ntind,acindvl)=(acpremean-initmean)/initsd;
        avls.acpstz(ntind,acindvl)=(acpstmean-initmean)/initsd;
        avls.acerracz(ntind,acindvl)=avls.acstderr(ntind,acindvl)/initsd;
        avls.acmean(ntind,acindvl)=acmean;
        avls.acpremean(ntind,acindvl)=acpremean;
        avls.acpstmean(ntind,acindvl)=acpstmean;
    if(muindvl)
   
        muvls=[avls.adjvls{ntind}{muindvl}(:,2)]
        
        avls.mustdv(ntind,muindvl)=std(muvls);
        avls.mustderr(ntind,muindvl)=std(muvls)/sqrt(length(muvls));
       
        mumean=mean(muvls);
        mun=length(muvls);
        avls.mun(ntind,muindvl)=mun;
        avls.mumean(ntind,muindvl)=mumean;
        mucv=std(muvls)/mean(muvls);
        avls.cvfac(ntind, muindvl)=accv/mucv
        
        avls.muz(ntind,muindvl)=(mumean-initmean)/initsd;
        
        avls.mustderrz(ntind,muindvl)=avls.mustderr(ntind,muindvl)/initsd;
        
    end
    end
end