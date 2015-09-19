function [avls]=calchistcomb(avls);
vls=avls.adjvls;
for ntind=1:length(avls.NT);
    muruns=find(avls.mulist>0)
    for ii=1:length(muruns)
        muind=muruns(ii);
        muindvl=avls.mulist(muind);
        acpreindvl=avls.aclist(muind,1);
        acpstindvl=avls.aclist(muind,2);
        acvls=[vls{ntind}{acpreindvl}(:,2) ;vls{ntind}{acpstindvl}(:,2)];
        muvls=[vls{ntind}{muindvl}(:,2)];
        edges=avls.edges{ntind};
        hstac=histc(acvls,edges);
        hstmu=histc(muvls,edges);
        hstac=hstac/sum(hstac);
        hstmu=hstmu/sum(hstmu);
        avls.hstac_comb{ntind}{muindvl}=hstac;
        avls.hstmu_comb{ntind}{muindvl}=hstmu
    end
end
