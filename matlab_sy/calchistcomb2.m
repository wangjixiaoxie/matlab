%calchistcomb2, modified 3/27/09, to make separate pre and post histograms.
function [avls]=calchistcomb2(avls);
vls=avls.adjvls;
for ntind=1:length(avls.NT);
    muruns=find(avls.mulist>0)
    for ii=1:length(muruns)
        muind=muruns(ii);
        muindvl=avls.mulist(muind);
        acpreindvl=avls.aclist(muind,1);
        acpstindvl=avls.aclist(muind,2);
        acvlspre=vls{ntind}{acpreindvl}(:,2);
        acvlspst=vls{ntind}{acpstindvl}(:,2);
        acvls=[acvlspre ;acvlspst];
        muvls=[vls{ntind}{muindvl}(:,2)];
        edges=avls.edges{ntind};
        hstac=histc(acvls,edges);
        hstacpre=histc(acvlspre,edges);
        hstacpst=histc(acvlspst,edges);
        hstacpre=hstacpre/sum(hstacpre);
        hstmu=histc(muvls,edges);
        hstac=hstac/sum(hstac);
        hstacpst=hstacpst/sum(hstacpst);
        hstmu=hstmu/sum(hstmu);
        avls.hstac_comb{ntind}{muindvl}=hstac;
        avls.hstacpre{ntind}{muindvl}=hstacpre;
        avls.hstmu_comb{ntind}{muindvl}=hstmu
        avls.hstacpst{ntind}{muindvl}=hstacpst;
    end
end
