function [stats]=calcsumplotstats(sumplot);
    crvls=sumplot(1)
stats.num_shinactiv=length(crvls.ac);
uniquebs=unique(crvls.bsnum);

stats.numbirds=length(uniquebs)-2;
stats.totshift=0;
for ii=1:length(uniquebs)
    bsvl=uniquebs(ii);
    crind=find(crvls.bsnum==bsvl);
    numshft=length(unique(crvls.shftnum(crind)));
    stats.totshift=stats.totshift+numshft; 
    
end

for ii=1:3
negind=find(sumplot(ii).ac<0)
offvls=sumplot(ii).off;
offvls(negind)=-offvls(negind);
stats.offmn(ii)=mean(offvls);
stats.offsd(ii)=std(offvls);
end


    