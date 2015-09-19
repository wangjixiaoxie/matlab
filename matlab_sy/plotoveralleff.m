function [ax]=plotoveralleff(sumbs,type)
 
col={'k' 'r'}


xcnt=0;
for ii=1:length(sumbs)
    if(ii>1)
        xcnt=xcnt+3;
    end
    ntlen=length(sumbs(ii).allnote);
    if ntlen>2
        ntlen=2;
    end
    for jj=1:ntlen
        
        ntvl=sumbs(ii).allnote(jj);
        effvls=sumbs(ii).combeff(ntvl,:);
   
    
        ind=find(effvls>0);
        logeffvls=log2(effvls(ind));
        mnvls=mean(logeffvls);
        stdvls=std(logeffvls);
        xvl=xcnt+(jj-1);
        
        xvls=xvl*ones(length(ind),1);
        plot(xvls,logeffvls,'o','Color',col{jj});
        hold on;
        bar(xvl,mnvls,'FaceColor','none','EdgeColor', col{jj});
        plot([xvl xvl], [mnvls-stdvls mnvls+stdvls],'Color',col{jj});
    end
end