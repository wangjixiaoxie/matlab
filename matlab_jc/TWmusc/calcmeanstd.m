function [avls]=calcmeanstd(avls);
vls=avls.adjvls;
% if(~isempty(avls.diron))
%     dirvls=avls.diradjvls;
% end
for ntind=1:length(avls.NT);
    avls.mnvl{ntind}=[];
    avls.stdv{ntind}=[];
    avls.dirmnvl{ntind}=[];
    avls.dirstdv{ntind}=[];
    
    for btind=1:length(vls{ntind})
        tmpmean=mean(vls{ntind}{btind}(:,2))
        tmpstd=std(vls{ntind}{btind}(:,2))
        avls.mnvl{ntind}=[avls.mnvl{ntind} tmpmean];
        avls.stdv{ntind}=[avls.stdv{ntind} tmpstd];
        
    end
    
%     if(~isempty(avls.diron))
%         for btind=1:length(dirvls{ntind})
%             tmpmean=mean(dirvls{ntind}{btind}(:,2))
%             tmpstd=std(dirvls{ntind}{btind}(:,2))
%             avls.dirmnvl{ntind}=[avls.dirmnvl{ntind} tmpmean];
%             avls.dirstdv{ntind}=[avls.dirstdv{ntind} tmpstd];
%         end
%     end
    


end
