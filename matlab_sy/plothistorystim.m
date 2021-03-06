%purpose of this code is to go through sumbs and calculate
%the offz for each initrun
%for each control run
%and for each baseline runs
%save to a masterstruct suminitshs
%with initacz, initmuz, initoff, init_bsvl, 
%modified 9.19.09 to return sumdynall and sumdynasymp

%this option tells whether to plot initind, selected as consecutive mu runs
%over some threshold.
function [shiftplot,mnout]=plothistorystim(sumbs,minpts,minshift,norm_asymp)

%first find appropriate stim runs.

for bsnm=1:length(sumbs)
   crbs=sumbs(bsnm);
    shiftruns=crbs.shiftruns;
    revruns=crbs.revruns;
    
    numshiftruns=length(shiftruns);
    numrevruns=length(revruns);
    
    if(numshiftruns<=numrevruns)
        numruns=numshiftruns 
    else
        numruns=numrevruns
    end
    
    for ii=1:numruns
        shiftplot(bsnm).drxn{ii}=crbs.drxn{ii};
        shiftruns=crbs.shiftruns{ii};
        tmpind=find(ismember(shiftruns,crbs.STANRUNS));
        shiftins=shiftruns(tmpind);
        tmpind=find(ismember(shiftins,crbs.asympruns{ii}))
        shiftins=shiftins(1:tmpind(1));
        
        revruns=crbs.revruns{ii};
        revind=find(ismember(revruns,crbs.STANRUNS));
        revind=revruns(revind);
        if(length(revind)>3)
            revind=revind(1:3)
        end
        shiftplot(bsnm).acshift{ii}=crbs.acz(shiftins);
        shiftplot(bsnm).mushift{ii}=crbs.muz(shiftins);
        
        
        
        analshiftinds=find(abs(crbs.acz(shiftins))>1);
        analshiftinds=shiftins(analshiftinds);
        analrevinds=find(abs(crbs.acz(revind))>1);
        analrevinds=revind(analrevinds);
        if(~isempty(analshiftinds))
            mnout(bsnm).acsh{ii}=mean(crbs.acz(analshiftinds));
            mnout(bsnm).mush{ii}=mean(crbs.muz(analshiftinds));
             mnout(bsnm).acrev{ii}=mean(crbs.acz(analrevinds));
            mnout(bsnm).murev{ii}=mean(crbs.muz(analrevinds));
        end
            
        shiftplot(bsnm).acrev{ii}=crbs.acz(revind);
        shiftplot(bsnm).murev{ii}=crbs.muz(revind);
    end
    
    ax(bsnm)=subplot(1,length(sumbs),bsnm);
    plotarrow(shiftplot(bsnm),mnout(bsnm),ax(bsnm));

end


function []=plotarrow(shiftplot,mnout,ax)
numdrxn=length(shiftplot.drxn);   
for ii=1:numdrxn
    drxn=shiftplot.drxn{ii};   
       yshvl1=shiftplot.acshift{ii};
       yshvl2=shiftplot.mushift{ii};
       yrvvl1=shiftplot.acrev{ii};
       yrvvl2=shiftplot.murev{ii};
       ymnrv1=mnout.acrev{ii}
       ymnrv2=mnout.murev{ii}
       ysh1=mnout.acsh{ii};
       ysh2=mnout.mush{ii};
       
       if(drxn=='up')
           xvl1=yshvl1;
           xvl2=yrvvl1
           mnx1=ysh1
           mnx2=ymnrv1
           
       else
           xvl1=-yshvl1;
           xvl2=-yrvvl1
           mnx1=-ysh1
           mnx2=-ymnrv1;
       end
           plot([xvl1;xvl1],[yshvl1;yshvl2],'k','Linewidth',2)
            hold on;
            plot([xvl1],[yshvl2],'ro','MarkerSize',4,'MarkerFaceColor','r')
            plot([mnx1;mnx1],[ysh1;ysh2],'r--','Linewidth',2)
           
            plot([xvl2;xvl2],[yrvvl1;yrvvl2],'b','Linewidth',2)
            hold on;
            plot([xvl2],[yrvvl2],'co','MarkerSize',4,'MarkerFaceColor','c')
             plot([mnx2;mnx2],[ymnrv1;ymnrv2],'c--','Linewidth',2)
            plot([0 8],[0 8],'k--');
            plot([0 8],[0 -8],'k--');
end
       
%        
%        for tgnum=1:length(tag)
%            crbs=tag(tgnum).bsnum;
%            crcol=tag(tgnum).col;
%            taginds=find(plotvls(ii).bsnum==crbs);
%            if(~isempty(taginds))
%                if(ii==1) 
%                     plot([xvl(taginds);xvl(taginds)],[yvl1(taginds);yvl2(taginds)],'Color',crcol,'Linewidth',3)
%                else
%                     plot([-xvl(taginds);-xvl(taginds)],[yvl1(taginds);yvl2(taginds)],'Color',crcol,'Linewidth',3)
%                end
%            end
%        end
%        
%             hold on;
%    end
%    axis([0 5 -5 5])
%    axis square;
%    plot([0 5],[0 5],'k','Linewidth',2)
%    hold on;
%    plot([0 5],[0 -5],'k','Linewidth',2)
%    


%get shiftlength
    
    %get revlength
    
    %now loop through each of the runs and assign 
%     %points, and direction
%     %output to something to plot
%     
%     
%     %calc a and b 
%     
%     
%     end
%         for shnm=1:length(runs)
%             shiftind=runs{shnm}
%             tmpind=find(ismember(shiftind,crbs.STANRUNS));
%             shiftind=shiftind(tmpind);
%             if(runvl==3)
%                 if(~isempty(shiftind))
%                     mulistind=find(crbs.STANRUNS==shiftind(1));
%                     shiftind=[crbs.STANRUNS(mulistind-1); shiftind];
%                 end
%           
%             end
% %             ind=find(crbs.mun(trgnote,shiftind)>minpts&abs(crbs.acz(trgnote,shiftind))>minshift);
% %             if(~isempty(ind))
%                 ct=ct+1;
%                 overallind=shiftind;
%                 ind=1:length(shiftind);
%                 if(~isempty(overallind))
% %                 if(exist('norm_asymp'))
% %                     sumdynasymp(ct).tms=crbs.asympshifttms{shnm}(ind)-crbs.asympshifttms{shnm}(ind(1))+1;
% %                 else
% %                     sumdynasymp(ct).tms=crbs.asympshifttms{shnm}(ind)-crbs.adjshifttms{shnm}(ind(1))+1;
% %                 end
% %                 sumout(ct).pct=crbs.pct(overallind);    
% %                 sumout(ct).off=crbs.offz(trgnote,overallind);
% %                 if(crbs.drxn{shnm}=='up')
% %                     sumout(ct).normoff=-crbs.offz(trgnote,overallind);
% %                 else
% %                     sumout(ct).normoff=crbs.offz(trgnote,overallind);
% %                 end
%                     sumout(ct).acz=crbs.acz(overallind);
%                     sumout(ct).muz=crbs.muz(overallind);
%                     sumout(ct).rawtms=crbs.tmvec(overallind,1)-crbs.baswntime{shnm};
% %                 sumout(ct).targeff=crbs.combeff(trgnote,overallind);
% % %                 sumout(ct).logtargeff=log2(sumdyn(ct).targeff);
% %                 sumout(ct).contreff=crbs.combeff(ctrlnote,overallind);
% %                 sumout(ct).logcontreff=log2(sumdyn(ct).contreff);
%                     sumout(ct).bsnum=bsnm
%                     if(runvl==1)
%                         sumout(ct).adjtms=crbs.adjshifttms{shnm}(ind);
%                         sumout(ct).muasympdist=crbs.muasympdist{shnm};
%                         sumout(ct).acasympdist=crbs.acasympdist{shnm};
% %                     sumout(ct).targeff=crbs.combeff(trgnote,overallind);
% % %                     sumout(ct).logtargeff=log2(sumdyn(ct).targeff);
% %                     sumout(ct).contreff=crbs.combeff(ctrlnote,overallind);
% %                     sumout(ct).logcontreff=log2(sumdyn(ct).contreff);
%                         sumout(ct).bsnum=bsnm
%                     elseif(runvl==3)
%                         if(~isempty(overallind))
%                             sumout(ct).adjtms=crbs.flrtmvec(overallind)-crbs.flrtmvec(overallind(1))+1;
%                         end
%                     else
%                         sumout(ct).adjtms=crbs.asympshifttms{shnm}(ind);
%                
%                     end
% %             end
%                 end
%         end
%         
% end
%     if(runvl==2)
%             sumdynasymp=sumout;
%     elseif(runvl==1)
%             sumdyn=sumout;
%     else
%         sumdynrev=sumout;
%     end
% end

        

