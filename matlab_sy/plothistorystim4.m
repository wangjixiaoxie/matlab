  %purpose of this code is to go through sumbs and calculate
%values for histograms of offsets.
%in this version I take out all paired calculations.
function [shiftplot,combvls]=plothistorystim(sumbs,minpts,minshift,norm_asymp, STIM)
figure
%first find appropriate stim runs.
colvec={'k' 'r' 'c' 'm'}
    for colind=1:length(colvec)
        combvls{colind}=[]
    end
for bsnm=1:length(sumbs)
   crbs=sumbs(bsnm);
   
    shiftruns=crbs.shiftruns;
    revruns=crbs.revruns;
    numshiftruns=length(shiftruns);
    numrevruns=length(revruns);
    
    %check to see whether runs are STANRUNS

    %THIS IS A BUG.
    for ii=1:numshiftruns
        shiftplot(bsnm).drxn{ii}=crbs.drxn{ii};
        shiftruns=crbs.shiftruns{ii};
        asympruns=crbs.asympruns{ii};
        if(STIM)
            tmpind=find(ismember(shiftruns,crbs.STANRUNS));
            shiftins1=shiftruns(tmpind);
        else
            shiftins1=shiftruns;
        end
        basruns=crbs.basruns;
        if(STIM)
            basind=find(ismember(basruns, crbs.STANRUNS));
            basind=basruns(basind);
        else
            basind=basruns
        end
        if(ii<=length(crbs.revruns))
            revruns=crbs.revruns{ii};
            if(STIM)
                revind=find(ismember(revruns,crbs.STANRUNS));
                revind=revruns(revind);
            else
                revind=revruns;
            end
            if(length(revind)>3)
                revind=revind(1:4)
            end
        else
            revind=[];
        end
        
        
        %include baseline runs in shiftins
        %limit revinds to the first three days after wn off...do this for
        %the plotting as well.
%         findpairs(revinds,shiftins);
    
        combins=[shiftins1' basind ]
%         cr_revind=revind(revct);
%         basins=basind;
%         craczrev=crbs.acz(cr_revind);
%         indshift=find(abs(crbs.acz(combins)-craczrev)<1);
%     %get average of all differences within one standard deviation
%     
%         if(indshift)
% %             combins=combins(indshift);
%             aczpr_vls{ii,revct}=crbs.acz(combins);
%             muzpr_vls{ii,revct}=crbs.muz(combins);
% %             aczbaspr_vls{ii,revct}=crbs.acz(basins);
% %             muzbaspr_vls{ii,revct}=crbs.muz(basins);
%             revacpr_vls{ii,revct}=craczrev;
%             revmupr_vls{ii,revct}=crbs.muz(cr_revind);
%         end
   
        if(STIM)
            aczvls=crbs.acz
            muzvls=crbs.muz
        else
            aczvls=crbs.acz(crbs.ntind,:)
            muzvls=crbs.muz(crbs.ntind,:)
        end
    
        %pick the inds for baseline measurements, and shift measurements
        %and reverse to baseline, and reverse shift.
        analshiftinds=find(abs(aczvls(combins))>1);
        analbasinds=setdiff(1:length(combins),analshiftinds);
        analshiftinds=combins(analshiftinds);
        analbasinds=combins(analbasinds);
        
        %now remove asympruns from this analysis
        if(asympruns>1)
           keepind=find(~ismember(analshiftinds,asympruns(2:end))) 
           analshiftinds=analshiftinds(keepind); 
        end
        
        if(~isempty(revind))
            analrevshiftinds=find(abs(aczvls(revind))>1);
            analrevbasinds=setdiff(1:length(revind), analrevshiftinds);
            analrevshiftinds=revind(analrevshiftinds);
            analrevbasinds=revind(analrevbasinds);
        end
%          mnout(bsnm).acsh{ii}=mean(crbs.acz(analshiftinds));
%         if(~isempty(analshiftinds))
%            
% %             mnout(bsnm).mush{ii}=mean(crbs.muz(analshiftinds));
% %              mnout(bsnm).acrev{ii}=mean(crbs.acz(analrevshiftinds));
% %             mnout(bsnm).murev{ii}=mean(crbs.muz(analrevinds));
%         end
        
        
        shiftplot(bsnm).acshift{ii}{1}=aczvls(analbasinds);
        shiftplot(bsnm).mushift{ii}{1}=muzvls(analbasinds);
        shiftplot(bsnm).acshift{ii}{2}=aczvls(analshiftinds);
        shiftplot(bsnm).mushift{ii}{2}=muzvls(analshiftinds);
        if(~isempty(revind))
            shiftplot(bsnm).acshift{ii}{4}=aczvls(analrevbasinds);
            shiftplot(bsnm).mushift{ii}{4}=muzvls(analrevbasinds);
            shiftplot(bsnm).acshift{ii}{3}=aczvls(analrevshiftinds);
            shiftplot(bsnm).mushift{ii}{3}=muzvls(analrevshiftinds);
        end
    end
%     if (~isempty('aczpr_vls'))
%         shiftplot(bsnm).aczpr_vls=aczpr_vls
%         shiftplot(bsnm).muzpr_vls=muzpr_vls;
%         shiftplot(bsnm).revacpr_vls=revacpr_vls;
%         shiftplot(bsnm).revmupr_vls=revmupr_vls;
%     end
        ax(bsnm)=subplot(2,length(sumbs),bsnm);
        axis([-1 7 -7 7]);
        axis square;
        box off;
    
    
    plotarrow(shiftplot(bsnm),ax(bsnm),colvec);
   
    axhor(bsnm)=subplot(2,length(sumbs),bsnm+length(sumbs));
    [combvls]=plothorizarrow(shiftplot(bsnm),axhor(bsnm),colvec,combvls);
    
end
linkaxes(ax);

function []=plotarrow(shiftplot,ax, colvec)
axes(ax);
hold on;
numdrxn=length(shiftplot.drxn);   
for ii=1:numdrxn
    for typevl=1:length(shiftplot.acshift{ii})
        drxn=shiftplot.drxn{ii};   
        
        y1=shiftplot.acshift{ii}{typevl};
        y2=shiftplot.mushift{ii}{typevl};
        
        
%         ymnrv1=mnout.acrev{ii}
%         ymnrv2=mnout.murev{ii}
%         ysh1=mnout.acsh{ii};
%         ysh2=mnout.mush{ii};
       
       if(drxn=='up')
           x=y1;
      
       else
           x=-y1;
       end
           plot([x;x],[y1;y2],'k','Linewidth',2)
            hold on;
            plot(x,y2,'o','MarkerEdgeColor',colvec{typevl},'MarkerSize',6,'MarkerFaceColor',colvec{typevl})
%             plot([mnx1;mnx1],[ysh1;ysh2],'r--','Linewidth',2)
           
%              plot([mnx2;mnx2],[ymnrv1;ymnrv2],'c--','Linewidth',2)
            plot([0 8],[0 8],'k--');
            plot([0 8],[0 -8],'k--');
    end
end
   
function [combvls]=plothorizarrow(shiftplot,ax, colvec,combvls)
axes(ax);
hold on;
numdrxn=length(shiftplot.drxn);   
ct=1;
for ii=1:numdrxn
    
    for typevl=1:length(shiftplot.acshift{ii})
        drxn=shiftplot.drxn{ii};   
        numvls=length(shiftplot.acshift{ii}{typevl});
        if(numvls)
            for vlind=1:numvls
                y1=shiftplot.acshift{ii}{typevl}(vlind);
                y2=shiftplot.mushift{ii}{typevl}(vlind);
                diff=y2-y1;
                if(drxn=='up')
                   combdiff=diff 
                else
                    combdiff=-diff;
                end
                if(typevl==1)&(ii==1)
                    combvls{1}=[combvls{1} combdiff]
                elseif(typevl==2)
                    combvls{2}=[combvls{2} combdiff]
                elseif(typevl==3)
                    combvls{3}=[combvls{3} combdiff]
                elseif(typevl==4)
                    combvls{4}=[combvls{4} combdiff]
                end
            combvls
                
                plot([0 diff],[ct ct],'k')
                hold on;
                plot(diff,ct,'o','MarkerEdgeColor',colvec{typevl},'MarkerSize',6,'MarkerFaceColor',colvec{typevl})
                ct=ct+1;
            end
        end
    end
    
end

        
%         ymnrv1=mnout.acrev{ii}
%         ymnrv2=mnout.murev{ii}
%         ysh1=mnout.acsh{ii};
%         ysh2=mnout.mush{ii};
       
     








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

        

