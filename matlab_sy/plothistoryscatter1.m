  %purpose of this code is to go through sumbs and calculate
%values for histograms of offsets.
%in this version I take out all paired calculations.
function [shiftplot,combvls]=plothistorystim5(sumbs,minpts,minshift,norm_asymp, STIM,axin,axbnds,axscatter,PLOTSCATTER,colin,MX_TM)

ctvl=1
 
 if(~iscell(colin))   
%first find appropriate stim runs.
% colvec{1}=[1 .6 .6]
% colvec{2}=[1 .6 .6]
% if(~exist(colin)||colin==0)
    colvec{1}{1}=[0.1 0.1 0.1]
    colvec{1}{2}=[0.5 0.5 0.5]
    
    colvec{2}{1}=[0.4 0.4 1]
    colvec{2}{2}=[0.7 0.7 1]
    colvec{3}=[.6 .6 1]
else
    colvec=colin;
end
    for colind=1:2
        for drxn=1:2
        combvls{colind}{drxn}=[]
        end
    end
for bsnm=1:length(sumbs)
   crbs=sumbs(bsnm);
    shiftruns=crbs.shiftruns;
    revruns=crbs.revruns;
    
    numshiftruns=length(shiftruns);
    numrevruns=length(revruns);
    for ii=1:numshiftruns
        asympruns=crbs.asympruns{ii};
        tms=crbs.adjshifttms{ii};
        shiftplot(bsnm).drxn{ii}=crbs.drxn{ii};
        
        for runtype=1:3
            if(runtype==1)
                inruns=crbs.shiftruns;
            elseif(runtype==2)
                inruns=crbs.basruns;
            else
                inruns=crbs.revruns;
            end
            if(ii<=length(inruns))
                if(runtype==2)
                    cr_runs{runtype}=inruns;
                else
                     cr_runs{runtype}=inruns{ii};
                end
                if(STIM)
                    tmpind=find(ismember(cr_runs{runtype},crbs.STANRUNS));
                    cr_runs{runtype}=cr_runs{runtype}(tmpind);
                end
            else
               cr_runs{runtype}=[]; 
            end
        end
    
        combins=[cr_runs{1}' cr_runs{2} ]
        shiftins=cr_runs{1};
        if(STIM)
            aczvls=crbs.acz
            muzvls=crbs.muz
        else
            aczvls=crbs.acz(crbs.ntind,:)
            muzvls=crbs.muz(crbs.ntind,:)
        end
        %DIVIDE DATA INTO THREE GROUPS
        %BASIND(1)  SHIFTIND(2)  REVIND(3)
        %pick the inds for baseline measurements, and shift measurements
        %and reverse to baseline, and reverse shift.
        analshiftinds=find(abs(aczvls(shiftins))>1);
         allinds=find(abs(aczvls(combins))>1);
        analbasinds=setdiff(1:length(combins),allinds);
        
        tminds=find(crbs.adjshifttms{ii}<=MX_TM);
        analshiftinds=intersect(tminds,analshiftinds);
        analshiftinds=shiftins(analshiftinds);
        outinds{1}=combins(analbasinds);
        
        %now remove asympruns from this analysis
        %and remove day max
        
        if(asympruns>1)
           keepind=find(~ismember(analshiftinds,asympruns(2:end))) 
         
           outinds{2}=analshiftinds(keepind); 
        else
            outinds{2}=analshiftinds;
        end
        
        if(~isempty(cr_runs{3}))
%             analrevshiftinds=find(abs(aczvls(revind))>1);
%             analrevbasinds=setdiff(1:length(revind), analrevshiftinds);
%             analrevshiftinds=revind(analrevshiftinds);
            outinds{3}=cr_runs{3};
        else
            outinds{3}=[];
        end
        
        for runvl=1:3
            shiftplot(bsnm).acshift{ii}{runvl}=aczvls(outinds{runvl});
            shiftplot(bsnm).mushift{ii}{runvl}=muzvls(outinds{runvl});
        end
    end
%     if (~isempty('aczpr_vls'))
%         shiftplot(bsnm).aczpr_vls=aczpr_vls
%         shiftplot(bsnm).muzpr_vls=muzpr_vls;
%         shiftplot(bsnm).revacpr_vls=revacpr_vls;
%         shiftplot(bsnm).revmupr_vls=revmupr_vls;
%     end
    if(PLOTSCATTER)

    
    
        plotscatter(shiftplot(bsnm),axscatter,colvec,axbnds);
    end
   
    if(bsnm==1)
%         figure
        axhor(1)=axin;
    end
    [combvls,ctvl]=plotmeanhorizarrow(shiftplot(bsnm),axhor(1),colvec,combvls,ctvl,axbnds);
    ctvl=ctvl+1;
end
% linkaxes(ax);


function []=plotscatter(shiftplot,ax, colvec,axbnds)
axes(ax);
hold on;
axis square
numdrxn=length(shiftplot.drxn);   
for ii=1:numdrxn
    for typevl=1:2
%         length(shiftplot.acshift{ii})
        drxn=shiftplot.drxn{ii};   
        
        y1=shiftplot.acshift{ii}{typevl};
        y2=shiftplot.mushift{ii}{typevl};
               
       if(drxn=='up')
           x=y2;
           yln=y1;
           ypt=y1;
           drxnvl=1;
       else
           x=y2;
           ypt=-y1
           yln=y1;
           drxnvl=2;
       end
%            plot([x;x],[y1;y2],'k','Linewidth',2)
           
            plot(x,ypt,'o','MarkerEdgeColor',colvec{typevl}{drxnvl},'MarkerSize',5,'MarkerFaceColor',colvec{typevl}{drxnvl})
            hold on;
            if(typevl==2)
            for jj=1:length(x)
            plot([yln(jj) x(jj)],[ypt(jj) ypt(jj)],'Color',colvec{typevl}{drxnvl},'Linewidth',2); 
            end
            end
            hold on;
%              axis(axbnds);
            %             plot([mnx1;mnx1],[ysh1;ysh2],'r--','Linewidth',2)
           
%              plot([mnx2;mnx2],[ymnrv1;ymnrv2],'c--','Linewidth',2)
            plot([-8 8],[-8 8],'k--');
            plot([-8 8],[8 -8],'k--');
            plot([-8 8],[0 0],'k--');
            plot([0 0],[-8 8],'k--');
    end
end






function []=plotarrow(shiftplot,ax, colvec,axbnds)
axes(ax);
hold on;
numdrxn=length(shiftplot.drxn);   
for ii=1:numdrxn
    for typevl=1:length(shiftplot.acshift{ii})
        drxn=shiftplot.drxn{ii};   
        
        y1=shiftplot.acshift{ii}{typevl};
        y2=shiftplot.mushift{ii}{typevl};
               
       if(drxn=='up')
           x=y1;
       else
           x=-y1;
       end
%            plot([x;x],[y1;y2],'k','Linewidth',2)
           
            plot(x,y2,'o','MarkerEdgeColor',colvec{typevl},'MarkerSize',3,'MarkerFaceColor',colvec{typevl})
             hold on;
             axis(axbnds);
            %             plot([mnx1;mnx1],[ysh1;ysh2],'r--','Linewidth',2)
           
%              plot([mnx2;mnx2],[ymnrv1;ymnrv2],'c--','Linewidth',2)
            plot([0 8],[0 8],'k--');
            plot([0 8],[0 -8],'k--');
    end
end
   
function [combvls]=plothorizarrow(shiftplot,ax, colvec,combvls,axbnds)
axes(ax);
hold on;
numdrxn=length(shiftplot.drxn);   
for ii=1:numdrxn
    
    for typevl=1:2
%         length(shiftplot.acshift{ii})
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
                if(typevl==1)
                    combvls{1}{ii}=[combvls{1}{ii} combdiff]
                elseif(typevl==2)
                    combvls{2}{ii}=[combvls{2}{ii} combdiff]
                elseif(typevl==3)
                    combvls{3}=[combvls{3} combdiff]
                end
                plot([0 diff],[ct ct],'k')
                hold on;
                plot(diff,ct,'o','MarkerEdgeColor',colvec{typevl},'MarkerSize',4,'MarkerFaceColor',colvec{typevl})
                
            end
        end
    end
    ctvl=ctvl+1;
end
function [combvls,ctvl]=plotmeanhorizarrow(shiftplot,ax, colvec,combvls,ctvl,axbnds)
axes(ax);
hold on;
numdrxn=length(shiftplot.drxn);  
% plot([-4 4],[ctvl-.3 ctvl-.3],'Color',[0.6 0.6 0.6])
for ii=1:numdrxn
    
    for typevl=1:length(shiftplot.acshift{ii})
        drxn=shiftplot.drxn{ii};   
        numvls=length(shiftplot.acshift{ii}{typevl});
        if(numvls)
                y1=shiftplot.acshift{ii}{typevl};
                y2=shiftplot.mushift{ii}{typevl};
                diff=y2-y1;
                if(drxn=='up')
                   combdiff=diff
                   drxnvl=1;
                else
                    combdiff=diff;
                    drxnvl=2;
                end
                nonaninds=find(~isnan(combdiff))
                if(typevl==1)
                    combvls{1}{drxnvl}=[combvls{1}{drxnvl} combdiff(nonaninds)]
                elseif(typevl==2)
                    combvls{2}{drxnvl}=[combvls{2}{drxnvl} combdiff(nonaninds)]
                elseif(typevl==3)
%                     combvls{3}=[combvls{3} combdiff(nonaninds)]
                end
   
                if(typevl>1)
%                     plot([0 mean(combdiff)],[ctvl ctvl],'k')
%                     hold on;
%                     plot(mean(combdiff),ctvl,'o','MarkerEdgeColor',colvec{typevl},'MarkerSize',6,'MarkerFaceColor',colvec{typevl})
                    hw=.35/abs(mean(combdiff(nonaninds)));
                    hh=.15/abs(mean(combdiff(nonaninds)))
           
                    axis(axbnds)
                    hold on;
%                     h=arrow([0 ctvl],[mean(combdiff(nonaninds)) ctvl],'Length',6,'width',1)
%                    set(h,'FaceColor',colvec{typevl});
%                     set(h,'EdgeColor','none');
%                 set(h,'Linewidth',2)
%                     xvec=combdiff(nonaninds)
%                     yvec=[(ctvl+.2)*ones(length(combdiff(nonaninds)),1)]
%                     plot(xvec,yvec,'o','Color',colvec{typevl},'MarkerSize',4)
                    ctvl=ctvl+.4
                end
            
        end
    end
    ctvl=ctvl+1;
end


