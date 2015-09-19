%10.9.09
%SCRIPT FOR THIRD VERSION OF STIM FIGURE

%start with sumbs script.

%first figure.

%plot upshift bk59w37(first), plot downshift bk57w35 (second)...with thin
%red lines


%second figure
%plot histograms.
%plot explanatory arrows

%third figure
%plot masterfiguree

%fourth figure
%plot a+b vs. pitch spread.
function [shiftplotstim,shiftplotpharm,combvls,phcombvls]=plothistoryscript2(sumbs,sts,phsumbs);

col.basrev=[1 .6 .6]
col.asyrev=[.6 .6 1]


figure
% figstoplot=[1 2 3 4 5 6]
figstoplot=[2 3 4 5 6]

%settings for first figure (time course)
ps(1).bsnum=[2 ]
ps(1).runstoplot{1}=[1 4 9 10 11 12 13 14 15 16 17]
ps(1).runstoplot{2}=[2 3 4 5 6 7 13  14 15 16 17]

ps(2).bsnum=[2 ]
ps(2).col=col;
ps(2).runstoplot{1}=[4 13]
% ps(2).runstoplot{2}=[3 14]
ps(2).totalplot=2;
ps(2).shiftdrxn=[1 0]
ps(2).arrowht=[0.4 0.45]
ps(2).distht=0.35
ps(2).shiftnum=1;
ps(2).ax(1)=subplot(434)
ps(2).ax(2)=subplot(437)
ps(2).axbnds=[2.25 2.6 0 0.5]

if(ismember(1,figstoplot))
    ps(1).ax=subplot(4,3,1);
    ps(1).axbnds=[-1 9 2.25 2.55]
    plotstim_tmcourse(sumbs,ps(1));
end

if(ismember(2,figstoplot))
    plothists(sts,sumbs,ps(2));
end

if(ismember(3,figstoplot))
   axin=subplot(4,3,[2 5])
   axbnds=[-4 4 0 39]
   PLOTEX=0
    [shiftplotstim,combvls]=plothistorystim5(sumbs,0,0,0,0,axin,axbnds,PLOTEX)
end

if(ismember(4,figstoplot))
   axin=subplot(4,3,[3 6])
   axbnds=[-4 4 0 39]
   PLOTEX=0
    [shiftplotpharm,phcombvls]=plothistorystim5(phsumbs,0,0,0,0,axin,axbnds,PLOTEX)
end


% if(ismember(4,figstoplot))
%     [shiftoffcomb,revoffcomb]=plotstim_revpairs(shiftplot);
%     figure
%     plot(shiftoffcomb,revoffcomb,'ko','MarkerFaceColor','k');
%     hold on;
%     axis([-2.5 2.5 -2.5 2.5])
%     box off
%     plot([-2.5 2.5],[0 0],'k--')
%     plot([0 0],[-2.5  2.5],'k--')
% end

if (ismember(5,figstoplot))
      clear ps
      
      ps.ax(1)=subplot(438)
      ps.ax(2)=subplot(4,3,11);
      ps.col=col;
      ps.edges=[-3.2:.4:3.2]
      ps.axbnds=[-4 4 0 0.5]
      
      plotgrouphsts(combvls,ps);
      box off;
    
end


if (ismember(6,figstoplot))
      clear ps
      
      ps.ax(1)=subplot(439)
      ps.ax(2)=subplot(4,3,12);
      ps.col=col;
      ps.edges=[-3.2:.4:3.2]
     
       ps.axbnds=[-4 4 0 0.5]
      plotgrouphsts(phcombvls,ps);
      box off;
    
end

function []=plotstim_tmcourse(sumbs,ps)
    bsln=length(ps.bsnum);
%      figure
    for ii=1:bsln
        crbs=ps.bsnum(ii);
        ps.STIM=1;
        ps.plotarrow=1;
%         axtmcourse(ii)=subplot(1,bsln,ii);
        plot_tmcourse2rev(sumbs(crbs),ps.runstoplot{ii},ps);
    end

function [axgrphst]=plotgrouphsts(combvls,ps)
%     colvec={'r' 'c' 'k' 'm'}
%     
    axgrphst=ps.ax;
    edges=ps.edges;
    col=ps.col;
    pair{1}=[2 3]
    pair{2}=[1 ]
    
    for ii=1:length(combvls)
        hst_tmp=histc(combvls{ii},edges)
        mnout{ii}=nanmean(combvls{ii})
        ster{ii}=nanstd(combvls{ii})./sqrt(length(combvls));
        hstout{ii}=hst_tmp./sum(hst_tmp);
        ct{ii}=sum(hst_tmp);
    end
    
    for figcnt=1:2
        axes(axgrphst(figcnt))
        crvls=pair{figcnt}
        colcnt=figcnt
        stairs(edges,hstout{crvls(1)},'Color',col.basrev,'Linewidth',3)
        axis(ps.axbnds);
        hold on;
        if(figcnt==1)
            stairs(edges,hstout{crvls(2)},'Color',col.asyrev,'Linewidth',3)
        end
            axis([-4 4 0 0.5])
      
        
        text(3, 0.35, ['n = ' num2str(ct{crvls(1)})],'Color',col.basrev,'Fontsize',16)  
        crmn=mnout{crvls(1)};
        crster=ster{crvls(1)};
        h1=arrow([0 0.35],[crmn 0.35],'Length',3,'Width',0.5)
        set(h1,'FaceColor',col.basrev)
        set(h1,'EdgeColor','none')
        
        
        
%         plot([0 crmn],[0.35 0.35],'Color',col.basrev,'Linewidth',3);
        plot([crmn-crster crmn+crster],[0.3 0.3],'Color',col.basrev,'Linewidth',3);
        if(figcnt==1)
        crmn=mnout{crvls(2)};
        crster=ster{crvls(2)};
         h1=arrow([0 0.35],[crmn 0.35],'Length',3,'Width',0.5)
        set(h1,'FaceColor',col.asyrev)
        set(h1,'EdgeColor','none')
        
%         arrowh([0 crmn],[0.35 0.35],col.asyrev,[200 200]);
%         plot([0 crmn],[0.35 0.35],'Color',col.asyrev,'Linewidth',3);
        plot([crmn-crster crmn+crster],[0.3 0.3],'Color',col.asyrev,'Linewidth',3);
        
        text(3, 0.3, ['n = ' num2str(ct{crvls(2)})],'Color',col.asyrev,'Fontsize',16)  
        end
        end

    
    
    
function []=plothists(sts,sumbs,ps)
    bsln=length(ps.bsnum);
    plotnum=1;
    figure
    
    dht=ps.distht;
    arht=ps.arrowht;
    
    for ii=1:bsln
        crbs=ps.bsnum(ii);
        crsumbs=sumbs(crbs);
       
        allinds=ps.runstoplot{ii};
        cmd=['cd ' sts(crbs).path]
         eval(cmd);
         cmd=['load datsum.mat']
         eval(cmd);
         
        for indnum=1:length(allinds)
            crind=allinds(indnum);
        %load avls 
            pathvl=avls.pvls{crind}
            btvl=avls.cvl{crind};
         
            if(isfield(avls,'baspath'))
                cmd=['cd ' avls.baspath pathvl]
                eval(cmd);
            else
                cmd=['cd ' pathvl]
                eval(cmd);
            end
            cmd=['load ' btvl '.mat']
            eval(cmd);

      %     load extradata.mat
       
    
        axes(ps.ax(plotnum));
        plotnum=plotnum+1;
        stairs(avls.HST_EDGES/3000,hsctnrm,'k','Linewidth',2)
        
        
        hold on;
        axis(ps.axbnds);
        if(ii==1)
            mnbas=ctmean/3;
            stdbas=stdct/3;
        end
       
%    plot([mnbas+stdbas mnbas+stdbas], [0 1], 'c--')
%    plot([mnbas-stdbas mnbas-stdbas], [0 1], 'c--')
     plot([crsumbs.initmean/3000 crsumbs.initmean/3000],[0 1],'Linestyle','--','Color',[0.6 0.6 0.6],'Linewidth',2)
     initpt=crsumbs.initmean+crsumbs.asympvl{ps.shiftnum}*crsumbs.initsd
     if((plotnum-1)==2)
        plot([initpt/3000 initpt/3000],[0 1],'Linestyle','--','Color',[0.6 0.6 0.6],'Linewidth',2)
     end
     tmvec(ii,:)=tms;
     mnvlsct(ii)=ctmean; 
     stdvlsct(ii)=stdct./sqrt(length(crctind));
    text(2.225,0.5,['catch=' num2str(length(crctind))],'Color','k');
   plot([(mnvlsct(ii)-stdvlsct(ii))/3000 (mnvlsct(ii)+stdvlsct(ii))/3000],[dht dht],'k','Linewidth',3)
   
   %this is for the stim trial
   if(~isempty(crfbind))
        mnvlsfb(ii)=fbmean;
        stdvlsfb(ii)=stdfb./sqrt(length(crfbind));
        plot([(mnvlsfb(ii)-stdvlsfb(ii))/3000 (mnvlsfb(ii)+stdvlsfb(ii))/3000],[dht dht],'r','Linewidth',3)
        text(2.225,0.2,['fb=' num2str(length(crfbind))],'Color','r');
        stairs(avls.HST_EDGES/3000,hsfbnrm,'r','Linewidth',2)
        box off;
   end
   if(ps.shiftdrxn(plotnum-1)==1)
       origarrow=[crsumbs.initmean/3000 arht(1) mnvlsct(ii)/3000  arht(1)]
       shiftarrow=[mnvlsct(ii)/3000 arht(2) mnvlsfb(ii)/3000  arht(2)]
       col=ps.col.basrev
   else
        
       origarrow=[initpt/3000 arht(1) mnvlsct(ii)/3000 arht(1)]
        shiftarrow=[mnvlsct(ii)/3000 arht(2) mnvlsfb(ii)/3000 arht(2)]
        col=ps.col.asyrev
   end
   
%    arrow(origarrow,shiftarrow, ps.shiftdrxn(plotnum-1),ps.col);
     h=arrow([origarrow(1) origarrow(2)],[origarrow(3)  origarrow(4)],'Length',8,'Width',1);
            set(h,'FaceColor',[0.6 0.6 0.6]);
            set(h,'EdgeColor','none');
            set(h,'Linewidth',2)
            
             h=arrow([shiftarrow(1) shiftarrow(2)],[shiftarrow(3)  shiftarrow(4)],'Length',8,'Width',1);
            set(h,'FaceColor',col);
            set(h,'EdgeColor','none');
            set(h,'Linewidth',2)
   
   end
        linkaxes(ps.ax);
    end

        
    
    function []=plotarrowh(origarrow,shiftarrow, shiftdrxn,col)
%         axes(axhist);
        if(shiftdrxn==0)
            arrow_col=col.asyrev
        else
            arrow_col=col.basrev
        end
        origlen=origarrow(3)-origarrow(1);
        shiftlen=shiftarrow(3)-shiftarrow(1);
%         headwidth(1)=.01/abs(origlen);
%         headwidth(2)=.01/abs(shiftlen);
%         headht(1)=.01/abs(origlen)
%         headht(2)=.01/abs(shiftlen)
          axis([2.25 2.55 0 0.6])
          axis square

        arrowh([origarrow(1) origarrow(3)],[origarrow(2) origarrow(2)],[0.6 0.6 0.6],[200 200]);
        plot([origarrow(1) origarrow(3)],[origarrow(2) origarrow(4)],'Color',[0.6 0.6 0.6],'Linewidth',2);
        arrowh([shiftarrow(1) shiftarrow(3)],[shiftarrow(2) shiftarrow(2)],arrow_col,[200 200])
        plot([shiftarrow(1) shiftarrow(3)],[shiftarrow(2) shiftarrow(4)],'Color',arrow_col,'Linewidth',2);
      
        
    
    
    
    
    
    
    
    function []=plotarrows(origarrow,shiftarrow, shiftdrxn,axhist)
%         axes(axhist);
        if(shiftdrxn==0)
            arrow_col='r'
        else
            arrow_col='r'
        end
        origlen=origarrow(3)-origarrow(1);
        shiftlen=shiftarrow(3)-shiftarrow(1);
        headwidth(1)=.01/abs(origlen);
        headwidth(2)=.01/abs(shiftlen);
        headht(1)=.01/abs(origlen)
        headht(2)=.01/abs(shiftlen)
        plot_arrow(origarrow(1),origarrow(2),origarrow(3),origarrow(4),'Color','k','Linewidth',3,'headwidth',headwidth(1),'headheight',headht(1),'facecolor','k');
        plot_arrow(shiftarrow(1),shiftarrow(2),shiftarrow(3),shiftarrow(4),'Color',arrow_col,'Linewidth',3,'headwidth',headwidth(2),'headheight',headht(2),'facecolor',arrow_col,'edgecolor',arrow_col)
        axis equal;
        
        
        
