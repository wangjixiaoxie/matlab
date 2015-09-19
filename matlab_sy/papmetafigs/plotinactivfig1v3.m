

%plotstimfig1.
%rewritten 12.13.09 to streamline
function [shiftplot,combvls,sumdata]=plotinactivfig1v3(sumbs,bs)
figstoplot=[ 2]

%settings for the summary data
GROUPHST=0;
BARPLOT=1;
%settings for first figure (time course)


  avlspath='/oriole/bk48w74/datasum'
  avlsmat='datsum.mat'
    
  cmd=['cd ' avlspath];
  eval(cmd);
  cmd=['load ' avlsmat];
  eval(cmd);

% ps(1).bsnum=[2 1]
% ps(1).runstoplot{1}=[1 4 9 10 11 12 13 14 15 16 17]
% ps(1).runstoplot{2}=[2 3 4 5 6 7 13  14 15 16 17]
% 
% ps(2).bsnum=[2 ]
% ps(2).runstoplot{1}=[4 13]
% % ps(2).runstoplot{2}=[3 14]
% ps(2).totalplot=2;
% ps(2).shiftdrxn=[1 0]
% ps(2).arrowht=[0.4 0.45]
% ps(2).distht=0.35
% ps(2).shiftnum=1;



if (ismember(1,figstoplot))
     clear ps
    
     ps.col{1}='k'
     ps.col{2}=[0.4 0.4 1]
     ps.col{3}=[0.6 0.6 0.6]
     
       ps.col{3}=[0.6 0.6 0.6]
       
       ps.shift_col=[0.6 0.6 0.6]
       ps.rev_col=[1 .6 .6]
     
     ps.arrowht=[0.4 0.45]
    ps.distht=0.65
    ps.ntind=1;

     indtoplot=[3 8  6  ]
    sumbsvl=[5  5 8 ]
    plotpre=[0 0 0 ]
    omitpost=[1 1  1]
    plot_triangles=[1 1 1]
    plotshiftarrow=[0 1 0 1 0]
    plotrevarrow=[0 0 1 0 1]
    
     for ii=1:length(indtoplot)
         avls=getavls(bs,sumbsvl(ii))
         ps.indtoplot=indtoplot(ii);
         ps.plotshiftarrow=plotshiftarrow(ii);
         ps.plotrevarrow=plotrevarrow(ii);
         ps.plotpre=plotpre(ii);
         ps.plot_triangles=plot_triangles(ii);
         ps.omitpost=omitpost(ii);
        subplot(5,2,(ii-1)*2+1)
         ps.initmean=sumbs(sumbsvl(ii)).initmean
         axhist(ii)=gca();
         ps.ax=gca();
         ps.axbnds= [2 2.7 0 0.7]
         
       plotexhistinactiv(avls,ps )
     end
     linkaxes(axhist);
    
end

if (ismember(2,figstoplot))
    subplot(5,2,[2 4]);
    ps.axscatter=gca();
    ps.axin=gca();
    ps.axbnds=[0 10 0 10]
    ps.MX_TM=4;
    ps.colin=0;
    ps.PLOTSCATTER=1;
    [shiftplot,combvls]=plothistoryscatter2(sumbs(1:9),0,0,0,0,axin,axbnds,axscatter,1,0,MX_TM) 
    sumdata.MX_TM=MX_TM;
    axis([-4 4 -1 4])
    axis square;
    
    
end

%SUMMARY PLOT
if (ismember(3,figstoplot))
    if(GROUPHST) 
      clear ps
     col{1}{1}=[0.1 0.1 0.1]
    col{1}{2}=[0.5 0.5 0.5]
    col{2}{1}=[0.4 0.4 1]
    col{2}{2}=[0.7 0.7 1]


      ps.ax(1)=subplot(526)
      ps.ax(2)=subplot(5,2,8);
      ps.col=col;
      ps.edges=[-3.3:.25:3.3]
     
       ps.axbnds=[-3 3 0 0.7]
      plotgrouphsts(combvls,ps);
      box off;
    elseif(BARPLOT)
        
        
        
        
    end
    
end


    if(ismember(4,figstoplot))
    ptsps.marksize=2.5
    ps.distht=19.5;
    ps.SCALEHIST=1;
    ps.SCALEGAIN=8;
    ps.omitzero=1;
    ps.SCALEOFFSET=15;
    ps.FLIP=1;
    ps.arrowht=[.65 .65];
    ptsps.rawvl='raw'
    ptsps.ntvl=1;
    ps.ntind=1;
    ps.plot_triangles=1;
    omitpost=[1 1  1]
   
    indtoplot=[3 8 6]
    sumbsvl=[5 5 8]
    xbndsrawpts{1}=[0 10]
    ptsps.STIM=0
    ptsps.plotextra=0
    plotpre=[0 0 0 1 0]
    
     ptsps.col{1}='k'
    ptsps.col{2}='k'
    ptsps.col{3}='r'
    ps.col=ptsps.col
    ps.col{2}='r';
    ptsps.initmean=sumbs(5).initmean;
    ps.initmean=sumbs(5).initmean;
    sfact=3000;
    figure
    
    
  
    
    for ii=1:length(indtoplot)
        
        avls=getavls(bs,sumbsvl(ii))
        ptsps.indtoplot=indtoplot(ii);
        ps.indtoplot=indtoplot(ii);
        ptsps.plotpre=plotpre(ii);
        ps.plotpre=plotpre(ii);
        ptsps.omitpost=omitpost(ii);
        ps.omitpost=omitpost(ii);
        axraw(ii)=subplot(1,3,ii)
        ptsps.ax=gca();
        ps.ax=gca();
        inactiv_rawpoints(avls,ptsps,sfact);
        axis square
        box off
        
         plotexhistinactiv(avls,ps);
    end
   
    linkaxes(axraw);
    axis([-4 24 2.05 2.75])
    
   
    
end

     
   


function []=plotextmcourse(sumbs,indlist,ax)
    plot_tmcourse2(sumbs,indlist,ax)



%TAKEN FROM PLOTHISTORYFIGSSCRIPT2, REWRITTEN IN GENERAL FORM
function []=plotexhist(avls,ps)

    dht=ps.distht;
    arht=ps.arrowht;
    
        for indnum=1:length(ps.indtoplot)
            crind=ps.indtoplot(indnum);
            axes(ps.ax);
            crct=avls.hsctnrm{crind};
            crfb=avls.hsfbnrm{crind};
            crctmn=avls.ctmean(crind)
            crfbmn=avls.fbmean(crind);
            stdct=avls.stdct(crind);
            stdfb=avls.stdfb(crind);
            crctind=avls.crctind{crind};
            crfbind=avls.crfbind{crind};
            
            crctste=stdct./sqrt(length(crctind));
            crfbste=stdfb./sqrt(length(crfbind));
            %plot catch value
            
            stairs(avls.HST_EDGES/3000,crct,'k','Linewidth',2)
            hold on;
       
            %    plot([mnbas+stdbas mnbas+stdbas], [0 1], 'c--')
            %    plot([mnbas-stdbas mnbas-stdbas], [0 1], 'c--')
            plot([ps.initmean/3000 ps.initmean/3000],[0 1],'k--','Linewidth',2)
            
            text(2.225,0.5,['catch=' num2str(length(crctind))],'Color','k');
            plot([(crctmn-crctste)/3000 (crctmn+crctste)/3000],[dht dht],'k','Linewidth',3)
   
             %this is for the stim trial
            if(~isempty(crfbind))
        
                plot([(crfbmn-crfbste)/3000 (crfbmn+crfbste)/3000],[dht dht],'r','Linewidth',3)
                text(2.225,0.2,['fb=' num2str(length(crfbind))],'Color','r');
                stairs(avls.HST_EDGES/3000,crfb,'r','Linewidth',2)
            box off;
            end
  
       origarrow=[ps.initmean/3000 arht(1) crctmn/3000  arht(1)]
       shiftarrow=[crctmn/3000 arht(2) crfbmn/3000  arht(2)]
   
   
   plotarrows(origarrow,shiftarrow);
   
   
   end
       
    
  


function [ax1,ax2]=plotexspect(exsong)
clear catchvl
clear ax

figure;
ax1(1)=subplot(12,1,2:6)
exsong(1).ax=ax1(1);
ax1(2)=subplot(12,1,1);

exsong(2)=exsong(1)
ax2(1)=subplot(12,1,8:12);
exsong(2).ax=ax2(1);
ax2(2)=subplot(12,1,7);

for ii=1:2
    [dat,fs]=plotcbin(exsong(ii));
    hold on;
    plot([0 2.5],[6600 6600],'c--')
end
% 
%     tms=0:1/fs:length(dat(:,2))/fs
%     tms2=tms(1:end-1);
%     ax2(2)=subplot(6,1,1)
%     plot(tms2,dat(:,2),'k','Lininactivtw5ewidth',2)

%plot stim component
for plotnum=1:2
    if plotnum==1
        axes(ax1(2));
    else
        axes(ax2(2));
    end
    %this is for first one. 
    cmd=['cd ' exsong(1).path]
    eval(cmd);
    cmd=['load ' exsong(1).cvl '.mat']
    eval(cmd);
    cmd=['cd ' exsong(1).stim_path]
    eval(cmd);
    cmd=['load ' exsong(1).stim_mat]
    eval(cmd);
    onstime(1)=fvst(2).ons(fvst(2).ind)
    onstime(2)=fvst(3).ons(fvst(3).ind)
    onstime(3)=fvst(4).ons(fvst(4).ind)
    % onstime(3)=fvst(9).ons(fvst(9).ind)
    stimofftime(1)=fvst(2).STIMTIME;
    catchvl(1)=fvst(2).STIMCATCH
    stimofftime(2)=fvst(3).STIMTIME;
    catchvl(2)=fvst(3).STIMCATCH
    stimofftime(3)=fvst(4).STIMTIME;
    catchvl(3)=fvst(3).STIMCATCH;

    dtime(1)=onstime(1)/1000+stimofftime(1)/1000-exsong(1).bnds(1);
    dtime(2)=onstime(2)/1000+stimofftime(2)/1000-exsong(1).bnds(1);
    dtime(3)=onstime(3)/1000+stimofftime(3)/1000-exsong(1).bnds(1);

    hold on;
    for ii=1:3
        plot([dtime(ii) dtime(ii)],[0 20000],'r','Linewidth',2)
    end

    for ii=1:3
        if catchvl(ii)<1
        plot(dtime(ii)+.065+stimtms,dat_stimnew,'k','Linewidth',2)
        end
    end
end
axes(ax1(2))
linkaxes(ax1,'x')
axis([0 2.63 0 2e4])

axes(ax2(2))
linkaxes(ax2,'x')
axis([1.34 1.5 0 2e4])

st_tm=onstime(2)/1000+.077-.016-exsong(1).bnds(1);
end_tm=st_tm+512/32000;
axes(ax2(1))
plotbox([st_tm end_tm 6000 8000],'c')


%this is for second one.
% axes(ax2(2));
% fvst=fvstex{2}
% %these happen to be 19 and 20.
% onstime(1)=fvst(40).ons(fvst(40).ind)
% onstime(2)=fvst(41).ons(fvst(41).ind)
% onstime(3)=fvst(42).ons(fvst(42).ind)
% % onstime(3)=fvst(9).ons(fvst(9).ind)
% stimofftime(1)=fvst(40).STIMTIME;
% catchvl(1)=fvst(40).STIMCATCH
% stimofftime(2)=fvst(41).STIMTIME;
% catchvl(2)=fvst(41).STIMCATCH;
% stimofftime(3)=fvst(42).STIMTIME;
% catchvl(3)=fvst(42).STIMCATCH
% dtime(1)=onstime(1)/1000+stimofftime(1)/1000-exsong(2).bnds(1);
% dtime(2)=onstime(2)/1000+stimofftime(2)/1000-exsong(2).bnds(1);
% dtime(3)=onstime(3)/1000+stimofftime(3)/1000-exsong(2).bnds(1);
% 
% hold on;
% for ii=1:3
%     plot([dtime(ii) dtime(ii)],[0 20000],'r','Linewidth',2)
% end
% 
% for ii=1:3
%     if catchvl(ii)<1
%     plot((dtime(ii)+.065+stimtms),stimblock,'k','Linewidth',2)
%     end
% end


title('plotstimexample.m')


function [axgrphst]=plotgrouphsts(combvls,ps)
%     colvec={'r' 'c' 'k' 'm'}
%     
    axgrphst=ps.ax;
    edges=ps.edges;
    col=ps.col;
    jjpair{1}=[1 2]
    jjpair{2}=[1 2]
    
    for ii=1:length(combvls)
        for jj=1:length(combvls{ii})
            hstvls=combvls{ii}{jj}.offz;
            hst_tmp=histc(hstvls,edges)
            mnout{ii}{jj}=nanmean(hstvls)
            ster{ii}{jj}=nanstd(hstvls)./sqrt(length(hstvls));
            hstout{ii}{jj}=hst_tmp./sum(hst_tmp);
            ct{ii}{jj}=sum(hst_tmp);
        end
    end
    
    for figcnt=1:2
        axes(axgrphst(figcnt))
        crvls=jjpair{figcnt}
        colcnt=figcnt
        stairs(edges,hstout{figcnt}{crvls(1)},'Color',col{figcnt}{crvls(1)},'Linewidth',3)
        axis(ps.axbnds);
        hold on;
%         if(figcnt==1)
            stairs(edges,hstout{figcnt}{crvls(2)},'Color',col{figcnt}{crvls(2)},'Linewidth',3)
%         end
            axis([-4 4 0 0.5])
    
        
        text(3, 0.35, ['n = ' num2str(ct{figcnt}{crvls(1)})],'Color',col{figcnt}{crvls(1)},'Fontsize',16)  
        for hstcount=1:2
%         
%         h1=arrow([0 0.35],[crmn 0.35],'Length',3,'Width',0.5)
%         set(h1,'FaceColor',col{figcnt}{crvls(hstcount)})
%         set(h1,'EdgeColor','none')
        for ii=1:2 
        crmn=mnout{ii}{crvls(hstcount)};
        crster=ster{ii}{crvls(hstcount)};
        plot([crmn-crster crmn+crster],[0.5 0.5],'Color',col{ii}{crvls(hstcount)},'Linewidth',3);
        plot(crmn,0.6,'v', 'Color',col{ii}{crvls(hstcount)},'MarkerFaceColor',col{ii}{crvls(hstcount)});
        plot([crmn-crster crmn+crster],[0.5 0.5],'Color',col{ii}{crvls(hstcount)},'Linewidth',3);
        plot(crmn,0.6,'v', 'Color',col{ii}{crvls(hstcount)},'MarkerFaceColor',col{ii}{crvls(hstcount)});
        end
        end
    
    
    end    
        
        
%         plot([0 crmn],[0.35 0.35],'Color',col.basrev,'Linewidth',3);
        plot([crmn-crster crmn+crster],[0.3 0.3],'Color',col{figcnt}{crvls(2)},'Linewidth',3);
        if(figcnt==1)
        crmn=mnout{crvls(2)};
        crster=ster{crvls(2)};
%          h1=arrow([0 0.35],[crmn 0.35],'Length',3,'Width',0.5)
%         set(h1,'FaceColor',col{figcnt}{crvls(2)})
%         set(h1,'EdgeColor','none')
        
%         arrowh([0 crmn],[0.35 0.35],col.asyrev,[200 200]);
% %         plot([0 crmn],[0.35 0.35],'Color',col.asyrev,'Linewidth',3);
%         plot([crmn-crster crmn+crster],[0.3 0.3],'Color',col{figcnt}{crvls(2)},'Linewidth',3);
%         
        text(3, 0.3, ['n = ' num2str(ct{crvls(2)})],'Color',col{figcnt}{crvls(2)},'Fontsize',16)  
        end
       

    
    




































































































