
%plotstimfig1.

function [outstr,combfvst,combstcr,combmatout,avgmod]=plotstimfigbkw(sps,sts,outstr,combstcr,avgmod)
% figstoplot=[ 1 2 3 4 5 6]
figstoplot=[ 5]
%settings for first figure (time course)e
if(~exist('outstr'))
    [outstr,combfvst,combstcr]= getresiduals7(sps,sts,1:length(sps));
    [residmod,avgmod,combmatout]=alignresids(outstr,[0 .09],1:8)    
end
  avlspath='/oriole/bk48w74/datasum'
  avlsmat='datsum.mat'
    
  cmd=['cd ' avlspath];
  eval(cmd);
  cmd=['load ' avlsmat];
  eval(cmd);

%just showing details of stim pulse
if(ismember(1,figstoplot))  
    stim_path='/oriole/bk48w74/'
    stim_mat='stimexample2.mat'
    cmd=['cd ' stim_path];
    eval(cmd);
    cmd=['load ' stim_mat]
    eval(cmd);
    subplot(541)
    plot(tms,datnew,'k','Linewidth',2)
    axis square
    axis([-.001 .011 -.1e4 2e4])
end

%histogram of stim times
if(ismember(2,figstoplot)) 
   stimcomb=combstcr*1000;
   pitchmeasure=96
   subplot(5,4,2)
    mkhist(stimcomb,pitchmeasure);
    axis([0 105 0 60])
    axis square 
    box off
end

if(ismember(3,figstoplot))
    clear ps
    osind=9;
    binind=2;
    outstrinds=[3 4 5 6  ]
    subplotinds=[5 9 13 17 ]
    stimdiff=-outstr(osind).tm{binind}(outstrinds)
    stimwidth=outstr(osind).tm{binind}(2)-outstr(osind).tm{binind}(1);
    for ii=1:length(outstrinds)
        cros=outstrinds(ii);
        ctcontours=outstr(osind).combctcon{binind}{cros};
        fbcontours=outstr(osind).combfbcon{binind}{cros};
        axcon(ii)=subplot(5,4,subplotinds(ii))
        ps.ax=gca();
        ps.basind=1
        ps.TIMEBNDS=[.07 .11]
        ps.OFFMEAN=[.079]
        ps.fillbnds=[.075 .083]
        
        ps.plotstim=1;
        ps.stimvals=datnew;  
        ps.stimdiff=stimdiff(ii);
        ps.STIMGAIN=.000001;
        ps.STIMOFFSET=2.3
        ps.stimwidth=stimwidth;
        plotcontour2(ctcontours,fbcontours,avls,ps)
%         axis square
    end
    linkaxes(axcon);
    axis([0.01 0.12 2.25 2.4])
    
end

%plot stimpriortime vs. offset
if(ismember(4,figstoplot))
    subplot(5,4,[6 7 10 11])
    binind=2;
    ps.frac=1;
    ps.plotinds=1:8;
    plotoffsets(outstr(9),binind,ps)
    axis square;
    box off;
     axis([-.09 .01 -.3 .4])
         plot([-.01 0.07],[0 0],'k--','Linewidth',2);
end


%plot stimpriortime vs. offset
if(ismember(10,figstoplot))
    offsetind=[2  5 6 7 9 12]
    tm_matrix=[0.01:.001:0.06]
    
    ps.frac=1;
           
    [avgoffsetstruct]=avgoffsets(outstr,offsetind,tm_matrix,ps);
    [outvec,outvecer]=calcmeanstder2(avgoffsetstruct);
    outvec=-outvec;
    subplot(5,4,[6 10])
    colvec={'k' 'r' 'c' }
    for ii=5
        crind=offsetind(ii);
        for jj=1:length(outstr(crind).mndiff)
            binind=2;
            ps.frac=1;
            ps.colvec=colvec{jj};
            plotoffsets(outstr(crind),jj,ps)
            axis square
            axis([-.01 .07 -.2 .5])
            plot([-.01 0.07],[0 0],'k--','Linewidth',2);
        end
    end
    axtst=subplot(5,4,[7 11])
    
    fillx1=tm_matrix;
    fillx2=tm_matrix(end:-1:1);
    filly1=outvec+outvecer;
    filly2=outvec-outvecer;
    yvls=[filly1 filly2(end:-1:1)]
    xvls=[fillx1 fillx2]
    hold on;
    
%     fill(xvls,yvls,[0.6 0.6 0.6],'edgecolor','none');
%     hold on;
    plot(tm_matrix,avgoffsetstruct,'Linewidth',2,'Color','k');
    
    
end



if(ismember(5,figstoplot))
    
    clear ps
    avgmodind=[1:8 ]
    avgmodindinit=1;
%     subplotinds=[5 8 11 14 ]
       
    for ii=1:length(avgmodindinit)  
            
            subplot(5,4,[15 16 19 20])
            ps.avgmodind=avgmodindinit(ii);
            
            ax(ii)=gca();
            ps.ax=ax(ii);
%           ps.basind=1
            ps.TIMEBNDS=[0 .09]
            ps.frac=0;
            ps.plotsigdiff=1;
            plotresid(avgmod,ps)
            plot([0 0.09],[0 0],'k--','Linewidth',2);
            axis square;
    end
     axis([-.01 0.09 -10 60] )

figure
for ii=1:length(avgmodind)  
            
    subplot(2,4,ii);
            ps.avgmodind=avgmodind(ii);
            
            ax(ii)=gca();
            ps.ax=ax(ii);
%           ps.basind=1
            ps.TIMEBNDS=[0 .09]
            ps.frac=0;
            ps.plotsigdiff=1;
            plotresid(avgmod,ps)
            plot([0 0.09],[0 0],'k--','Linewidth',2);
            axis square;
       end
   
    linkaxes(ax);
  
        axis([-.01 0.09 -10 60] )
    
    
    

title('plotstimfig2.m')
end

if(ismember(6,figstoplot))
    
    clear ps
   
        subplot(5,4,[15 16 19 20])
        ps.ax=gca();
        ps.frac=1;
        ps.plotall=1;
%         ps.basind=1
%         ps.TIMEBNDS=[0 .07]
        plotresidgen(combmatout,ps)
        plot([0 0.08],[0 0],'k--','Linewidth',2);
        axis square;
   
    linkaxes(ps.ax);
    axis([-.01 0.09 -.5 .2])
    

title('plotstimfig2.m')
end



     
 %rewritten to subtract fb from matched catch mean.  
 
 
    function []=mkhist(combstcr,pitchmeasure)
        delayvls=pitchmeasure-combstcr
        edges=min(1.05*delayvls):4:max(1.05*delayvls);
        hstout=histc(delayvls,edges);
        stairs(edges,hstout,'k','Linewidth',2)
 
 
function []=plotresid(avgmod,ps)
crmod=avgmod{ps.avgmodind};
ind=find(crmod.tms>=ps.TIMEBNDS(1)&crmod.tms<=ps.TIMEBNDS(2)); 
tms=crmod.tms(ind);
mnvls=crmod.mnvls(ind);
stevls=crmod.stevls(ind);

if(isfield(ps,'frac'))
    if(ps.frac==1)
        mnvls=crmod.mnfrac(ind)
        stevls=crmod.stefrac(ind);
    end
end

     fillx=[tms tms(end:-1:1)]
     yvls=[mnvls+stevls]
        yvls2=[mnvls-stevls]
        filly=[yvls yvls2(end:-1:1)]
        
        acfillcol=[0.6 0.6 0.6]
         fill(fillx,filly,acfillcol,'edgecolor','w');
         hold on;
          plot(tms,mnvls,'k','Linewidth',2);
         if(ps.plotsigdiff)
             if(isfield(crmod,'initsigtime'))
                init_tm=crmod.initsigtime;
                if(crmod.consecsigind(end)<=length(tms));
                    fin_tm=tms(crmod.consecsigind(end));
                else
                    fin_tm=tms(end);
                end
%                 plot([init_tm fin_tm],[50 50],'c','Linewidth',2);
%                 plot([init_tm init_tm],[0 50],'r--');
                strout=[num2str(floor(init_tm*1000)) ' ms'];
             text(0,40,strout);
             end
         end
             
             
             
        
   
    
   function []=plotresidgen(plotstruct,ps)
frac=0;
mnvls=plotstruct.mn;
stevls=plotstruct.ste;

mnvlsall=plotstruct.mnall;
if(isfield(ps,'frac'))
    if(ps.frac==1)
        mnvls=plotstruct.frac
        mnvlsall=plotstruct.fracall;
        stevls=plotstruct.stefrac;
        frac=1;
    end
end

tms=plotstruct.tms;
 acfillcol=[0.6 0.6 0.6]
 rawcol=[0.4 0.4 0.4];
tmsind=1:length(mnvls);
tms=tms(tmsind);
     fillx=[tms tms(end:-1:1)]
     
     yvls=[mnvls+stevls]
        yvls2=[mnvls-stevls]
        filly=[yvls yvls2(end:-1:1)]
        
       
         fill(fillx,filly,acfillcol,'edgecolor','w');
         hold on;

         plot(tms,mnvls,'k','Linewidth',3);
         
        hold on;
        if(ps.plotall)
          for ii=1:length(mnvlsall(:,1))
            plot(tms,mnvlsall(ii,:),'Color',rawcol,'Linewidth',1);
          end
        end
function []=plotoffsets(outstr,binind,ps)
    diffvls=[];
    ervls=[]
    tmvls=[];
    difflist=[];
    crdiff=outstr.mndiff{binind};
    if(ps.frac)
        crdiff=outstr.mnfrac{binind};
    end
    if(isfield(ps,'plotinds'))
        plotinds=ps.plotinds
    else
        plotinds=1:length(crdiff)
    end
        
    
        for ii=plotinds
        if(~isempty(crdiff{ii}))
            difflist=[difflist ii];
        end
    end
%     crdiffind=find(~isempty(crdiff));
    for ii=difflist
        diffvls=[diffvls -crdiff{ii}]
        if(~ps.frac)
            ervls=[ervls outstr.mner{binind}{ii}]
        else
            ervls=[ervls outstr.mnerfrac{binind}{ii}]
        end
        tmvls=[outstr.tm{binind}(plotinds)]
    end
    plot(tmvls,diffvls,'Color','k','Linewidth',2);
    hold on;
    plot([-.08 0],[0 0],'k--');
    
    for ii=difflist
        plot([tmvls(ii) tmvls(ii)],[diffvls(ii)-ervls(ii) diffvls(ii)+ervls(ii)],'k','Linewidth',2)
    end

 function [outvec]=avgoffsets(outstr,osind,tm_matrix,ps)
     %average individual traces into tm_matrix
     %initialize tm_matrix
     ln=length(tm_matrix);
     for ii=1:length(osind)
         crind=osind(ii);
         cros=outstr(crind);
         rowlen=length(cros.mndiff);
         outvls=zeros(rowlen,ln);
         
         
         for jj=1:length(cros.mndiff)
             tmvls=[];
             diffvls=[];
             crdiff=cros.mndiff{jj};
            if(ps.frac)
                crdiff=cros.mnfrac{jj};
            end
            for tmind=1:length(cros.tm{jj})
                diffvls=[diffvls crdiff{tmind}]
                tmvls=[tmvls cros.tm{jj}{tmind}]
            end
            if(length(tmvls)>1)
                [outvls(jj,:)]=interp1(tmvls,diffvls,tm_matrix);
            end
         end

    outvec(ii,:)=calcmeanstder2(outvls);
     end
         
     
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




































































































