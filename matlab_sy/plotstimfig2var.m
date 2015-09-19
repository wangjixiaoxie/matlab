
%plotstimfig1.

function [outstr,combfvst,combstcr]=plotstimfig2(sumbs,sts,sps,outstr,combstcr)
figstoplot=[ 5]
%settings for first figure (time course)
if(~exist('outstr'))
    [outstr,combfvst,combstcr]= getresiduals4(sps,sts,12);
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
    subplot(531)
    plot(tms,datnew,'k','Linewidth',2)
    axis square
    axis([-.001 .011 -.1e4 2e4])
end

%histogram of stim times
if(ismember(2,figstoplot)) 
   stimcomb=combstcr*1000;
   pitchmeasure=96
   subplot(5,3,2)
    mkhist(stimcomb,pitchmeasure);
    axis([0 105 0 60])
    axis square 
    box off
end

if(ismember(3,figstoplot))
    clear ps
    outstrinds=[3 4 5 6 ]
    subplotinds=[4 7 10 13 16 ]
    stimdiff=[.045 .035 .025 .015]
    for ii=1:length(outstrinds)
        ctcontours=outstr.combctcon{outstrinds(ii)};
        fbcontours=outstr.combfbcon{outstrinds(ii)};
        axcon(ii)=subplot(5,3,subplotinds(ii))
        ps.ax=gca();
        ps.basind=1
        ps.TIMEBNDS=[.07 .11]
        
        ps.plotstim=1;
        ps.stimvals=datnew;
        ps.stimdiff=stimdiff(ii);
        ps.STIMGAIN=.000001;
        ps.STIMOFFSET=2.3
        plotcontour2(ctcontours,fbcontours,avls,ps)
        axis square
    end
    linkaxes(axcon);
    axis([0.02 0.12 2.25 2.65])
    
end

%plot stimpriortime vs. offset
if(ismember(4,figstoplot))
    subplot(5,3,6)
    plotoffsets(outstr)
    axis square
     axis([-.01 .07 -50 30])
         plot([-.01 0.07],[0 0],'k--','Linewidth',2);
end

if(ismember(5,figstoplot))
    
    clear ps
    outstrinds=[3 4 5 6]
    subplotinds=[5 8 11 14 17 ]
    for ii=1:length(outstrinds)
        outresidfb=outstr.outmatfb{outstrinds(ii)};
        out_tmsfb=outstr.fb_tms{outstrinds(ii)};
        
        outresidct=outstr.outmatct{outstrinds(ii)};
        out_tmsct=outstr.ct_tms{outstrinds(ii)};
        
        subplot(5,3,subplotinds(ii))
        ps.ax(ii)=gca();
        ps.basind=1
        ps.TIMEBNDS=[.07 .12]
        plotresid(out_tmsfb,outresidfb,out_tmsct,outresidct)
        plot([0 0.08],[0 0],'k--','Linewidth',2);
        axis square;
    end
    linkaxes(ps.ax);
    axis([0 0.09 -30 50])
    
    
    
end
title('plotstimfig2.m')

     
 %rewritten to subtract fb from matched catch mean.  
 
 
    function []=mkhist(combstcr,pitchmeasure)
        delayvls=pitchmeasure-combstcr
        edges=min(1.05*delayvls):4:max(1.05*delayvls);
        hstout=histc(delayvls,edges);
        stairs(edges,hstout,'k','Linewidth',2)
 
 
function []=plotresid(fbtms,residfb,ct_tms,residct)
%first align catch and fb
lnct=length(ct_tms);
lnfb=length(fbtms);
if(lnct<lnfb)
    ln=lnct
    tms=ct_tms
else
    ln=lnfb
    tms=fbtms
end

mnct=mean(residct,2);

for ii=1:ln
    
    indfb=find(residfb(ii,:)~=0);
    indct=find(residct(ii,:)~=0);
        lnindfb=length(indfb);
        lnindct=length(indct);
        if(lnindfb>10&lnindct>10)
            mnvls(ii)=mean(residfb(ii,indfb))-mean(residct(ii,indct));
            stevls(ii)=std(residfb(ii,indfb))./sqrt(lnindfb);
        else
            mnvls(ii)=0;
            stevls(ii)=0;
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
   
    
   


function []=plotoffsets(outstr)
    diffvls=[];
    ervls=[]
    tmvls=[];
    for ii=1:length(outstr.mndiff)
        diffvls=[diffvls -outstr.mndiff{ii}]
        ervls=[ervls outstr.mner{ii}]
        tmvls=[tmvls outstr.tm{ii}]
    end
    plot(tmvls,diffvls,'k','Linewidth',2);
    hold on;
    for ii=1:length(outstr.mndiff)
        plot([tmvls(ii) tmvls(ii)],[diffvls(ii)-ervls(ii) diffvls(ii)+ervls(ii)],'k','Linewidth',2)
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




































































































