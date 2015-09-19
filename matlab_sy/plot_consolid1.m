
%plotstimfig1.

function [combvls]=plot_consolid1(phsumbs,bs,sumbs)
figure

% figstoplot=[1 2 3 4 5  ]
% figstoplot=[5 6 7 ];
figstoplot=[7 8 10 11 12 13]
%settings for first figure (time course)

% 
%   avlspath='/oriole/bk48w74/datasum'
%   avlsmat='datsum.mat'
%     
%   cmd=['cd ' avlspath];
%   eval(cmd);
%   cmd=['load ' avlsmat];
%   eval(cmd);
% 
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

%Example asymptote plot.
% figure
%   [sumdyn,sumdynasymp,sumdynrev]=selectdynamicruns3(sumbs,0,0,0,0)
if(ismember(1,figstoplot))
    ax(1)=subplot(3,2,3)
    
    ps.minx=6
    ps.maxx=6
    ps.ac_col='k'
    ps.mu_col=[0.4 0.4 1]
    ps.flip=1
    ps.plotraw=0
    ps.addzero=0
    ps.plotsep=1;
    ps.type='nor';
    ps.insert=1
    ps.aligntimes=1
    ps.usepct=1
    ps.subdayinterp=0
    ps.plotarrow=1;
      ps.arrowind=[1 2 4 6]
    [outvlaczcomb,outvlmuzcomb,sdind]=plotasympdynamics(sumdynasymp(3),ps)
      axis([0 7 -.2 1.2])
      plot([0 7],[0 0],'k--','Linewidth',2)
end 
if(ismember(2,figstoplot))
    ax(2)=subplot(324);
   ps.minx=3
    ps.maxx=4
    ps.ac_col='k'
    ps.mu_col=[0.4 0.4 1]

    ps.flip=1
    ps.plotraw=0
    ps.addzero=0
    ps.plotsep=1;
    ps.type='nor';
    ps.insert=1
    ps.aligntimes=1
    ps.plotpct=1
    
    ps.subdayinterp=1

    [outvlaczcomb,outvlmuzcomb]=plotrevdynamics(sumdynrev(2),ps)
     plot([0 7],[0 0],'k--','Linewidth',2)
    axis([0 5 -.2 1.2])

end

if(ismember(3,figstoplot))
   ax(1)=subplot(3,2,5)
    
    ps.minx=6
    ps.maxx=6
    ps.ac_col='k'
    ps.mu_col=[0.4 0.4 1]
    ps.flip=1
    ps.plotraw=0
    ps.addzero=0
    ps.plotsep=1;
    ps.type='nor';
    ps.insert=1
    ps.aligntimes=1
    ps.usepct=1
    ps.subdayinterp=0
     ps.plotarrow=1;
      ps.arrowind=[1 2 3 4 5 6]
    [outvlaczcomb,outvlmuzcomb,sdind]=plotasympdynamics(sumdynasymp,ps)
    axis([0 7 -.2 1.2])
    plot([0 7],[0 0],'k--','Linewidth',2)
    
end

if(ismember(4,figstoplot))
     ax(2)=subplot(326);
   ps.minx=3
    ps.maxx=4
    ps.ac_col='k'
    ps.mu_col=[0.4 0.4 1]

    ps.flip=1
    ps.plotraw=0
    ps.addzero=0
    ps.plotsep=1;
    ps.type='nor';
    ps.insert=1
    ps.aligntimes=1
    ps.plotpct=1
    
    ps.subdayinterp=1
    ps.plotarrow=1;
    [outvlaczcomb,outvlmuzcomb]=plotrevdynamics(sumdynrev,ps)
     plot([0 7],[0 0],'k--','Linewidth',2)
    axis([0 5 -.2 1.2])
end

if (ismember(5,figstoplot))
    crbs=5;
%     ps.STIM=0;
        ps.plotarrow=1;
%         ps.runstoplot=[5 7 8 9 10 11 12 13 14 15]
        ps.runstoplot=[2 3 6 8 10 11 12  ]
        crax=subplot(3,3,1)
        ps.ax=crax;
        ps.adjx=1;
%         axtmcourse(ii)=subplot(1,bsln,ii);
        plot_tmcourse2rev(phsumbs(crbs),ps.runstoplot,ps);
        axis square
        axis([ -2 16 2.15 2.45])
        
end

if (ismember(6,figstoplot))
    crbs=1;
%     ps.STIM=0;
        ps.plotarrow=1;
%         ps.runstoplot=[5 7 8 9 10 11 12 13 14 15]
        ps.runstoplot=[ 1     2     3     4     5     6     7 13  ]
        crax=subplot(3,3,2)
        ps.ax=crax;
        ps.adjx=1;
        ps.STIM=1
%         axtmcourse(ii)=subplot(1,bsln,ii);
        plot_tmcourse2rev(sumbs(crbs),ps.runstoplot,ps);
        axis square
        axis([ -2 16 2.15 2.45])      
end



% 
% if (ismember(6,figstoplot))
%     crbs=2;
% %     ps.STIM=0;
%         ps.plotarrow=1;
%         ps.runstoplot=[5 7 8 9 10 11]
%         crax=subplot(322)
%         ps.ax=crax;
% %         axtmcourse(ii)=subplot(1,bsln,ii);
%         plot_tmcourse2rev(sumbs(crbs),ps.runstoplot,ps);
% end

if (ismember(7,figstoplot))
     clear ps
     ps.ax=subplot(337);
%      [sumdyn]=selectdynamicruns3(sumbs,0,0,0,0)
    
    ps.minx=1.4
ps.maxx=4
ps.ac_col='k'
ps.mu_col=['r']
ps.mufillcol=[1 0.7 0.7]
ps.acfillcol=[0.82 0.82 0.82]


ps.exmufillcol=[1 0.7 0.7]
ps.exacfillcol=[0.82 0.82 0.82]
ps.flip=1
ps.plotraw=1
ps.addzero=1
ps.plotsep=0;
ps.ploter=0;
ps.type='nor';
ps.insert=1
ps.aligntimes=1
ps.usepct=0
ps.plotbnds=[-1 5  -1 5]
ps.analbnds=[0 4]
% figure
[sumdyn,sumdynasymp,sumdynrev]=selectdynamicrunsstim2(sumbs,0,0,0)
[outvlaczcomb,outvlmuzcomb]=plotcombdynamics(sumdyn,sumbs,ps)
hold on;
plot([-1 5],[0 0],'k--');
plot([-1 5],[3.5 3.5],'k--');

end

if (ismember(8,figstoplot))
     clear ps
     ps.ax=subplot(338);
%      [sumdyn]=selectdynamicruns3(sumbs,0,0,0,0)
    
    ps.minx=5
ps.maxx=5
ps.ac_col='k'
ps.mu_col=['r']
ps.mufillcol=[1 0.7 0.7]
ps.acfillcol=[0.82 0.82 0.82]
ps.exmufillcol=ps.mufillcol
ps.exacfillcol=ps.acfillcol;

ps.ploter=0;
ps.flip=1
ps.plotraw=1
ps.addzero=1
ps.plotsep=0;
ps.type='nor';
ps.insert=1
ps.aligntimes=1
ps.usepct=1
ps.plotbnds=[-1 5 -1 5]
ps.analbnds=[0 6]
% figure
[sumdyn,sumdynasymp,sumdynrev]=selectdynamicrunsstim2(sumbs,0,0,0)
[outvlaczcomb,outvlmuzcomb]=plotcombdynamics(sumdynasymp,sumbs,ps)
hold on;
 plot([-1 5],[0 0],'k--');
plot([-1 5],[3.5 3.5],'k--');

end
if (ismember(9,figstoplot))
     clear ps
     ps.ax=subplot(339);
%      [sumdyn]=selectdynamicruns3(sumbs,0,0,0,0)
    
    ps.minx=1
ps.maxx=4
ps.ac_col='k'
ps.mu_col=['r']
ps.mufillcol=[1 0.7 0.7]
ps.acfillcol=[0.82 0.82 0.82]

ps.flip=1
ps.plotraw=0
ps.addzero=1
ps.plotsep=0;
ps.type='nor';
ps.insert=1
ps.aligntimes=1
ps.usepct=1
ps.plotbnds=[-1 5 -1 5]
ps.analbnds=[0 4]
% figure
[sumdyn,sumdynasymp,sumdynrev]=selectdynamicrunsstim2(sumbs,0,0,0)
[outvlaczcomb,outvlmuzcomb]=plotcombdynamics(sumdynrev,sumbs,ps)
hold on;
 plot([-1 5],[0 0],'k--');
plot([-1 5],[3.5 3.5],'k--');

end
if (ismember(10,figstoplot))
     clear ps
     ps.ax=subplot(335);
%      [sumdyn]=selectdynamicruns3(sumbs,0,0,0,0)
    
    ps.minx=5
ps.maxx=5
ps.ac_col='k'
ps.mu_col=['r']
ps.mufillcol=[1 0.7 0.7]
ps.acfillcol=[0.82 0.82 0.82]
ps.exmufillcol=ps.mufillcol
ps.exacfillcol=ps.acfillcol;

ps.flip=1
ps.plotraw=1
ps.ploter=0;
ps.addzero=1
ps.plotsep=0;
ps.type='nor';
ps.insert=1
ps.aligntimes=1
ps.usepct=1
ps.plotbnds=[-1 5 -1 5]
ps.analbnds=[0 6]
% figure
[phsumdyn,phsumdynasymp,phsumdynrev]=selectdynamicruns4(phsumbs,0,0,0)
[outvlaczcomb,outvlmuzcomb]=plotcombdynamics(phsumdynasymp,phsumbs,ps)
hold on;
 plot([-1 5],[0 0],'k--');
plot([-1 5],[3.5 3.5],'k--');

end

if (ismember(11,figstoplot))
     clear ps
     ps.ax=subplot(334);
%      [sumdyn]=selectdynamicruns3(sumbs,0,0,0,0)
    
    ps.minx=1
ps.maxx=4
ps.ac_col='k'
ps.mu_col=['r']
ps.mufillcol=[1 0.7 0.7]
ps.acfillcol=[0.82 0.82 0.82]
ps.exmufillcol=ps.mufillcol
ps.exacfillcol=ps.acfillcol;

ps.flip=1
ps.plotraw=1
ps.ploter=0;
ps.addzero=1
ps.plotsep=0;
ps.type='nor';
ps.insert=1
ps.aligntimes=1
ps.usepct=0
ps.plotbnds=[-1 5 -1 5]
ps.analbnds=[0 4]
% figure
[phsumdyn,phsumdynasymp,phsumdynrev]=selectdynamicruns4(phsumbs,0,0,0)
[outvlaczcomb,outvlmuzcomb]=plotcombdynamics(phsumdyn,phsumbs,ps)
hold on;
 plot([-1 5],[0 0],'k--');
plot([-1 5],[3.5 3.5],'k--');

end

if (ismember(12,figstoplot))
     clear ps
     ps.ax=subplot(336);
%      [sumdyn]=selectdynamicruns3(sumbs,0,0,0,0)
    
    ps.minx=6
ps.maxx=7
ps.ac_col='k'
ps.mu_col=['r']
ps.mufillcol=[1 0.7 0.7]
ps.acfillcol=[0.82 0.82 0.82]
ps.exmufillcol=ps.mufillcol
ps.exacfillcol=ps.acfillcol;

ps.flip=1
ps.plotraw=1
ps.ploter=0;
ps.addzero=1
ps.plotsep=0;
ps.type='nor';
ps.insert=1
ps.aligntimes=1
ps.usepct=1
ps.plotbnds=[-1 5 -1 5]
ps.analbnds=[0 6]
% figure
[phsumdyn,phsumdynasymp,phsumdynrev]=selectdynamicruns4(phsumbs,0,0,0)
[outvlaczcomb,outvlmuzcomb]=plotcombdynamics(phsumdynasymp,phsumbs,ps)
hold on;
 plot([-1 5],[0 0],'k--');
plot([-1 5],[3.5 3.5],'k--');

end

if (ismember(13,figstoplot))
     clear ps
     ps.ax=subplot(339);
%      [sumdyn]=selectdynamicruns3(sumbs,0,0,0,0)
    
    ps.minx=6
ps.maxx=7
ps.ac_col='k'
ps.mu_col=['r']
ps.mufillcol=[1 0.7 0.7]
ps.acfillcol=[0.82 0.82 0.82]
ps.exmufillcol=ps.mufillcol
ps.exacfillcol=ps.acfillcol;

ps.ploter=0;
ps.flip=1
ps.plotraw=1
ps.addzero=1
ps.plotsep=0;
ps.type='nor';
ps.insert=1
ps.aligntimes=1
ps.usepct=1
ps.plotbnds=[-1 5 -1 5]
ps.analbnds=[0 6]
% figure
[sumdyn,sumdynasymp,sumdynrev]=selectdynamicrunsstim2(sumbs,0,0,0)
[outvlaczcomb,outvlmuzcomb]=plotcombdynamics(sumdynasymp,sumbs,ps)
hold on;
 plot([-1 5],[0 0],'k--');
plot([-1 5],[3.5 3.5],'k--');

end



function avls=getavls(bs,sumbsvl)
    pth=bs(sumbsvl).path;
    mtname=bs(sumbsvl).matfilename;
    cmd=['cd ' pth 'datasum'];
    eval(cmd);
    cmd=['load ' mtname];
    eval(cmd);
    test=1;
    
     
   


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




































































































