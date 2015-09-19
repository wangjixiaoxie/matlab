
%plotstimfig1.

function [combvls]=plotstimfig1(sumbs,sts)
% figstoplot=[1 2 3 4 5 7]
 figstoplot=[1 2 3  5 6 7]
%settings for first figure (time course)


  avlspath='/oriole/bk48w74/datasum'
  avlsmat='datsum.mat'
    
  cmd=['cd ' avlspath];
  eval(cmd);
  cmd=['load ' avlsmat];
  eval(cmd);

ps(1).bsnum=[2 1]
ps(1).runstoplot{1}=[1 4 9 10 11 12 13 14 15 16 17]
ps(1).runstoplot{2}=[2 3 4 5 6 7 13  14 15 16 17]

ps(2).bsnum=[2 ]
ps(2).runstoplot{1}=[4 13]
% ps(2).runstoplot{2}=[3 14]
ps(2).totalplot=2;
ps(2).shiftdrxn=[1 0]
ps(2).arrowht=[0.4 0.45]
ps(2).distht=0.35
ps(2).shiftnum=1;


if(ismember(1,figstoplot))
  
    
    exsong(1).path=[avls.baspath avls.pvls{2}]
    exsong(1).cvl=avls.cvl{2}
    exsong(1).fn='bk48w74_210809_105358.168.cbin';
    exsong(1).bnds=[4 6.7]
    exsong(1).stim_path='/oriole/bk48w74/stimexample'
    exsong(1).stim_mat='stimexample.mat'
    
    
%     exsong(1).fvst=fvstex{1};

    exsong(2).avlspath='oriole/bk48w74/datasum'
    exsong(2).avlsmat='datsum.mat'
    exsong(2).path=[avls.baspath avls.pvls{4}]
    exsong(2).cvl=avls.cvl{4}
    exsong(2).fn='bk48w74_230809_180859.15.cbin';
    exsong(2).bnds=[7.3 10]
%     exsong(2).fvst=fvstex{2}

    [axspects1,axspects2]=plotexspect(exsong);
end

if(ismember(2,figstoplot))
    ptsps.marksize=3
    ptsps.rawvl='raw'
    ptsps.ntvl=1;
    indtoplot=[1 4 68]
    ptsps.STIM=1
    ptsps.plotextra=0
    sfact=3000;
    figure
    for ii=1:length(indtoplot)
        ptsps.indtoplot=indtoplot(ii);
        axraw(ii)=subplot(length(indtoplot),4,(ii-1)*4+2)
        ptsps.ax=gca();
        inactiv_rawpoints(avls,ptsps,sfact);
        axis square
        box off
    end
    linkaxes(axraw)
    axis([-.1 .25 2.1 2.7])
    
end

if(ismember(3,figstoplot))
    clear ps
    indtoplot=[1 4 68]
    ps.basind=1;
    ps.TIMEBNDS=[.07 .11]
    for ii=1:length(indtoplot)
        ps.indtoplot=indtoplot(ii);
        axcon(ii)=subplot(length(indtoplot),4,(ii-1)*4+1)
        ps.ax=gca();
        ps.boxbnds=[.079 .095]
        plotcontour(avls,ps)
        axis square
        box off
    end
    linkaxes(axcon)
    axis([.06 .12 2.2 2.55])
end

if(ismember(4,figstoplot))
    clear ps
    indlist{1}=[2 3 4 5 ]
    indlist{2}=[67 68 70 75 ]
    bsvl=1
    for ii=1:length(indlist)
         subplot(3,4,(ii-1)*4+4)
         ps.ax=gca();
        plotextmcourse(sumbs(1),indlist{ii},ps)
    end
end

if (ismember(5,figstoplot))
     clear ps
     indtoplot=[1 4 68]
     bsnum=1;
     ps.arrowht=[0.4 0.45]
    ps.distht=0.35

     for ii=1:length(indtoplot)
         ps.indtoplot=indtoplot(ii);
         ps.initmean=sumbs(bsnum).initmean
           subplot(length(indtoplot),4,(ii-1)*4+3)
           axhist(ii)=gca();
           ps.ax=gca();
           ps.shift_col='k'
           ps.rev_col='r'
       plotexhist(avls,ps )
     axis square;
     end
     
     linkaxes(axhist);
     axis([2.1 2.7 0 0.5])
end



if (ismember(6,figstoplot))
    subplot(3,4,12);
    axscatter=gca();
    axin=gca();
    axbnds=[0 10 0 10]
    colin{1}{1}=[0.1 0.1 0.1]
    colin{1}{2}=[0.5 0.5 0.5]
    
    colin{2}{1}=[1 0 0]
    colin{2}{2}=[1 .6 .6]
    colin{3}=[.6 .6 1]
    [shiftplot,combvls]=plothistoryscatter1(sumbs,0,0,0,0,axin,axbnds,axscatter,1,colin) 
    axis([-2 4 -2 4])
    axis square;   
end
if (ismember(7,figstoplot))
     
      clear ps
     col{1}{1}=[0.1 0.1 0.1]
    col{1}{2}=[0.5 0.5 0.5]
    col{2}{1}=[1 0 0]
    col{2}{2}=[1 0.7 .7]


      ps.ax(1)=subplot(346)
      ps.ax(2)=subplot(3,4,9);
      ps.col=col;
      ps.edges=[-3.2:.4:3.2]
     
       ps.axbnds=[-4 4 0 0.5]
      plotgrouphsts(combvls,ps);
      box off;
    
end

% 
% if (ismember(7,figstoplot))
%      clear ps
%      figure
%      [sumdyn]=selectdynamicrunsstim(sumbs,0,0,0)
%      
%      ps.minx=4
%     ps.maxx=4
%     ps.ac_col='k'
%     ps.mu_col=[0.4 0.4 1]
% 
%     ps.flip=1
%     ps.plotraw=0
%     ps.addzero=1
%     ps.plotsep=0;
%     ps.type='nor';
%     ps.insert=1
%     ps.aligntimes=1
%     ps.usepct=1
%     [outvlaczcomb,outvlmuzcomb]=plotinitstimdynamics(sumdyn,sumbs,ps)
%     axis([-1 5 -.5 4])
%     plot([-1 5],[ 0 0],'k--','Linewidth',2);
% end


     
     
   


function []=plotextmcourse(sumbs,indlist,ps)
    
%      figure
   
        
        ps.STIM=1;
        ps.plotarrow=1;
%         axtmcourse(ii)=subplot(1,bsln,ii);
        plot_tmcourse2rev(sumbs,indlist,ps);
   








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
        ps.plotshiftarrow=1;
        ps.plotrevarrow=1;
   plotarrows(origarrow,shiftarrow,ps);
   
   
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




































































































