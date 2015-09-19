
%plotstimfig1.

function [plotvls]=plotap5fig(sumbs,bs)
figstoplot=[4]
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

if (ismember(4,figstoplot))
     clear ps
    ps.marksize=2.5;
    ps.rawvl='raw'
    ps.ntvl=1;
    ps.STIM=0;
    plotnextday=[0 1 0 1];
    nextind=[0 25 0 20];
    ps.plotextra=1;
    omitpost=[1 1 1 1];
    indtoplot=[2 5 2  4]
    sumbsvl=[10 10 11 11];
     ps.col{1}='k'
     ps.col{2}=['k']
     ps.col{3}=['r']
      
        ps.plot_triangles=1;
        ps.triangle_xvl=9.5;
    
    ps.distht=19.5
    
    ps.ntind=1;
    
    sfact=3000;
    axinds=[1 2 3 4]
     for ii=1:length(indtoplot)
         avls=getavls(bs,sumbsvl(ii))
         ps.indtoplot=indtoplot(ii);
        subplot(2,4,axinds(ii))
         ps.initmean=sumbs(sumbsvl(ii)).initmean
         ps.PLOTNEXTDAY=plotnextday(ii);
         ps.nextdayind=nextind(ii)
        
         ps.plotshiftarrow=[1 1 1 1 ]
         ps.plotrevarrow=[1 1 1 1 ]
         axhist(ii)=gca();
         ps.ax=gca();
         ps.initmean=sumbs(sumbsvl(ii)).initmean
         ps.plotpre=0;
         ps.plot_triangles=1;
         ps.omitpost=1;
         ps.SERIAL=0;
        sumdataex(ii)=inactiv_rawpoints(avls,ps,sfact )
        box off;
     end
     linkaxes(axhist(1:2));
     axes(axhist(1));
      axis ([-2 28 2.2 2.8])
     linkaxes(axhist(3:4));
     axes(axhist(3));
      axis ([-2 28 1.9 2.65])
   
%   
%    linkaxes(axhist([2 4]));
%    
%     axis ([-2 14 2 2.8])
   

end



if (ismember(5,figstoplot))
     clear ps
    
     ps.col{1}='k'
     ps.col{2}='r'
     ps.col{3}=[0.6 0.6 0.6]
     ps.arrowht=[0.4 0.45]
    ps.distht=0.35
    ps.ntind=1;
    
    plotpre=[0 0 0 0]
    omitpost=[1 1  1 1]

     indtoplot=[2 5 2 4]
    sumbsvl=[10 10 11 11]
    axinds=[1 2 3 4]
     for ii=1:length(indtoplot)
         avls=getavls(bs,sumbsvl(ii))
         ps.indtoplot=indtoplot(ii);
        subplot(2,4,axinds(ii)+4)
         ps.initmean=sumbs(sumbsvl(ii)).initmean
         ps.plotshiftarrow=[1 1 1 1 ]
         ps.plotrevarrow=[1 1 1 1 ]
         axhist(ii)=gca();
         ps.ax=gca();
         ps.plotpre=0;
         ps.plot_triangles=1;
         ps.omitpost=1;
        plotexhistinactiv2(avls,ps )
       axis square
     end
     linkaxes(axhist([1 2]));
    axes(axhist(1));
    axis ([2.1 2.8 0 0.4])
  
   linkaxes(axhist([3 4]));
   axes(axhist(3));
    axis ([1.9 2.6 0 0.4])
   

end

if (ismember(6,figstoplot))
figure
clear ps
ps.axsum=gca();
ps.axpct=gca();
ps.axbas=gca();
%     clear ps
    %     subplot(3,2,[4 6]);
%     ps.axscatter=gca();
%     ps.axin=gca();
%     ps.axbnds=[0 10 0 10]
%     ps.colin=0;
    ps.PLOTSCATTER=1;
    ps.calcz=0;
    ps.USEPRE=1;
    ps.minpts=0;
    ps.minshift=0;
    ps.DIFF=0;
    ps.combduplicates=1;
    ps.norm_asymp=0;
    ps.plotcv=0;
    ps.STIM=0;
    ps.NTANAL=[1 2];
   ps.MX_TM=4;
   plotin(1).type=1;
         plotin(1).nt=1;
         plotin(1).xvl=[1 2 3]
         plotin(1).comb=1;
         plotin(1).norm=1;
         plotin(1).drxn=1;
         plotin(1).CTRL=0;
          plotin(2).type=2;
         plotin(2).nt=1;
         plotin(2).xvl=[1 2 3]
         plotin(2).comb=1;
         plotin(2).norm=1;
         plotin(2).drxn=1;
         plotin(2).CTRL=0;
         
         
         
   ps.TYPE='plotsum'
   
    [shiftplot,combvls]=plothistoryscatterB(sumbs(10:11),ps) 
    [plotvls,stats]=plotgroupbarinit(combvls,ps,plotin);
    sumdata.MX_TM=ps.MX_TM;
%     axis([-1 4 -1 4])
%     axis square;
    
    
end



function avls=getavls(bs,sumbsvl)
    pth=bs(sumbsvl).path;
    mtname=bs(sumbsvl).matfilename;
%     cmd=['cd ' pth 'datasum'];
cmd=['cd ' pth ];
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




































































































