 %plotdynamics script,
 
 %4.15.09, calls plotdynamics2

 
 %1 is asymptote example of pu34 (first panel)
 %2 is example of measured dynamics for pu34
 %3 
 %4 is slope calculation, plot.
 %5 is example of reversal.
 %6 is summary plot of muscimol offsets.
 
 figstoplot=[6]
 
 
 if(find(figstoplot==1))
     figure
     clear ps
     ps.ax=gca();
     ps.notevl=1;
%      ps.indin=find(sumbs(2).mulist~=0)
ps.indin=[1:7]
     ps.plotasymp=1;
     ps.addx=1;
     %plot_tmcourse3, 
     plot_tmcourse3(sumbs(9), ps)
     clear ps
 end
 
 
%takes sumdyn,
%and a ps (plotstruct).
%ps.norm - whether to normalize or not
%ps.minx - minimum x value to include
%ps.maxx - maximum x value to include.
%ps.colors,  col={ 'k' 'r','c'}
%ps.colvec
%ps.type, (pct or off)
%ps.addx

if(find(figstoplot==2))
    ps.normeff=1
    ps.comb=0;
    ps.norm=0;
    ps.minx=6;
    ps.maxx=6;
    ps.col={'k' 'r' 'c'}
    ps.colvec=[ones(length(sumdyn),1)]
    ps.colvec(3)=2;
    ps.colvec(end)=3;
    ps.colvec(end-1)=3;
    ps.type='off'
    ps.plotavg=0
    ps.eff='both'
    ps.excludeind=11;
    [sumdynmean]=plotdynamics3(sumdyn(3),ps);
end

%this is for slope calculation
if(find(figstoplot==4))
    ps.normeff=1
    ps.comb=1;
    ps.norm=0;
    ps.minx=3;
    ps.maxx=6;
    ps.col={'k' 'r' 'c'}
    ps.colvec=[ones(length(sumdyn),1)]
    ps.colvec(3)=2;
    ps.colvec(end)=3;
    ps.colvec(end-1)=3;
    ps.type='dis'
    ps.plotavg=1
    ps.eff='both'
    ps.excludeind=11;
    [sumdynmean]=plotdynamics3(sumdyn,ps);
    figure;
    plot(sumdynmean.effslope,-sumdynmean.yslope,'o','MarkerSize',12,'Color','k')
    hold on;
    plot([-3 3],[0 0],'k--')
    plot([0 0],[-3 3],'k--')
end


if(find(figstoplot==5))
     figure
     clear ps
     ps.ax=gca();
     ps.notevl=1;
%      ps.indin=find(sumbs(2).mulist~=0)
ps.indin=[5 7:15]
     ps.plotasymp=1;
     ps.addx=1;
     %plot_tmcourse3, 
     plot_tmcourse3(sumbs(2), ps)
     clear ps
end
 
if(find(figstoplot==6))
     figure
    plotbarinitind5(sumplotrevout, 'bar',1,1)
    
 end