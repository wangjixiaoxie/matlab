%rewritten 10.7.08

% plot_timecourse.m
%assumes proper avl plot is loaded
%to plot against real time, do it like

%plotspecs
timepts=1:11
lw=3
msize=9
errlw=2
mucol=[0.4 0.4 1]
mufillcol=[0.82 0.82 0.82]
acfillcol=[0.92 0.96 0.98]
accol='k'
meanlw=2
errht=50
acmark='o'
mumark='o'
clear acmean mumean acstderr mustderr
notevl=1;
subval=floor(avls.rawtimes(1,1));
indvl=find(avls.mulist~=0)
 for ii=1:length(timepts)
     ind=timepts(ii)
        indvl2=avls.mulist(indvl(ind));
        acmean(ind)=avls.acmean(notevl,indvl2);
        acstderr(ind)=avls.acstderr(notevl,indvl2)
        mumean(ind)=avls.mumean(notevl,indvl2);
        mustderr(ind)=avls.mustderr(notevl,indvl2)
 end
    rwtms=avls.rawtimes(avls.mulist(indvl(timepts)),1)
    figure
    
     %now try to fill in, first the direct motor pathway, blue and then the
    %motor pathway, black.
    xpts=rwtms-subval
    xvec=[1:12 12:-1:1]
    mupts=[0 0 -1 -1.5 -1.9 -2  -2 -1.8 -1 -0.2 0 0 ]
    acpts=[0 0 -1.5 -2 -2 -2  -2 -1 -.2 0 0 0]
    yvec=[avls.initmean{1}*ones(length(mupts),1);mupts(end:-1:1)']
    yvec2=[mupts acpts(end:-1:1)]
hold on;
%     fill(xvec,yvec,acfillcol);
    fill(xvec,yvec2,mufillcol);
    
 
    
    
    plot(xvec(1:12),acpts,'Color',accol,'Linewidth',3);
    hold on;
     plot(xvec(1:12),acpts,'Marker',acmark,'Color',accol,'Markersize',msize,'MarkerFaceColor',accol); 
    plot(xvec(1:12),mupts,'Color',mucol,'Linewidth',3);
    plot(xvec(1:12),mupts,'Marker',mumark,'Color',mucol,'Markersize',msize,'MarkerFaceColor',mucol);
    
%     %plot vertical lines as error bars
%     plot([xvec(1,:);xvec(1,:)],[(acmean-acstderr); (acmean+acstderr)],'k','Linewidth',errlw);
%     plot([(rwtms-subval)' ;(rwtms-subval)'],[(mumean-mustderr); (mumean+mustderr)],'Color',mucol,'Linewidth', errlw);
    
%     hold on;
%     plot([0 30],[avls.initmean{1},avls.initmean{1}],'k--','Linewidth',2)
  
   
    
    title('/pu34/figs/timecourse.m, created with plot_timecourse.m')