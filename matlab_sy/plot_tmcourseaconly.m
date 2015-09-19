%rewritten 10.7.08

% plot_timecourse.m
%assumes proper avl plot is loaded
%to plot against real time, do it like
function [ax]=plot_tmcourse_aconly(avls, nt,indin,syn_flag,ax)
if(exist('ax'))
    axes(ax(1))
else
    figure
end

%plotspecs
if(exist('indin'))
    muind=find(avls.mulist(indin)>0);
   
    
else
    muind=find(avls.mulist>0)
end

timepts=1:length(muind)
lw=2
msize=9
errlw=2
mucol=[0.4 0.4 1]
mufillcol=[0.82 0.82 0.82]
acfillcol=[0.92 0.96 0.98]
accol='k'
col{1}=mucol;
col{2}=accol;
meanlw=2
errht=50
acmark='o'
mumark='o'
clear acmean mumean acstderr mustderr
notevl=nt;
subval=floor(avls.rawtimes(1,1));
indvl=indin;
 for ii=1:length(indvl)
     ind=indvl(ii)
        if(~isempty(avls.mulist(ind)))
%             indvl2=avls.mulist(ind);
            acmean(ii)=avls.acmean(notevl,ind);
            acstdv(ii)=avls.acstdv(notevl,ind);
            acstderr(ii)=avls.acstderr(notevl,ind)
%             mumean(ii)=avls.mumean(notevl,indvl2);
%             mustderr(ii)=avls.mustderr(notevl,indvl2)
        end
 end
    rwtms=floor(avls.rawtimes(indvl,1))
%     figure
    
     %now try to fill in, first the direct motor pathway, blue and then the
    %motor pathway, black.
    xpts=rwtms-subval
    xvec=[(rwtms-subval);xpts(end:-1:1)]
%     mupts=mumean(1:length(indvl))
    acpts=acmean(1:length(indvl))
%     yvec=[avls.initmean{notevl}*ones(length(mupts),1);mupts(end:-1:1)']
%     yvec2=[mupts';acpts(end:-1:1)']
    hold on;
%     fill(xvec,yvec,acfillcol);
%     fill(xvec,yvec2,mufillcol);
    
%     plot(rwtms-subval,acmean,'Color',accol,'Linewidth',3);
    hold on;
     plot(rwtms-subval,acmean,'Marker',acmark,'Color',accol,'Markersize',msize,'MarkerFaceColor',accol); 
%     plot(rwtms-subval,mumean,'Color',mucol,'Linewidth',3);
%     plot(rwtms-subval,mumean,'Ma2rker',mumark,'Color',mucol,'Markersize',msize,'MarkerFaceColor',mucol);
    
    %plot vertical lines as error bars
    plot([(rwtms-subval)'; (rwtms-subval)'],[(acmean-acstdv); (acmean+acstdv)],'k','Linewidth',errlw);

    
    if(syn_flag)
    axes(ax(2));
        for ii=1:2
            plot([(rwtms-subval)'],[avls.ntrat(:,ii)'],'Color', col{ii},'Linestyle','-')
            hold on;
            plot([(rwtms-subval)'; (rwtms-subval)'],[avls.ntrat(:,ii)';avls.ntrat(:,ii)'],'o','Color',col{ii})
        end
    end
%      
%     
%     
%     
%     plot([(rwtms-subval)' ;(rwtms-subval)'],[(mumean-mustderr); (mumean+mustderr)],'Color',mucol,'Linewidth', errlw);
    
    hold on;
    plot([-15 30],[3456 3456],'k--','Linewidth',2)
    hold on;
%     plot
    
    title('/pu34/figs/timecourse.m, created with plot_timecourse.m')