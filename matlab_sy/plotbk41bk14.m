%rewritten 10.7.08

% plot_timecourse.m
%assumes proper avl plot is loaded
%to plot against real time, do it like
function [ax]=plot_tmcourse_aconly(avls, nts,ax,indin,syn_flag)
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
    indin=muind;
end

xbnds{1}={'2009-5-05 7' '2009-5-10 7'}
xbnds{2}={'2009-5-19 7' '2009-5-27 7'}
xbnds{3}={'2009-5-23 7' '2009-5-27 7'}
xbnds{4}={'2009-6-08 7' '2009-6-14 7'}


for ii=1:length(xbnds)
    xline{ii}=datenum(xbnds{ii},'yyyy-mm-dd HH');
end

timepts=1:length(muind)
lw=2
msize=4
errlw=1
mucol=[0.4 0.4 1]
mufillcol=[0.82 0.82 0.82]
acfillcol=[0.92 0.96 0.98]
accol='k'
col{1}='r';
col{2}=accol;
meanlw=2
errht=50
acmark='o'
mumark='o'
clear acmean mumean acstderr mustderr

for ntind=1:length(nts);
    crnt=nts(ntind);
    axes(ax(ntind));
    
   subval=datenum([2009,5,28])
    indvl=indin;
    for ii=1:length(indvl)
        ind=indvl(ii)
            if(~isempty(avls.mulist(ind)))
%             indvl2=avls.mulist(ind);
                acmean(ind)=avls.acmean(crnt,ind);
                acstdv(ind)=avls.acstdv(crnt,ind);
                acstderr(ind)=avls.acstderr(crnt,ind)
%             mumean(ii)=avls.mumean(notevl,indvl2);
%             mustderr(ii)=avls.mustderr(notevl,indvl2)
            end
    end
    rwtms=floor(avls.rawtimes(:,1))
%     figure
    [out,sortind]=sort(rwtms(indvl));
     %now try to fill in, first the direct motor pathway, blue and then the
    %motor pathway, black.
    sortind=indvl(sortind);
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
     plot(rwtms(sortind)-subval,acmean(sortind)./1000,'Linestyle','none','Marker',acmark,'Color',accol,'Markersize',msize,'MarkerFaceColor',accol); 
%     plot(rwtms-subval,mumean,'Color',mucol,'Linewidth',3);
%     plot(rwtms-subval,mumean,'Ma2rker',mumark,'Color',mucol,'Markersize',msize,'MarkerFaceColor',mucol);
    
    %plot vertical lines as error bars
    plot([(rwtms(sortind)-subval)'; (rwtms(sortind)-subval)'],[((acmean(sortind)-acstdv(sortind))./1000); ((acmean(sortind)+acstdv(sortind)))./1000],'k','Linewidth',errlw);
    
    for ii=1:length(xline)
        plot([(xline{ii}(1)-subval) xline{ii}(2)-subval],[3500 3500],'k--')
    end
end
    if(syn_flag)
    axes(ax(3));
        for ii=1:2
            plot([(rwtms(sortind)-subval)'],[avls.ntrat(sortind,ii)'],'Color', col{ii},'Linestyle','-','Linewidth',2)
            hold on;
            plot([(rwtms(sortind)-subval)'; (rwtms(sortind)-subval)'],[avls.ntrat(sortind,ii)';avls.ntrat(sortind,ii)'],'o','Color',col{ii},'MarkerFaceColor',col{ii})
        end
        
%         lowind=find(avls.ntrans(sortind)<50);
%         for kk=1:length(lowind)
%             plot([(rwtms(lowind(kk))-subval) (rwtms(lowind(kk))-subval)],[.95 1],'r-')
%         end
        for ii=1:length(xline)
%         plot([(xline{ii}(1)-subval) xline{ii}(2)-subval],[0.95 0.95],'k--')
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
    axes(ax(1))
    axis([-15 35 3200 3800])
    axes(ax(2))
    axis([-15 35 4650 5200])
    
    axes(ax(3))
    axis([-15 75 0 1])
    title('/pu34/figs/timecourse.m, created with plot_timecourse.m')