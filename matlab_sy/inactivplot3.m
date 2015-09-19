%plotalldays_inactivation.m
%written 8.20 in order to show each day's raw data which goes into overall
%stats.

function [ax]=inactivplot3(avls,graphvals,noteval,rawvl,ppf, indtoplot,ax)
%control of the subplots,ppf is plotsperfigure


figcount=1;
%1st column is raw data
%2nd column is histogram.
if rawvl=='raw'
    datvls=avls.rawvls{noteval}
else
    datvls=avls.adjvls{noteval}
end

%this is in avls.mulist coordinates...i.e. the 10th muscimol run is indexed
%with 10
if ~exist('indtoplot')
    indtoplot=avls.aclist(:,1);
end


for ii=1:length(indtoplot)
    %calculate the mod using plotsperfigure 
    mulstind=indtoplot(ii);
    pn=mod(ii,ppf)
    %if this is the last plot
    if(pn==0)
        pn=ppf;
    end
    
    if(pn==1)
        figure;
    end
    %this code will go through each all the acsf start_times and find the unique
%days to create
    ax1(pn)=subplot(ppf,2,pn*2-1)
            %plot acsf
            
      %need to create a zero reference time. this should be start of the mutimevec-the offset. if no timevec, 
      %this can be halfway between the acsf start_time and end_time.
      if(avls.mulist(mulstind)) 
        reftime=avls.tmvc(avls.mulist(mulstind))-avls.muoffset; 
        
      else
          reftime=mean([avls.tmvc(avls.aclist(mulstind),1) avls.tmvc(avls.aclist(mulstind),2)])
         
      end
           for jj=1:2
                    ind=avls.aclist(mulstind,jj)
                    if(ind)
                        plot(datvls{ind}(:,1)-reftime,datvls{ind}(:,2),'.','Color',graphvals.col(2),'MarkerSize',9)
                        hold on;
                        plot([-.15 .3], [avls.initmean{noteval} avls.initmean{noteval}],'k--');
                        xlim([-.15 .4])
%                         ylim([avls.edges{noteval}(1)*.95 avls.edges{noteval}(end)*1.05])
                    end
                    if(jj==1&pn==1)
                        titlestr=[avls.bname '-plot' num2str(figcount) '.eps']
                        title([avls.pvls{ind} avls.cvl{ind} ])
                    end
                    end

           %plot muscimol
            ind=avls.mulist(mulstind)
                    if(ind)
                    plot(datvls{ind}(:,1)-reftime,datvls{ind}(:,2),'.','Color',graphvals.col(1),'MarkerSize',9)
                    
                    hold on;
                    plot([avls.rawtimes(ind,1)-reftime avls.rawtimes(ind,2)-reftime],[avls.edges{1}(3) avls.edges{1}(3)],'r')
                    end
        %now work on the other subplot to make histograms.
        ax2(pn)=subplot(ppf,2,pn*2)
            for jj=1:2
                ind=avls.aclist(mulstind,jj);
                if(ind)
%                     stairs(avls.edges{noteval},avls.hst{noteval}{ind},'k')
                    hold on;
                    if (jj==1)
                        plotgaussian(ax2(pn),avls.edges{noteval},avls.mnvl{noteval}(ind),avls.stdv{noteval}(ind),'k','--',1);
                    else
                        plotgaussian(ax2(pn),avls.edges{noteval},avls.mnvl{noteval}(ind),avls.stdv{noteval}(ind),'k','--',1);
                    end
                    hold on;
                    nvl=length(avls.adjvls{noteval}{ind}(:,2));
                    axes(ax2(pn))
                    text(avls.edges{noteval}(end-3), .1*jj, ['n=' num2str(nvl)],'Color','k');
                    text(avls.edges{noteval}(3),0.15,['day ' num2str(ii)]);
                    plot([avls.initmean{noteval} avls.initmean{noteval}], [0 .2], 'k--');
                end
                muind=avls.mulist(mulstind);
                if(muind)
%                     stairs(avls.edges{noteval},avls.hst{noteval}{muind},'r')
                   
                    plotgaussian(ax2(pn),avls.edges{noteval},avls.mnvl{noteval}(muind),avls.stdv{noteval}(muind),'r','--',1);
                    nvl=length(avls.adjvls{noteval}{muind}(:,2));
                    axes(ax2(pn))
                    text(avls.edges{noteval}(end-3), 0.15, ['n=' num2str(nvl)],'Color','r');
                    
                    
                    end
            end
                axes(ax1(pn))
                box off;
                %set x dimensions to standard
                
%                 scale_ln=[0 avls.mnvl{noteval}]
                set(ax1(pn),'Xcolor','w')
                                %call fxn modplot
             axes(ax2(pn))
                box off;
                %set ydimensions to standard value
                %turn y_off
                %set xdimensions to min_max
                %scale_bar
        if (pn==ppf)||length(avls.aclist)
            axis(ax1)
            xlabel('time')
            ylabel('frequency')
            box off;
            axis(ax2)
            ylabel('probability')
            xlabel('frequency')
            box off;
            
                axes(ax2(1))        
                titlestr=[avls.bname '-plot' num2str(figcount) 'note' num2str(noteval) '.eps']
                title(titlestr);
                figcount=figcount+1;
            
            strcmd=['print -depsc ../figs/' titlestr];
            eval(strcmd);
        end
end