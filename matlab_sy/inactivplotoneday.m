%plotalldays_inactivation.m
%written 8.20 in order to show each day's raw data which goes into overall
%stats.

function [ax]=inactivplotoneday(avls,graphvals,noteval,rawvl,indtoplot,ax)
%control of the subplots,ppf is plotsperfigure


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
    indtoplot=avls.aclist;
end
ax(1)=subplot(2,1,1);
ax(2)=subplot(212);
for ii=1:length(avls.aclist(indtoplot))
    %calculate the mod using plotsperfigure 
      %if there was a muscimol run
      if(avls.mulist(indtoplot)) 
        reftime=avls.tmvc(avls.mulist(indtoplot))-avls.muoffset; 
      %if just acsf
      else
          reftime=mean([avls.tmvc(avls.aclist(indtoplot),1) avls.tmvc(avls.aclist(indtoplot),2)])
      end
      for jj=1:2
           ind=avls.aclist(indtoplot,jj)
           if(ind)
                plot(datvls{ind}(:,1)-reftime,datvls{ind}(:,2),'.','Color',graphvals.col(2),'MarkerSize',12)
                hold on;
                plot([-.15 .3], [avls.initmean{noteval} avls.initmean{noteval}],'k--');
                xlim([-.15 .4])
            end
            if(jj==1)
               titlestr=[avls.bname '-plot.eps']
               title([avls.pvls{ind} avls.cvl{ind} ])
            end
       end

           %plot muscimol
       ind=avls.mulist(indtoplot)
       if(ind)
          plot(datvls{ind}(:,1)-reftime,datvls{ind}(:,2),'.','Color',graphvals.col(1),'MarkerSize',12)
          hold on;
          plot([avls.rawtimes(ind,1)-reftime avls.rawtimes(ind,2)-reftime],[avls.edges{1}(3) avls.edges{1}(3)],'r')
       end
        %now work on the other subplot to make histograms.
%         ax2(pn)=subplot(ppf,2,pn*2)
%             for jj=1:2
%                 ind=avls.aclist(mulstind,jj);
%                 if(ind)
% %                     stairs(avls.edges{noteval},avls.hst{noteval}{ind},'k')
%                     hold on;
%                     if (jj==1)
%                         plotgaussian(ax2(pn),avls.edges{noteval},avls.mnvl{noteval}(ind),avls.stdv{noteval}(ind),'k');
%                     else
%                         plotgaussian(ax2(pn),avls.edges{noteval},avls.mnvl{noteval}(ind),avls.stdv{noteval}(ind),'k','--');
%                     end
%                     hold on;
%                     nvl=length(avls.adjvls{noteval}{ind}(:,2));
%                     axes(ax2(pn))
%                     text(avls.edges{noteval}(end-3), .1*jj, ['n=' num2str(nvl)],'Color','k');
%                     text(avls.edges{noteval}(3),0.15,['day ' num2str(ii)]);
%                     plot([avls.initmean{noteval} avls.initmean{noteval}], [0 .2], 'k--');
%                 end
%                 muind=avls.mulist(mulstind);
%                 if(muind)
% %                     stairs(avls.edges{noteval},avls.hst{noteval}{muind},'r')
%                    
%                     plotgaussian(ax2(pn),avls.edges{noteval},avls.mnvl{noteval}(muind),avls.stdv{noteval}(muind),'r');
%                     nvl=length(avls.adjvls{noteval}{muind}(:,2));
%                     axes(ax2(pn))
%                     text(avls.edges{noteval}(end-3), 0.15, ['n=' num2str(nvl)],'Color','r');
%                     
%                     
%                     end
%             end
                axes(ax(1))
                box off;
                %set x dimensions to standard
                
%                 scale_ln=[0 avls.mnvl{noteval}]
                set(ax(1),'Xcolor','w')
                                %call fxn modplot
%              axes(ax2(pn))
%                 box off;
                %set ydimensions to standard value
                %turn y_off
                %set xdimensions to min_max
                %scale_bar
        
            axis(ax(1))
            xlabel('time')
            ylabel('frequency')
            box off;
            axis(ax(2))
            ylabel('probability')
            xlabel('frequency')
            box off;
            
                axes(ax(1))        
                titlestr=[avls.bname '-plotnote' num2str(noteval) '.eps']
                title(titlestr);
                
            strcmd=['print -depsc ../figs/' titlestr];
            eval(strcmd);
        
end