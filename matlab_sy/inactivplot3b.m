%plotalldays_inactivation.m
%written 8.20 in order to show each day's raw data which goes into overall
%stats.

%changed so as only to plot inactivation days.

function [ax]=inactivplot3b(bs,bsind,rawvl,ppf)
    [avls,graphvals]=loadsumdata(bs,bsind);

figcount=1;
ntvl=bs(bsind).ntind
%1st column is raw data
%2nd column is histogram.
if rawvl=='raw'
    datvls=avls.rawvls{ntvl}
else
    datvls=avls.adjvls{ntvl}
end

%this is in avls.mulist coordinates...i.e. the 10th muscimol run is indexed
%with 10
if ~exist('indtoplot')
    indtoplot=find(avls.mulist>0);
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
    
    ax1(pn)=subplot(ppf,2,pn*2-1)
            %plot acsf
            
      %need to create a zero reference time. this should be start of the mutimevec-the offset. if no timevec, 
      %this can be halfway between the acsf start_time and end_time.
      if(avls.mulist(mulstind)) 
        reftime=avls.tmvc(avls.mulist(mulstind))-avls.muoffset; 
        
      else
          reftime=mean([avls.tmvc(avls.aclist(mulstind),1) avls.tmvc(avls.aclist(mulstind),2)])
         
      end
      
      %this plots rawdata
      %ps struct requires .ax, .col, .notevl, .indtoplot,
      %.rawvl,marksz;
      ps.ax=ax1(pn);
      col{1}=[0 0 0]
      col{2}=[0 0 0]
      col{3}=[0.4 0.4 1]
      ps.col=col;
      ps.ntvl=ntvl;
      ps.marksize=11;
      ps.indtoplot=mulstind;
      ps.rawvl='raw';
      
      inactiv_rawpoints(avls,graphvals,ps);
% 
%                     end
%                     if(jj==1&pn==1)
%                         titlestr=[avls.bname '-plot' num2str(figcount) '.eps']
%                         title([avls.pvls{ind} avls.cvl{ind} ])
%                     end
%                     end

           
        %now work on the other subplot to make histograms.
        ax2(pn)=subplot(ppf,2,pn*2)
        
        for ii=1:3
            hs(ii).wid=[2];
            hs(ii).orient=['norm'];
            hs(ii).ntvl=ps.ntvl;
            hs(ii).lw=2;
            hs(ii).col=col{ii}
            hs(ii).ax=ax2(pn);
            hs(ii).ls='-'
        end
        hs(2).lw=1;
        hs(2).ls='--'
        
        hs=gethsts(hs, mulstind,avls);
        hs=getedges(hs,mulstind,avls);
        hs(ii).ls     
        hs(1).marklnht=.65
        hs(2).marklnht=.75
        hs(3).marklnht=.7
        hs(1).markht=.65
        hs(2).markht=.75
        hs(3).markht=.7
        plothists(hs(1:3));
        
%         
%         for jj=1:2
%                 ind=avls.aclist(mulstind,jj);
%                 if(ind)
% %                     stairs(avls.edges{noteval},avls.hst{noteval}{ind},'k')
%                     hold on;
%                    
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
%                     text(avls.edges{noteval}(end-3), 0.15, ['n=' num2str(nvl)],'Color','r')            
%                     end
%             end
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
                titlestr=[avls.bname '-plot' num2str(figcount) 'note' num2str(ntvl) '.eps']
                title(titlestr);
                figcount=figcount+1;
            
            strcmd=['print -depsc ../figs/' titlestr];
            eval(strcmd);
        end
end