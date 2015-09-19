function [ax]=plotmeanstdev2(avls, datsource,graphvals,noteval,plotdir,axraw,axsum)
%rewritten extensively 6.24.08

%assume one plot for raw data, one plot for meanstd
noinput=0;
%    axes(axsum) 
    subval=floor(avls.rawtimes(1,1))

    if datsource=='raw'
        datvls=avls.rawvls{noteval}
    else
        datvls=avls.adjvls{noteval}
    end
    if(plotdir)
        dirdatvls=avls.dirvls{noteval}
    end

    ax=[];
    
if(~exist('axsum')&~exist('axraw'))
    noinput=1;
    ax(1)=subplot(2,1,1)
    ax(2)=subplot(212)
else
    if(~isempty(axraw))    
        ax(1)=axraw;
    else
        ax(2)=axsum;
    end
end
if (noinput==1||~isempty(axraw))
        axes(ax(1));
        if(exist('graphvals.wn'))
            fillind=[1 .9 .9]
       %this for loop plots white noise in blocks.
            plotwn(ax, fillind, graphvals,subval,noteval)
        end 
        for ii=1:graphvals.numcol
           
                for jj=1:length(graphvals.colvals{ii})
                    ind=graphvals.colvals{ii}(jj);
                    jj
                    ind
                    indtmp=find(datvls{ind}(:,1)-subval>20&datvls{ind}(:,1)-subval<23)
          if(indtmp)
              jjj=4;
          end
                    plot(datvls{ind}(:,1)-subval,datvls{ind}(:,2),'.','Color',graphvals.col(ii),'MarkerSize',9)
                    hold on;
                    
                    
                end
%                 if(plotdir)
%                     for jj=1:length(avls.diron)
%                         indvl=avls.diron(jj)
%                         plot(dirdatvls{ind}(:,1)-subval,dirdatvls{ind}(:,2),'.','Color',graphvals.col(3),'MarkerSize',9)
%                     end  
%                 end
                
        end
    end
if(exist('axsum')||noinput==1)
   if(noinput==1)
    ax(2)=subplot(212);
   else
    ax(2)=axsum;
    end

        axes(ax(2));
   
      if(exist('graphvals.wn'))
        fillind=[1 .9 .9]
        %this for loop plots white noise in blocks.
        plotwn(ax, fillind, graphvals,subval,noteval)
      end

           %get parameters for mean plot
          x=mean(avls.adjtimes,2)
          x=x-subval
          y=avls.mnvl{noteval};
          e=avls.stdv{noteval};
    
       %this line plots the mean, std.   
       for ii=2:2%graphvals.numcol
              colind=graphvals.colvals{ii}
              xvec=[x(colind) x(colind)]
              xvec=xvec'
              yvec=[y(colind)+e(colind);y(colind)-e(colind)]
              yvec2=[y(colind);y(colind)]
              
              plot(xvec, yvec, 'Color',graphvals.col(ii),'Linewidth',3);
              hold on;
              plot(xvec(1,:), yvec2(1,:), '-','Color',graphvals.col(ii))
       end
       
       if(~isempty(avls.diron))

          x=mean(avls.dirdays,2)
          x=x-subval
          y=avls.dirmnvl{noteval};
          e=avls.dirstdv{noteval};
          colind=graphvals.colvals{3};
          xvec=[x x]'
          yvec=[y+e;y-e]
          yvec2=[y;y]
          plot(xvec,yvec,'Color',graphvals.col(3),'Linewidth',3);
          hold on;
          plot(xvec(1,:),yvec2(1,:),'o','Color',graphvals.col(ii))
       end




end
%     if(graphvals.plttext)
%         for ii=1:graphvals.numcol
%             for jj=1:length(graphvals.colvals{ii})
%                 ind=graphvals.colvals{ii}(jj)
%                 numval=length(chkarray{ind}(:,1))
%             %offset x value by half a day
% %                 text(x(ind)+.5,y(ind),[ num2str(numval)],'Color',graphvals.col(ii));
%         end
%         end
%     end
    
%     ylabel('frequency (hz)', 'Fontsize',18)
%     xlabel('day', 'Fontsize',18)

% if exist('plotxlim')
%     for ii=1:length(ax)
%         axes(ax(ii));
%         xlim([plotxlim(1) plotxlim(2)])
%     end
% end
 end


    