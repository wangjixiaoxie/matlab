function [ax]=plotmeanstdev(chkarray, stimes, etimes, graphvals, plotbounds,mnval,stdev, chunkvals)

 
   
   %NEED TO WRITE SOMETHING IN CASE GRAPHVALS IS EMPTY
   %this is loop for every category of color
   %if(chkarray)
   
   %fillbounds should be a 2-D vector with start times and end times of
   %muscimol on.
  
   %used to start first day at 1
    if (exist('mnval'))
        plotnum=2
    else
        plotnum=1
    end
    
    
   subplot(plotnum,1,1)
   
   subval=floor(chkarray{1}(1,1))
   
   if(plotbounds)
       fillind=[1 .9 .9]
       
%        wnon=graphvals.wnon;
%        wnoff=graphvals.wnoff;
%        minbnds=graphvals.minbnds;
%        maxbnds=graphvals.maxbnds;
%        
%        for ii=1:length(wnon)
%            x1=datenum([wnon{ii}],'yyyy-mm-dd  HH:MM:SS')-subval;
%            x2=datenum([wnoff{ii}],'yyyy-mm-dd  HH:MM:SS')-subval;
%            y1=minbnds{ii};
%            y2=maxbnds{ii};
%            ax=fill([x1 x2 x2 x1], [y1 y1 y2 y2], fillind,'EdgeColor','none')
%            
%            hold on; 
%         
%        end
       
       
       for ii=1:length(graphvals.timoff)
%            if(~isempty(graphvals.timoff{ii}))
%                x1=datenum([graphvals.date{ii} ' ' graphvals.timon{ii}],'yyyy-mm-dd  HH:MM:SS')-subval
%                x2=datenum([graphvals.date{ii} ' ' graphvals.timoff{ii}], 'yyyy-mm-dd  HH:MM:SS')-subval
%                
%            
%                 y1=graphvals.tickht
%                 
%                 plot([x1 x2], [y1 y1],'r','Linewidth',7)
%            end
       end
       
       
   end
   
   
   
   for ii=1:graphvals.numcol
       if(graphvals.pltpnts)
       for jj=1:length(graphvals.colvals{ii})
            ind=graphvals.colvals{ii}(jj);
            plot(chkarray{ind}(:,1)-subval,chkarray{ind}(:,2),'.','Color',graphvals.col(ii),'MarkerSize',9)
            hold on;
        
       end
       end
   end
   
   
   if exist('chunkvals')
       for ii=1:length(graphvals.timon)
          
           if(~isempty(graphvals.timon{ii}))
               plotind=find(chunkvals.mean{ii}~=0);
               errorbar(chunkvals.tmvl{ii}(plotind)-subval,chunkvals.mean{ii}(plotind),chunkvals.std{ii}(plotind),'r+');
           end
       end
   end
           if (exist('mnval'))
          subplot(2,1,2)
%            for ii=1:length(wnon)
%            x1=datenum([wnon{ii}],'yyyy-mm-dd  HH:MM:SS')-subval;
%            x2=datenum([wnoff{ii}],'yyyy-mm-dd  HH:MM:SS')-subval;
%            y1=minbnds{ii};
%            y2=maxbnds{ii};
%            ax=fill([x1 x2 x2 x1], [y1 y1 y2 y2], fillind,'EdgeColor','none')
%            
%            hold on; 
%         
%        end    
          x=mean([stimes ;etimes],1)
          x=x-x(1)+1
          y=mnval;
          e=stdev;
%           
%           for ii=1:length(wnon)
%            x1=datenum([wnon{ii}],'yyyy-mm-dd  HH:MM:SS')-subval;
%            x2=datenum([wnoff{ii}],'yyyy-mm-dd  HH:MM:SS')-subval;
%            y1=minbnds{ii};
%            y2=maxbnds{ii};
%            ax=fill([x1 x2 x2 x1], [y1 y1 y2 y2], fillind,'EdgeColor','none')
%            
%            hold on; 
%         
%        end
          
          
          
          for ii=1:graphvals.numcol
              colind=graphvals.colvals{ii}
              xvec=[x(colind);x(colind)]
              yvec=[y(colind)+e(colind);y(colind)-e(colind)]
              yvec2=[y(colind);y(colind)]
             
%                 ofst=graphvals.offset(ii)
             
              
             
                plot(xvec, yvec, 'Color',graphvals.col(ii),'Linewidth',3);
              hold on;
             plot(xvec(1,:), yvec2(1,:), 'o','Color',graphvals.col(ii))
          end
       end
    if(graphvals.plttext)
        for ii=1:graphvals.numcol
            for jj=1:length(graphvals.colvals{ii})
                ind=graphvals.colvals{ii}(jj)
                numval=length(chkarray{ind}(:,1))
            %offset x value by half a day
%                 text(x(ind)+.5,y(ind),[ num2str(numval)],'Color',graphvals.col(ii));
        end
        end
    end
    
    ylabel('frequency (hz)', 'Fontsize',18)
    xlabel('day', 'Fontsize',18)
    box off;
    