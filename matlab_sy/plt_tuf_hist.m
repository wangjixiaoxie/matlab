function [TH]=plt_tuf_hist(yvals)


        
        
        ln=length(yvals);
        x1=[1:ln;1:ln]
        y1=[zeros(1,ln);yvals.*ones(1,ln)];
        plot(x1,y1,'k','Linewidth',2)
        hold on;
        x2=[(1:ln)-.5;1:ln]
        y2=[yvals;yvals]
        plot(x2,y2,'k','Linewidth',2)
      
%             numspikes=length(inarray{kk});
%             spikes=inarray{kk}(:,1);
%             spikes=spikes';
%             x=[spikes; spikes];
%             y=[ones(1,numspikes)+kk-0.4;ones(1,numspikes)+kk+0.4];
%         
%             if(kk==length(inarray)&nargin~=1)
%                 
%                 [axt,h1,h2]=plotyy(x,y,edges,hist)
%                 linkaxes(axt,'x');
%                 flag=1;
%             elseif (nargin==1)
%                 plot(x,y,'r');
%                 hold on;
%                 flag=1;
%             
%             else
%                 
%                 plot(x,y,'r')
%                 
%             hold on;
%         end
%     end