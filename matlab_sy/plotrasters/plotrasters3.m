%This function takes as input a 2 column vector with spike times and with
%trials, which don't need to be unique.  It plots rasters and outputs a
%vector which maps each trial number to a trial number in the matrix.


function [axt,h1,h2]=plotrasters3(inarray,edges,hist,datalength)

flag=0;
for kk=1:length(inarray)
    if(isempty(inarray{kk})==0)
      
            numspikes=length(inarray{kk});
            spikes=inarray{kk}(:,1);
            spikes=spikes';
            x=[spikes; spikes];
            y=[ones(1,numspikes)+kk-0.4;ones(1,numspikes)+kk+0.4];
        
            if(kk==length(inarray)&nargin~=1)
                
                [axt,h1,h2]=plotyy(x,y,edges,hist)
                linkaxes(axt,'x');
                flag=1;
            elseif (nargin==1)
                plot(x,y,'k','Linewidth',1.6);
                hold on;
                flag=1;
            
            else
                
                plot(x,y,'r')
                
            hold on;
        end
    end
   
box off;
end
if (flag==0)
    [axt,h1,h2]=plotyy(x,y,edges,hist)
end