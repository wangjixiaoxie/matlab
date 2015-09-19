%plotcompar
%originally written 9.12.09 for plotting stim data
%takes a plotstruct with 
%ps.ax
%ps.col
%ps.xvecs - each column (matches columns of yvecs) is a series of xvalues
%to plot y values against
%ps.yvecs - each column (matches xvecs) is yvalues
%ps.col

function []=plotcompar(ps)

axes(ps.ax)
numvec=length(ps.xvecs)

for colind=1:numvec
    plot(ps.xvecs{colind}, ps.yvecs{colind},'.','Color',ps.col);
    hold on;
end

%now plot lines
if(numvec>1)
    for colind=2:2:numvec
       plot([ps.xvecs{colind-1}'; ps.xvecs{colind}'],[ps.yvecs{colind-1};ps.yvecs{colind}],'Color',ps.col) 
        
    end
    
end

if(ps.plotbar)
   for(colind=1:numvec)
       mnout=mean(ps.yvecs{colind})
       steout=std(ps.yvecs{colind})/length(ps.yvecs{colind})
       xvl=ps.xvecs{colind}(1);
       bar(xvl,mnout,'FaceColor','none');
       plot([xvl xvl], [mnout-steout mnout+steout],'k');
   end
end