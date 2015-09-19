
function [ax]=plotwn(ax, fillind, graphvals, subval,ntvl)
           
         wnon=graphvals.wn
           x1=datenum([wnon{ntvl}.tmon],'yyyy-mm-dd  HH:MM')-subval;
           x2=datenum([wnon{ntvl}.tmoff],'yyyy-mm-dd  HH:MM')-subval;
           y1=wnon{ntvl}.freqlo;
           y2=wnon{ntvl}.freqhi;
           for ii=1:length(x1)
            ax=fill([x1(ii) x2(ii) x2(ii) x1(ii)], [y1(ii) y1(ii) y2(ii) y2(ii)], fillind,'EdgeColor','none')
            hold on; 
        end   