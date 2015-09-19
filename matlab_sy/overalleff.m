function [sumbs]=overalleff(sumbs)  
       for ii=1:length(sumbs)
         crbs=sumbs(ii);
         
             ind=find(crbs.combeff>0&isnan(crbs.combeff)==0&isinf(crbs.combeff)==0)
             sumbs(ii).mneff=mean(crbs.combeff(ind));
             sumbs(ii).stdeff=std(crbs.combeff(ind))
             
        end