%12.10.08
%plotbas2, plots baseline cv reduction as xy scatter.
%combine days.
%combdays 1, calculate mncv, stdcv across all bas runs, combdays=0 plot
%ever days individually

%same 
function []=plotbas2(sumbs, combdays);
figure;
for ii=1:length(sumbs)
   %acmean
   acmean=sumbs(ii).acmean;
   %mumean
   mumean=sumbs(ii).mumean;
   acstdv=sumbs(ii).acstdz
   mustdv=sumbs(ii).mustdz
   aczer=sumbs(ii).acerrz
   muzer=sumbs(ii).muerrz
   basruns=sumbs(ii).basruns
   %switch this to all inds, not just targetind
   allnote=sumbs(ii).allnote;
   
   [ac_cv]=calc_cv(acmean(allnote(1),basruns), acstdv(allnote(1), basruns));
   [mu_cv]=calc_cv(mumean(allnote(1),basruns), mustdv(allnote(1),basruns));
   
   
       if(combdays)
            plot(mean(ac_cv)*100, mean(mu_cv)*100,'o');
       else
           plot(ac_cv*100, mu_cv*100,'o');
       end
       hold on;
  
   plot([-1 100],[-1 100],'k-')
end
function [cv_vec]=calc_cv(meanmat, stdvmat);
    cv_vec=stdvmat./meanmat;
    
    
