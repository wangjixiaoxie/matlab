%plotscatterxy3
function [sumstats]=plotscatterbas(sumbs);
figure

colvec={'b' 'k' 'c'}

acvlcomb=[];
mcvlcomb=[];
acmncomb=[];
mumncomb=[];
bsvls=[];
for ii=1:length(sumbs)
   [ac_cvls, mu_cvls,acmn,mumn]=calc_cv(sumbs(ii));
   acvlcomb=[acvlcomb ac_cvls]
   mcvlcomb=[mcvlcomb mu_cvls]
   bsvls=[bsvls ii*ones(1,length(ac_cvls))]
   %i.e. example
   if ii==2 
    plot(ac_cvls, mu_cvls, 'ro');
   else
       plot(ac_cvls,mu_cvls,'ko');
   end
          hold on;
            
plot([-10 10],[-10 10],'k-')

sumstats.acvl=acvlcomb;
sumstats.mcvlcomb=mcvlcomb;
sumstats.bsvls=bsvls
end
figure
for ii=1:length(sumbs)
   [ac_cvls, mu_cvls,acmn,mumn]=calc_cv(sumbs(ii));
   acmncomb=[acmncomb acmn]
   mumncomb=[mumncomb mumn]
   bsvls=[bsvls ii*ones(1,length(acmn))]
   %i.e. example
   if ii==2 
    plot(acmn, mumn, 'ro');
   else
       plot(acmn,mumn,'ko');
   end
          hold on;
            
plot([-10 10],[-10 10],'k-')

sumstats.acmn=acmncomb;
sumstats.mmumn=mumncomb;
sumstats.bsvls=bsvls
end
function [ac_cvls, mu_cvls,sc_acmean,sc_mumean]=calc_cv(sumbs);
    ntind=sumbs.ntind;
    acmean=sumbs.acmean(ntind,sumbs.basruns)
    sc_acmean=acmean./sumbs.sfact(ntind);
    mumean=sumbs.mumean(ntind,sumbs.basruns)
    sc_mumean=mumean./sumbs.sfact(ntind);
    acstdv=sumbs.acstdv(ntind,sumbs.basruns);
    mustdv=sumbs.mustdv(ntind,sumbs.basruns);
    
    ac_cvls=acstdv./acmean;
    mu_cvls=mustdv./mumean;








    
