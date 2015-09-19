%plotbaselinestandard_dev.m
%switch to STANDARD DEVIATION
for ii=1:length(sumbs)
   acstdv=sumbs(ii).acstdz
   mustdv=sumbs(ii).mustdz
   aczer=sumbs(ii).acerrz
   muzer=sumbs(ii).muerrz
   basruns=sumbs(ii).basruns
   %switch this to all inds, not just targetind
   allnote=sumbs(ii).allnote;
   for jj=1:length(allnote)
       ntvl=allnote(jj);
       plot(acstdv(ntvl, basruns), mustdv(ntvl,basruns),'o');
       hold on;
   end
   plot([-1 100],[-1 100],'k-')
end
   
    
