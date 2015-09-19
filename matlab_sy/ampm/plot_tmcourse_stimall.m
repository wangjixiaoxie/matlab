%plot_tmcourse_stimall.m
function [stimvls,catchvls] = plot_tmcourse_stimall(sumbs)
%set bsind
bsind=1;
%num2average=20;
endptvls=15;

%shiftind=4;
shiftindvl=2;
crbs=sumbs(1);

shiftruns=crbs.shiftruns{shiftindvl};
shiftstanruns=intersect(shiftruns,crbs.STANRUNS);

indstim=find(crbs.NOSTIM(shiftstanruns)==0);
indnostim=find(crbs.NOSTIM(shiftstanruns));
shiftstim=shiftstanruns(indstim);
shiftnotstim=shiftstanruns(indnostim);

ybnds=[7000 8000];
   
   %search for a nonstim run on same day.
   plotedge(ybnds,1:50); 
   hold on;
   [xvls,stimvls, catchvls]=plot_stim(crbs,shiftindvl,shiftstim);
   [inxvls,endxvls]=plot_notstim(crbs,shiftindvl,shiftnotstim,endptvls);
   
   
  function [xvls,stimvls,catchvls] = plot_stim(crbs,shiftindvl,shiftstim)
    
    xind=find(ismember(crbs.shiftruns{shiftindvl}, shiftstim))
    xvls=crbs.exactshifttms{shiftindvl}(xind);
    
    stimvls=crbs.mumean(shiftstim);
   catchvls=crbs.acmean(shiftstim);
    plot(xvls,stimvls,'ro','Markersize',6,'MarkerFaceColor','r');
    hold on;
     plot(xvls,catchvls,'ko','Markersize',6,'MarkerFaceColor','k');
    
  function [initxvls,endxvls] = plot_notstim(crbs,shiftindvl,shiftnotstim,endptvls)
      shiftvl=shiftindvl;
      runvls=shiftnotstim;
      edgevl=endptvls;
      
      [initxvls, inityvls]=getedgevls(crbs,shiftvl,runvls,edgevl,'ini');
      [endxvls, endyvls]=getedgevls(crbs,shiftvl,runvls,edgevl,'fin');
       plot(initxvls,inityvls,'bo','Markersize',6,'MarkerFaceColor','b');
    hold on;
     plot(endxvls,endyvls,'bo','Markersize',6,'MarkerFaceColor','b');
     
   function [xvls,yvls] = getedgevls(crbs,shiftvl,runvls,edgevl,type)
       for ii=1:length(runvls)
           runind=runvls(ii);
           crptvls=crbs.ptvls{runind}
           ln=length(crptvls(:,1));
           bastime=crbs.baswntime{shiftvl};
           initind=1:edgevl;
           finind=ln:-1:ln-edgevl
           inxtmsinit=crptvls(initind,1)-bastime;
           inyvlsinit=crptvls(initind,2);
           finxtms=crptvls(finind,1)-bastime;
           finyvls=crptvls(finind,2);
           if(type=='ini')
               yvls(ii)=mean(crptvls(1:edgevl,2));
               xvls(ii)=mean(crptvls(1:edgevl,1))-crbs.baswntime{shiftvl};
               
               plot(inxtmsinit,inyvlsinit,'b.');
           else
               yvls(ii)=mean(crptvls(end-edgevl:end,2));
               xvls(ii)=mean(crptvls(end:-1:end-edgevl,1))-crbs.baswntime{shiftvl};
               plot(finxtms,finyvls,'b.');
           end
       end
       function []=plotedge(ybnds,xvecin)
%            minx=min([makerow(initxvls) makerow(xvls) makerow(enxvls)]);
%            maxx=max([makerow(initxvls) makerow(xvls) makerow(enxvls)]);
%            st_x=floor(minx);
%            max_x=ceil(maxx);
           
           for ii=[xvecin]
               crx=xvecin(ii);
              xbnds(1)=crx-3/24;
              xbnds(2)=crx+6/24;
              xvec=[xbnds xbnds(end:-1:1)]
              yvec=[ybnds(1) ybnds(1) ybnds(2) ybnds(2)]
              fill(xvec,yvec,[0.7 0.7 0.7],'edgecolor','none');
              hold on;
               
               
           end
               
       
    

