%plot_tmcourse_stimall.m
function [stimvls,catchvls] = plot_tmcourse_pharmall(phsumbs,avls)
%set bsind
bsind=2;
%num2average=20;
avls.endptvls=15;

%shiftind=4;
shiftindvl=8;
crbs=phsumbs(bsind);

%load avls.

indtoplot=[1:43]
ybnds=[6000 8000]


% shiftruns=crbs.shiftruns{shiftindvl};
% shiftstanruns=intersect(shiftruns,crbs.STANRUNS);

%plot edges


for ii=1:length(indtoplot);
    crind=indtoplot(ii);
    if(ii==1)
       avls.bastime=crbs.flrtmvec(crind); 
       plotedge(ybnds,1:30);
    end
    %acsf only.
    if(crbs.aclist(crind,2)==0)
        %just plot ends
        ind=crbs.aclist(crind,1);
        plotends(crbs,avls,ind,[1 2]);
    else
        acinitind=crbs.aclist(crind,1);
        acfinind=crbs.aclist(crind,2);
        muind=crbs.mulist(crind);
%         avls.endptvls=endptvls;
        plotends(crbs,avls,acinitind,[1 2]);
        plotends(crbs,avls,acfinind,2);
        if(muind)
        plotmu(crbs,avls,muind); 
        end
        end
  tst=1;      
end
    

   
   
  function [xvls,stimvls,catchvls] = plotends(crbs,avls,plotind,edgevls)
    bastime=avls.bastime;
    edgenum=avls.endptvls;
    %plotbeginning
    if(ismember(1,edgevls))
        ln=avls.rawvls{1}{plotind}(:,2);
        if(ln>edgenum)
            ind=1:edgenum;
        else
            ind=1:ln
        end
        yvls=avls.rawvls{1}{plotind}(ind,2);
        xvls=avls.rawvls{1}{plotind}(ind,1)-bastime;
        plot(xvls,yvls,'b.');
        hold on;
        plot(mean(xvls),mean(yvls),'ko','Markersize',8,'MarkerFaceColor','k');
    end
      
    
    if(ismember(2,edgevls))
        numvls=length(avls.rawvls{1}{plotind}(:,1));
        if(numvls>edgenum)
            ind=numvls:-1:numvls-edgenum;
        else
            ind=numvls:-1:1;
        end
            yvls=avls.rawvls{1}{plotind}(ind,2);
            xvls=avls.rawvls{1}{plotind}(ind,1)-bastime;
            plot(xvls,yvls,'b.');
            hold on;
            plot(mean(xvls),mean(yvls),'ko','Markersize',8,'MarkerFaceColor','k');
        end
  
  
    
  function [initxvls,endxvls] = plotmu(crbs,avls,muind)
      bastime=avls.bastime;
      edgenum=avls.endptvls;
         yvls=avls.adjvls{1}{muind}(:,2);
        xvls=avls.adjvls{1}{muind}(:,1)-bastime;
      
      
       plot(mean(xvls),mean(yvls),'ro','Markersize',8,'MarkerFaceColor','r');
    
     
   
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
               
       
    

