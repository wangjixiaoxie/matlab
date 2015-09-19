


function []=plotrawlearning(avls,startind,endind,ps)
if (isfield(ps,'runstoplot'))
    nozeroind=ps.runstoplot
else
    nozeroind=find(avls.aclist(:,1)~=0);
end
start_time=[avls.rawtimes(startind,1)] ;

if(isfield(ps,'plotedge'))
  if(ps.plotedge)
    plotedge(ps.plotedgeybnds,ps.plotedgexbnds);
  end
end





for ii=1:length(nozeroind)
    crvls=avls.allvls{1}{avls.aclist(nozeroind(ii),1)}
    hold on;
   [tms,avgvls,ervls]=calcavgvls(crvls)
   if(~isempty(tms)) 
   
   hold on;
   plot(crvls(:,1)-start_time,crvls(:,2)/3000,'o','MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerSize',1.5)
   plot(tms-start_time,avgvls/3000,'k','Linewidth',2)
   end
    %calcrunningaverage

end


% nozeroind=find(avls.aclist(:,2)~=0);

for ii=1:length(nozeroind)
    if((avls.aclist(nozeroind(ii),2)))
        crvls=avls.allvls{1}{avls.aclist(nozeroind(ii),2)}
        hold on;
        [tms,avgvls,ervls]=calcavgvls(crvls)
        if(~isempty(tms))
            plot(crvls(:,1)-start_time,crvls(:,2)/3000,'o','MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerSize',1.5)
            plot(tms-start_time,avgvls/3000,'k','Linewidth',2)
        end
    end
    %calcrunningaverage
end

% nozeroind=find(avls.mulist~=0);
for ii=1:length(nozeroind)
   if(avls.mulist(nozeroind(ii)))
    crvls=avls.allvls{1}{avls.mulist(nozeroind(ii))}
    hold on;
    
    [tms,avgvls,ervls]=calcavgvls(crvls)
    if(~isempty(tms))
        plot(crvls(:,1)-start_time,crvls(:,2)/3000,'o','MarkerFaceColor',[.7 .7 1],'MarkerEdgeColor',[.7 .7 1],'MarkerSize',1.5)
    plot(tms-start_time,avgvls/3000,'Color',[.3 .3 1],'Linewidth',3)
    end
    %calcrunningaverage
   end
end

    
    
    
    

   

function [tms,avgvls,ervls]=calcavgvls(crvls)
    num2avg=15;
    if(length(crvls(:,1))<=num2avg)
        num2avg=length(crvls(:,1));
    end
        for ii=num2avg:length(crvls(:,1))
            avgvls(ii-num2avg+1)=mean(crvls((ii-num2avg+1):ii,2));
            tms(ii-num2avg+1)=mean(crvls((ii-num2avg+1):ii,1))
            ervls(ii-num2avg+1)=std(crvls((ii-num2avg+1):ii,2))./sqrt(num2avg);
        end

        
   function []=plotedge(ybnds,xvecin)
%            minx=min([makerow(initxvls) makerow(xvls) makerow(enxvls)]);
%            maxx=max([makerow(initxvls) makerow(xvls) makerow(enxvls)]);
%            st_x=floor(minx);
%            max_x=ceil(maxx);
           
           for ii=1:length(xvecin)
               crx=xvecin(ii);
              xbnds(1)=crx-9/24;
              xbnds(2)=crx;
              xvec=[xbnds xbnds(end:-1:1)]
              yvec=[ybnds(1) ybnds(1) ybnds(2) ybnds(2)]
              fill(xvec,yvec,[0.7 0.7 0.7],'edgecolor','none');
              hold on;
               
               
           end
        
        
% %     else
%         avgvls=[];
%         tms=[];
%         ervls=[];
%     end
    

