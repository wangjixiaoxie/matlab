function []=plotexhistinactiv(avls,ps)
    if(isfield(ps,'FLIP'))
        if(ps.FLIP)
            FLIP=1;
        else
            FLIP=0
        end
    else
        FLIP=0;
    end
    
    if(isfield(ps,'SCALEHIST'))
        if(ps.SCALEHIST)
            SCALEHIST=1;
        else
            SCALEHIST=0;
        end
    else
        SCALEHIST=0;
    end
    
    
        
    
    dht=ps.distht;
    arht=ps.arrowht;
    ntind=ps.ntind;
        for indnum=1:length(ps.indtoplot)
            crind=ps.indtoplot(indnum);
            crind=avls.mulist(crind);
            axes(ps.ax);
            
            crvls{1}=avls.hstacpre{ntind}{crind};
            crvls{2}=avls.hstmu_comb{ntind}{crind};
            crvls{3}=avls.hstacpst{ntind}{crind};
            
            %throw out last value
            crvls{1}=crvls{1}(1:end-1);
            crvls{2}=crvls{2}(1:end-1);
            crvls{3}=crvls{3}(1:end-1);
            
            
            indvl=find(avls.mulist==crind);
            acind=avls.aclist(indvl,1);
            muind=avls.mulist(indvl);
            
            cracvls=avls.adjvls{ntind}{acind}(:,2);
            crmuvls=avls.adjvls{ntind}{muind}(:,2);
            crmn{1}=mean(cracvls);
            crmn{2}=mean(crmuvls);
            %throw out last value
            credges=avls.edges{ntind}(1:end-1);
            credgediff=credges(2)-credges(1);
            credges=credges+.5*credgediff
            
            crste{1}=std(cracvls)./sqrt(length(cracvls));
            crste{2}=std(crmuvls)./sqrt(length(crmuvls));
%             crctind=avls.crctind{crind};
%             crfbind=avls.crfbind{crind};
%             
%             crctste=stdct./sqrt(length(crctind));
%             crfbste=stdfb./sqrt(length(crfbind));
            %plot catch value
            if(ps.plotpre)
             
                numvls=1;
            else
                if(isfield(ps,'omitpost'))
                    if(ps.omitpost)
                        numvls=[1 2];
                    else
                        numvls=[1 2 3];
                    end
                else
                    numvls=[1 2 3]
                end
            end
          
            
            
            for ii=numvls
                if(ps.omitzero)
                    ind=1:length(credges);
                   indout=find(crvls{ii}>0);
                   if(indout(1)>1)
                       instart=indout(1)-1;
                   else
                       instart=1;
                   end
                   if(indout(end)<ind(end))
                       indend=indout(end)+1;
                   else
                       indend=ind(end);
                   end
                       
                   ind=ind(instart):ind(indend);
                else
                    ind=1:length(credges);
                    
                end
                  if(SCALEHIST)
                       crvlsout{ii}=crvls{ii}*ps.SCALEGAIN;
                    crvlsout{ii}=crvlsout{ii}+ps.SCALEOFFSET;
                    
%                     crmnout{ii}=crmn{ii}*ps.SCALEGAIN;
%                     crmnout{ii}=crmnout{ii}+ps.SCALEOFFSET;
                   
                  end
               
                
                  
                  
                
                if(FLIP)
                   xvls=crvlsout{ii}(ind);
                   yvls=credges(ind)/3000;
                else
                    xvls=credges(ind)/3000;
                    yvls=crvlsout{ii}(ind);
                end
                stairs(xvls,yvls,'Color',ps.col{ii},'Linewidth',2)
                hold on;
            %    plot([mnbas+stdbas mnbas+stdbas], [0 1], 'c--')
            %    plot([mnbas-stdbas mnbas-stdbas], [0 1], 'c--')
                plot([ps.initmean/3000 ps.initmean/3000],[0 1],'k--','Linewidth',2)
                
                if(ps.plot_triangles)
                    if(FLIP)
                        x=dht
                        y=crmn{ii}/3000;
                        marktype='<'
                    else
                        x=crmn{ii}/3000;
                        y=dht
                        marktype='^'
                    end
                    plot(x,y,'Marker',marktype,'Color',ps.col{ii},'MarkerFaceColor',ps.col{ii})
                end
            
%                 text(2.225,0.5,['catch=' num2str(length())],'Color','k');
%                 if(ii<3)
%                     plot([(crmn{ii}-crste{ii})/3000 (crmn{ii}+crste{ii})/3000],[dht dht],'k','Linewidth',3)
%                     
%                 end
            end
             %this is for the stim trial
            
           
  
       origarrow=[ps.initmean/3000 arht(1) crmn{1}/3000  arht(1)]
       shiftarrow=[crmn{1}/3000 arht(2) crmn{2}/3000  arht(2)]
   
%    ps.plotshiftarrow=1;
%    ps.plotrevarrow=1;

if(isfield(ps,'axbnds'))
    axis(ps.axbnds);
end
% plotarrows(origarrow,shiftarrow,ps);
   
   
   end