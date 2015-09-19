%written 4.1.10 to plot group data histograms
%input combvls takes the following form.
%combvls{1} - baseline runs
%combvls{2} - shift runs

% combvls{X}{1}-upshift
% combvls{X}{2}-downshift

%combvls{X}{X}{1}-shift
%combvls{X}{X}{2}-control.

%try to write compatible with stim and shift

%rewritten to output mean value if multiple shifts.

function [plotvls,stats,ctrlvls]=plotgroupbarinit(combvls,ps,plotin)   
    
%     ax=ps.ax;
    ps.col={'k' 'r' 'k'}
    xvls=1:3
    if(ps.TYPE=='plotall')
        [sumdata,sumdatarev]=plotall(combvls,xvls,ps);
    elseif(ps.TYPE=='plotsum')
        axes(ps.axsum)
        [plotvls,ctrlvls]=plotsum(combvls,xvls,ps,plotin);
        axes(ps.axpct)
        plotpct(plotvls,ps);
        axes(ps.axbas)
        if(ps.plotcv)
        plotcvbas(plotvls,ps);
        end
        [stats]=calcsumstats(plotvls)
    end
           function []=plotcvbas(plotvls,ps)
            figure
            %type 1 is baseline, type 2 is shift
            for type=1:2
                plotxvl=type;
            sumdatash=plotvls{type}{1}
            % this determines whether to combine across upshift and
            % downshift
            if(type==1)
                COMB=1;
            else
                COMB=1
            end
             
            indup=find(sumdatash.drxn==1);
            inddn=find(sumdatash.drxn==2);
           lidind=sumdatash.lid;
           muind=setdiff(1:length(sumdatash.drxn),lidind);
               yvl=sumdatash.cvred -1
               yplotind=find(~isnan(yvl));
                     lnvls=length(find(~isnan(yvl)));
                     if(COMB)
                        xvl=type*ones(lnvls,1);   
                     else
                        xvl1=type*ones(length(indup),1);
                        xvl2=(type+1)*ones(length(inddn),1);
                     end
                        mnvls=nanmean(yvl)
                        dnmn=nanmean(yvl(inddn));
                        upmn=nanmean(yvl(indup));
                        upster=nanstd(yvl(indup))./sqrt(length(indup));
                        dnster=nanstd(yvl(inddn))./sqrt(length(inddn));
                        stervls=nanstd(yvl)./sqrt(lnvls);
                        if(COMB)
                            bar(plotxvl,mnvls,0.6,'FaceColor','none');
                            
                            hold on;
                            plot(plotxvl+.2,nanmean(yvl(inddn)),'c<')
                            plot(plotxvl+.2,nanmean(yvl(indup)),'r<')
                            plot(plotxvl+.4,nanmean(yvl(lidind)),'g<')
                            plot(plotxvl+.4,nanmean(yvl(muind)),'k<')
                             plot([plotxvl plotxvl],[mnvls-stervls mnvls+stervls],'k')
                        else
                            bar(plotxvl,upmn,0.6,'FaceColor','none');
                            plot([plotxvl plotxvl],[upmn-upster mnvls+upster],'k')
                            hold on;
%                              plot([plotxvl plotxvl],[mnvls-stervls mnvls+stervls],'k')
                             bar(plotxvl+1,dnmn,0.6,'FaceColor','none');
                             plot([plotxvl+1 plotxvl+1],[dnmn-dnster mnvls+dnster],'k')
                        end
                       
                        if(COMB)
%                         plot(xvl,yvl(yplotind),'ko','MarkerSize',3);
                        else
%                          plot(xvl1,yvl(indup),'ko','MarkerSize',3);
%                          plot(xvl2,yvl(inddn),'ko','MarkerSize',3);
                        end
                        
                        if(~ps.STIM)
                            %plot lid
                            lidind=sumdatash.lid;
                            if(COMB)
%                             plot((type)*ones(1,length(lidind)),yvl(lidind),'co','MarkerSize',3);
                           
                            else
                               liddn=intersect(lidind,inddn);
                               lidup=intersect(lidind,indup);
%                                plot(type*ones(1,length(lidup)),yvl(lidup),'co','MarkerSize',3);
                               hold on;
                               plot((type+1)*ones(1,length(liddn)),yvl(liddn),'co','MarkerSize',3);
                            end    
%                             lidmn=nanmean(yvl(sumdatash.lid));
%                             plot(1.7, lidmn,'c<');
%                             text(2,lidmn,['n=' num2str(length(sumdatash.lid))])
%                             plot(1.7, mnvls,'k<')
%                             
%                             text(2,mnvls,['n=' num2str(lnvls)]);
%                             plot(1.7,dnmn,'b<')
%                             text(2,dnmn,'down','Color','b')
%                             plot(1.7, upmn,'r<')
%                             text(2, upmn,'up','Color','r');
                        end   
               end
%              if(isfield(ps,'aspectsumrever'))
%                     axis([ps.bnds_cv]);
%                     daspect(ps.aspect_cv)
%                     box off;
%                     end
    
    function [plotvls,ctrlvls]=plotsum(combvls,xvls,ps,plotin)
       %ps.runtypetoplot
       %1 is baseline runs, 
       %2 is shiftruns
       %3 is shiftruns rev
    for pnum=1:length(plotin);
        crtype=plotin(pnum).type;
        ntvl=plotin(pnum).nt;
        crxvl=plotin(pnum).xvl;
        ps.COMB=plotin(pnum).comb;
        ps.FLIP=plotin(pnum).FLIP
        ps.NORM=plotin(pnum).norm;
        ps.drxn=plotin(pnum).drxn;
        if(plotin(pnum).CTRL)
            CTRL=1;
        else
            CTRL=0;
            ctrlvls=[];
        end
                   crvlsup=combvls{crtype}{1}{ntvl};
                   crvlsdn=combvls{crtype}{2}{ntvl};
                   
                   if(~ps.STIM)
                       
                       crvlsup.acanal=crvlsup.acpreshift;
                       crvlsdn.acanal=crvlsdn.acpreshift;
                       crvlsup.acpstanal=crvlsup.acpostshift;
                       crvlsdn.acpstanal=crvlsdn.acpostshift;
                   else
                       crvlsup.acanal=crvlsup.acshift;
                       crvlsdn.acanal=crvlsdn.acshift;
                   end
                   
                   if(~isempty(crxvl))
                       %determine net difference in crvls between ac and mu
                      
                           analvls.ac=[crvlsup.acanal crvlsdn.acanal];
                           analvls.mu=[crvlsup.mushift crvlsdn.mushift];
                           if(~ps.STIM)
                               analvls.acpst=[crvlsup.acpstanal crvlsdn.acpstanal];
                           end
                           analvls.bsnm=[crvlsup.bsnm crvlsdn.bsnm]
                           analvls.cvred=[crvlsup.cvred crvlsdn.cvred];
                           analvls.shnm=[crvlsup.shnum crvlsdn.shnum]
                           analvls.pct=[crvlsup.pct crvlsdn.pct]
                           lnup=length(crvlsup.acanal);
                           lndn=length(crvlsdn.acanal);
                           for ii=1:lnup
                           analvls.drxnind(ii)=1;
                           end
                           for ii=lnup+1:lnup+lndn
                            analvls.drxnind(ii)=2;
                           end
                           if(isfield(crvlsup,'lidflag'))
                               analvls.lidind=find(crvlsup.lidflag==1)
                               analvls.lidind=[analvls.lidind lnup+find(crvlsdn.lidflag==1)];
                           else
                                analvls.lidind=[];
                           end
                       
                           
                                [duplicate_inds]=find_duplicates(analvls);
                           
                                 [plotvls{crtype}{ntvl},ps.lidindout]=average_duplicates(analvls,duplicate_inds,ps);
                   end
          
                 if(CTRL)
                    ctrvlsup=combvls{crtype}{1}{2};
                   ctrvlsdn=combvls{crtype}{2}{2};
                   
                   if(~ps.STIM)
                       
                       ctrvlsup.acanal=ctrvlsup.acpreshift;
                       ctrvlsdn.acanal=ctrvlsdn.acpreshift;
                       ctrvlsup.acpstanal=ctrvlsup.acpostshift;
                       ctrvlsdn.acpstanal=ctrvlsdn.acpostshift;
                   else
                       ctrvlsup.acanal=ctrvlsup.acshift;
                       ctrvlsdn.acanal=ctrvlsdn.acshift;
                   end
                   
                   if(~isempty(crxvl))
                       %determine net difference in crvls between ac and mu
                      
                           ctanalvls.ac=[ctrvlsup.acanal ctrvlsdn.acanal];
                           ctanalvls.mu=[ctrvlsup.mushift ctrvlsdn.mushift];
                           if(~ps.STIM)
                               ctanalvls.acpst=[ctrvlsup.acpstanal ctrvlsdn.acpstanal];
                           end
                           ctanalvls.bsnm=[ctrvlsup.bsnm ctrvlsdn.bsnm]
                           ctanalvls.cvred=[ctrvlsup.cvred ctrvlsdn.cvred];
                           ctanalvls.shnm=[ctrvlsup.shnum ctrvlsdn.shnum]
                           ctanalvls.pct=[ctrvlsup.pct ctrvlsdn.pct]
                           lnup=length(ctrvlsup.acanal);
                           lndn=length(ctrvlsdn.acanal);
                           for ii=1:lnup
                           ctanalvls.drxnind(ii)=1;
                           end
                           for ii=lnup+1:lnup+lndn
                            ctanalvls.drxnind(ii)=2;
                           end
                           if(isfield(crvlsup,'lidflag'))
                               ctanalvls.lidind=find(ctrvlsup.lidflag==1)
                               ctanalvls.lidind=[ctanalvls.lidind lnup+find(ctrvlsdn.lidflag==1)];
                           else
                                ctnalvls.lidind=[];
                           end
                       
                           
                                [duplicate_inds]=find_duplicates(ctanalvls);
                           
                                 [ctrlvls{crtype}{ntvl},ps.lidindout]=average_duplicates(ctanalvls,duplicate_inds,ps);
                   end
                 end
                   %get notNAN bs/sh combinations.
                   [indlist]=find(~isnan(plotvls{crtype}{ntvl}.ac));
                   bslist=plotvls{crtype}{ntvl}.bs(indlist);
                   shlist=plotvls{crtype}{ntvl}.sh(indlist);
                        if(ps.STIM)
                        ps.offz=plotrawpts(crxvl,plotvls{crtype}{ntvl},ps,bslist,shlist);
                        end
                    [mnvls,stervls]=plotsumbar(crxvl,plotvls{crtype}{ntvl},ps);
                    if(isfield(ps,'aspectsumrever'))
                    axis([ps.bndssumrever]);
                    daspect(ps.aspectsumrever)
                    box off;
                    end
    end
                  
        function []=plotpct(plotvls,ps)
            figure
            sumdatash=plotvls{2}{1}
            plotxvl=1; 
            indup=find(sumdatash.drxn==1);
            inddn=find(sumdatash.drxn==2);
           if(~ps.STIM)
               lidind=sumdatash.lid;
           lnlid=length(lidind);
           end    
           yvl=sumdatash.pct 

                     lnvlsdn=length(find(~isnan(yvl(inddn))));
                     lnvlsup=length(find(~isnan(yvl(indup))));
                     xvldn=plotxvl*ones(lnvlsdn,1);   
                      xvlup=plotxvl*ones(lnvlsup,1);
                      if(~ps.STIM)
                      xlid=plotxvl*ones(lnlid,1);
                        
                      end
                      mnvls=nanmean(yvl)
                        dnmn=nanmean(yvl(inddn));
                        upmn=nanmean(yvl(indup));
                        stervls=nanstd(yvl)./sqrt(lnvlsdn+lnvlsup);
                        bar(plotxvl,mnvls*100,0.6,'FaceColor','none');
                        hold on;
                        plot([plotxvl plotxvl],[mnvls*100-stervls*100 mnvls*100+stervls*100],'k')
                        plot(xvlup,yvl(indup)*100,'ko','MarkerSize',5);
                        
                         plot(xvldn,yvl(inddn)*100,'ko','MarkerSize',5);
                         
                         
                       if(isfield(ps,'aspectsumrever'))
                        axis([ps.bnds_pct]);
                        daspect(ps.aspect_pct)
                    box off;
                    end  
                        if(~ps.STIM)
%                             plot(xlid,yvl(lidind),'co','MarkerSize',5);
                            lidmn=nanmean(yvl(sumdatash.lid));
                            muind=find(~ismember(1:length(sumdatash.sh),sumdatash.lid))
                            mumn=nanmean(yvl(muind));
                            plot(0.6, lidmn*100,'c>','MarkerSize',7);
                            text(2,lidmn,['n=' num2str(length(sumdatash.lid))])
                            plot(0.6, mumn*100,'k>','MarkerSize',7)
                            
%                             text(2,mnvls,['n=' num2str(lnvls)]);
                            plot(1.2,dnmn*100,'b<','MarkerSize',7)
                            text(2,dnmn,'down','Color','b')
                            plot(1.2, upmn*100,'r<','MarkerSize',7)
                            text(2, upmn*100,'up','Color','r');
                        else
                            text(2,mnvls,['n=' num2str(lnvlsdn+lnvlsup)]);
                            plot(1.2,dnmn*100,'b<','MarkerSize',7)
                            text(2,dnmn,'down','Color','b')
                            plot(1.2, upmn*100,'r<','MarkerSize',7)
                            text(2, upmn,'up','Color','r');
                               
                         
                       if(isfield(ps,'aspectsumrever'))
                        axis([ps.bnds_pct]);
                        daspect(ps.aspect_pct)
                    box off;
                       end 
                        end
                 
        function [sumdata]=calcstats(analvls,plotvls,ps)
             offz=plotvls.ac-plotvls.mu
            inddn=find(plotvls.drxn==2)
            indup=find(plotvls.drxn==1);  
            if(~ps.splitupdown)
            offz(indup)=-offz(indup);
            end
            sumdata.offzall=offz;
            sumdata.offzup=offz(indup);
            sumdata.pctvlsall=plotvls.pct;
            sumdata.bs=plotvls.bs;
            sumdata.sh=plotvls.sh;
            sumdata.pctvlsup=plotvls.pct(indup);
            sumdata.pctvlsdn=plotvls.pct(inddn);

        function [mnvls,stervls]=plotsumbar(xvls,plotvls,ps)
       %calc sum stats
        
            if(ps.COMB)
                dnind=find(plotvls.drxn==2)
                upind=find(plotvls.drxn==1);
                lidind=plotvls.lid;
                muind=setdiff(1:length(plotvls.drxn),lidind);
                for ii=1:length(xvls)
                    crxvl=xvls(ii)
                    if(ii==1)
                        vls=plotvls.ac;
                       
                        
                        
                        if(ps.FLIP)
                            
                            vls(dnind)=-vls(dnind);
                        end
                         vlsinit=vls;
                    elseif(ii==2)
                        vls=plotvls.mu;
                        if(ps.FLIP)
                            vls(dnind)=-vls(dnind);
                        end
                    else
                        vls=plotvls.acpst;
                        if(ps.FLIP)
                            vls(dnind)=-vls(dnind)
                        end
                    end
                    if(ps.NORM)
                        vls=vls-vlsinit
                    end
                    lnvls=length(find(~isnan(vls)));
                    mnvls=nanmean(vls)
                    stervls=nanstd(vls)./sqrt(lnvls);
                    bar(crxvl,mnvls,0.6,'FaceColor','none');
                    hold on;
                    plot(crxvl+.2, nanmean(vls(dnind)),'c<')
                   text(crxvl+.2, nanmean(vls(dnind)),'dn','Color','c')
                    plot(crxvl+.2, nanmean(vls(upind)),'r<')
                    plot(crxvl+.4,nanmean(vls(lidind)),'g<');
                    text(crxvl+.2, nanmean(vls(lidind)),'lid','Color','g')
                     plot(crxvl+.4,nanmean(vls(muind)),'k<');
                    plot([crxvl crxvl],[mnvls-stervls mnvls+stervls],'k')
                    if(ii==1)
                    text(crxvl,mnvls*1.5,['n=' num2str(lnvls)]);
                    end
                end
            else
               for ii=1:length(xvls)
                    
                        drxn=ps.drxn
                            crxvl=xvls(ii);
                          if(ii==1)
                            vls=plotvls.ac;
                          elseif(ii==2)
                            vls=plotvls.mu;
                          else
                            vls=plotvls.acpst;
                        end 
                           
                        crind=find(plotvls.drxn==drxn)  
                        lnvls=length(find(~isnan(vls(crind))));
                        
                        mnvls=nanmean(vls(crind))
                        stervls=nanstd(vls(crind))./sqrt(lnvls);
                        bar(crxvl,mnvls,0.6,'FaceColor','none');
                        hold on;
                        plot([crxvl crxvl],[mnvls-stervls mnvls+stervls],'k')
                        text(crxvl,mnvls*1.5,['n=' num2str(lnvls)]);
                        end
                   
                
            end
       
       
        
                   
        function  [offzout]=plotrawpts(crxvl,inplotvls,ps,bslist,shlist)
            %since this is sumplot...
            %plot all points, but plot different directions in different
            %color.
         for ii=1:length(crxvl)
              crx(ii)=crxvl(ii)
              crx(ii)=crxvl(ii);
         end
              plotvls=inplotvls
            
               
           
            %get matching ind
            matchind=[];
            for bsvl=1:length(bslist)
                crbs=bslist(bsvl)
                crsh=shlist(bsvl);
                outind=find(plotvls.bs==crbs&plotvls.sh==crsh);    
                if(~isempty(outind))
                    matchind=[matchind outind]
                else
                    matchind=[matchind NaN]
                end
            end
            
            offz=plotvls.mu(matchind)-plotvls.ac(matchind)
            inddn=find(plotvls.drxn(matchind)==2)
            indup=find(plotvls.drxn(matchind)==1);  
            if(~ps.COMB)
            offz(inddn)=-offz(inddn);
            end
            crxdn=crx'*ones(1,length(inddn));
            
                crxup=crx'*ones(1,length(indup));
                if(ps.NORM)
                    plotvls.ac=plotvls.ac-plotvls.ac;
                    plotvls.mu=plotvls.mu-plotvls.ac;
                   if(~ps.STIM)
                       plotvls.acpst=plotvls.acpst-plotvls.ac;
                   end
                end
                
                
                
                if(ps.STIM==0)
                plot(crxdn,[plotvls.ac(inddn);plotvls.mu(inddn);plotvls.acpst(inddn)],'Color',[0.5 0.5 0.5],'Linewidth',1);               
                hold on;
                plot(crxup,[plotvls.ac(indup);plotvls.mu(indup);plotvls.acpst(indup)],'Color',[0.5 0.5 0.5],'Linewidth',1);
                else
                    if(ps.COMB)
                        
                        plot(crxdn,[plotvls.ac(inddn);plotvls.mu(inddn)],'Color',[0.5 0.5 0.5],'Linewidth',1);               
                        hold on;
                        plot(crxup,[plotvls.ac(indup);plotvls.mu(indup)],'Color',[0.5 0.5 0.5],'Linewidth',1);
                    else
                        if(ps.drxn==2)
                            plot(crxdn,[plotvls.ac(inddn);plotvls.mu(inddn)],'Color',[0.5 0.5 0.5],'Linewidth',1);               
                        else
                            plot(crxup,[plotvls.ac(indup);plotvls.mu(indup)],'Color',[0.5 0.5 0.5],'Linewidth',1);
                        end
                    end
                end
              offzout=offz;  
          
        
    function [sumdata,sumdatarev]=plotall(combvls,xvls,ps) 

    for typevl=1:length(combvls)
        for drxnvl=1:length(combvls{typevl})
%             for ntvl=1:length(combvls{typevl}{drxnvl})
              for ntvl=1
                crvls=combvls{typevl}{drxnvl}{ntvl}
                crxvl=xvls{typevl}{drxnvl}{ntvl}
                
                if(~isempty(crxvl))
                
                
                axes(ax)
                for ii=1:length(crxvl)
                    if(ii==1)
                        if(~ps.STIM)
                            axvls{1}=crvls.acpreshift;
                     
                             lidind=find(crvls.lidflag==1)
                        else
                            axvls{1}=crvls.acshift;
                            lidind=[];
                        end
                    elseif(ii==2)
                        axvls{2}=crvls.mushift;
                    else
                        axvls{3}=crvls.acpostshift;
                    end


                    [duplicate_inds]=find_duplicates(crvls);
                    [plotvls{ii},lidindout,pctvls_ave]=average_duplicates(crvls,axvls{ii},duplicate_inds,lidind)
                    
                        crx{ii}=crxvl(ii)*ones(1,length(plotvls{ii}));
                    
                    plot(crx{ii},plotvls{ii},'o','Color',ps.col{ii},'MarkerSize',2);
                    hold on;
                    if(ii>1)
                        plot([crx{ii-1} ;crx{ii}],[plotvls{ii-1}; plotvls{ii}],'Linewidth',0.8,'Color','k');
                        
                        if(~isempty(lidindout))
                            plot([makerow(crx{ii-1}(lidindout)) ;makerow(crx{ii}(lidindout))],[plotvls{ii-1}(lidindout); plotvls{ii}(lidindout)],'Linewidth',0.8,'Color','c');
                        end
                        if (ii==2)
                            text(crx{1}(1)-1,1,['n=' num2str(length(find(~isnan(plotvls{ii}))))]);
                        end
                    end
              
                %calcmean and ster of each value for bar plots
                    
                            mnoutpre_all=nanmean(axvls{ii})
                            sterpre_all=nanstd(axvls{ii})./sqrt(length(axvls{ii}));
                            mnoutpre_ave=nanmean(plotvls{ii});
                            sterpre_ave=nanstd(plotvls{ii})./sqrt(length(axvls{ii}));
                            bar(crxvl(ii),mnoutpre_ave,0.6,'FaceColor','none');
                            plot([crxvl(ii) crxvl(ii)],[mnoutpre_ave-sterpre_ave mnoutpre_ave+sterpre_ave],'k')
                    %shift_target notes
                     if(typevl==2)
                        if(ntvl==1)
                            if(ii==2)
                            sumdata(drxnvl).mnoutpre_all=mnoutpre_all;
                            sumdata(drxnvl).sterpre_all=sterpre_all;
                            sumdata(drxnvl).mnoutpre_ave=mnoutpre_ave;
                            sumdata(drxnvl).sterpre_ave=sterpre_ave;
                            sumdata(drxnvl).pctall=crvls.pct
                            sumdata(drxnvl).pctave=pctvls_ave;
                            sumdata(drxnvl).aveplotvls=plotvls{ii};
                            sumdata(drxnvl).aven=length(find(~isnan(plotvls{ii})));
                            sumdata(drxnvl).alln=length(mnoutpre_all);
                            sumdata(drxnvl).lidvls=length(lidindout);
                            sumdata(drxnvl).pctlid=pctvls_ave(lidindout);
                            end
                            end
                     end
                    if(typevl==3)
                        if(ntvl==1)
                            if(ii==2)
                            sumdatarev(drxnvl).mnoutpre_all=mnoutpre_all;
                            sumdatarev(drxnvl).sterpre_all=sterpre_all;
                            sumdatarev(drxnvl).mnoutpre_ave=mnoutpre_ave;
                            sumdatarev(drxnvl).sterpre_ave=sterpre_ave;
                            sumdatarev(drxnvl).pctall=crvls.pct
                            sumdatarev(drxnvl).pctave=pctvls_ave;
                            sumdatarev(drxnvl).aveplotvls=plotvls{ii};
                            sumdatarev(drxnvl).aven=length(find(~isnan(plotvls{ii})));
                            sumdatarev(drxnvl).alln=length(mnoutpre_all);
                            sumdatarev(drxnvl).lidvls=length(lidindout);
                            sumdatarev(drxnvl).pctlid=pctvls_ave(lidindout);
                            end
                            end
                     end
                
                end
                
                
                
            end
        end
        end
    end







    function [xvls]=setxvls(ps)
    
    if(~ps.STIM)
        
        if(ps.TYPE=='plotsum')
    
    %baseline runs  target notes are going to 12 13 14
    xvls{1}{1}{1}=[1]
    xvls{1}{1}{2}=[3]
%     
    %shift runs, target notes
    xvls{2}{1}{1}=[5]
%     xvls{2}{2}{1}=[0 1 2]
    
    %shift runs, control notes
    xvls{2}{1}{2}=[7]
%     xvls{2}{2}{2}=[8 9 10]
    
    %baseline runs, control notes
%     xvls{1}{1}{2}=[]
%     xvls{1}{2}{2}=[]
    
    xvls{3}{1}{1}=[9]
    xvls{3}{1}{2}=[11]
    
%     xvls{3}{2}{1}=[24:26]
%     
%     xvls{3}{1}{2}=[];
%     xvls{3}{2}{2}=[];
    
    
    
    xvls{4}{1}{1}=[13]
    xvls{4}{1}{2}=[15]
    
%     xvls{4}{2}{1}=[32:34]
%     
%     xvls{4}{1}{2}=[];
%     xvls{4}{2}{2}=[];
    
        end
    else
    %shiftruns, targetnotes
    xvls{2}{1}{1}=[5]
    xvls{2}{1}{2}=[7]
    
    %baseline runs, target notes
    xvls{1}{1}{1}=[1]
    xvls{1}{1}{2}=[3]
    
    %shiftruns,rev
    xvls{3}{1}{1}=[9 ]
    xvls{3}{1}{2}=[11]
    
    
%      xvls{4}{1}{1}=[15 16]
%     xvls{4}{2}{1}=[18 19]
    
    
    end
    
    
    function[plotvls,lidindout]=average_duplicates(analvls,duplicate_inds,ps)
        %loop through each of the axvls values
        %confirm that value has not been included yet
        %check if ind is duplicate_ind
        %
        
        anal_completelist=[];
        pctvls_ave=[];
        lidindout=[];
        plotvls.ac=[];
        plotvls.mu=[];
        plotvls.acpst=[];
        plotvls.cvred=[];
        plotvls.pct=[];
        plotvls.drxn=[];
        plotvls.bs=[];
        plotvls.sh=[];
        if(~ps.STIM)
            plotvls.lid=[];
        end
        crctr=0;
        for ii=1:length(analvls.ac)
           if ~ismember(ii,anal_completelist)
               duplicateflag=0;
               for jj=1:length(duplicate_inds)
                  crduplicates=duplicate_inds{jj}
                  if(ismember(ii,crduplicates))
                     if(ps.combduplicates)

                         plotvls.ac=[plotvls.ac nanmean(analvls.ac(crduplicates))];
                         if(~ps.STIM) 
                         plotvls.acpst=[plotvls.acpst nanmean(analvls.acpst(crduplicates))];
                         end
                         plotvls.mu=[plotvls.mu nanmean(analvls.mu(crduplicates))];
                      plotvls.bs=[plotvls.bs analvls.bsnm(crduplicates(1))];
                      plotvls.cvred=[plotvls.cvred nanmean(analvls.cvred(crduplicates))];
                      plotvls.sh=[plotvls.sh analvls.shnm(crduplicates(1))];
                      plotvls.pct=[plotvls.pct nanmean(analvls.pct(crduplicates))];
                      plotvls.drxn=[plotvls.drxn nanmean(analvls.drxnind(crduplicates))];
                      duplicateflag=1;
                      anal_completelist=[anal_completelist crduplicates];
                      crctr=crctr+1;
                                      
                     
                      if(ismember(ii,analvls.lidind))
                          plotvls.lid=[plotvls.lid crctr];
                      end
                     end
                  end
               end
               if(duplicateflag==0)
                      plotvls.ac=[plotvls.ac analvls.ac(ii)];
                      plotvls.mu=[plotvls.mu analvls.mu(ii)];
                      plotvls.pct=[plotvls.pct analvls.pct(ii)];
                       plotvls.bs=[plotvls.bs analvls.bsnm(ii)];
                       plotvls.cvred=[plotvls.cvred analvls.cvred(ii)];
                      plotvls.sh=[plotvls.sh analvls.shnm(ii)];
                    if(~ps.STIM)
                      plotvls.acpst=[plotvls.acpst analvls.acpst(ii)];
                    end 
                      plotvls.drxn=[plotvls.drxn analvls.drxnind(ii)];
                      anal_completelist=[anal_completelist ii];
                      crctr=crctr+1;
                      if(ismember(ii,analvls.lidind))
                          plotvls.lid=[plotvls.lid crctr];
                      end
               end
           end
            
            
            
        end
    tst=1;
    
    
    function [duplicate_inds]=find_duplicates(analvls)
        redinds=[];
        bsind=unique(analvls.bsnm);
        shind=unique(analvls.shnm);
        ct=1;
        for ii=1:length(bsind)
            crbs=bsind(ii);
            for jj=1:length(shind);
                crsh=shind(jj);
                bsinds=find(analvls.bsnm==crbs)
                shinds=find(analvls.shnm==crsh);
                inds=intersect(bsinds,shinds);
                
                if(~isempty(inds))
                    if(length(inds)>=2)
                        duplicate_inds{ct}=inds;
                        ct=ct+1;
                    end
                end
            
            end
        end
        if(~exist('duplicate_inds'))
            duplicate_inds=[];
        end
    
    
    %this function takes the maximum value if there are multiple runs on
    %the same shift.
    function [redinds]=reduce_vls(crvls)
        redinds=[];
        bsind=unique(crvls.bsnm);
        shind=unique(crvls.shnum);
        for ii=1:length(bsind)
            crbs=bsind(ii);
            for jj=1:length(shind);
                crsh=shind(jj);
                bsinds=find(crvls.bsnm==crbs)
                shinds=find(crvls.shnum==crsh);
                inds=intersect(bsinds,shinds);
                
                if(~isempty(inds))
                    redinds=[redinds max(inds)];
                end
            
            end
        end
     function [stats]=calcsumstats(plotvls)   
         bas=plotvls{1}{1};
         shift=plotvls{2}{1};
         basdn=find(bas.drxn==2);
         basup=find(bas.drxn==1);
         shiftdn=find(shift.drxn==2);
         shiftup=find(shift.drxn==1);
         
         bascomb.ac=[bas.ac(basdn) bas.ac(basup)];
         bascomb.mu=[bas.mu(basdn) bas.mu(basup)];
         
         shdn.ac=[shift.ac(shiftdn)];
         shdn.mu=[shift.mu(shiftdn)];
        
         
         shup.ac=[shift.ac(shiftup) ];
         shup.mu=[shift.mu(shiftup) ];
         
         cvdn=[shift.cvred(shiftdn)];
         cvup=[shift.cvred(shiftup)];
         cvbas=[bas.cvred(basdn) bas.cvred(basup)];
         
         
         
         
         [stats.hbas,stats.pbas]=ttest(bascomb.ac,bascomb.mu);
         [stats.hshdn,stats.pshdn]=ttest(shdn.ac,shdn.mu);
         [stats.hshup,stats.pshup]=ttest(shup.ac,shup.mu);
         
           [stats.hcvbas,stats.pcvbas]=ttest(cvbas,1,.05,'left');
         [stats.hcvdn,stats.pcvdn]=ttest(cvdn,1,.05,'left');
         [stats.hcvup,stats.pcvup]=ttest(cvup,1,.05,'left');
         [stats.hcvcomb,stats.pcvcomb]=ttest([cvbas cvdn cvup],1,.05,'left');
         stats.mncvcomb=mean([cvbas cvdn cvup])
         
         
         
         
