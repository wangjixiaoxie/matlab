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

%NEED TO MODIFY TO USE PREINDS...
%%%which correspond to acvalues and muvalues
%from combvls{5}, check that shift is greater than 1sd!

function [plotvls]=plotbarpre(sumbs,combvls,ps)   
    ps.TYPE='plotsum'
%     figure
%     ps.axsum=gca();
%     ax=ps.ax;
    ps.col={'k' 'r' 'k'}
    [xvls]=setxvls(ps);
    if(ps.TYPE=='plotall')
        [sumdata,sumdatarev]=plotall(combvls,xvls,ps);
    elseif(ps.TYPE=='plotsum')
        axes(ps.axsum)
        [plotvls]=plotsum(sumbs,combvls,xvls,ps);
%         axes(ps.axpct)
%         plotpct(sumdatash(1));
    end
    
    function [plotvls]=plotsum(sumbs,combvls,xvls,ps)
       %only look at runtype 4, asymp
       crtype=2;
       ntvl=1;
                shvlsup=combvls{2}{1}{1};
                shvlsdn=combvls{2}{2}{1};
       
           crtype=4;
           ntvl=1;
           
          
                crvlsup=combvls{4}{1}{ntvl};
                crvlsdn=combvls{4}{2}{ntvl};
                
%                 allvlsup=combvls{5}{1}{1};
%                 allvlsdn=combvls{5}{2}{1};
                
                    crxvls=xvls{crtype}{1}{ntvl}
                   
                   crvlsup.acanal=crvlsup.acpreshift;
                   crvlsdn.acanal=crvlsdn.acpreshift;
                  
                  
                   
                   if(~isempty(crxvls))
                       %determine net difference in crvls between ac and mu
                      
                           analvls.ac=[crvlsup.acanal -crvlsdn.acanal];
                           analvls.mu=[crvlsup.mushift -crvlsdn.mushift];
                           analvls.bsnm=[crvlsup.bsnm crvlsdn.bsnm]
                           analvls.shnm=[crvlsup.shnum crvlsdn.shnum]
                           analvls.pct=[crvlsup.pct crvlsdn.pct]
                           analvls.cvred=[crvlsup.cvred crvlsdn.cvred]
                           analvls.asymptms=[crvlsup.asymptms crvlsdn.asymptms];
                           allvls.adjtimes=[shvlsup.adjtimes shvlsdn.adjtimes];
                           allvls.ac=[shvlsup.acshift -shvlsdn.acshift];
                           allvls.pct=[shvlsup.pct shvlsdn.pct];
                           allvls.cvred=[shvlsup.cvred shvlsdn.cvred];
                           allvls.mu=[shvlsup.mushift -shvlsdn.mushift];
                           allvls.bsnm=[shvlsup.bsnm shvlsdn.bsnm];
                           allvls.shnm=[shvlsup.shnum shvlsdn.shnum];
                           analvls.pretms=[crvlsup.pretms crvlsdn.pretms];
                           analvls.preinds=[crvlsup.preinds crvlsdn.preinds];
                           analvls.preindbs=[crvlsup.preindbs crvlsdn.preindbs];
                           analvls.preindsh=[crvlsup.preindsh crvlsdn.preindsh];
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
                       
                           
                                [duplicate_inds]=find_duplicates(allvls);
                           
                                 [plotvls]=calc_duplicates(sumbs,analvls,allvls,duplicate_inds,ps);
                         
                    %THEN MAKE PLOT OF THIS PARTICULAR VALUE, WRITE VALUE
                    %OUT TO OUTVALUE   
                   %start with plotting in zscore units.  %then I can
                   %switch to other values.
                    
%                     ps.offz=plotrawz(crxvl,plotvls,ps);
                   
                    [mnvls,stervls]=plotsumbar(crxvls,plotvls,ps,allvls);
                    
%                     if(crtype==2)
%                         [sumdatash(ntvl)]=calcstats(analvls,plotvls,ps);
%                     elseif(crtype==3)
%                         [sumdatarev(ntvl)]=calcstats(analvls,plotvls,ps);
%                     elseif(crtype==1)
%                         [sumdatabas(ntvl)]=calcstats(analvls,plotvls,ps);
%                         end
                   
                   end
        
       
        function []=plotpct(sumdatash)
            crxvl=1;  
            for ii=1:2
                     if(ii==2)
                        plotxvl=crxvl
                        yvl=sumdatash.pctvlsup
                        
                     else
                        plotxvl=crxvl+1;
                        yvl=sumdatash.pctvlsdn;
                     end
                     
                     lnvls=length(find(~isnan(yvl)));
                     xvl=plotxvl*ones(lnvls,1);   
                        mnvls=nanmean(yvl)
                        stervls=nanstd(yvl)./sqrt(lnvls);
                        bar(plotxvl,mnvls,0.6,'FaceColor','none');
                        hold on;
                        plot([plotxvl plotxvl],[mnvls-stervls mnvls+stervls],'k')
                        plot(xvl,yvl,'k.','MarkerSize',1);
                        text(plotxvl,mnvls*1.5,['n=' num2str(lnvls)]);
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
            sumdata.pctvlsup=plotvls.pct(indup);
            sumdata.pctvlsdn=plotvls.pct(inddn);

        function [mnvls,stervls]=plotsumbar(xvl,plotvls,ps,allvls)
       %calc sum stats
            axes(ps.axsum)
            if(ps.MKPLOT)
                for ii=1:length(plotvls.pct(1,:))
                    crxvl=xvl(ii);
                    crvls=plotvls.pct(:,ii);
                    lnvls=length(crvls);
                    mnvls=nanmean(crvls)
                    stervls=nanstd(crvls)./sqrt(lnvls);
                 
                    bar(crxvl,mnvls,0.6,'FaceColor','none');
                    hold on;
                    plot([crxvl crxvl],[mnvls-stervls mnvls+stervls],'k')
                    text(crxvl,mnvls*1.5,['n=' num2str(lnvls)]);
                end
           
            axes(ps.axz)
            
                   
               for ii=1:length(plotvls.pct(1,:))
                    crxvl=xvl(ii);
                    crvls=plotvls.cvred(:,ii);
                    lnvls=length(crvls);
                    mnvls=nanmean(crvls-1)
                    stervls=nanstd(crvls-1)./sqrt(lnvls);
                 
                    bar(crxvl,mnvls,0.6,'FaceColor','none');
                    hold on;
                    plot([crxvl crxvl],[mnvls-stervls mnvls+stervls],'k')
                    text(crxvl,mnvls*1.5,['n=' num2str(lnvls)]);
                end
            else
                mnvls=[];
                stervls=[];
            end
                %addpre
%                 lnvls=length(find(~isnan(plotvls.meanprecvred)));
%                 mnpre=nanmean(plotvls.meanprecvred-1);
%                 sterpre=nanstd(plotvls.meanprecvred-1)./sqrt(lnvls);
%                 bar(0,mnpre,0.6,'FaceColor','none');
%                 plot([0 0],[mnpre+sterpre mnpre-sterpre],'k');
%                 text(0,mnpre*1.5,['n=' num2str(lnvls)]);
%                 for ii=1:length(plotvls.ac(:,1))
%                     crxvl=plotvls.avex(ii,:);
%                     crvlsac=plotvls.ac(ii,:);
%                     crvlsmu=plotvls.mu(ii,:);
%                     if(plotvls.drxn(ii)==2)
%                         crvlsac=-crvlsac;
%                         crvlsmu=-crvlsmu;
%                     end
%                     
%                     plot(crxvl, crvlsac,'ko');
%                     hold on;
%                     plot(crxvl, crvlsmu,'ro');
%                     plot(crxvl,crvlsac,'k');
%                     plot(crxvl,crvlsmu,'r');
%             
%                 end
%              
               
                
                   
       
       
        
                   
        function  [offz]=plotrawpts(crxvl,plotvls,ps)
            %since this is sumplot...
            %plot all points, but plot different directions in different
            %color.
            
            if(~ps.splitupdown)
                crxvlup=crxvl
                
            else
                crxvlup=crxvl+1;
            end
                
            offz=plotvls.mu-plotvls.ac
            inddn=find(plotvls.drxn==2)
            indup=find(plotvls.drxn==1);  
            if(~ps.splitupdown)
            offz(inddn)=-offz(inddn);
            end
            crxdn=crxvl*ones(1,length(inddn));
            
                crxup=crxvlup*ones(1,length(indup));
                
           
                plot(crxdn,offz(inddn),'.','Color','k','MarkerSize',ps.marksize);
                crlidind=find((ismember(indup,ps.lidindout)))
                plot(crxup(crlidind),offz(indup(crlidind)),'.','Color','c','MarkerSize',ps.marksize);
                
                hold on;
                plot(crxup,offz(indup),'.','Color','k','MarkerSize',ps.marksize);
                crlidind=find((ismember(inddn,ps.lidindout)))
                plot(crxdn(crlidind),offz(inddn(crlidind)),'.','Color','c','MarkerSize',ps.marksize);
                
            
        
    
    function [xvls]=setxvls(ps)
    
    
        
        if(ps.TYPE=='plotsum')
    
    %baseline runs  target notes are going to 12 13 14
    xvls{4}{1}{1}=[1 2 3]
%     xvls{4}{2}{1}=[3]
% %     
%     %shift runs, target notes
%     xvls{2}{1}{1}=[5]
% %     xvls{2}{2}{1}=[0 1 2]
%     
%     %shift runs, control notes
%     xvls{2}{1}{2}=[7]
%     xvls{2}{2}{2}=[8 9 10]
    
    %baseline runs, control notes
%     xvls{1}{1}{2}=[]
%     xvls{1}{2}{2}=[]
%     
%     xvls{3}{1}{1}=[9]
%     xvls{3}{1}{2}=[11]
%     
% %     xvls{3}{2}{1}=[24:26]
% %     
% %     xvls{3}{1}{2}=[];
% %     xvls{3}{2}{2}=[];
%     
%     
%     
%     xvls{4}{1}{1}=[13]
%     xvls{4}{1}{2}=[15]
%     
% %     xvls{4}{2}{1}=[32:34]
% %     
%     xvls{4}{1}{2}=[];
%     xvls{4}{2}{2}=[];
%     
%         end
%     else
% %     %shiftruns, targetnotes
%     xvls{2}{1}{1}=[5]
%     xvls{2}{1}{2}=[7]
%     
%     %baseline runs, target notes
%     xvls{1}{1}{1}=[1]
%     xvls{1}{1}{2}=[3]
%     
%     %shiftruns,rev
%     xvls{3}{1}{1}=[9 ]
%     xvls{3}{1}{2}=[11]
%     
%     
%      xvls{4}{1}{1}=[15 16]
%     xvls{4}{2}{1}=[18 19]
    
    
    end
    
    
    function[plotvls,lidindout]=calc_duplicates(sumbs,analvls,allvls,duplicate_inds,ps)
        %THIS FUNCTION NEEDS TO BE CHANGED..
        %algorithm:
%         1. loop through duplcate inds.
%         2. check analvls.tms versus existing time inds.
%         3. if matchinds.
%             4. add to plotvls in appropriate bins  (DONE)
        plotvls.pct=[];
        plotvls.mu=[];
        plotvls.ac=[];
        
        plotvls.avex=[];
        plotvls.drxn=[];
        plotvls.bsnm=[];
        plotvls.cvred=[];
        plotvls.meanprecvred=[];
        plotvls.shnm=[];
        plotvls.meanpre=[];
        plotvls.asymp=[];
        plotvls.shiftpct=[];
        plotvls.mupct=[];
        fixct=0;  
        for ii=1:length(duplicate_inds)  
               crinds=duplicate_inds{ii}    
               adjtms=allvls.adjtimes(crinds);
               fixct=fixct+1; 
               for crtmind=1:length(ps.tmwins)
                   
                   crtms=ps.tmwins{crtmind}
                    [vls,matchinds]=intersect(adjtms,crtms);
                    
                    if(~isempty(matchinds))
                        inds=crinds(matchinds);
%                          [meanpre(fixct),meanprecvred(fixct)]=getmeanpre(allvls,analvls,inds);
                        meanpct(fixct,crtmind)=mean(calcmeanstder2(allvls.pct(inds)));
                        meancvred(fixct,crtmind)=mean(calcmeanstder2(allvls.cvred(inds)));
                        meanacz(fixct,crtmind)=mean(calcmeanstder2(allvls.ac(inds)));
                        meanx(fixct,crtmind)=mean(calcmeanstder2(allvls.adjtimes(inds)));
%                         meanptvl(fixct,crtmind)=mean(calcmeanstder2(ptvl(matchinds)));
                        meanmuz(fixct,crtmind)=mean(calcmeanstder2(allvls.mu(inds)));
%                         meaninactvl(fixct,crtmwin)=mean(calcmeanstder2(inactvl(matchinds))); 
            
                    else
                       
                        meanpct(fixct,crtmind)=NaN;
                        meanacz(fixct,crtmind)=NaN;
                        meanmuz(fixct,crtmind)=NaN;
                       
                    end
               end
                   
                    ps.matchreq=1:length(ps.tmwins)
                    notnaind=find(~isnan(meanpct(fixct,:)));

                    nomatchvls=find(~ismember(ps.matchreq, notnaind)); 
                    mnvls=abs(meanacz(fixct,:));
                    if(isempty(nomatchvls)&isempty(find(mnvls<1)))
                          plotvls.pct=[plotvls.pct;meanpct(fixct,:)];
%                           plotvls.meanpre=[plotvls.meanpre;meanpre(fixct)];
%                           plotvls.meanprecvred=[plotvls.meanprecvred;meanprecvred(fixct)];
                          plotvls.cvred=[plotvls.cvred;meancvred(fixct,:)];
                          plotvls.ac=[plotvls.ac;meanacz(fixct,:)];
                          plotvls.mu=[plotvls.mu;meanmuz(fixct,:)];
                          plotvls.avex=[plotvls.avex;meanx(fixct,:)];
%                           plotvls.drxn=[plotvls.drxn allvls.drxnind(inds(1))]; 
                          plotvls.bsnm=[plotvls.bsnm allvls.bsnm(inds(1))];
                          plotvls.shnm=[plotvls.shnm allvls.shnm(inds(1))];
                           plotvls.asymp=getasymp(sumbs,allvls.bsnm(inds(1)),allvls.shnm(inds(1)));
                           if(~isempty(plotvls.asymp))
                                plotvls.shiftpct=[plotvls.shiftpct;meanacz(fixct,:)./plotvls.asymp];
                                plotvls.mupct=[plotvls.mupct;meanmuz(fixct,:)./plotvls.asymp]
                           else
                               plotvls.shiftpct=[plotvls.shiftpct;meanacz(fixct,:)];
                               plotvls.mupct=[plotvls.mupct;meanmuz(fixct,:)];
                           
                           end
            
               
                    end
        end
    
    
        function [meanpct,meancvred]=getmeanpre(allvls,analvls,inds)
            %this tells us which bs/sh we are working with.
            bsnm=analvls.bsnm(inds(1));
            shnm=analvls.shnm(inds(1));
            %which of the PREINDS have this bsnm and shnum
           newinds=find(analvls.preindbs==bsnm&analvls.preindsh==shnm);
           if(~isempty(newinds))
               pretms=analvls.pretms(newinds);
               preinds=analvls.preinds(newinds);
               allinds=find(allvls.bsnm==bsnm&allvls.ind==shnm);
               
               acvls=allvls.ac(allinds(preinds));
               cv_vls=allvls.cvred(allinds(preinds));
               muvls=allvls.mu(allinds(preinds));
               shiftinds=find((abs(acvls)>1)&pretms>-3);
               if(~isempty(shiftinds))
                    meanpct=mean(muvls(shiftinds)./acvls(shiftinds));
                    meancvred=mean(cv_vls(shiftinds));
               else
                   meanpct=NaN;
                   meancvred=NaN;
               end
           else
               meanpct=NaN;
               meancvred=NaN;
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
       function [asympvl]=getasymp(sumbs,bsnm,shnm)
            asympvl=sumbs(bsnm).asympvl{shnm};
       
