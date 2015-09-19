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

function [sumdata,sumdatarev]=plotgroupbar(combvls,ps)   
    ax=ps.ax;
    ps.col={'k' 'r' 'k'}
%     edges=ps.edges;
%     col=ps.col;
% %     jjpair{1}=[1 2]
%     jjpair{2}=[1 2]
    
    %XVLS
    %baseline runs are going to x-7, 8, 9
    
    if(~ps.STIM)
    
    %baseline runs  target notes are going to 12 13 14
    xvls{1}{1}{1}=[12 13 14]
    xvls{1}{2}{1}=[16 17 18]
    
    %shift runs, target notes
    xvls{2}{1}{1}=[0 1 2]
    xvls{2}{2}{1}=[0 1 2]
    
    %shift runs, control notes
    xvls{2}{1}{2}=[4 5 6]
    xvls{2}{2}{2}=[8 9 10]
    
    %baseline runs, control notes
    xvls{1}{1}{2}=[]
    xvls{1}{2}{2}=[]
    
    xvls{3}{1}{1}=[20:22]
    xvls{3}{2}{1}=[24:26]
    
    xvls{3}{1}{2}=[];
    xvls{3}{2}{2}=[];
    
    
    
    xvls{4}{1}{1}=[28:30]
    xvls{4}{2}{1}=[32:34]
    
    xvls{4}{1}{2}=[];
    xvls{4}{2}{2}=[];
    
    
    else
    %shiftruns, targetnotes
    xvls{2}{1}{1}=[0 1]
    xvls{2}{2}{1}=[0 1]
    
    %baseline runs, target notes
    xvls{1}{1}{1}=[3 4]
    xvls{1}{2}{1}=[6 7]
    
    %shiftruns,rev
    xvls{3}{1}{1}=[9 10]
    xvls{3}{2}{1}=[12 13]
    
    
     xvls{4}{1}{1}=[15 16]
    xvls{4}{2}{1}=[18 19]
    
    
    end
    
    

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

    function[plotvls,lidindout,pctvls_ave]=average_duplicates(crvls,axvls,duplicate_inds,lidind)
        %loop through each of the axvls values
        %confirm that value has not been included yet
        %check if ind is duplicate_ind
        %
        
        anal_completelist=[];
        pctvls_ave=[];
        lidindout=[];
        plotvls=[];
        crctr=0;
        for ii=1:length(axvls)
           if ~ismember(ii,anal_completelist)
               duplicateflag=0;
               for jj=1:length(duplicate_inds)
                  crduplicates=duplicate_inds{jj}
                  if(ismember(ii,crduplicates))
                      plotvls=[plotvls nanmean(axvls(crduplicates))];
                      pctvls_ave=[pctvls_ave nanmean(crvls.pct(crduplicates))];
                      ind=find(isnan(axvls(crduplicates)))
                      if(~isempty(ind))
                          tst=1;
                      end
                      duplicateflag=1;
                      anal_completelist=[anal_completelist crduplicates];
                      crctr=crctr+1;
                      if(ismember(ii,lidind))
                          lidindout=[lidindout crctr];
                      end
                  end
               end
               if(duplicateflag==0)
                      plotvls=[plotvls axvls(ii)];
                      if(isempty(axvls(ii)))
                          tst=1;
                      end
                      pctvls_ave=[pctvls_ave crvls.pct(ii)];
                      anal_completelist=[anal_completelist ii];
                      crctr=crctr+1;
                      if(ismember(ii,lidind))
                          lidindout=[lidindout crctr];
                      end
               end
           end
            
            
            
        end
    
    
    
    function [duplicate_inds]=find_duplicates(crvls)
        redinds=[];
        bsind=unique(crvls.bsnm);
        shind=unique(crvls.shnum);
        ct=1;
        for ii=1:length(bsind)
            crbs=bsind(ii);
            for jj=1:length(shind);
                crsh=shind(jj);
                bsinds=find(crvls.bsnm==crbs)
                shinds=find(crvls.shnum==crsh);
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
        
