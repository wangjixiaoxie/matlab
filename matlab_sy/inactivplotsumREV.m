%written 9.16.08.  
%This will call inactivlplotmu in an orderly way across birds
%variables which are set IN THIS SCRIPT
%are ps.ht (an array of heights)
% ps.indtoplot
%ps.ntind
%ps.style
%pss.ax
ps.ax=figure;
ps.style='mu';
%determines whether or not to plot summary vector.
ps.vec=1;
bascolor=[0.31 1 0.25];
concolor=[0.29 .45 1]
%every run gets offset in height
%have a separate negative counter and a posititive counter. 
%     uphtcnt=.15;
%     downhtcnt=0;
    zrat=1;
 for birdind=1:length(bs)      
    ntind=bs(birdind).ntind;
    ps.ntind=ntind;
    cmd=['cd ' bs(birdind).path 'datasum']
    eval(cmd);
    cmd=['load ' bs(birdind).matfilename]
        eval(cmd);
    for blocind=1:length(avls.revanal{ntind})
            %set the ind to plot
                ps.indtoplot=avls.revanal{ntind}{blocind}   
                ind=ps.indtoplot
                if(ind>2)
                    ind=ind(1:2)
                end
                %need to specifgy THE ROW OF AVLS.ACZ
                if(~isempty(ps.indtoplot))
        subplot(211)
                    scatter(avls.acz(ntind,ind),avls.muz(ntind,ind),15,'k','filled')
                    
                    hold on;
                    for kk=1:length(ind)
                        test=2
                    end
                        
                        
                    
                    plot([avls.acz(ntind,ind)-avls.acerracz(ntind,ind) ;avls.acz(ntind,ind)+avls.acerracz(ntind,ind)],[avls.muz(ntind,ind); avls.muz(ntind,ind)],'k','Linewidth',1);
                    hold on;
                    plot([avls.acz(ntind,ind);avls.acz(ntind,ind)],[avls.muz(ntind,ind)-avls.mustderrz(ntind,ind) ;avls.muz(ntind,ind)+avls.mustderrz(ntind,ind)], 'k','Linewidth',1);
                    %calculate average vector.
                    mnacz=mean(avls.acz(ntind,ind));
                    mnmuz=mean(avls.muz(ntind,ind));
                    subplot(212)
                    plot([mnacz mnacz],[mnacz mnmuz],'k-','Linewidth',2);
                    hold on;
                end
    end

    for blocind=1:length(avls.basind{ntind})
            %set the ind to plot
                ps.indtoplot=avls.basind{ntind}{blocind}   
                ind=ps.indtoplot
                
                %need to specifgy THE ROW OF AVLS.ACZ
                if(~isempty(ps.indtoplot))
   subplot(211)    
                 scatter(avls.acz(ntind,ind),avls.muz(ntind,ind),15,bascolor,'filled')
                    hold on;
                    
                        
                    subplot(211)
                    plot([avls.acz(ntind,ind)-avls.acerracz(ntind,ind) ;avls.acz(ntind,ind)+avls.acerracz(ntind,ind)],[avls.muz(ntind,ind); avls.muz(ntind,ind)],'Color',bascolor,'Linewidth',1);
                    plot([avls.acz(ntind,ind);avls.acz(ntind,ind)],[avls.muz(ntind,ind)-avls.mustderrz(ntind,ind) ;avls.muz(ntind,ind)+avls.mustderrz(ntind,ind)], 'Color',bascolor,'Linewidth',1);
                    %calculate average vector.
                    mnacz=mean(avls.acz(ntind,ind));
                    mnmuz=mean(avls.muz(ntind,ind));
                    subplot(212)
                    plot([mnacz mnacz],[mnacz mnmuz],'Color',bascolor,'Linewidth',2);
                    hold on;
                end
    end
    
    
    if ~isempty(bs(birdind).contrind)
        clind=bs(birdind).contrind
        for blocind=1:length(avls.revanal{clind})
            in=avls.revanal{clind}{blocind};
              if(~isempty(in))
                    subplot(211)
                    scatter(avls.acz(clind,in),avls.muz(clind,in),15,concolor,'filled')
                    hold on;
                    test=1;
                    subplot(211)
                    plot([avls.acz(clind,in)-avls.acerracz(clind,in) ;avls.acz(clind,in)+avls.acerracz(clind,in)],[avls.muz(clind,in); avls.muz(clind,in)],'Color',concolor,'Linewidth',1);
                    plot([avls.acz(clind,in);avls.acz(clind,in)],[avls.muz(clind,in)-avls.mustderrz(clind,in) ;avls.muz(clind,in)+avls.mustderrz(clind,in)], 'Color',concolor,'Linewidth',1 );
                    %calculate average vector.
                    mnacz=mean(avls.acz(clind,in));
                    mnmuz=mean(avls.muz(clind,in));
                    subplot(212)
                    plot([mnacz ;mnacz],[mnacz ;mnmuz],'Color',concolor,'Linewidth',2);
                    hold on;
              end
        end
    end
    
            
        

 end
 
    
 
% text(-5,htcnt,avls.bname)
 subplot(211)
plot([-10 10],[-10 10],'k--');
subplot(212)
plot([-10 10],[-10 10],'k--');

xlabel('pitch z-score');
box off;