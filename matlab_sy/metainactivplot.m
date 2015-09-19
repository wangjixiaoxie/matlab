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
%every run gets offset in height
%have a separate negative counter and a posititive counter. 
    uphtcnt=.15;
    downhtcnt=0;
    zrat=1;
 for birdind=1:length(bs)      
    
    ntind=bs(birdind).ntind;
    ps.ntind=ntind;
    cmd=['cd ' bs(birdind).path 'datasum']
    eval(cmd);
    cmd=['load ' bs(birdind).matfilename]
        eval(cmd);
        
    for blocind=1:length(avls.asympind{ntind})
            %set the ind to plot
                
                ps.indtoplot=avls.asympind{ntind}{blocind}   
                if(~isempty(ps.indtoplot))
                htvls=[.33:.33:length(ps.indtoplot)/3;];
                %set the heights
                
                %downshift
                if(avls.acmaxz(ntind,blocind)<0)
                    ps.ht=downhtcnt+htvls;
                    ps.vecht=downhtcnt+htvls(end)+.5;
                    downhtcnt=downhtcnt+htvls(end)+1.5;
                else

                    ps.ht=uphtcnt+htvls;;
                    ps.vecht=uphtcnt+htvls(end)+.5;
                    uphtcnt=uphtcnt+htvls(end)+1.5;
                end
                inactivplotmu(avls,ps);  
            end
    end
   
 end
% text(-5,htcnt,avls.bname)
 
plot([0 0],[0 50],'k--');
xlabel('pitch z-score');
box off;