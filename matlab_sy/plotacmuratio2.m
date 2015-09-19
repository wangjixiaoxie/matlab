function [ax]=plotacmuratio(birdstruct,ratvls,zratio)
ax=figure    
for birdind=1:length(birdstruct)
    ntind=birdstruct(birdind).ntind;    
    cmd=['cd ' birdstruct(birdind).path 'datasum']
        eval(cmd);
        cmd=['load ' birdstruct(birdind).matfilename]
        eval(cmd);
        
        hold on;
        for ii=1:length(avls.muanal)
            muind=avls.muanal(ii);
            %for each of these points, plot a vertical line going from
            %(acmean-mean) to (mumean-acmean)
        
            acmean=ratvls(birdind).acmean{ntind}{muind};
            initmean=avls.initmean{ntind};
            mumean=ratvls(birdind).mumean{ntind}{muind};
            xvl=acmean-initmean;
            yvl=mumean-acmean;
            cvfac=ratvls(birdind).cvfac{ntind}{muind};
            if(exist('zratio'))
                xvl=(xvl/avls.initsd{ntind});
                yvl=(yvl/avls.initsd{ntind});
            end
            plot([xvl xvl],[0 yvl],'Linestyle','-','Color',birdstruct(birdind).color,'Linewidth',2);
            plot(xvl,yvl,'o','MarkerSize',cvfac*5,'Color',birdstruct(birdind).color);
            
        end
    end
        plotline([-500 500] ,[500 -500],'--','k')
        plot([0 0], [-500 500],'k','Linewidth',1)
        plot([-500 500], [0 0],'k','Linewidth',1)

        if(~exist('zratio'))
        xlabel('Current Pitchshift (Hz)','Fontsize',16)
        ylabel('Acute Inactivation Effect (Hz)', 'Fontsize', 16)
    else
        xlabel('Current Pitchshift (zscore)','Fontsize',16)
        ylabel('Acute Inactivation Effect (zscore)', 'Fontsize', 16)
    end
    
    