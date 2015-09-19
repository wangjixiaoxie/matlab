function []=plotacmuratio(avls,ratvls,ntind, zratio)
    ax=figure
    hold on;
    for ii=1:length(avls.muanal)
        muind=avls.muanal(ii);
        %for each of these points, plot a vertical line going from
        %(acmean-mean) to (mumean-acmean)
        
        acmean=ratvls(muind).acmean{ntind};
        initmean=avls.initmean{ntind};
        mumean=ratvls(muind).mumean{ntind};
        xvl=acmean-initmean;
        yvl=mumean-acmean;
        cvfac=ratvls(muind).cvfac;
        if(exist('zratio'))
            xvl=(xvl/avls.initsd{ntind});
            yvl=(yvl/avls.initsd{ntind});
        end
            plotvertline([xvl 0],[xvl yvl],'-');
            plot(xvl,yvl,'ko','MarkerSize',cvfac*5);
            plotline([-500 500] ,[500 -500],'--','k')
            plot([0 0], [-500 500],'k','Linewidth',2)
            plot([-500 500], [0 0],'k','Linewidth',2)
    end
    if(~exist('zratio'))
        xlabel('Current Pitchshift (Hz)','Fontsize',16)
        ylabel('Acute Inactivation Effect (Hz)', 'Fontsize', 16)
    else
        xlabel('Current Pitchshift (zscore)','Fontsize',16)
        ylabel('Acute Inactivation Effect (zscore)', 'Fontsize', 16)
    end
    
    