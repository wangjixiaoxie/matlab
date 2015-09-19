%this program called by inactivplotmeta, calls plotline, also calls
%getstyle to get the style for plot.
%used to make inactivation figures.



function [ax] = inactivplotmu(avls,psin)

    ax=psin.ax;
%psin has the following fields
% psin.ht is the list of heights for each of the indtoplot
% psin.indtoplot is the list of muind to plot(following original btind convention
% psin.ntind is the nt to plot
%psin.ax

%acmean is a list of acmeans for indtoplot
%mumean is a series of mumeans

acz=avls.acz(psin.ntind,psin.indtoplot);
muz=avls.muz(psin.ntind,psin.indtoplot);
[psout,msout]=getstyle(psin.style)
%loop through and plot the lines

% psout.ax=ax;

for ii=1:length(acz)
    %just set the pts and the markerpts.
    psout(1).pts=[acz(ii) psin.ht(ii) muz(ii) psin.ht(ii)]
    psout(2).pts=[acz(ii) psin.ht(ii) 0 psin.ht(ii)]
    msout(1).markpts(1)=acz(ii);
    msout(1).markpts(2)=psin.ht(ii);
    msout(2).markpts(1)=muz(ii);
    msout(2).markpts(2)=psin.ht(ii);
    msdm=[];
    lineplot(psout(2),msdm)
   if((abs(muz(ii))-abs(acz(ii)))>0)
        psout(1).col='r';
    else
        psout(1).col='k';
   end
    lineplot(psout(1),msout);
    if(psin.vec)
        [pts]=calcvec(avls,psin.ntind, psin.indtoplot,psin.vecht);
        ln=abs(pts(3)-pts(1));
        plot_arrow(pts(1), pts(2) ,pts(3) ,pts(4),'Linewidth',3,'Headwidth',.1/ln,'Headheight',0.1/ln);
    end
end

