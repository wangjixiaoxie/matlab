%inactivplot5.m
%this program designed to plot for one specific set of runs, as called by
%inactivplotmeta.



function [ax] = inactivplotmu(avls,ps)

%
%acmean is a list of acmeans for indtoplot
%mumean is a series of mumeans

acmeans=avls.





runind=[];
if (isempty(ps.ax))
   figure
   ax=gcf;
else
    ax=ps.ax;
end

if(isempty(ps.hist))
    hist=0;
else
    hist=ps.hist;
end
if(isempty(ps.ht))
    ht=[0 1];
else
    ht=ps.ht;
end

dyind=ps.dyind;
notevl=ps.notevl;

    

if(hist)
    %first create acsf part of struct
    for ii=1:length(hist)
        hstvl=hist(ii)
        %acsf
        if(hstvl<3)
            ind=avls.aclist(dyind,hstvl);
            hs(hstvl).col='k';
            
        else
            ind=avls.mulist(dyind)
            hs(hstvl).col='r';
            
        end
        
        hs(hstvl).edges=avls.edges{notevl};
        hs(hstvl).hst=avls.hst{notevl}{ind};
        
        if(hstvl==2)
            hs(hstvl).ls='--'
            hs(hstvl).wid=1;
        else
            hs(hstvl).ls='-'
            hs(hstvl).wid=2;
        end
    end
    plothists(hs,ax);
end

%plot_arrow(ax, arrowstruct), plot arbitrary number of arrows.
    %for reverse direction one blackrow goes from maxmean to
    %the current acmean
    %one red arrow above goes from last reported pitch to acute target
    %pitch.
     indvl=avls.mulist(dyind);
    %that is for every one of the up-down runs. 
    %start off by assuming this is a reverse run.
    %treat each up/down run as unique.
    if(ps.plt=='all')
        muanal=avls.muanal;
    elseif(ps.plt=='rev')
        muanal=avls.revanal;
    elseif(ps.plt=='asymp'
        muanal=avl
    end
    for runbloc=1:length(muanal{notevl})
        runind=find(muanal{notevl}{runbloc}==indvl);
            if~isempty(runind)
                runblocvl=runbloc;
            end
    end
        %get the maxval.
        if(ps.mu=='re')
            mxvl=avls.maxmeans(runblocvl);
        end
        acvl=avls.acmean(notevl,indvl);
        muvl=avls.mumean(notevl,indvl);
        mnvl=avls.initmean{notevl} 
        if(ps.zrat)
                 if(ps.mu=='re')   
                    mxvl=(mxvl-mnvl)/avls.initsd{notevl};
                 end
                    acvl=(acvl-mnvl)/avls.initsd{notevl};
                    muvl=(muvl-mnvl)/avls.initsd{notevl};
                    mnvl=0;
        end
        if(ps.mu=='mu')
            x1=[mnvl acvl];
            x2=[acvl muvl];
        else
            x1=[mxvl acvl] 
            x2=[acvl muvl]
        end
            y=ht
        col(1)='k';
        col(2)='r';
        
            for ii=1:2
                axis(axis);
                if(ps.line)
                    plot([x1(ii) x2(ii)],[y(ii) y(ii)],'Color',col(ii),'Linewidth',4-ii); 
                end
                hold on;
                plot([x2(ii) x2(ii)], [y(ii)-ps.tickht y(ii)+ps.tickht], 'Color', col(ii)','Linewidth',1);
            end
        end
