figure
for ii=1:4
    %plot initup
    
    if (ii==1)
        xvals=ii*ones(length(initup),1);
        plot(ii,initup,'ko');
        mninitup=mean(initup);
        sterinitup=std(initup)/sqrt(length(initup));
        hold on;
        bar(ii,mninitup,'FaceColor','none');
          plot([ii ii], [mninitup-sterinitup mninitup+sterinitup],'k');
    end
    %plot initdown
    if (ii==2)
        xvals=ii*ones(length(initdown),1);
        plot(ii,initdown,'ko');
        mninitdown=mean(initdown);
        sterinitdown=std(initdown)/sqrt(length(initdown));
        bar(ii,mninitdown,'FaceColor','none');
         plot([ii ii], [mninitdown-sterinitdown mninitdown+sterinitdown],'k');
    end
    %plot bas
    if (ii==3)
        xvals=ii*ones(length(bas.off),1);
        plot(ii,bas.off,'ko');
        mnbasoff=mean(bas.off);
        sterbasoff=std(bas.off)/sqrt(length(bas.off));
        bar(ii,mnbasoff,'FaceColor','none');
        
        plot([ii ii], [mnbasoff-sterbasoff mnbasoff+sterbasoff],'k');
    end
    %plot control
    if(ii==4)
        xvals=ii*ones(length(ctrl.off),1);
        plot(ii,ctrl.off,'ko');
        mnctrloff=mean(ctrl.off);
        sterctrloff=std(ctrl.off)/sqrt(length(ctrl.off));
        bar(ii,mnctrloff,'FaceColor','none');
        plot([ii ii], [mnctrloff-sterctrloff mnctrloff+sterctrloff],'k');
    end
end