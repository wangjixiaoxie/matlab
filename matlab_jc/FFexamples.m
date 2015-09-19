load /cardinal/b4o59/Data0124.mat
% /cardinal/LearningExample.eps
    [a,b]=hist(valsPre(:,2),[2200:10:2500]);
    [c,d]=hist(valsWN123(end-100:end,2),[2200:10:2500]);
    figure;hold on;
    stairs(a/sum(a),b,'b','Linewidth',2)
    stairs(c/sum(c),d,'r','Linewidth',2)
    plot([0 0.3],[2350 2350],'k')
    plot(0.25,mean(valsPre(:,2)),'b.')
    plot(0.25,mean(valsWN123(end-100:end,2)),'r.')
    ylim([2220 2440])
    xlim([0 0.3])
% /cardinal/b4o59Example.eps
    cd /cardinal/b4o59/0119_wnoff
    [avY,t,f]=get_avn('batchY','y',0.2,3,'','','obs0');
    %
    % bring in WN from another b4o59 folder
    avY(:,905:930)=avK(:,110:135);
    %
    figure
    imagesc(t,f,log(avY));syn;ylim([0,1e4]);
    hold on;plot([0.4 0.5],[1000 1000],'g','Linewidth',4)
    xlim([-0.05 2.98])
    % colormap 'hot' min=7 max=15
