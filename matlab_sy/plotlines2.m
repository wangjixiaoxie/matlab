function []=plotlines2(lns)
for ii=1:length(lns)
    axes(lns(ii).ax);
    x=lns(ii).xpts;
    y=lns(ii).ypts;
    ls=lns(ii).ls;
    lw=lns(ii).lw;
    col=lns(ii).col;
    plot(x, y,'Color',col,'LineStyle',ls, 'Linewidth',lw);
    if(lns(ii).plotmarker)
        mrkpts=lns(ii).markloc
        plot([mrkpts(1) mrkpts(2)],[mrkpts(3) mrkpts(4)],'Color', col,'Linewidth',2*lw);
    end
end