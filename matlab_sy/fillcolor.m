function [aclist,mulist]=fillcolor(ps)
    for ii=1:length(ps)
        fill(ps(ii).fillpts(:,1),ps(ii).fillpts(:,2),ps(ii).fillcol);
        hold on;
    end