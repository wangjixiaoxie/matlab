function [hs]=getedges(hs,muind,avls);

for ii=1:3
    hs(ii).edges=avls.edges{hs(ii).ntvl}
end