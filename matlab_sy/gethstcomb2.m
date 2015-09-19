%hs struct returned has mnvl and marked.
function [hs]=gethstcomb2(hs,avls);



for ii=1:length(hs)
  
    
    if(hs(ii).acvec==1)
        hs(ii).hst=avls.hstacpre{hs(ii).ntvec}{hs(ii).muvec}
    elseif(hs(ii).acvec==2)
        hs(ii).hst=avls.hstmu_comb{hs(ii).ntvec}{hs(ii).muvec}
    else
        hs(ii).hst=avls.hstacpst{hs(ii).ntvec}{hs(ii).muvec}
    end
    hs(ii).edges=avls.edges{hs(ii).ntvec};
end
    
