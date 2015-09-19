%hs struct returned has mnvl and marked.
%modified 4.14.09 in order to allow for avls, with multiple fields.
function [hs]=gethstcomb3(hs,avlsin, avls_fields);



for ii=1:length(hs)
  if(exist('fields'))
      avls=avlsin(hs(ii).avlsind);
  end
    
    if(hs(ii).acvec==1)
        hs(ii).hst=avls.hstacpre{hs(ii).ntvec}{hs(ii).muvec}
    elseif(hs(ii).acvec==2)
        hs(ii).hst=avls.hstmu_comb{hs(ii).ntvec}{hs(ii).muvec}
    else
        hs(ii).hst=avls.hstacpst{hs(ii).ntvec}{hs(ii).muvec}
    end
    hs(ii).edges=avls.edges{hs(ii).ntvec};
end
    
