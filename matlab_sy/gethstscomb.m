%hs struct returned has mnvl and marked.
function [hs]=gethstcomb(hs,avls);
muvec=hs.muvec;
acmuvec=hs.acmuvec
ntvec=hs.ntvec
for ii=1:length(muvec)
    muind=muvec(ii); 
    if(ac(ii))
        hs(ii).hst=avls.hstac_comb{ntvec(ii)}{muind}
    else
        hs(ii).hst=avls.hstmu_comb{ntvec(ii)}{muind}
    end
    hs(ii).edges=avls.edges{ntvec(ii)};
end
    
