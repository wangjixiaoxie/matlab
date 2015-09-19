%hs struct returned has mnvl and marked.
function [hs]=gethsts(hs,muind,avls);


for ii=1:3
        if(ii<3)
            ind=avls.aclist(muind,ii);
        else
            ind=avls.mulist(muind);
        end
        hs(ii).hst=avls.hst{hs(ii).ntvl}{ind}
        mnvl=avls.mnvl{hs(ii).ntvl}(muind);
        stdv=avls.stdv{hs(ii).ntvl}(muind);
        hs(ii).plotline=1;
        hs(ii).marked=[mnvl-stdv mnvl+stdv]; 
        hs(ii).markx=[mnvl];    
        hs(ii).mark='o';
end
    
