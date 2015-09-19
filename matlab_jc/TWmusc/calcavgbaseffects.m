
function [shiftplot,sumdata]=calcavgbaseffects(combvls,shiftplot,bsgrp)
    for bsnm=1:length(bsgrp)
        bsvls=bsgrp{bsnm};
        offzcomb=[];
        for ii=1:2
            crind=find(ismember(combvls{1}{ii}.bsnm,bsvls));
            offzcomb=[offzcomb combvls{1}{ii}.offz(crind);]
            
        end
        for kk=1:length(bsvls)
                crvl=bsvls(kk);
                shiftplot(crvl).offzcomb=offzcomb
        end
    end