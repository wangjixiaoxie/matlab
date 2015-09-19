%purpose of this code is to go through sumbs and calculate
%the offz for each initrun
%for each control run
%and for each baseline runs
%save to a masterstruct suminitshs
%with initacz, initmuz, initoff, init_bsvl, 


%this option tells whether to plot initind, selected as consecutive mu runs
%over some threshold.
function [sumdyn]=selectdynamicruns2(sumbs,minpts,minshift,norm_asymp)

ct=0;
for bsnm=1:length(sumbs)
    crbs=sumbs(bsnm)
    trgnote=crbs.allnote(1);
    if(length(crbs.allnote)>1)
        ctrlnote=crbs.allnote(2);
    else
        ctrlnote=trgnote;
    end
    for shnm=1:length(crbs.asympruns)
        
        shftind=crbs.asympruns{shnm};
        ind=find(crbs.mun(trgnote,shftind)>minpts&abs(crbs.acz(trgnote,shftind))>minshift);
        if(~isempty(ind))
        ct=ct+1;
        overallind=shftind(ind);
        if(exist('norm_asymp'))
            sumdyn(ct).tms=crbs.asympshifttms{shnm}(ind)-crbs.asympshifttms{shnm}(ind(1))+1;
        else
            sumdyn(ct).tms=crbs.asympshifttms{shnm}(ind)-crbs.adjshifttms{shnm}(ind(1))+1;
        end
            sumdyn(ct).pct=crbs.pct(overallind);    
        sumdyn(ct).off=crbs.offz(trgnote,overallind);
        if(crbs.drxn{shnm}=='up')
            sumdyn(ct).normoff=-crbs.offz(trgnote,overallind);
        else
            sumdyn(ct).normoff=crbs.offz(trgnote,overallind);
        end
        sumdyn(ct).acz=crbs.acz(trgnote,overallind);
        sumdyn(ct).muasympdist=crbs.muasympdist{shnm};
        sumdyn(ct).acasympdist=crbs.acasympdist{shnm};
        sumdyn(ct).targeff=crbs.combeff(trgnote,overallind);
        sumdyn(ct).logtargeff=log2(sumdyn(ct).targeff);
        sumdyn(ct).contreff=crbs.combeff(ctrlnote,overallind);
        sumdyn(ct).logcontreff=log2(sumdyn(ct).contreff);
        sumdyn(ct).bsnum=bsnm
        end
    end
end

        

