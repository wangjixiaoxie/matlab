%purpose of this code is to go through sumbs and calculate
%the offz for each initrun
%for each control run
%and for each baseline runs
%save to a masterstruct suminitshs
%with initacz, initmuz, initoff, init_bsvl, 
%modified 9.19.09 to return sumdynall and sumdynasymp

%this option tells whether to plot initind, selected as consecutive mu runs
%over some threshold.
function [sumdyn,sumdynasymp,sumdynrev]=selectdynamicruns2(sumbs,minpts,minshift,norm_asymp)
for runvl=1:3
    ct=0;
    clear sumout
for bsnm=1:length(sumbs)
    crbs=sumbs(bsnm)
%     trgnote=crbs.allnote(1);
%     if(length(crbs.allnote)>1)
%         ctrlnote=crbs.allnote(2);
%     else
%         ctrlnote=trgnote;
%     end
    if(runvl==1)
        runs=crbs.shiftruns;
        
    elseif(runvl==2)
        runs=crbs.asympruns;
    %for revruns, add one before as well
    else
       runs=crbs.revruns;
      
    end
   
        
       
        for shnm=1:length(runs)
            shiftind=runs{shnm}
            tmpind=find(ismember(shiftind,crbs.STANRUNS));
            shiftind=shiftind(tmpind);
            if(runvl==3)
                if(~isempty(shiftind))
                    [out,sortind]=sort(crbs.tmvec(crbs.STANRUNS));
                    sortruns=crbs.STANRUNS(sortind);
                    mulistind=find(sortruns==shiftind(1));
                    shiftind=[sortruns(mulistind-1); shiftind];
                end
                

            end
%             ind=find(crbs.mun(trgnote,shiftind)>minpts&abs(crbs.acz(trgnote,shiftind))>minshift);
%             if(~isempty(ind))
                ct=ct+1;
                overallind=shiftind;
                ind=1:length(shiftind);
                if(~isempty(overallind))
%                 if(exist('norm_asymp'))
%                     sumdynasymp(ct).tms=crbs.asympshifttms{shnm}(ind)-crbs.asympshifttms{shnm}(ind(1))+1;
%                 else
%                     sumdynasymp(ct).tms=crbs.asympshifttms{shnm}(ind)-crbs.adjshifttms{shnm}(ind(1))+1;
%                 end
%                 sumout(ct).pct=crbs.pct(overallind);    
%                 sumout(ct).off=crbs.offz(trgnote,overallind);
%                 if(crbs.drxn{shnm}=='up')
%                     sumout(ct).normoff=-crbs.offz(trgnote,overallind);
%                 else
%                     sumout(ct).normoff=crbs.offz(trgnote,overallind);
%                 end
                    sumout(ct).acz=crbs.acz(overallind);
                    sumout(ct).muz=crbs.muz(overallind);
                    sumout(ct).runs=overallind;
                    sumout(ct).rawtms=crbs.tmvec(overallind,1)-crbs.baswntime{shnm};
                    
                     sumout(ct).asympvl=crbs.asympvl{shnm};
                     sumout(ct).asympacpct=sumout(ct).acz/sumout(ct).asympvl
                    sumout(ct).asympmupct=sumout(ct).muz/sumout(ct).asympvl;
                    %                 sumout(ct).targeff=crbs.combeff(trgnote,overallind);
% %                 sumout(ct).logtargeff=log2(sumdyn(ct).targeff);
%                 sumout(ct).contreff=crbs.combeff(ctrlnote,overallind);
%                 sumout(ct).logcontreff=log2(sumdyn(ct).contreff);
                    sumout(ct).bsnum=bsnm
                    if(runvl==1)
                        sumout(ct).adjtms=crbs.adjshifttms{shnm}(tmpind);
                        sumout(ct).exadjtms=crbs.exactshifttms{shnm}(tmpind);
                        sumout(ct).muasympdist=crbs.muasympdist{shnm};
                        sumout(ct).acasympdist=crbs.acasympdist{shnm};
%                     sumout(ct).targeff=crbs.combeff(trgnote,overallind);
% %                     sumout(ct).logtargeff=log2(sumdyn(ct).targeff);
%                     sumout(ct).contreff=crbs.combeff(ctrlnote,overallind);
%                     sumout(ct).logcontreff=log2(sumdyn(ct).contreff);
                        sumout(ct).bsnum=bsnm
                    elseif(runvl==3)
                        if(~isempty(overallind))
                            sumout(ct).adjtms=crbs.tmvec(overallind)-crbs.flrtmvec(overallind(1))+1;
                            sumout(ct).adjtms=sumout(ct).adjtms-sumout(ct).adjtms(1)+1;
                           
                        end
                    else
                        sumout(ct).adjtms=crbs.asympshifttms{shnm}(tmpind);
                        sumout(ct).exadjtms=crbs.exactasymptms{shnm}(tmpind);
               
                    end
%             end
                end
        end
        
end
    if(runvl==2)
            sumdynasymp=sumout;
    elseif(runvl==1)
            sumdyn=sumout;
    else
        sumdynrev=sumout;
    end
end

        

