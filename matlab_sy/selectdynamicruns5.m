%purpose of this code is to go through sumbs and calculate
%the offz for each initrun
%for each control run
%and for each baseline runs
%save to a masterstruct suminitshs
%with initacz, initmuz, initoff, init_bsvl, 
%modified 9.19.09 to return sumdynall and sumdynasymp

%this option tells whether to plot initind, selected as consecutive mu runs
%over some threshold.
function [sumdyn,sumdynasymp,sumdynrev]=selectdynamicruns2(sumbs,excludebs,minpts,minshift,norm_asymp,TYPE)
bslist=[];
for ii=1:length(sumbs)
    if(ii~=excludebs)
       bslist=[bslist ii]; 
    end
end
for runvl=1:3
    ct=0;
    clear sumout
    for bsnm=bslist
        crbs=sumbs(bsnm)
        
        if(TYPE=='phar')
            trgnote=crbs.allnote(1);
            if(length(crbs.allnote)>1)
                ctrlnote=crbs.allnote(2);
            else
                ctrlnote=trgnote;
            end
        end
    if(runvl==1)
        runs=crbs.shiftruns;
    elseif(runvl==2)
        runs=crbs.allasympruns;
    %for revruns, add one before as well
    else
       runs=crbs.revruns;
      
    end
  
        for shnm=1:length(runs)
            shiftind=runs{shnm}
            if(TYPE=='stim')
                tmpind=find(ismember(shiftind,crbs.STANRUNS));
                shiftind=shiftind(tmpind);
            end
            if(runvl==3)  
                if(~isempty(shiftind))
                nonzeroind=find(crbs.mulist~=0);
                mulistnonzero=crbs.mulist(nonzeroind);
                mulistind=find(ismember(mulistnonzero,shiftind));
                shiftind=[mulistnonzero(mulistind(1)-1) mulistnonzero(mulistind)];
                 overallind=shiftind;
                else
                    overallind=[];
                end
            end
            
            if(runvl==2)  
%                 mulistind=find(crbs.mulist==shiftind(1));
                shiftind=[crbs.asympruns{shnm}];
                 overallind=shiftind;
            end
            if(runvl==1)
                 cr_shiftruns=crbs.shiftruns{shnm};
                basruns=crbs.basruns;
                nonzeroind=find(crbs.mulist~=0);
                mulistnonzero=crbs.mulist(nonzeroind);
                shind=find(ismember(mulistnonzero,cr_shiftruns'));
%             shindloc=find(nonzeroind==shind(1));
                basind=find(ismember(mulistnonzero(shind(1)-1),basruns))
                if(~isempty(basind))
               
                    shiftind=nonzeroind([shind(1)-1 shind]);
                      
                else
                       shiftind=nonzeroind(shind);
                end
                overallind=crbs.mulist(shiftind);
            end
            
            
%             ind=find(crbs.mun(trgnote,shiftind)>minpts&abs(crbs.acz(trgnote,shiftind))>minshift);
     
%                 if(~isempty(ind))
                ct=ct+1;
               
              
            
%                 if(exist('norm_asymp'))
%                     sumdynasymp(ct).tms=crbs.asympshifttms{shnm}(ind)-crbs.asympshifttms{shnm}(ind(1))+1;
%                 else
%                     sumdynasymp(ct).tms=crbs.asympshifttms{shnm}(ind)-crbs.adjshifttms{shnm}(ind(1))+1;
%                 end
                
                sumout(ct).pct=crbs.pct(overallind);
                sumout(ct).off=crbs.offz(trgnote,overallind);
                sumout(ct).drxn=crbs.drxn{shnm}
                if(crbs.drxn{shnm}=='up')
                        sumout(ct).normoff=-crbs.offz(trgnote,overallind);
                else
                        sumout(ct).normoff=crbs.offz(trgnote,overallind);
                end
%                     
%                      sumout(ct).targeff=crbs.combeff(trgnote,overallind);
% %                   sumout(ct).logtargeff=log2(sumdyn(ct).targeff);
%                     sumout(ct).contreff=crbs.combeff(ctrlnote,overallind);
%                 
                
                sumout(ct).asympvl=crbs.asympvl{shnm};
                sumout(ct).shnum=shnm;
                
                
                sumout(ct).acz=crbs.acprez(trgnote,overallind);
               
                sumout(ct).muz=crbs.muz(trgnote,overallind);
                if(runvl==3)
                    shsize=sumout(ct).asympvl-sumout(ct).acz;
                    sumout(ct).motsize=sumout(ct).asympvl-sumout(ct).muz;
                    sumout(ct).lmansize=sumout(ct).muz-sumout(ct).acz;
                    sumout(ct).ac_pct=sumout(ct).motsize./shsize;
                    sumout(ct).mu_pct=sumout(ct).lmansize./shsize;
                end
%                 sumout(ct).logcontreff=log2(sumdyn(ct).contreff);
                sumout(ct).bsnum=bsnm
                if(runvl==1)
                    sumout(ct).exadjtms=crbs.flrtmvec(overallind)-crbs.flrtmvec(crbs.shiftruns{shnm}(1))+crbs.adjshifttms{shnm}(1);
                    ind=find(sumout(ct).exadjtms<0);
                    sumout(ct).exadjtms(ind)=0;
%                     sumout(ct).muasympdist=crbs.muasympdist{shnm};
%                     sumout(ct).acasympdist=crbs.acasympdist{shnm};
                    
%                     sumout(ct).logcontreff=log2(sumdyn(ct).contreff);
                    sumout(ct).bsnum=bsnm
                elseif(runvl==3)
                    if(~isempty(overallind))
                        sumout(ct).exadjtms=crbs.flrtmvec(overallind,1)-crbs.flrtmvec(overallind(1),1)+1;
                    else
                        sumout(ct).exadjtms=[];
                    end
                elseif(runvl==2)
                    
                    sumout(ct).asympmotpct=sumout(ct).muz/sumout(ct).asympvl;
                    if(~isempty(crbs.asympshifttms{shnm}))
                    sumout(ct).exadjtms=crbs.flrtmvec(overallind)-crbs.flrtmvec(crbs.allasympruns{shnm}(1))+1;
                    end
                end
%             end
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

        

