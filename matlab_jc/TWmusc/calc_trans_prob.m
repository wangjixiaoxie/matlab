%written by tw to find indices of a given sequence transition.
%currently written to evaluate three cases: 1. a specified transition from
%a note, 2. ending on a note, and 3. all other occurrences.
function [transdata]=calc_trans_prob(fvdata,targetnote,postnote);
    
    transdata=struct('endnotes',[],'match',[],'nomatch',[]);

    for ii=1:length(fvdata)
        lblstr=fvdata(ii).lbl;
        %see if the bird ended on the target note
        if (fvdata(ii).ind)==length(lblstr)
            transdata.endnotes=[transdata.endnotes ii];
        %targetnote-postnote
        elseif(lblstr(fvdata(ii).ind+1)==postnote)
            transdata.match=[transdata.match ii];
        %all other transitions
        else
            transdata.nomatch=[transdata.nomatch ii];
            
        end
    end
    %transdata=transdata2;