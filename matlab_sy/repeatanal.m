%how to quantify number of repeats

%for each of the first a's which I've picked out, 
%pick out the fv.lbl string and look for how many as till not an a.

%algorithm can be simply to create a truncated string with the rest of the
%characters and look for index for not an a

%take fvcomb 

function [fvout]=repeatanal(fv)
    fvout=fv
    for ii=1:length(fv)
        NT=fv(ii).lbl(fv(ii).ind)
        startind=fv(ii).ind;
        rptstr=fv(ii).lbl(fv(ii).ind+1:end);
        rptind=find(rptstr~=NT)
        %what about the diabolical case that bird ends with repeat, here if
        %rptind is 0, then just call the size of the rptstr
        if (isempty(rptind))
            fvout(ii).nrpeat=length(rptstr);
        else
            fvout(ii).nrpeat=rptind(1);
        end

    end