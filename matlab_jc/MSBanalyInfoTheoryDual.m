% FRONT & REAR
    for i=1:8
        if isequal(DShifts(i).dirA,'up')
            valu=51;
        else
            valu=49;
        end
        [ITfirst(i) ITsecond(i)]=ContingSimMSBDual(DShifts(i).toffset-DShifts(i).onset,jc_residuals(DShifts(i).pitchBaseline(DShifts(i).onset:DShifts(i).onset+400,:)),valu);
        Actualfirst(i)=round(median(DShifts(i).toffset-DShifts(i).onset));
        Actualsecond(i)=round(Actualfirst(i)+192);
    end


