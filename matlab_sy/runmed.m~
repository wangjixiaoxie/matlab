%assumes data is already separated into continuous sample
function [datout]=runave2(datin,win, REMOVE_OUTLIERS)
   
            medvl=median(datin);
            stdvl=std(datin);
            if(exist('REMOVE_OUTLIERS'))
            if(REMOVE_OUTLIERS)
                removeind=find((outvls{jj}.pt<(medvl-outlier_factor*stdvl)))
                removeind2=find((outvls{jj}.pt>(medvl+outlier_factor*stdvl)));
                origind=1:length(outvls{jj}.pt)
                [cmbremoveind]=[removeind; removeind2]
                finalind=setdiff(origind,cmbremoveind);
            end
            end

    startind=ceil(win/2);
    endind=length(datin)-ceil(win/2);
    inds=startind:endind;
    for ii=1:length(inds)
        crind=inds(ii);
        datout(ii)=median(datin(crind)-floor(win/2):datin(crind)+floor(win/2)-1);
    end