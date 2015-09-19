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
            
           dataout=filtfilt(ones(1,win)/win,1,datin);
           datout=dataout(ceil(win/2):(end-ceil(win/2)));
          
        