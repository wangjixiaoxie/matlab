
%this function used to timestamp the first and last file in a batch file.
function [starttimes endtimes] = batch2time(pathvl,catchvl)

    for ii=1:length(pathvl)
        if(~isempty(pathvl{ii})) 
            strcmd=['cd ' pathvl{ii} ]
            eval(strcmd);
            ff=load_batchf(catchvl{ii});
        
%now go through and find the hyphens

            for jj=[1 length(ff)]
                ind=find(ff(jj).name=='_');
                strdate{ii}=ff(jj).name(ind(1)+1:ind(1)+6)
                strtime{ii}=ff(jj).name(ind(2)+1:ind(2)+4)
                strcom=[strdate{ii} strtime{ii}]
                if jj==1
                    timesvec{ii}=datevec(strcom,'ddmmyyHHMM');
                    starttimes(ii)=datenum(timesvec{ii});
                else
                    timesvec{ii}=datevec(strcom,'ddmmyyHHMM');
                    endtimes(ii)=datenum(timesvec{ii});
                end
        end
        end
    end
