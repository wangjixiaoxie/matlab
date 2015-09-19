function [longer,middle,allcontinga]=Fig3ContingSimulator(Alldata,targregion,conting,ACSF,mm)
% targregion --- 1=beginning, 2=middle, 3=end
% conting --- 20 is 20 above and below
% ACSF --- 1 if test ACSFs, 0 if test INAs
% above =1 if percentile tested is above 50 (e.g. 80th percentile)
%%% This loop takes in a structure and selects the curves that meet the contingency
    middle=[];
    allcontinga=[];
    for k=1:length(Alldata)
        longer(k)=Alldata(k).basevaroffset-Alldata(k).basevaronset; % length of segment
        if targregion==1;  targpoint=1; end % beginning
        if targregion==2;                   % middle
            targpoint=round((longer(k))/2);
            middle(k)=targpoint;
        end
        if targregion==3; targpoint=(longer(k)); end   % end
        
        % Determines the LMAN presence condition to be tested
        if ACSF==1; pitch=Alldata(k).baselineAC(Alldata(k).basevaronset:Alldata(k).basevaroffset,:);
        else if ACSF==0; pitch=Alldata(k).baselineINA(Alldata(k).basevaronset:Alldata(k).basevaroffset,:);
            else if ACSF==2; pitch=Alldata(k).pitchUDpre(Alldata(k).basevaronset:Alldata(k).basevaroffset,:);
        else if ACSF==3; pitch=Alldata(k).pitchDpre(Alldata(k).basevaronset:Alldata(k).basevaroffset,:);
            else if ACSF==4; pitch=Alldata(k).pitchUDpost(Alldata(k).basevaronset:Alldata(k).basevaroffset,:);
                end
            end
                end
            end
        end
        mnbase=mean(pitch');
        clear pointvalues
        % Find the values at the contingency point
        for j=1:size(pitch,2) % for each note
            pointvalues(j)=pitch(targpoint,j);
        end
        % Determine the contingency based on the percentile of these values
            contingency=prctile(pointvalues,conting);
        % Choose the pitch curves that clear this contingency
            count=1;
            selectedcurves=[];
        for i=1:size(pitch,2)
            if mm==2
                if pitch(targpoint,i)<=contingency
                    selectedcurves(:,count)=pitch(:,i);
                    count=count+1;
                end
            else
                if pitch(targpoint,i)>contingency
                    selectedcurves(:,count)=pitch(:,i);
                    count=count+1;
                end
            end
        end
        % Normalize these curves and put them in a matrix
            allcontinga(k,1:length(mnbase))=(mean(selectedcurves')-mnbase);
            allcontinga(k,1:length(mnbase))=allcontinga(k,1:length(mnbase))/allcontinga(k,targpoint);
    end
