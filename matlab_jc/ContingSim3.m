function aavfin=ContingSim3(toffsets,pitch,cont)

% filter residuals

% make predictions as before
indTrack=[];
count=0;
toffsets=toffsets;
    for i=1:length(toffsets)
        offset=round(toffsets(i)); %round(mean(toffsets)+0.8*std(toffsets)*randn);
        ffs=(pitch(offset,:)); % estimate of pitch within the window
        L=prctile(ffs,cont);
        if cont>50
            ind=find(ffs>L);
        else
            ind=find(ffs<L);
        end
        if ~isempty(ind)
            count=count+1;
            aavTot(count,:)=mean(pitch(:,ind)');
        end
        indTrack=[indTrack ind];
    end
aavfin=aavTot;