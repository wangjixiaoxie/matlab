function aavfin=ContingSim2B(toffsets,pitch,cont,dir)



count=0;
    for i=1:length(toffsets)
        offset=round(toffsets(i)); %round(mean(toffsets)+0.8*std(toffsets)*randn);
        ffs=(pitch(offset,:)); % estimate of pitch within the window
        L=prctile(ffs,cont);
        if isequal(dir,'up')
            ind=find(ffs>L);
        else
            ind=find(ffs<L);
        end
        if ~isempty(ind)
            count=count+1;
            aavTot(count,:)=mean(pitch(:,ind)');
        end
    end
aavfin=aavTot;