function aavfin=ContingSim(toffsets,pitch,cont)



%for k=1:100
    for i=1:size(pitch,2)
        rr=round(rand*length(toffsets));
        if rr==0
            rr=1;
        end
        offset(k)=round(toffsets(rr)); %round(mean(toffsets)+0.8*std(toffsets)*randn);
        ffs(i)=median(pitch(offset(k)-64:offset(k),i)); % estimate of pitch within the window
    end
    L=prctile(ffs,cont);
    if cont>50
        ind=find(ffs>L);
    else
        ind=find(ffs<L);
    end

    aavTot(k,:)=mean(pitch(:,ind)');
%end
aavfin=mean(aavTot);