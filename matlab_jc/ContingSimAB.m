function aavfin=ContingSimAB2(toffsets,pitch,dirA,dirB)

between=192; % 24ms - 3 windows

count=0;
    for i=3:length(toffsets)
        offset=round(toffsets(i)); %round(mean(toffsets)+0.8*std(toffsets)*randn);
        ffsA=(pitch(offset,:)); % estimate of pitch within the window
        ffsB=(pitch(offset+between,:));
        LA=prctile(ffsA,50);
        LB=prctile(ffsB,50);
        if isequal(dirA,'up')
            indA=find(ffsA>LA);
        else
            indA=find(ffsA<LA);
        end
        if isequal(dirB,'up')
            indB=find(ffsB>LB);
        else
            indB=find(ffsB<LB);
        end
        
        indC=find(indB==indA);
        if ~isempty(indC)
            count=count+1;
            if length(indC)==1
                aavTot(count,:)=(pitch(:,indC));
            else
                aavTot(count,:)=mean(pitch(:,indC)');
            end
        end
    end
aavfin=aavTot;