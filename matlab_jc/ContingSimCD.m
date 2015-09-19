function aavfin=ContingSimCD(toffsets,pitch,dirA,dirB)
toffsets=[400];
between=300; % 24ms - 3 windows

count=0;
    for i=1:length(toffsets)
        offset=round(toffsets(i)); %round(mean(toffsets)+0.8*std(toffsets)*randn);
        ffsA=(pitch(offset,:)); % estimate of pitch within the window
        ffsB=(pitch(offset+between,:));
        LA=prctile(ffsA,80);
        LB=prctile(ffsB,80);
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
        indC=[];
        for j=1:size(pitch,2)
            if sum(indA==j)*sum(indB==j)>0;
                indC=[indC j];
            end
        end
        if ~isempty(indC)
            count=count+1;
            if length(indC)==1
                aavTot(count,:)=(pitch(:,indC));
            else
                aavTot(count,:)=mean(pitch(:,indC)');
            end
        end
        atot(i).data=pitch(:,indC);
    end
aavfin=aavTot;