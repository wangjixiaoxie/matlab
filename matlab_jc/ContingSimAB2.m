function CSpred=ContingSimAB2(toffsets,pitch,dirA,dirB,onset,offset)

between=192; % 24ms - 3 windows

count=0;
indALL=[];

for i=1:length(toffsets)
        offset2=round((toffsets(i))); %round(mean(toffsets)+0.8*std(toffsets)*randn);
        ffsA=(pitch(offset2,:)); % estimate of pitch within the window
        ffsB=(pitch(offset2+between,:));
        LA=prctile(ffsA,40);
        LB=prctile(ffsB,60);

        if isequal(dirA,'up')
            indA=find(ffsA>LA);
        else
            indA=find(ffsA<LB);
        end
        if isequal(dirB,'up')
            indB=find(ffsB(indA)>LA);
        else
            indB=find(ffsB(indA)<LB);
        end
        indC=indA(indB);
        indALL=[indALL indC];
        

%         if ~isempty(indC)
%             
%             count=count+1;
%             if length(indC)==1
%                 aavTot(count,:)=(pitch(:,indC));
%             else
%                 aavTot(count,:)=mean(pitch(:,indC)');
%             end
%         end
end
% aavfin=aavTot;

            for j=1:max(indALL)
                nums(j)=length(find(indALL==j));
            end
CSpred=nums*pitch(:,1:length(nums))';
%CSpred=CSpreda-mean(CSpreda(onset:offset));
