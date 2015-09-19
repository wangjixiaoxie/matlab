function CSpred=ContingSimAB3(toffsets,pitch,dirA,dirB,onset,offset)

between=192; % 24ms - 3 windows

count=0;
indALL=[];

for i=1:length(toffsets)
        offset2=round((toffsets(i))); %round(mean(toffsets)+0.8*std(toffsets)*randn);
        ffsA=(pitch(offset2,:)); % estimate of pitch within the window
        ffsB=(pitch(offset2+between,:));
        LA1=prctile(ffsA,60);
        LA2=prctile(ffsA,40);
        LB1=prctile(ffsB,60);
        LB2=prctile(ffsB,40);

        if isequal(dirA,'up')
            indA=find(ffsA>LA2);
        else
            indA=find(ffsA<LA1);
        end
        if isequal(dirB,'up')
            indB=find(ffsB(indA)>LB2);
        else
            indB=find(ffsB(indA)<LB1);
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
CSpred=nums*jc_residuals(pitch(:,1:length(nums)))';
%CSpred=CSpreda-mean(CSpreda(onset:offset));
