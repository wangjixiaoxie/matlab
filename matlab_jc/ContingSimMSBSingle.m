function [shiftup]=ContingSimMSBSingle(toffsets,pitch,cont)
%aavfin=ContingSimMSB(Predict(1).Targeting-Predict(1).onset,Predict(1).ResidAC(Predict(1).onset:Predict(1).onset+350,:),60);
% pitch - just the residuals within the note (between onset and offset)
pitch=jc_residuals(pitch+2500);
pitch=pitch(:,floor(rand(1,200)*size(pitch,2)+1));
if cont<50
    pitch=pitch*-1;
end

% only consider targeting within the syllable
ind1=find(toffsets>1);
ind2=find(toffsets(ind1)<size(pitch,1));
toffsets=toffsets(ind1(ind2));
count=0;
indices=zeros(1,size(pitch,2));
% for each targeting point, find pitch curves that are compliant 
    for i=1:length(toffsets)
        offset=round(toffsets(i)); %round(mean(toffsets)+0.8*std(toffsets)*randn);
        ffs1=(pitch(offset,:)); % estimate of pitch within the window
        L1=prctile(ffs1,60);
        ind=(ffs1>L1);
        if sum(ind)>0
            count=count+1;
            indices=[indices+ind];
        end
    end
histog=indices/sum(indices); % probability that each curve will escape
for i=1:size(pitch,1)
    correlatedval(i)=pitch(i,:)*histog';
end

% Decision: just do most likely 5ms 
wideness=40; % 5ms
    for i=1:size(pitch,1)-wideness
        thr5ms(i)=median(correlatedval(i:i+wideness-1));
    end
    [p,siteup]=max(thr5ms);
    shiftup=round(siteup+wideness/2);
