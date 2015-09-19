% Note - sigwind is different from window to be CONSERVATIVE
    % results are same either way
    % sigwind is in the middle of the note vs. edges where sigma may be
        % quite large because of noise.
        


% RawData3
for i=[3:11] % For 1 it never gets there
    windowcenter=ceil(median(Predict(indrevised3(i)).Targeting));
    window=[windowcenter-64:windowcenter];
    sigwind=Predict(indrevised3(i)).onset+round((Predict(indrevised3(i)).offset-Predict(indrevised3(i)).onset)/2);
    sigwind=[sigwind-20:sigwind+20];
    meanBase=mean(mean(RawData3(i).pitchBase(window,:)'));
    sigmaBase=mean(std(RawData3(i).pitchBase(sigwind,:)'));
    if i==8
        sigmaBase=mean(std(RawData3(i).pitchBase(300,:)'));
    end
    j=0;
    slidingmean=meanBase;
    coef=1-2*isequal(Predict(indrevised3(i)).direction,'down');

        while coef*(slidingmean-meanBase)<2*sigmaBase
            j=j+1;
            slidingmean=mean(mean(RawData3(i).pitchShift(window,j:j+50)));
        end
    for k=0:50
    Learned3(i).data(k+1,:)=coef*(RawData3(i).pitchShift(:,j+k)'-mean(RawData3(i).pitchBase'));
    end
end
    
for i=[1:4] 
    windowcenter=ceil(median(Predict(indrevised4(i)).Targeting));
    window=[windowcenter-64:windowcenter];
    sigwind=Predict(indrevised4(i)).onset+round((Predict(indrevised4(i)).offset-Predict(indrevised4(i)).onset)/2);
    sigwind=[sigwind-20:sigwind+20];
    meanBase=mean(mean(RawData4(i).pitchBase(window,:)'));
    sigmaBase=mean(std(RawData4(i).pitchBase(sigwind,:)'));
    j=0;
    slidingmean=meanBase;
    coef=1-2*isequal(Predict(indrevised4(i)).direction,'down');
    if i==1
        coef=coef*-1;
        RDS=RawData4(i).pitchShift(:,find(std(RawData4(i).pitchShift(400:700,:))<100));
    else
    RDS=RawData4(i).pitchShift;
    end
    if i==3
        coef=-1;
    end
        while coef*(slidingmean-meanBase)<2*sigmaBase
            j=j+1;
            slidingmean=mean(mean(RDS(window,j:j+50)));
        end
    for k=0:50
    Learned4(i).data(k+1,:)=coef*(RDS(1:1500,j+k)'-mean(RawData4(i).pitchBase(1:1500,:)'));
    end
end
% RawData5 is #1 and #2 from RawData2 - Evren's
for i=[1:2] 
    windowcenter=ceil(median(Predict(indrevised5(i)).Targeting));
    window=[windowcenter-64:windowcenter];
    sigwind=Predict(indrevised5(i)).onset+round((Predict(indrevised5(i)).offset-Predict(indrevised5(i)).onset)/2);
    sigwind=[sigwind-20:sigwind+20];

    meanBase=mean(mean(RawData5(i).pitchBase(window,:)'));
    sigmaBase=mean(std(RawData5(i).pitchBase(sigwind,:)'));
    j=0;
    slidingmean=meanBase;
    coef=1-2*isequal(Predict(indrevised5(i)).direction,'down');
    RDS=RawData5(i).pitchShift;
        while coef*(slidingmean-meanBase)<2*sigmaBase
            j=j+1;
            slidingmean=mean(mean(RDS(window,j:j+50)));
        end
    for k=0:50
    Learned5(i).data(k+1,:)=coef*(RDS(1:1500,j+k)'-mean(RawData5(i).pitchBase(1:1500,:)'));
    end
end
% RawData1 is my data
for i=[1 2] 
    
    windowcenter=ceil(median(Predict(indrevised1(i)).Targeting));
    window=[windowcenter-64:windowcenter];
    sigwind=Predict(indrevised1(i)).onset+round((Predict(indrevised1(i)).offset-Predict(indrevised1(i)).onset)/2);
    sigwind=[sigwind-20:sigwind+20];

    meanBase=mean(mean(RawData1(i).pitchBase(window,:)'));
    sigmaBase=mean(std(RawData1(i).pitchBase(sigwind,:)'));
    j=0;
    slidingmean=meanBase;
    coef=1-2*isequal(Predict(indrevised1(i)).direction,'down');
    RDS=RawData1(i).pitchShift;
        while coef*(slidingmean-meanBase)<2*sigmaBase
            j=j+1;
            slidingmean=mean(mean(RDS(window,j:j+50)));
        end
    for k=0:50
    Learned1(i).data(k+1,:)=coef*(RDS(1:1500,j+k)'-mean(RawData1(i).pitchBase(1:1500,:)'));
    end
end

clear Learned
Learned=[Learned3(3:11) Learned4 Learned5 Learned1];
indrevised=[indrevised3(3:11) indrevised4 indrevised5 indrevised1];
 for i=1:length(Learned)
    Learned(i).data=mean(Learned(i).data);
end
aligned=zeros(28,1250);
for i=1:length(Learned)
    maximum=max(Learned(i).data(Predict(indrevised(i)).onset:Predict(indrevised(i)).offset-100));
    pre=Learned(i).data(Predict(indrevised(i)).onset:ceil(median(Predict(indrevised(i)).Targeting)));
    post=Learned(i).data(ceil(median(Predict(indrevised(i)).Targeting))+1:Predict(indrevised(i)).offset);
    aligned(i,700-length(pre)+1:700)=pre/maximum;
    aligned(i,701:700+length(post))=post/maximum;
end
for i=[1 5 8:9 12 16:18 20:22]
    aligned(i,:)=zeros(size(aligned(i,:)));
    pre=Predict(i).LearnedNorm(Predict(i).onset:ceil(median(Predict(i).Targeting)));
    post=Predict(i).LearnedNorm(ceil(median(Predict(i).Targeting))+1:Predict(i).offset);
    aligned(i,700-length(pre)+1:700)=pre;
    aligned(i,701:700+length(post))=post;
end
clear mnaligned
for i=1:1250
    if ~isempty(find(aligned(:,i)~=0))
        mnaligned(i)=mean(aligned(find(aligned(:,i)~=0),i));
    else
        mnaligned(i)=0;
    end
end
figure;plot(mnaligned(700-32-200:700-32+200)/0.91)
ylim([0 1])

clear aligned2
for i=[1:28]
    pre=Predict(i).LearnedNorm(Predict(i).onset:ceil(median(Predict(i).Targeting)));
    post=Predict(i).LearnedNorm(ceil(median(Predict(i).Targeting))+1:Predict(i).offset);
    aligned2(i,700-length(pre)+1:700)=pre;
    aligned2(i,701:700+length(post))=post;
end
clear mnaligned2
for i=1:1250
    if ~isempty(find(aligned2(:,i)~=0))
        mnaligned2(i)=mean(aligned2(find(aligned2(:,i)~=0),i));
    else
        mnaligned2(i)=0;
    end
end
hold on;plot(mnaligned2(700-32-200:700-32+200)/0.94,'r')
ylim([0 1])
