fs=8000;
clear psdvaluesAC
clear psdvaluesINA
clear signal
exp=13;
firstpoint=400;
lastpoint=750;
sizer=lastpoint-firstpoint+1;

normalized=Predict(exp).ResidAC(firstpoint:lastpoint,:);
for i=1:size(normalized,2)
    pitchN=normalized(:,i);
    [psdvaluesAC(i,:),f]=periodogram(pitchN,[],[],fs);
end
normalized=Predict(exp).ResidINA(firstpoint:lastpoint,:);
for i=1:size(normalized,2)
    pitchN=normalized(:,i);
    [psdvaluesINA(i,:),f]=periodogram(pitchN,[],[],fs);
end
thous=1000*ones(1,length(f));
thous=thous';
fmod=thous./f;
diffs=mean(psdvaluesAC)-mean(psdvaluesINA);
cvdiffs=std(psdvaluesAC)./mean(psdvaluesAC);
diffs=diffs/max(diffs);
signal=zeros(sizer,50);
clear acorrsignal
for j=1:200
    t=1:1:sizer;
    % DC component
        phase=(sizer/round(rand*sizer))*2*pi;
        magnitude=diffs(1)*(rand*2-1)*cvdiffs(1);
        period=8500;
        signal(:,j)=(magnitude*cos(1/period*2*pi*t+phase))';
    for i=2:100
        phase=(sizer/round(rand*sizer))*2*pi;
        magnitude=diffs(i)*(rand*2-1)*cvdiffs(i);
        period=fmod(i)*8;
        signal(:,j)=signal(:,j)+(magnitude*cos(1/period*2*pi*t+phase))';
        acorrsignal(:,j)=xcorr(signal(:,1))/max(xcorr(signal(:,1)));
    end
end
figure;plot(signal)



            notelength=sizer;
            numms=floor(notelength/8);
            clear crosscoAC
            clear crosscoINA
            clear mnccAC
            clear mnccINA
            clear mnccLMAN
            clear crosscoLMAN
            for ii=1:numms
                first=ii*8;
                middle=firstpoint+first;
                init=500-first+1;
                for j=1:notelength
                    ab=corrcoef(Predict(exp).ResidAC(j+firstpoint,:),Predict(exp).ResidAC(middle,:));
                    crosscoAC(init+j,ii)=ab(2);
                    ai=corrcoef(Predict(exp).ResidINA(j+firstpoint,:),Predict(exp).ResidINA(middle,:));
                    crosscoINA(init+j,ii)=ai(2);
                    index1=find(~isnan(mean(signal)));
                    aL=corrcoef(signal(j,index1),signal(middle-firstpoint,index1));
                    crosscoLMAN(init+j,ii)=aL(2);
                end
            end
            for i=1:size(crosscoAC,1)
                ind1=find(crosscoAC(i,:)>0);
                mnccAC(i)=mean(crosscoAC(i,ind1));
                ind2=find(crosscoINA(i,:)>0);
                mnccINA(i)=mean(crosscoINA(i,ind2));
                ind3=find(crosscoLMAN(i,:)>0);
                mnccLMAN(i)=mean(crosscoLMAN(i,ind3));
            end

            figure;plot(mnccAC)
            hold on;plot(mnccINA,'r')
            plot(mnccLMAN,'g')
%%%%%%%%%%%
