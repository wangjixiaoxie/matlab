% ......
for i=1:14
    clear indShift
    clear indBase
    tvalsBase=timing3(RawDataLONG(i).fvBase);
    tvalsShift=timing3(RawDataLONG(i).fvShift);
    MeanpitchBase=mean(RawDataLONG(i).pitchBase(RawDataLONG(i).start:RawDataLONG(i).end,:));
    MeanpitchShift=mean(RawDataLONG(i).pitchShift(RawDataLONG(i).start:RawDataLONG(i).end,:));
    % find days
        count=1;
        indBase(1)=1;
        for j=2:length(tvalsBase)
            if tvalsBase(j)>tvalsBase(j-1)+8
                count=count+1;
                indBase(count)=j;
            end
        end
    % for each day, get the CV
        count=0;
        for j=2:length(indBase)
            pitchDay=MeanpitchBase(indBase(j-1):indBase(j)-1);
            if length(pitchDay)>50
                residDay=detrend(pitchDay,'linear');
                Baseline(i).day(j-1)=std(residDay)/mean(pitchDay);
            end
        end
    % find days
        count=1;
        indShift(1)=1;
        for j=2:length(tvalsShift)
            if tvalsShift(j)>tvalsShift(j-1)+8
                count=count+1;
                indShift(count)=j;
            end
        end
    % for each day, get the CV
        count=0;
        for j=2:length(indShift)
            pitchDay=MeanpitchShift(indShift(j-1):indShift(j)-1);
            if length(pitchDay)>50
                count=count+1;
                residDay=detrend(pitchDay,'linear');
                Shifting(i).day(count)=std(residDay)/mean(pitchDay);
            end
        end
          
end

Shift2=[];
for i=1:14
    Shift2=[Shift2 Shifting(i).day];
end
Base2=[];
for i=1:14
    Base2=[Base2 Baseline(i).day];
end
mean(Shift2(find(Shift2>0)))
mean(Base2(find(Base2>0)))
FFrangeShift=[Shift1 Shift2];
FFrangeBase=[Base1 Base2];

figure;hold on;
bar([SummaryData.meanCVbase SummaryData.meanCVshift],'w')
errorbar(1,SummaryData.meanCVbase,SummaryData.semCVbase,'k')
errorbar(2,SummaryData.meanCVshift,SummaryData.semCVshift,'k')
xlim([0.4 2.6])