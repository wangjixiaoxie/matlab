for i=1:length(ExperimentPC)
    window1=ExperimentPC(i).time-10:ExperimentPC(i).time+10;
    coef=1-2*isequal(ExperimentPC(i).DIR,'down');
    LearnEnd(i)=coef*(mean(mean(ExperimentPC(i).pitchWN(window1,end-15:end)))-mpitchprePC(i));
    PostFF(i).data=coef*(mean(ExperimentPC(i).pitchPost(window1,:))-mpitchprePC(i));
    TimeFF(i).data=ExperimentPC(i).timePost-ExperimentPC(i).timeWN(end);
end
figure;hold on;
for i=[1 3:14]
    plot(0,LearnEnd(i),'r+')
    plot(TimeFF(i).data,PostFF(i).data,'b.')
end
for i=1:14
    a(i)=max(TimeFF(i).data)
end
indSameDay=[1 7 11 13 14];   
figure;hold on;
for i=indSameDay
       plot(0,LearnEnd(i),'r+')
    plot(TimeFF(i).data,PostFF(i).data,'b.')
end
figure;hold on
for i=3:14
    m=polyfit(runningaverage(TimeFF(i).data,20),runningmedian(PostFF(i).data,20),1);
    plot(runningaverage(TimeFF(i).data,10),runningmedian(PostFF(i).data,10),'.')
    slope1(i)=m(1); % per hour
    yint1(i)=m(2);
end

indNext=[1 3:6 8:10 12 14];
for i=indNext
    indexNext=find(TimeFF(i).data>9 & TimeFF(i).data<35);
    mnNext(i)=mean(PostFF(i).data(indexNext));
end