for i=ind
    window1=FFmwin(i).data;
    coef=1-2*isequal(Experiment(i).DIR,'down');
    PostFF(i).data=coef*(mean(Experiment(i).pitchACpost(window1,:))-mpitchpreTarg(i));
    TimeFF(i).data=Experiment(i).timeACpost-Experiment(i).timeAPVwn(end);
end
figure;hold on;
for i=ind
    plot(0,LearnEnd(i),'r+')
    plot(TimeFF(i).data,PostFF(i).data,'b.')
end
for i=1:14
    a(i)=max(TimeFF(i).data)
end
indSameDay=[1 7 11 13 14];   
figure;hold on;
for i=ind%SameDay
      % plot(0,LearnEnd(i),'r+')
    plot(TimeFF(i).data,PostFF(i).data,'b.')
end
figure;hold on
for i=ind
    m=polyfit(runningaverage(TimeFF(i).data,20),runningmedian(PostFF(i).data,20),1);
   plot(runningaverage(TimeFF(i).data,10),runningmedian(PostFF(i).data,10))
    datsize1(i)=max(TimeFF(i).data)-min(TimeFF(i).data);
    slope1(i)=m(1); % per hour
    yint1(i)=m(2);
    [b,a]=sort(TimeFF(i).data);
%     starters(i)=mean(PostFF(i).data(a(1:1*round(end/10))));
%     enders(i)=mean(PostFF(i).data(a(9*round(end/10)):end));
end

indNext=[1 3:6 8:10 12 14];
for i=indNext
    indexNext=find(TimeFF(i).data>9 & TimeFF(i).data<35);
    mnNext(i)=mean(PostFF(i).data(indexNext));
end