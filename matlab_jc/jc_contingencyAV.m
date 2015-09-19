function [avgnn]=jc_contingencyAV(norm,chspec,BarHeight,BarDirection,chunk)
%%%Calibrate 10ms period
if chspec==0
    width=round(10*(44100/(4*1000)));
else
    width=(10*(32000/(4*1000)));
end

%%%%Assign the contingency to this region:
center=round(mean(chunk));
for i=1:size(norm,2)
    mm(i)=mean(norm(chunk(1):chunk(2),i));
end
ss=std(mm);

%Contingency
if strcmp(BarDirection,'above')
    count=1;
    %%%Above contingency
    for i=1:size(norm,2)
        if mm(i)>BarHeight*ss
            nn(:,count)=norm(:,i);
            count=count+1;
        end
    end
else
    count=1;
    %%%Below contingency
    for i=1:size(norm,2)
        if mm(i)<-1*BarHeight*ss
            nn(:,count)=norm(:,i);
            count=count+1;
        end
    end
end
num=count-1;
for i=1:length(nn)
    avgnn(i)=mean(nn(i,:));
end
%if chspec==0
%    avgnn=resample(avgnn,32,44);
%end
%%%%%%%%%%Previous Output Attempts%%%%%%%%%%%%
%avgchop=avgnn(150:length(avgnn));
%hhTime=60+find(avgchop>(min(avgnn(60:100))+0.01),1,'first');
%hhTime=150+find(avgchop<(avgnn(150)-0.01),1,'first');

%Output
%falloff=(abs(avgnn(center))-abs(avgnn(center+width)+avgnn(center-width))/2)/abs(avgnn(center));
%one side
%falloff=(abs(avgnn(center))-abs(avgnn(center+width)))/abs(avgnn(center));

%Absolute
%falloff=(abs(avgnn(center))-abs(avgnn(center+width)+avgnn(center-width))/2);



    