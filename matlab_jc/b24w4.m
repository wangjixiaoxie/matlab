% b24w4.m
% 
cleandir4('batch',100000,500,5,5);
% segment at 1e5

% Previously implanted with cannulae in RA and used by Mark Miller

% 06.22.11 - first baseline day
% 06.25.11 (saturday) - wn on 8 and above 
% 06.28 - last day of wn

% really unstable

figure;hold on;
plot(runningaverage(timevalsPRE,20),runningaverage(distntPRE,20),'b')
plot(runningaverage(timevalsWN,20),runningaverage(distntWN,20),'r')
plot(runningaverage(timevalsPost,20),runningaverage(distntPost,20),'b')
figure;hold on;plot(hist(distntWN(173:end),[1:1:8])/length(distntWN(173:end)),'r');plot(hist(distntPre,[1:1:8])/length(distntPre))

length(find(distntPRE>5))/length(find(distntPRE>4))
length(find(distntPRE>6))/length(find(distntPRE>5))
length(find(distntPRE>7))/length(find(distntPRE>6))
length(find(distntPRE>8))/length(find(distntPRE>7))

length(find(distntWN>5))/length(find(distntWN>4))
length(find(distntWN>6))/length(find(distntWN>5))
length(find(distntWN>7))/length(find(distntWN>6))
length(find(distntWN>8))/length(find(distntWN>7))

