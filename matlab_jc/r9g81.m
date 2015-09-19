% r9g81.m
% screening on launchpad
% segment at 1e5

% 6.10 - 6.12 - screened
% computer died
% 6.14 - started screening again
    % good distn (mean=4.3, s.d.=1.03)
    % 1-(100)->2-(96)->3-(83)->4-(51)->5-(23)->6-(16)->...
    % Questions...
    % If I hit 4-->5, does 3-->4 reduce?

% 6.15 - load template - hit 5 and up, refrac 0.6
% 6.
% 6.19 - perfect targeting & already lots of learning - but singing rate has reduced
% 6.21 (pm) - wn off
% 6.24 - adjusted template
% 6.25 - WN on at dawn - hit 4 and up, refrac 0.6
figure;hold on;
plot(runningaverage(timevalsPRE,20),runningaverage(distntPRE,20),'b')
plot(runningaverage(timevalsWN,20),runningaverage(distntWN,20),'r')

[timevalsWN2,distntWN2] = jcrepeatdist('a','batchnotes'); % 2.87
figure;hold on;
plot(runningaverage(timevalsPRE2,20),runningaverage(distntPRE2,20),'b')
plot(runningaverage(timevalsWN2,20),runningaverage(distntWN2,20),'r')
plot(runningaverage(timevalsPOST2,20),runningaverage(distntPOST2,20),'b')

length(find(distntPRE2>1))/length(find(distntPRE2>0)) % 1
length(find(distntPRE2>2))/length(find(distntPRE2>1)) % 0.96
length(find(distntPRE2>3))/length(find(distntPRE2>2)) % 0.67 HIT
length(find(distntPRE2>4))/length(find(distntPRE2>3)) % 0.36 HIT
length(find(distntPRE2>5))/length(find(distntPRE2>4)) % 0.13 HIT
length(find(distntPRE2>6))/length(find(distntPRE2>5)) % 0 HIT

length(find(distntWN2>1))/length(find(distntWN2>0)) % 0.98
length(find(distntWN2>2))/length(find(distntWN2>1)) % 0.93
length(find(distntWN2>3))/length(find(distntWN2>2)) % 0.31 HIT
length(find(distntWN2>4))/length(find(distntWN2>3)) % 0.11 HIT
length(find(distntWN2>5))/length(find(distntWN2>4)) % 0.25 HIT
length(find(distntWN2>6))/length(find(distntWN2>5)) % 0 HIT


length(find(distntPRE>1))/length(find(distntPRE>0)) 
length(find(distntPRE>2))/length(find(distntPRE>1)) 
length(find(distntPRE>3))/length(find(distntPRE>2)) 
length(find(distntPRE>4))/length(find(distntPRE>3))
length(find(distntPRE>5))/length(find(distntPRE>4)) 
length(find(distntPRE>6))/length(find(distntPRE>5)) 

length(find(distntWN>1))/length(find(distntWN>0))
length(find(distntWN>2))/length(find(distntWN>1)) 
length(find(distntWN>3))/length(find(distntWN>2)) 
length(find(distntWN>4))/length(find(distntWN>3))
length(find(distntWN>5))/length(find(distntWN>4)) 
length(find(distntWN>6))/length(find(distntWN>5)) 



plot(runningaverage(timevalsPost,20),runningaverage(distntPost,20),'b')
figure;hold on;plot(hist(distntWN2(173:end),[1:1:8])/length(distntWN2(173:end)),'r');plot(hist(distntPre,[1:1:8])/length(distntPre))