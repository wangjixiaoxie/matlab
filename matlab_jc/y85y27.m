% y85y27

% 1.23.12
% REPEATS
% range 6-17, lots of 6,7,8,9,10
% made template - good refrac is 0.08
dirf('*.cbin.not.mat','batchnotes')
 [timevalsPre,distntPre] = jcrepeatdist('a','batchnotes');
 min(distntPre) % 5
 length(find(distntPre>=4))/length(find(distntPre>=3)) % 1
  length(find(distntPre>=5))/length(find(distntPre>=4)) % 1
 length(find(distntPre>=6))/length(find(distntPre>=5)) % 0.9811
 length(find(distntPre>=7))/length(find(distntPre>=6)) % 0.9423 - reduced?
 length(find(distntPre>=8))/length(find(distntPre>=7)) % 0.9184 -- assuming there are 7, does 8 happen
length(find(distntPre>=9))/length(find(distntPre>=8)) % 0.7333 - reduced?
 length(find(distntPre>=10))/length(find(distntPre>=9)) % 0.6667
 length(find(distntPre>=11))/length(find(distntPre>=10)) % 0.5909

 % Friday 1.27.12 at dawn - hit transition from 7 to 8
    % will he sing through?  fairly low volume wn
        % 3pm - reduced wn volume further
 % Saturday
 % Sunday - took until midday today to really sing through consistently 
 % Monday
 
  length(find(distntWN>=4))/length(find(distntWN>=3)) % 1
  length(find(distntWN>=5))/length(find(distntWN>=4)) % 1
 length(find(distntWN>=6))/length(find(distntWN>=5)) % 0.94
 length(find(distntWN>=7))/length(find(distntWN>=6)) % 0.93 - reduced?
 length(find(distntWN>=8))/length(find(distntWN>=7)) % 0.97 -- assuming there are 7, does 8 happen
length(find(distntWN>=9))/length(find(distntWN>=8)) % 0.86 - reduced?
 length(find(distntWN>=10))/length(find(distntWN>=9)) % 0.71
 length(find(distntWN>=11))/length(find(distntWN>=10)) % 0.71

figure;plot(runningaverage([timevalsPre timevalsPre1],20),runningaverage([distntPre distntPre1],20),'b')
hold on;plot(runningaverage(timevalsWN,20),runningaverage(distntWN,20),'r')
hold on;plot(runningaverage(timevalsPost,20),runningaverage(distntPost,20),'k')

% Feb 8 at dawn - Feb 11 at dusk - WN on 11 and up - if anything, this
% caused an increase in repeat length



 [a,b]=hist([distntPre distntPre1],[1:1:22]);
[c,d]=hist([distntWN(100:end)],[1:1:22]);
[e,f]=hist([distntPost],[1:1:22]);

figure;hold on;
stairs(b,a/sum(a),'Linewidth',2)
stairs(d,c/sum(c),'color','r','Linewidth',2)
 stairs(f,e/sum(e),'color','k','Linewidth',2)
