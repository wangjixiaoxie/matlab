% b65o7.m
% 

% MMAN lesion on 8.26.11 
% started singing soon but sequence (initially) and acoustics (after) were
% degraded

% Song looked good by 8.31.11
% Increased stereotypy at the branch point!!!
% 9.1.11-9.2.11 --> 87% A-->B, 13% of either A-->C or A-->X (other)
    % note that A-->B %age may be higher in the mornings
    % needs some template tweaking

% 9.5.11 - a bit of instability - decrease in p(A-->B)
% 9.6.11 - wn on at dawn
% 9.8.11 - adjusted template at 10:30am
% 9.9.11 - wn off at end of day (4days of wn)

% 9.13.11 - new template for repeat learning
%   - mean = 3.1
%   - hit 2 and up (as before MMAN lesion)
% 9.14.11 - wn on at dawn for repeat learning
% 9.18.11 - wn off at dawn (actually 8am)

% 9.20.11 - accidental wn on at 7pm
% 9.21.11 - wn off at 10am

% 9.27.11 - FF - wn on at 2:30pm
% 9.28.11 - ---- wn off at 2:30pm - about 30-40Hz learning AM to AM


load DataSeq912.mat
    figure;hold on;
    plot(runningaverage(tvsorted902,60),runningaverage(numsorted902,60),'Color','k')
    plot(runningaverage(tvsorted905,60),runningaverage(numsorted905,60),'Color','k')
    plot(runningaverage(tvsorted906,60),runningaverage(numsorted906,60),'Color','r')
    plot(runningaverage(tvsorted909,60),runningaverage(numsorted909,60),'Color','r')
    plot(runningaverage(tvsorted911,60),runningaverage(numsorted911,60),'Color','k')

load DataRepeats0916.mat
    figure;hold on;
    plot(runningaverage(timevalsPre,30),runningaverage(distntPre,30),'Color','k')
    plot(runningaverage(timevalsWN,30),runningaverage(distntWN,30),'Color','r')
    plot(runningaverage(timevalsPost2,30),runningaverage(distntPost2,30),'Color','k')
     plot(runningaverage(timevalsPost3,30),runningaverage(distntPost3,30),'Color','k')
   

    
% Repeat experiments (see below)
dirf('*.cbin.not.mat','batchnotes')
fvalsB=findwnoteJC('batchnotes','b','a','',0,[2000 2700],8500,1,'obs0',1)
fvalsC=findwnoteJC('batchnotes','c','a','',0,[2000 2700],8500,1,'obs0',1)
fvalsX=findwnoteJC('batchnotes','x','a','',0,[2000 2700],8500,1,'obs0',1)
tv912B=timing4(fvalsB);
tv912C=timing4(fvalsC);
tv912X=timing4(fvalsX);
tvall912=[tv912B tv912C tv912X];
numsall2=[ones(1,length([tv912B])) zeros(1,length([tv912C])) zeros(1,length([tv912X]))];
[b2,ind2]=sort(tvall912);
numsorted912=numsall2(ind2);
tvsorted912=tvall912(ind2);
% probability of B as a function of time
figure;hold on;
plot(runningaverage(tvsorted902,60),runningaverage(numsorted902,60),'Color','k')
plot(runningaverage(tvsorted905,60),runningaverage(numsorted905,60),'Color','k')
plot(runningaverage(tvsorted906,60),runningaverage(numsorted906,60),'Color','r')
plot(runningaverage(tvsorted909,60),runningaverage(numsorted909,60),'Color','r')
plot(runningaverage(tvsorted911,60),runningaverage(numsorted911,60),'Color','k')
plot(runningaverage(tvsorted912,60),runningaverage(numsorted912,60),'Color','k')


% segment at 5e4

% previously used for syntax by Tim

% 06.28 at 8:27am - wn on, hit repeats 5 and above using template that also
    % recognizes note pvs to repeats (so tmp hits 6 and above)
    % (supposed by go on at 7am but time switch problem)
    % 10:08am - slight template adjustment (added amplitude threshold)
% 06.28 at 6pm(?) - changed soundboxes and he stopped singing
    % also made a mistake with the amplitude threshold
% 06.29 at 12:30pm - corrected the problem with the amplitude threshold
% 06.30 - conclusion - only targeted ones decreased

% 07.05 - wn on again - hit repeats 2 and above using same template above
    % big, rapid learning
% 07.06 - wn off at end of day

    
 [timevalsWN2,distntWN2] = jcrepeatdist('a','batchnotes'); % 2.87   
 [timevalsPost2,distntPost2] = jcrepeatdist('a','batchnotes'); % 2.87   

figure;hold on;
plot(runningaverage(timevalsPre1,20),runningaverage(distntPre1,20),'b')
plot(runningaverage(timevalsWN1,20),runningaverage(distntWN1,20),'r')
plot(runningaverage(timevalsPost1,20),runningaverage(distntPost1,20),'b')

plot(runningaverage(timevalsPre2,20),runningaverage(distntPre2,20),'b')
plot(runningaverage(timevalsWN2,20),runningaverage(distntWN2,20),'r')
plot(runningaverage(timevalsPost2,20),runningaverage(distntPost2,20),'b')


length(find(distntPRE>1))/length(find(distntPRE>0)) % 0.98
length(find(distntPRE>2))/length(find(distntPRE>1)) % 0.92
length(find(distntPRE>3))/length(find(distntPRE>2)) % 0.38 change?
length(find(distntPRE>4))/length(find(distntPRE>3)) % 0.23 HIT
length(find(distntPRE>5))/length(find(distntPRE>4)) % 0.21 HIT


length(find(distntWN630>1))/length(find(distntWN630>0)) % 0.99
length(find(distntWN630>2))/length(find(distntWN630>1)) % 0.96
length(find(distntWN630>3))/length(find(distntWN630>2)) % 0.43 change? NO
length(find(distntWN630>4))/length(find(distntWN630>3)) % 0.14 HIT - decrease
length(find(distntWN630>5))/length(find(distntWN630>4)) % 0 HIT - decrease