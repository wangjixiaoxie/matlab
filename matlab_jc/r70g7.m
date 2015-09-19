% r70g7

% REDUCTION
    % 11.16.10 - templates for repeat experiment
    % 11.17.10 - DAY ONE - dawn - hit 3 and above to reduce
    % 11.18.10 - learns great - robust and fast
[timevalsPre1015,distntPre1015] = jcrepeatdist('c','batchnotes1015'); % 2.87
[timevalsPre1016,distntPre1016] = jcrepeatdist('c','batchnotes1016'); % 2.86
[timevalsWN1017,distntWN1017] = jcrepeatdist('c','batchnotes'); % rapid - down to 2

figure;hold on;
plot(runningaverage(timevalsPre1015,30),runningaverage(distntPre1015,30))
plot(runningaverage(timevalsPre1016,30),runningaverage(distntPre1016,30))
plot(runningaverage(timevalsWN1017,30),runningaverage(distntWN1017,30))
plot(runningaverage(timevalsPost1019,30),runningaverage(distntPost1019,30))
plot(runningaverage(timevalsPost1020,30),runningaverage(distntPost1020,30))
% EXPANSION
    % 11.19.10 - templates - hit 0-2
    % 11.20.10 
    %   10:45am - test WN - doesn't cause infinite loop
    %   12:01pm - change BTAF MAX to 40 to eliminate false positives
    %   9pm - WN on
    % 11.21.10 - DAY ONE - hit 0-2
    % 11.22.10 -7am - problem - WN is missing first note - change tmp
    % 11.23.10 -dawn - WN back on - looks good - "TemplatesRepExpand2.mat"
    

[timevalsPre1022,distntPre1022] = jcrepeatdist('c','batchnotes'); % 
[timevalsWN1023,distntWN1023] = jcrepeatdist('c','batchnotes'); % 
[timevalsPost1025,distntPost1025] = jcrepeatdist('c','batchnotes'); %
figure;hold on;plot(runningaverage(timevalsPre1022,30),runningaverage(distntPre1022,30))
plot(runningaverage(timevalsWN1023,30),runningaverage(distntWN1023,30))
plot(runningaverage(timevalsPost1025,30),runningaverage(distntPost1025,30))


%%%%%%%%%%%%%%%
%%%%% POSTLESION
%%%%%%%%%%%%%%%
% REDUCTION
    % 2.10.10 - first preday (I think Tim listed a non-existent folder as
        % output so we lost the previous few days)
        % mean=2.93
        % 12:30pm - load templates - same as before - hit 3 and above
        %
    % 2.11.10
        % looks good - one concern is that third syllable sometimes avoids
        % getting recognized but this is rare - another concern is wn
        % refrac but I think it should be okay
    % 2.12.11 - DAY ONE
[timevalsPre0211,distntPre0211] = jcrepeatdist('c','batchnotes'); % 
    % 2.14.11 - it doesn't appear that he sang enough for much learning.
[timevalsWN0212,distntWN0212] = jcrepeatdist('c','batchnotes'); % 
[timevalsPre0215,distntPre0215] = jcrepeatdist('c','batchnotes'); % 
    % 2.16.11 - DAY ONE - Repeat reduction post-lesion
        % 12:04pm (computer time) - adjusted template - there was a problem
        % that the long flat stack note (used for BTAF AND) had changed so that it
        % was occassionally not recognized
        % note to self - include catch trials at end of second day (2.17 - thursday)
    % 2.17.11
        % 1:37pm computer time - catch trials = 20%
        % 7pm - switch back to 10%
    % 2.18.11
        % 2:30pm - catch trials=20%
        
    % 2.19.11 - DAY FOUR - wn off at dark - good learning (repeat reduction)
    
    % 2.20.11
    
[timevalsWN0216,distntWN0216] = jcrepeatdist('c','batchnotes'); % 
[timevalsPost0220,distntPost0220] = jcrepeatdist('c','batchnotes'); %     
    
figure;hold on;
plot(runningaverage(timevalsPre0211,30),runningaverage(distntPre0211,30))
plot(runningaverage(timevalsWN0212,15),runningaverage(distntWN0212,15),'r')
plot(runningaverage(timevalsPre0215,30),runningaverage(distntPre0215,30))
plot(runningaverage(timevalsWN0216,30),runningaverage(distntWN0216,30),'r')
plot(runningaverage(timevalsPost0220,15),runningaverage(distntPost0220,15),'k')

[timevalsPre0302,distntPre0302] = jcrepeatdist('a','batchnotes'); %   
[timevalsWN0303,distntWN0303] = jcrepeatdist('a','batchnotes'); %   
[timevalsPost0308,distntPost0308] = jcrepeatdist('a','batchnotes'); %   
figure;hold on;
plot(runningaverage(timevalsPre0302,15),runningaverage(distntPre0302,15),'k')
plot(runningaverage(timevalsWN0303,20),runningaverage(distntWN0303,20),'r')
% 03.02.11 - mean pre is around 2.1
    %   3:10pm - change BTAF MAX to 40 to eliminate false positives
% 03.03.11 - WN ON HIT first two repeats MIN=1, MAX=2
% looked like learning on the first day but then a decrease
% 03.05.11 - wn accidentally interrupted from 10:30am to 12:30pm computer
% time
% 03.06.11 - as of dark, no learning
