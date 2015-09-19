% r29g89

% tried repeat learning pre-lesion and didn't learn, but Tim did LMAN lesion anyway (for branchpts)

% REDUCTION
    % 11.16.10 - templates for repeat experiment
    % 11.17.10 - dawn - hit 3 and above to reduce
    % 11.18.10 - doesn't learn much at all

[timevalsPre1015,distntPre1015] = jcrepeatdist('c','batchnotes1015'); % 2.61
[timevalsPre1016,distntPre1016] = jcrepeatdist('c','batchnotes1016'); % 2.66
[timevalsWN1017,distntWN1017] = jcrepeatdist('c','batchnotes'); % slow - down to 2.1

figure;hold on;
plot(runningaverage(timevalsPre1015,30),runningaverage(distntPre1015,30))
plot(runningaverage(timevalsPre1016,30),runningaverage(distntPre1016,30))
plot(runningaverage(timevalsWN1017,30),runningaverage(distntWN1017,30))

% EXPANSION
    % 11.19.10 - templates 
    % 11.20.10 - 
    %   % 10:45am - test WN - doesn't cause infinite loop
    %   % 12:21pm - decrease BTAF MAX to 40 and change template to hit 0-3
    %   9pm - WN on
    % 11.21.10 - DAY ONE - hit 0-3
    % Targeting was perfect but change wasn't robust
 %figure;hold on;
  plot(runningaverage(timevalsPre1020,30),runningaverage(distntPre1020,30))   
 plot(runningaverage(timevalsWN1021,30),runningaverage(distntWN1021,30))   
 [timevalsPre1020,distntPre1020] = jcrepeatdist('c','batchnotes'); % 2.87        
[timevalsWN1021,distntWN1021] = jcrepeatdist('c','batchnotes'); % 2.87        
    
% Post-lesion
[timevalsPre0117,distntPre0117] = jcrepeatdist('c','batchnotes'); % slow - down to 2.1
 figure; plot(runningaverage(timevalsPre0117,30),runningaverage(distntPre0117,30))   
% 1.18.11 - hit 3 and above to reduce
    %   2:45pm - adjusted acquisition because of speaker noise
% 1.19.11 - 7:30am - error in acq caused wn to stop
    % no learning prior to that point - not worth repeating the experiment