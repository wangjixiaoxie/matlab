% bk5bk25

% REDUCTION
    % 11.30.10 - templates for repeat experiment ('A' reduction)
        % 3:20pm - wn test
    % 12.01.10 - DAY ONE - dawn - hit 5 and above (5-25) to reduce (38% escape)
    % doesn't seem to sing much...
    % 12.02.10 - DAY TWO
  % No learning, not much singing -- never managed to sing through wn
    % Returned to TW CR cage
    
    
    [timevalsPre1130a,distntPre1130a] = jcrepeatdist('a','batchnotes'); % 5.36 
%[timevalsPre1130b,distntPre1130b] = jcrepeatdist('b','batchnotes'); % range: 3-6 (possible)
[timevalsWN1201a,distntWN1201a] = jcrepeatdist('a','batchnotes'); % range: 3-6 (possible)

figure;hold on;
plot(runningaverage(timevalsPre1130a,20),runningaverage(distntPre1130a,20))
plot(runningaverage(timevalsWN1201a,20),runningaverage(distntWN1201a,20))



