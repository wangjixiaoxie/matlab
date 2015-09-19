% w1o21
% MMAN LESION
    % electrolytic bilateral MMAN lesion - JDC - 05.05.11
    % stopped singing for a while after lesion
        % note that TW anesthetized him again ~05.07 to epoxy skull cap

% REPEATS
    % 05.14.11 - dawn - hit repeats 5 to end
    % initially reduces singing rate and causes premature song stoppages
    % 05.16.11 - 9am corrected
    %          - 6:20pm - hit repeats 8 to end???? NO I DIDN'T???
    % 05.18.11 - no change thus far
load /bulbul3/SyntaxBirds/w1o21_MMAN/repeatData.mat    
    [timevalsPrelesion,distntPrelesion] = jcrepeatdist('f','batchnotes'); 
    [timevalsWN,distntWN] = jcrepeatdist('f','batchnotes');
    
    figure;hold on;
    plot(runningaverage(timevalsPrelesion,20),runningaverage(distntPrelesion,20))
    plot(runningaverage(timevalsPostlesion,20),runningaverage(distntPostlesion,20))
    plot(runningaverage(timevalsWN,20),runningaverage(distntWN,20),'r') 
    plot(runningaverage(timevalsPostWN,20),runningaverage(distntPostWN,20),'b')    
    
% PITCH
    % 05.19.11 - dawn -switch to pitch shift
    %          - hit below 3325Hz, note "A", TH=3, MIN=2
    % good template but no learning 
    % 05.23.11 - by now learning is large and obvious
load /cardinal9/SyntaxMMAN/w1o21/0519_wnon_pitchshift/pitchData.mat
% pitch variability
    figure;hold on;
    plot(std(pitchprelesionA'),'b')
    plot(std(pitchpostlesionA'),'r')
    mean(std(pitchprelesionA(150:250,:)')) % 74.3
   mean(std(pitchpostlesionA(150:250,:)')) % 69.2    
% pitch learning
    figure;hold on;
    plot(runningaverage(timing3(fvalsprelesionA),20),runningaverage(mean(pitchprelesionA(150:250,:)),20))
    plot(runningaverage(timing3(fvalspostlesionA),20),runningaverage(mean(pitchpostlesionA(150:250,:)),20))
    plot(runningaverage(timing3(fvalswnA),20),runningaverage(mean(pitchwnA(150:250,:)),20),'r')
