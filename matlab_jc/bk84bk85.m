% bk84bk85

% syntax/lesion candidate

% TW collected pre data

% 10.18.11 - lesion MMAN & LMAN

% POSTLESION
    % SYNTAX SHIFT - TW

    % FF SHIFT (DataFF.mat)
        % 11.01 - wn on for FF - hit below 3215Hz (~70%)
            % 2:22pm - hit below 3225hz
        % 11.02 - wn off at dark
        % no FF change
        
load DataFF.mat
figure;hold on;
plot(timing3(fvPre),mean(pitchPre(210:250,:)),'.')
plot(timing3(fvWN),mean(pitchWN(210:250,:)),'r.')
plot(timing3(fvPOST),mean(pitchPOST(210:250,:)),'.')

plot(runningaverage(timing3(fvPre),20),runningaverage(mean(pitchPre(210:250,:)),20),'Linewidth',3)
plot(runningaverage(timing3(fvWN),20),runningaverage(mean(pitchWN(210:250,:)),20),'Linewidth',3,'Color','r')
plot(runningaverage(timing3(fvPOST),20),runningaverage(mean(pitchPOST(210:250,:)),20),'Linewidth',3)
