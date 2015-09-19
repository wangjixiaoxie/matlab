
% r3g74
% stim candidate
% 7.11.11 - started screening
% 7.12.11 - made very good stim template - use REFRAC=0.1
    % BT not off low stack
% 7.12.11 - made good wn template - REFRAC=0.2
    % BT not off low stack (to miss high stacks)
    % BT not off note beforehand that gets stim (to miss that note)
% 7.12.11
    % made new wn template 
    % BT not off low stack (to miss high stacks)
    % NOT note beforehand that gets stim (to miss that note)
        % not a BT, just a straight-up not 
    % BT off that NOT note too

% 7.14.11 - wn on at 11:40am - hit below 3250Hz (~70th prctile)
    % wn off at 5:30pm - learned a lot
% good example figure of rapid learning and gradual recovery
   

% 7.21.11 - bilateral LMAN stim array implantation

% 7.25.11 - leads in
% 7.26.11 - made new templates/adjusted old templates

% stim params - R23, L23, stim 20%, 100uA both sides, 20ms delay, 80ms dur
% 7.29.11 - WN on at 1pm - hit below 3310Hz 
            % 5:25pm - hit below 3375Hz (may have set orig threshold too
            % low)
% 7.30.11 - 12:30pm - hit below 3425Hz, stim R34 L34 same params

    
figure;hold on;
plot(runningmedian(timing3(fvalsPRE),20),runningmedian(pitchPRE(250,:),20))
plot(runningmedian(timing3(fvalsWN),20),runningmedian(pitchWN(250,:),20),'r') 
plot(runningmedian(timing3(fvalsPOST),20),runningmedian(pitchPOST(250,:),20))


