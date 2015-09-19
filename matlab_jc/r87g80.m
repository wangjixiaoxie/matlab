% r87g80

% September 2011 - use for Directed Reversion test? 
    % YES - good directed singer
    
    % 10.03.11 - tested with several females - sang to both
        % 1004_screen --> sang [205,225,227,229] around 3pm
            % variability reduction in low stack note (b) 
            % unclear in other notes
        % 4:15pm - started template test
     % 


r87g80stimanaly.m

% 9.21.2010 - baseline data
% 9.22 - WN on at 9pm
    % Two contexts for low stack
    %   1. 'AB' - after intro-ish downsweep - more common
    %   2. 'CB' - shortly after a high stack note - less common
    % FF for AB > FF for CB
    % Hit all below median FF for AB (which is ~2630Hz)
        % Hits all CBs and half of ABs
% 9.23
    % Learning - hit below 2670Hz at 11:46am
    
% 9.30
    % SURGERY - implanted bilateral LMAN stim arrays
% 10.06 - leads on - singing - but not singing very much
    % template at 1:10pm
% 10.07 - stim on at 10:30am
    % R14, L23 - 70uA - caused song stoppages
    % 10:55am - lowered to 50uA
    % 12:33pm - lowered to 30uA
        % GREAT
% 10.08 - wn on at dawn - hit C above 3230Hz to shift down
% R14, L23, 160ms delay, 80ms duration - 70uA
        % 11am - adjusted WN template (it was pretty good already)
        % 11:30am - slightly adjusted Stim template
        % 3:10pm - switch to 120ms delay, 90ms duration
% 10.09 - 10:00am - lowered to hit above 3220Hz (around 30th prctile)
% 10.10 - 9pm - wn off - allow recovery

% Plan to increase stim duration to allow 100% hits
% 10.12 - recovering well...
% 10.13 - may be still recovering - ready to go
        % hit below 3290Hz (?)
        
% 10.14 
% 9:30am - 100% stim test - 30uA
% 10am - wn on - hit below 3300Hz with 100% stim
% noon - 10% wn catch trials (i.e. notes) - 0% catch songs
% 12:38pm - 3320Hz
% ~3pm - turned wn off and stim to 50%- catch trials
    % no change in pitch of non-stimmed syllables
% 5:20pm - turned stim to 50uA
% 10pm - wn on, 0% catch

% 10.15 - DAY ONE - 100% stim to interfere with acquisition
% 10am - raised threshold to 3370Hz
% 11:20am - lowered threshold to 3360Hz

% 10.16 - DAY TWO - 100% stim to interfere with acquisition
% 11:18am - eliminated catch trials
% 11:32am - lowered threshold to 3330Hz (70% for the morning)

% 10.17 - DAY THREE - 100% stim to interfere with acquisition
% ~3:30pm - wn and stim off as a 'catch trial' - looks like learning
% 4:20pm - wn off - stim 50%

% 10.18
% 3:36pm - test - stim on 2nd low stack (B)
% 7pm - stim on C
% 8:30pm - stim on B at 100%

% 10.19 - DAY ONE - 100% stim + WN
% 7-9am - big facilitation of stim effect (shows it isn't wn-dependent)
% 9:50am - wn on - hit below 2680Hz
    % 10:35am - moved to 2670Hz
    % 11:26am - back to 2680Hz

% 10.20 - DAY TWO - 100% stim + WN
% 7pm - WN off, 0% stim - test for learning - YES
% 8pm - 30% stim

% 10.21
% 4pm - test 120ms stim
% 7pm - flipped stimulators off/on because error lights - test 40ms stim
% 7:30pm - stim off

% 10.22- 8:20pm - no wn, 100% stim
% 10.23 - 100% stim
% 10.24 - 100% stim
%   - 7:30pm - 0% stim
%   - 8pm - 50% stim
%   - 8:30pm - 0% stim

% 10.25
%    - 2:54pm - handled bird to clean cage, check connections, etc.
%   - 50% stim - plan to do downshift
% - 100% stim at ~5:25pm
% - wn on at ~7:30pm - hit above 2600Hz

% 10.26
% - 10:30am - hit above 2585Hz
% 2pm - tested stim - looked to be good
% 4:30pm - wn off

% 10.27
% - wn on, stim 20%, hit above 2600Hz
% - 9pm - hit above 2583Hz
% Test other stim configs midday tomorrow

% 10.28
% AM data - same as PM for comparison
% noon data - checking stim configs
% PM data - same as AM for comparison
% 9pm - changed template and threshold (to 2515Hz)


% 10.29 - wn off at noon

% 10.31 (all times computer time)
%   10:45am - stim 100%
%   11:50am - wn on - hit above ???
%   3:40pm - wn off

% 11.01 (recovery has occurred, ready to reverse)
%   10:32am - stim 100%
%   11:50am - wn on, stim 100% - hit below 2620Hz
%   12:15pm - hit below 2610Hz
%   3:45 to 3:55 - catch - no learning yet?
%   4:15pm - back to 2620Hz
%   7:30pm - wn off

% 11.02 
%   10:22am - stim 100%
%   11:35am - hit above 2560Hz
%   1:35pm - 2550
%   3:03pm - 2540
%   3:55pm - wn off
% NOTE THAT STIM VALUE REMAINS 20 Volts

% 11.03
%   10:21am - stim 100%
%   11:40am - hit below 2625Hz
%   6:45pm - wn off

% 11.04 (control - chronic stim) -- looks like no learning
%   10:32am - stim 100%
%   3:40pm - stim off

% 11.05 (control - chronic stim) -- looks like no learning
%   10:40am - stim 100%

% 11.07 
%   7:30pm - delay 110ms, stim 90ms - hits 'C' - stim 20%

% 11.08
%   11:10am (computer time) - stim 100%
%   ~12:30pm (computer time) - wn on - hit below 3260Hz
%   3:40pm - hit below 3280Hz
%   5:47pm - wn off, stim 20%




% EXPERIMENT 35 - 43
figure;hold on;
subplot(141);hold on;
plot(mean(experiment(36).contours(:,experiment(36).crctind(32:end))'),'b')
plot(mean(experiment(37).contours(:,experiment(37).crctind(306:end))'),'r')
plot(mean(experiment(38).contours(:,experiment(38).crctind(1:50))'),'g')
plot(mean(experiment(38).contours(:,experiment(38).crctind(1:110))'),'k') % entire afternoon
xlim([300 700])
subplot(142);hold on;
plot(mean(experiment(38).contours(:,experiment(38).crctind(236:end))'),'b')
plot(mean(experiment(39).contours'),'y')
plot(mean(experiment(40).contours(:,539:end)'),'r')
plot(mean(experiment(41).contours(:,experiment(41).crctind(1:50))'),'g')
plot(mean(experiment(41).contours(:,experiment(41).crctind(1:357))'),'k')
xlim([300 700])



% STIM B2
% 11.10
%   stim on 20% around 7pm - NO - Was off due to testing PIC chip
% 11.11
%   11:45am (computer time) - stim on 20% - wn on - hit below 2645Hz (mean+20Hz)
% 11.12
%   11:45am - wn off

% 11.14
%   10:36am - stim on 100%
%   11:45am - hit below 2665Hz - mean+20Hz
% 11.15
%   12:30pm - stim off, wn off

% 11.16
%   10:30am - stim on 100% (no wn)

% 11.17
%   11:45am - stim back to 20%
%***********************************************%
% 3:49pm - changed to 0ms delay stim triggerring off WN template - 20% stim

% 11.18
% 10:30am - 100% stim (triggering off WN template)
% 11:45am - wn on - hit below 2645-2655Hz (varied a bit over day)

% 11.19
% 11:45am - wn off, stim 20% (note that stim is still acausal window)
%***********************************************%
% 12:53pm - stim 20% in premotor window

% 11.20
% Decide to wait until 11.21 for full recovery

% 11.21
% 1:15pm computer time - wn on hit above 2616Hz (20Hz below median=2636Hz)

% 11.22
% 1:15pm computer time - wn off - learning is ~40Hz

% 11.23 
% Decide to wait until 11.24 for full recovery
% 7:23pm - NOTICED WIRE 2 pulled out - corresponds to R3 (not being used)

% 11.24
% 9:30am computer time - leads were getting loose, so replugged them
% Note that stim effect was still very good as of AM
% 10:30am - started singing again
% 11:20am - stim on 100%
% 12:50am computer time (noon real time) - wn on - hit below 2624Hz
% 6:40pm - hit above 2615Hz

% 11.25
% stim failed to block expression
% wn off around noon

% 11.26 
% 100% stim on at 11:37am computer time
% wn on at 12:09pm computer time

% 11.27  ---- PROBLEM - STIM DID NOT BLOCK LEARNING - STIM EFFECT NEGLIGIBLE
% wn off at 12:09pm computer time
% 11.27 - 7pm - stim off
% STIM OFF 
% 11.28 - 3pm - stim on - good offset effect

% 11.29
% 10:51am - 100% stim on at 10:51am computer time
% 12:30pm - wn on, hit above 2600Hz 

% Looks like stim isn't blocking learning
% Switch stim configs
% 12R - definite improvement, 34L

% 11.30
% 11am - stim on 100% - soundbox door accidentally open

% 12.01
% 10:20am - stim on 100%
% 12:30pm - wn on, hit above 2615Hz

% 12.02
% 12:30pm - wn off
% 7pm - switch to acausal stim

% 12.03 
% 10:20am - acausal stim - (wn template, 0ms delay)
    % surprisingly this caused an acute effect
% 12:30pm - wn on, hit above 2580Hz (median-20Hz)

% 12.04
% 12:30pm - wn off
% 12:45pm - switched to causal/premotor stim



%%%%% 
load experimentREV.mat
%%%%%
% 12.19 - TW WN on at dawn - hit below 5325Hz
% noon - hit below 5345Hz
% 8pm - hit below 5400Hz

% 12.20
% pretty stable for most of the day
% 7:45pm - hit below 5430Hz

% 12.21
% 10:45am (computer time) - hit below to 5490Hz (solid learning)


figure;hold on;
plot(runningmedian(mean(experiment(37).contours(350:400,:)),60),'b')
plot(runningmedian(mean(experiment(46).contours(350:400,:)),60),'r')
plot(runningmedian(mean(experiment(40).contours(350:400,:)),60),'k')