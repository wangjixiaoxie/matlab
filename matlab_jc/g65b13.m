% g65b13

% Possible syntax candidate

% 05.10.2010 - transferred to Maxwell and began screening

% A (branch point) ---> the last of multiple downsweeps
% B ---> high stack note (generally repeated multiple times)
% C ---> intro-type note generally followed by high stack notes
% Initial probability is around 50/50

% 05.13 ---> hit B at dawn (actually began at 8pm on 05.12)
% 05.13 --- 11am - changed wn to normal file (started with 100ms delay)
%           7:30pm - changed TH from 2.5 to 2 for first three templates to
%                   reduce false positives on C
% 05.14 --- noon - up to 85%


% 08.06 - bilateral RA cannulae implantation - went perfectly
    % used 40degs and 0.13mm caudal of y0 - down ~2500um to central RA, implanted 1200um down
% 08.07 - in maxwell soundbox recording
%   song looks good - both transitions
% 08.23 - probes implanted in afternoon
% 08.24 - morning - pulled left probe out
% 08.24 -new probes implanted at 4pm
    % Measure FF in note D -- 2nd stack note after note C
    % Choose this note because 
        % never hit by WN when hitting B or C (WN may go into 1st note after
            % note C)
% 08.25 - singing starting at 9am
%           %   70.15% C for the afternoon
%      - 3:45pm - 5mMapv on at 1.5uL/min
%           % 35percent variability reduction
%           %   68.79% C
%      - 4:45pm - acsf on
%           % variability recovers by 6:30 (last 20 notes)
%       - 6:41pm - 5mMapv on at 1.5uL/min
%           % similar variability reduction
%           %   70.40% C
%       - 8:27pm  - acsf back on
%       -  - wn on C transition
% 8.26 - probes clogged  
% 8.27 - new probes at noon
% 8.28 - singing at dawn - 74% C
%      - WN on at ~8am
%      - template adjusted at ~11am
%      - avg of 55% C from 3pm to 6pm - learning
%      - 5mM apv at ~6pm - avg of 58% C - no reversion
%           - approximately a 35% variability reduction
% 8.29 - 3pm - clog - removed probes, cleaned and reinserted them at 3:33pm
% 8.30 - morning - singing - 55% C
%      - apv - 46 catch transitions - 54% C
%      - probes out because of leak/clog around 2pm
%      - wn off at 5:05pm
% 8.31 - new probes at 3pm
%      - tmp test - FF reversion experiment
%      - WN on at 10pm - hit below 3500Hz on all first high stacks - around  75% 
%           pretty stable - within 10Hz of previous day
% 9.01 - TW template adjustment at 9am
%      - hit below 3580Hz at 11:30am
%      - hit below 3610Hz at 3:45pm
%      - 5mM apv on at 5:30pm or so at 2.0uL/min
%      - apv off at 7:30pm or so at 3.0uL/min
%      - hit below 3625Hz

% --- Switched AP5 --- these experiments had bad inactivations
    % 9.02 
    %      - hit below 3700Hz at 12:30pm
    %      - made new DL-AP5
    %
    
    
    
% Tracking the variation - load Data907.mat  
    figure;plot(timing3(fv828apvECX),mean(pitch828apvECX(200:250,:)),'*','Color','r')
    hold on;plot(timing3(fv828acsfECX),mean(pitch828acsfECX(200:250,:)),'*','Color','b')
    hold on;plot(timing3(fv825apvECX),mean(pitch825apvECX(200:250,:)),'*','Color','r')
    hold on;plot(timing3(fv825acsfECX),mean(pitch825acsfECX(200:250,:)),'*','Color','b')
    hold on;plot(timing3(fv830apvECX),mean(pitch830apvECX(200:250,:)),'*','Color','r')
    hold on;plot(timing3(fv830acsfECX),mean(pitch830acsfECX(200:250,:)),'*','Color','b')
    hold on;plot(timing3(fv901acsfECX),mean(pitch901acsfECX(200:250,:)),'*','Color','b')
    hold on;plot(timing3(fv901apvECX),mean(pitch901apvECX(200:250,:)),'*','Color','r')
    hold on;plot(timing3(fv902acsfECX),mean(pitch902acsfECX(200:250,:)),'*','Color','b')    
    hold on;plot(timing3(fv902apvECX),mean(pitch902apvECX(200:250,:)),'*','Color','r')    
    hold on;plot(timing3(fv903apvECX),mean(pitch903apvECX(200:250,:)),'*','Color','r')

figure;hold on;
plot(timing3(fv901wnEC),pitch901wnEC(160,:),'*')
plot(timing3(fv902acsfEC),pitch902acsfEC(160,:),'*')

% FF analysis - note 'e'
% 830acsf - pre1
% 830apv - apv baseline effect
% 831acsf - pre2
% 901acsf - wn on
% 901apv - apv for reversion



% Templates for reversal --- C hits C, D hits the high stack after C
    % Neither one does that great, but should be enough to speed up
    % reversal.  Start out with C, use D if it doesn't work.
    
tvallWN=[tv513B tv513C];
numsall2=[zeros(1,length([tv513B])) ones(1,length([tv513C]))];
[b2,ind2]=sort(tvallWN);
numsortedWN=numsall2(ind2);
tvsortedWN=tvallWN(ind2);

tvallWN=[tv513B tv513C];
numsall2=[zeros(1,length([tv513B])) ones(1,length([tv513C]))];
[b2,ind2]=sort(tvallWN);
numsortedWN=numsall2(ind2);
tvsortedWN=tvallWN(ind2);


% probability of B as a function of time
figure;hold on;
plot(runningaverage(tvsortedWN,40),runningaverage(numsortedWN,40),'*','Color','r')
plot(runningaverage(tvsortedPRE,40),runningaverage(numsortedPRE,40),'*','Color','b')


