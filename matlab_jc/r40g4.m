% r40g4
% LMAN implant
% DECEASED?





% Definite covert candidate --- But not implanted yet

% Possible syntax candidate--- NO
    % Bad for syntax - 100% WN experiments show that he won't change syntax
    % significantly for multiple days.
% Multiple variable repeats
% good singer

% Short stack repeats
% 2xlong stack followed by sweep repeats

% 03.10.10 @ 3pm - began screening on launchpad
% 5pm - tmptest for note A
    %     800ms refractory period is essential
    
% 03.11
    % new template -- BTaf AND off the note that typically precedes A
        % 1:10pm - WN on --- hit below 3285 (lowest 70th percent)
        % 2:50pm - hit below 3270 (hit rate was really high)\
% 03.12 - ampoff at 12:15pm
    % Significant learning appears during the morning
    
% Do a better + control - start Friday at noon
    % hit below mean + 0.5sigma
% 03.18 - 4pm - adjusted cntrng to improve hit rate.

% 03.19 - 11:20 - wn on A
    % median: 3220Hz
    % 0.5*sigma: 24Hz
    % WN hit below 3244Hz
    
    % 2:40pm - looks like learning
% 03.22 - Data322.mat has all learning and post data
    % 11:20 - wn off, lights off
    % 2:20 - lights on
    % Lots of learning, substantial but incomplete recovery

% 03.23 - 11:20 - wn on A
    % median (of pvs 24hrs): 3250Hz
    % median (of that morning): 3240Hz
    % 0.5*sigma=26Hz
    % hit above 3224 Hz
    % 6:30pm (change) -- hit above 3214Hz
% 03.24 - 11:20 - wn off A - lights out
%       -  2:50 - lights on (3.5hrs - slightly more than planned)
% ROBUST LEARNING, PARTIAL RECOVERY
    % as of 3.26 - recovery remains partial
% WN on syntax starting Saturday at dawn 
    % HIT the A's.
% 03.26 - WN on at 9:15pm --- sang <10 songs with WN before dark
% A - high stacks (WN)
% B - low stack often followed by another low stack and then downsweeps
% C - downsweeps
% 03.27
    % NO LEARNING AS OF NOON
    % NO LEARNING AS OF 4pm
% 03.28 --- clearly not a branch point, so it's reasonable to call this a
%           100% hit rate experiment (3.28 is day #2, do two more days)
% 03.29 --- sti ll on learning
% 03.30 at 9pm --- allow only the lowest 10% to escape --- hit above 3165
% 04.02 ---- SEE CODE BELOW - LOOKS GOOD

figure;plot(tvals100prct,pitch100prct(250,:),'*')
hold on;plot(tvals90prct,pitch90prct(250,:),'*','Color','r')
hold on;plot(tvalspre,pitchpre(250,:),'*','Color','k')
hold on;plot(timing3(fv325postA),pitch325postA(250,:),'*','Color','k')
hold on;plot(timing3(fv326postA),pitch326postA(250,:),'*','Color','k')







tvallPRE=[tv326A tv326B tv326C];
numsall2=[zeros(1,length([tv326A])) ones(1,length([tv326B])) ones(1,length([tv326C]))];
[b2,ind2]=sort(tvallPRE);
numsortedPRE=numsall2(ind2);
tvsortedPRE=tvallPRE(ind2);

tvallWN=[tv327wnA tv327wnB tv327wnC];
numsall2=[zeros(1,length([tv327wnA])) ones(1,length([tv327wnB])) ones(1,length([tv327wnC]))];
[b2,ind2]=sort(tvallWN);
numsortedWN=numsall2(ind2);
tvsortedWN=tvallWN(ind2);


% probability of not A
figure;hold on;
plot(runningaverage(tvsortedPRE,20),runningaverage(numsortedPRE,20),'*','Color','b')
plot(runningaverage(tvsortedWN,20),runningaverage(numsortedWN,20),'*','Color','r')
    
figure;hold on;
plot(tvals310preA,pitch310preA(260,:),'*')
plot(tvals311preA,pitch311preA(260,:),'*')
plot(tvals311wnA,pitch311wnA(260,:),'*','Color','r')
plot(tvals312postA,pitch312postA(260,:),'*','Color','k')

figure;hold on;
plot(tv323preA,pitch323preA(260,:),'*')
