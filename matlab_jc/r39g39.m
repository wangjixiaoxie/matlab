% r39g39
% 100d old on 05.07.10

% screened 06.01.10

% syntax is complicated and there is some history dependence
% TRY THE FOLLOWING:
    % branch point as short sweepy harmonic note
        % A - intro + at least one long low stack - 25% 
        % B - one high stack + downsweep          - 56%
        % C - more than one high stack            - 19%

% 06.02.10 at dusk - hit B
% Definitely  learning but not a lot


% 10.04.10 at dusk - hit B/C - loud wn
% Problem with WN background noise

% 10.06.10 at dawn - hit B/C with loud WN
% Doesn't look like possible to get learning in one day...
% 7pm - catch trial ratio increased to 0.5
    % There might be learning
    % BETTER STRATEGY - only hit the A in one context (e.g. allow multiple
    % high stacks to escape but hit the high stack--> weird sweep)



tvallWN=[tv603A tv603B tv603C];
numsall2=[ones(1,length([tv603A])) zeros(1,length([tv603B])) ones(1,length([tv603C]))];
[b2,ind2]=sort(tvallWN);
numsortedWN=numsall2(ind2);
tvsortedWN=tvallWN(ind2);

tvallPRE=[tv602preA tv602preB tv602preC];
numsall2=[ones(1,length([tv602preA])) zeros(1,length([tv602preB])) ones(1,length([tv602preC]))];
[b2,ind2]=sort(tvallPRE);
numsortedPRE=numsall2(ind2);
tvsortedPRE=tvallPRE(ind2);

% probability of A&C as a function of time (i.e. not B)
figure;hold on;
plot(runningaverage(tvsortedWN,60),runningaverage(numsortedWN,60),'*','Color','r')
plot(runningaverage(tvsortedPRE,60),runningaverage(numsortedPRE,60),'*','Color','b')

