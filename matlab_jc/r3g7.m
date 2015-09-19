% r3g7

% Possible syntax candidate

% 03.10.10 @ 3pm - began screening on launchpad
% 03.11 --- stable at ~50% A (high stacks) and 50% B (low stacks)
%       5:20pm - tested template - hits A decently
        % adjusted it to templaA2 and cntrngA2
        % But...much easier to hit B's - templaB and cntrngB
         % So... decided to hit B's instead
%       5:50pm - tested template - hits B perfectly

% 03.12 --- 10:05am - wn on B's (syntax)
%         -- 5:10pm - adjust template (templaB2) so hit both B's and hit
%                       earlier in the note
% 03.13 --- 4pm - 27% B's and 73% A's --- great success
%           --- WN off at 4:55pm
%
% 03.15 - complete recovery

% 03.16 - SURGERY - bilateral RA cannulae implantation
% 03.17 - singing again

% 03.19 - Covert positive control - C up
    % 11:20am - WN on
    % median: 2510
    % 0.5*sigma: 20
    % Hit below 2530Hz
    
    % 2:40pm - hitting properly

% 03.22 - robust learning, no recovery

tvallPRE=[tv311A tv311B];
numsall2=[zeros(1,length([tv311A])) ones(1,length([tv311B]))];
[b2,ind2]=sort(tvallPRE);
numsortedPRE=numsall2(ind2);
tvsortedPRE=tvallPRE(ind2);

tvallPRE1=[tv312preA tv312preB];
numsall2=[zeros(1,length([tv312preA])) ones(1,length([tv312preB]))];
[b2,ind2]=sort(tvallPRE1);
numsortedPRE1=numsall2(ind2);
tvsortedPRE1=tvallPRE1(ind2);

tvallWN=[tv312wnA tv312wnB];
numsall2=[zeros(1,length([tv312wnA])) ones(1,length([tv312wnB]))];
[b2,ind2]=sort(tvallWN);
numsortedWN=numsall2(ind2);
tvsortedWN=tvallWN(ind2);

tvallPOST=[tv315postA tv315postB];
numsall2=[zeros(1,length([tv315postA])) ones(1,length([tv315postB]))];
[b2,ind2]=sort(tvallPOST);
numsortedPOST=numsall2(ind2);
tvsortedPOST=tvallPOST(ind2);




% plot probability of B as a function of time
figure;hold on;
plot(runningaverage(tvsortedPRE(1:55),30),runningaverage(numsortedPRE(1:55),30),'*','Color','b')
plot(runningaverage(tvsortedPRE(56:end),40),runningaverage(numsortedPRE(56:end),40),'*','Color','b')
plot(runningaverage(tvsortedPRE1,40),runningaverage(numsortedPRE1,40),'*','Color','b')
plot(runningaverage(tvsortedWN(1:29),20),runningaverage(numsortedWN(1:29),20),'*','Color','r')
plot(runningaverage(tvsortedWN(30:end),40),runningaverage(numsortedWN(30:end),40),'*','Color','r')
plot(runningaverage(tvsortedPOST,30),runningaverage(numsortedPOST,30),'*','Color','k')
