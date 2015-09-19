% pu16bk67
% Screened 05.15.11 - data on /cardinal4

% Segment at 30,10,50,2
cleandir4('batch',100,500,6,10);

% syntax variability
% MMAN lesion candidate
% AB-->CD or AB-->E
% 5.15.2011
    % analyze recordings from 5.15 on wav file machine
    % 54% C, 46% E
    % LOOKS GREAT IN TERMS OF VARIABILITY AND TARGETABILITY
% 5.21.11
    % 59% C, 41% E
    % template to hit C - evtaf on C, birdtaf on B
    % template test at 1pm
% 5.22.11 - wn on at dawn to reduce p(C)
% 5.23.11
% 5.24.11
    % noon-end of day - 39% C, 61% E (i.e. lots of learning)
% 5.31.11
    % all day - 65% C, 35% E (i.e. complete reversion)
    
% plan to lesion Friday 6.3.11


tv520preC=timing4(fvalspreC);
tv520preE=timing4(fvalspreE);
tvall520pre=[tv520preC tv520preE];
numsall2=[zeros(1,length([tv520preC])) ones(1,length([tv520preE]))];
[b2,ind2]=sort(tvall520pre);
numsorted520pre=numsall2(ind2);
tvsorted520pre=tvall520pre(ind2);
%
tv522wnC=timing4(fvalswnC);
tv522wnE=timing4(fvalswnE);
tvall522wn=[tv522wnC tv522wnE];
numsall2=[zeros(1,length([tv522wnC])) ones(1,length([tv522wnE]))];
[b2,ind2]=sort(tvall522wn);
numsorted522wn=numsall2(ind2);
tvsorted522wn=tvall522wn(ind2);
%
tv531C=timing4(fvalswnC);
tv531E=timing4(fvalswnE);
tvall531=[tv531C tv531E];
numsall2=[zeros(1,length([tv531C])) ones(1,length([tv531E]))];
[b2,ind2]=sort(tvall531);
numsorted531=numsall2(ind2);
tvsorted531=tvall531(ind2);

% probability of B as a function of time
figure;hold on;
plot(runningaverage(tvsorted520pre,50),runningaverage(numsorted520pre,50),'*','Color','k')
plot(runningaverage(tvsorted522wn,50),runningaverage(numsorted522wn,50),'*','Color','r')    
plot(runningaverage(tvsorted531,50),runningaverage(numsorted531,50),'*','Color','k')        