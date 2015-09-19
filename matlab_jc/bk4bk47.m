% bk4bk47

% Claimed by TW 
% Potential syntax birds (also has high stacks for pitch shifting)
% 
% Started screening on 03.05.10
% high stack repeats (A) 
% low stacks (B) followed by high stacks (C)e
tvallPRE1=[tv307A tv307B];
numsall2=[zeros(1,length(tv307A)) ones(1,length(tv307B))];
[b2,ind2]=sort(tvallPRE1);
numsortedPRE1=numsall2(ind2);
tvsortedPRE1=tvallPRE1(ind2);

tvallPRE=[tv309preA tv309preB];
numsall2=[zeros(1,length(tv309preA)) ones(1,length(tv309preB))];
[b2,ind2]=sort(tvallPRE);
numsortedPRE=numsall2(ind2);
tvsortedPRE=tvallPRE(ind2);

tvallWN=[tv309wnA tv309wnB];
numsall2=[zeros(1,length(tv309wnA)) ones(1,length(tv309wnB))];
[b2,ind2]=sort(tvallWN);
numsortedWN=numsall2(ind2);
tvsortedWN=tvallWN(ind2);

tvallPOST=[tv310postA tv310postB];
numsall2=[zeros(1,length(tv310postA)) ones(1,length(tv310postB))];
[b2,ind2]=sort(tvallPOST);
numsortedPOST=numsall2(ind2);
tvsortedPOST=tvallPOST(ind2);


figure;plot(runningaverage(tvsortedPRE(1:213),50),runningaverage(numsortedPRE(1:213),50),'*','Color','b')
hold on;plot(runningaverage(tvsortedPRE(214:end),50),runningaverage(numsortedPRE(214:end),50),'*','Color','b')
plot(runningaverage(tvsortedPRE1,50),runningaverage(numsortedPRE1,50),'*','Color','b')
plot(runningaverage(tvsortedWN,30),runningaverage(numsortedWN,30),'*','Color','r')
plot(runningaverage(tvsortedPOST,30),runningaverage(numsortedPOST,30),'*','Color','k')
% 03.07 --
    % 40% B and 60% A


%% 03.09.10 - Experiment 1
    % Made template that does a birdtaf AND off of the note before A and a birdtaf NOT
    % off of B - Template308.mat
    % Hit A---> assay for increased p(B-C)
    % 3.09 --- Start WN at 9:45am
        % by 3:30pm it appears to have shifted, but interpretation is difficult
            % because of drifting baseline beforehand
        % WN OFF at 7:50pm - obvious and considerable learning
        % No singing between then and lights out (unfortunately)
    % 3.10 --- rapid forgetting
 %%
%  *** worth tweaking template to avoid hitting after B's ***

% 3.12 - cannulae implantation surgery
% 3.13 - began singing crappy song
% 3.15 - began singing good song
% 3.18 - plan a covert + control --- hit high stack after low stacks ('C')
    % template on at 5pm
% 3.19 - WN on at 11:20am
    % median = 3317 for notes during the morning
    % sigma*0.5=31
    % Hit below 3348Hz
    
    % 2:30pm - looking like learning
    
% 3.22 - lots of learning, gradual recovery


% DECEASED - I think this is the one that died from not being able to
% breathe while inserting probes


j=3;
figure;hold on;
plot(PC24(j).tvPRE,PC24(j).pitchPRE(round(median(PC24(j).targeting-32)),:),'*')
plot(PC24(j).tvWN,PC24(j).pitchWN(round(median(PC24(j).targeting-32)),:),'*','Color','r')
plot(PC24(j).tvPOST,PC24(j).pitchPOST(round(median(PC24(j).targeting-32)),:),'*','Color','k')    
    
    
    