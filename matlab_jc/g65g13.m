
% r37g7

% Refrac=0.2

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
j=4;
figure;hold on;
plot(PC24(j).tvPRE,PC24(j).pitchPRE(round(median(PC24(j).targeting-32)),:),'*')
plot(PC24(j).tvWN,PC24(j).pitchWN(round(median(PC24(j).targeting-32)),:),'*','Color','r')
plot(PC24(j).tvPOST,PC24(j).pitchPOST(round(median(PC24(j).targeting-32)),:),'*','Color','k')    
    
    
    

% 03.22 - probes in
% 03.23 - clog, probes out
%       pm - star
ted singing again - song looks good
% 03.24 - probes in
% 03.25 - started singing in the morning
    % branch point intact
% 03.26 - made new apv batch
    % 5mM apv on at 3:50pm at 1.0uL/min
    % made new templaA3 - does a good job
  % WN on at 9:15 --- sang 3 songs and received 3 bursts of WN before dark
% 03.27 - WN on at dawn (also 45min the night before)
%   12:30pm - adjusted template (missing occassionally becuase of BT NOT)
            % cntrngA3 --> cntrngA4
            % Reduced NOT MAX from 20 to 9
            % Reduced AND TH from 2 to 2.5
%   3:30pm - wn off because not singing enough to learn
%  03.28 WN on at dawn
%   noon - looking good
%   not much learning
%   5:45 - apv on
%   6pm - clog, removed probes

cleandir4('batch',10000,500,6,10);
% 04.12 - probes in
% 04.13 - singing
% probabilities somewhat off - more A's - wait until probabilities either
    % return to normal or are stable
% 04.15 - start at dawn
% 04.15 - 8:40pm - 55% B (from 45% B) adjust config to hit A's more reliably
    % (5) BTAF AND - TH from 2.5 to 3 (more permissive)
    % (4) BTAF NOT - TH from 2 to 1.3 (more permissive)
    % (4) BTAF NOT - MAX from 30 to 50 (more strict - avoid hitting 'D' after BC)
% 04.16 - noon 
    % (5) BTAF AND - MIN from 5 to 1
    % (4) BTAF NOT - TH from 1.3 to 1
% 04.16 - 5:12pm chick time (4:50pm real time) - apv on at 0.5uL/min
    % 7:00 (5) BTAF AND - TH from 3 to 2 to avoid restarting on A
% 04.17 - 11am - reversed
    % 2pm - new template on B
    % 5:56pm chick time - apv on at 1.5uL/min (made folder 10min later)
% 04.18 - noon - adjusted template
    % 6:08pm - 55% A, 45% B - apv on at 1.5uL/min
% 04.19
    % 06:05pm chick time - apv on at 1.5uL/min
    % REVERSE - contingency reversed at 9pm - no singing afterwards
% 04.20
    % 06:16pm chick time - apv on at 1.5uL/min
% 04.21
    % 06:08pm chick time - apv on at 1.5uL/min
% 04.22 - ATTEMPT TO BLOCK REVERSION
    % 07:45am real time - apv on at 1.5uL/min
    % 08:00am real time - lights on
    % 08:45am real time - switch template to hit B's
    % 02:22am real time - switch template to hit A's
% 04.24 - Tim - attempt to block pitch learning - didn't sing much
    % NOTE C
    % 07:45am - 5mMapv on at 1.5uL/min
    % 09:40am - wn on - hit below median+0.5*(acsf sigma)
    % Bird didn't sing much - only 93 hits - stopped singing entirely @ 4pm
  % It appears that there was still apv from 4.22 or something
% 04.26 - 8:45am - wn on note C - reversion experiment - hit below 2220Hz
    % 12:15pm - hit below 2250Hz
    % definite learning but not a lot - 4:50pm - for 1-4, MIN=2 (from 1), TH=2.5 (from 2)
    % 6:05pm real time - 5mM apv on at 1.5uL/min
% 04.27 - adjusted template upward throughout the day
    % - lots of learning by 4pm
    % apv on around 5:20 at 1.5uL/min
    % apv off around 7:45
% 04.28
    % sing 20 songs
    % apv 5mM on at 1.5uL/min around 9am
    % acsf on around 11am
    % apv 5mM on at 1.5uL/min at 5:58pm
    %
% 04.30 - probes out
% 05.04 - probes in
% 05.06 - 5mM apv on at 1.5uL/min at 2:55 chick time to test baseline effect
    % ACSF back on at 1.5uL/min at 5:40 chick time
 % WN on at 8pm
% 05.07 - 4pm - adjusted template (effectively removed 3rd synshift) to
    % reduce triggering on note 'D' - also increase WN volume
    % cntrng507A
% 05.08 - noon - not learning enough - trigger with MIN=1 instead of MIN=2
    % cntrng508A
% 05.09 - more adjustments - increased wn volume

figure;hold on;
plot(tv506preB(62:end),(pitch506preB(220,62:end)),'*','Color','b')
plot(tv506apvB,(pitch506apvB(220,:)),'*','Color','r')
plot(tv509preB,(pitch509preB(220,:)),'*','Color','b')
plot(tv510wnB,(pitch510wnB(220,:)),'*','Color','k')
plot(tv511B,(pitch511B(220,:)),'*','Color','b')
% 05.10 - no learning ---> use bird for covert
    % 11:30 - wn off
    % noon chick time - wn on B hit above 2590Hz = median-0.5*sigma (shift down)
    % 1:40pm chick time - wn on first B (error had it on second B) - above 2600Hz
    % 1:57pm - better template
% 05.11 - no learning ---> try another note
    % 11:15 chick time - wn on - hit A below 2660Hz = median+0.5*sigma
    
    
    
    
    
    
figure;hold on;
plot(tv424apvC,pitch424apvC(210,:),'*','Color','r')
plot(tv422apvC3,pitch422apvC3(210,:),'*','Color','r')    
plot(tv422apvC1,pitch422apvC1(210,:),'*','Color','r')    
plot(tv422apvC2,pitch422apvC2(210,:),'*','Color','r')    
plot(tv421apvC,pitch421apvC(210,:),'*','Color','r')    
plot(tv425Catch,pitch425Catch(210,:),'*','Color','k')
plot(tv426Catch,pitch426Catch(210,:),'*','Color','b')
plot(tv427apvC,pitch427apvC(210,:),'*','Color','r')
plot(tv427wnC,pitch427wnC(210,:),'*','Color','b')
plot(tv428AMapvC(30:end),pitch428AMapvC(210,30:end),'*','Color','r')
plot(tv428AMacsfC,pitch428AMacsfC(210,:),'*','Color','b')
plot(tv423Cx,pitch423Cx(210,:),'*','Color','k')
plot(tv423C1,pitch423C1(210,:),'*','Color','k')
plot(tv423C2,pitch423C2(210,:),'*','Color','k')
plot(tv423C3,pitch423C3(210,:),'*','Color','k')
plot(tv423C4,pitch423C4(210,:),'*','Color','k')
plot(tv428C2,pitch428C2(210,:),'*','Color','b')
plot(tv426apvC,pitch426apvC(210,:),'*','Color','r')
plot(tv427apv2,pitch427apv2(210,:),'*','Color','k')

figure;plot(median(pitch427wnC(:,100:end)')-median(pitch425Catch'))
hold on;plot(median(pitch427apvC')-median(pitch425Catch'),'r')
hold on;plot(median(pitch424apvCatch')-median(pitch425Catch'),'r')
hold on;plot([0 500],[0 0],'k')
xlim([180 300])
hold on;plot(median(pitch427apvC')-median(pitch424apvCatch'),'g')    
hold on;plot(median(pitch428AMapvC')-median(pitch424apvCatch'),'g')  
hold on;plot(median(pitch428AMacsfC')-median(pitch425Catch'),'k')    
%%%%%%%
%%%%%%%%5
%%%%%%%
tvall413am=[tv413amA tv413amB];
numsall2=[zeros(1,length([tv413amA])) ones(1,length([tv413amB]))];
[b2,ind2]=sort(tvall413am);
numsorted413am=numsall2(ind2);
tvsorted413am=tvall413am(ind2);
    
    
tvall413=[tv413A tv413B];
numsall2=[zeros(1,length([tv413A])) ones(1,length([tv413B]))];
[b2,ind2]=sort(tvall413);
numsorted413=numsall2(ind2);
tvsorted413=tvall413(ind2);

tvall414=[tv414A tv414pmA tv414B tv414pmB];
numsall2=[zeros(1,length([tv414A tv414pmA])) ones(1,length([tv414B tv414pmB]))];
[b2,ind2]=sort(tvall414);
numsorted414=numsall2(ind2);
tvsorted414=tvall414(ind2);

tvall416wn=[tv416wnA tv416wnB];
numsall2=[zeros(1,length([tv416wnA])) ones(1,length([tv416wnB]))];
[b2,ind2]=sort(tvall416wn);
numsorted416wn=numsall2(ind2);
tvsorted416wn=tvall416wn(ind2);

tvall415wn4=[tv415wn4A tv415wn4B];
numsall2=[zeros(1,length([tv415wn4A])) ones(1,length([tv415wn4B]))];
[b2,ind2]=sort(tvall415wn4);
numsorted415wn4=numsall2(ind2);
tvsorted415wn4=tvall415wn4(ind2);

tvall416apv=[tv416apvA tv416apvB];
numsall2=[zeros(1,length([tv416apvA])) ones(1,length([tv416apvB]))];
[b2,ind2]=sort(tvall416apv);
numsorted416apv=numsall2(ind2);
tvsorted416apv=tvall416apv(ind2);

tvall417apv=[tv417apvA tv417apvB];
numsall2=[zeros(1,length([tv417apvA])) ones(1,length([tv417apvB]))];
[b2,ind2]=sort(tvall417apv);
numsorted417apv=numsall2(ind2);
tvsorted417apv=tvall417apv(ind2);


tvall418apv=[tv418apvA tv418apvB];
numsall2=[zeros(1,length([tv418apvA])) ones(1,length([tv418apvB]))];
[b2,ind2]=sort(tvall418apv);
numsorted418apv=numsall2(ind2);
tvsorted418apv=tvall418apv(ind2);

tvall416wn=[tv416wnA tv416wnB];
numsall2=[zeros(1,length([tv416wnA])) ones(1,length([tv416wnB]))];
[b2,ind2]=sort(tvall416wn);
numsorted416wn=numsall2(ind2);
tvsorted416wn=tvall416wn(ind2);

tvall417=[tv417A tv417pmA tv417B tv417pmB];
numsall2=[zeros(1,length([tv417A tv417pmA])) ones(1,length([tv417B tv417pmB]))];
[b2,ind2]=sort(tvall417);
numsorted417=numsall2(ind2);
tvsorted417=tvall417(ind2);
% 
% tvall417pm=[tv417pmA tv417pmB];
% numsall2=[zeros(1,length([tv417pmA])) ones(1,length([tv417pmB]))];
% [b2,ind2]=sort(tvall417pm);
% numsorted417pm=numsall2(ind2);
% tvsorted417pm=tvall417pm(ind2);


tvall418am=[tv418amA tv418amB];
numsall2=[zeros(1,length([tv418amA])) ones(1,length([tv418amB]))];
[b2,ind2]=sort(tvall418am);
numsorted418am=numsall2(ind2);
tvsorted418am=tvall418am(ind2);

tvall418pm=[tv418pmA tv418pmB];
numsall2=[zeros(1,length([tv418pmA])) ones(1,length([tv418pmB]))];
[b2,ind2]=sort(tvall418pm);
numsorted418pm=numsall2(ind2);
tvsorted418pm=tvall418pm(ind2);

tvall419am=[tv419amA tv419amB];
numsall2=[zeros(1,length([tv419amA])) ones(1,length([tv419amB]))];
[b2,ind2]=sort(tvall419am);
numsorted419am=numsall2(ind2);
tvsorted419am=tvall419am(ind2);

tvall419pm=[tv419pmA tv419pmB];
numsall2=[zeros(1,length([tv419pmA])) ones(1,length([tv419pmB]))];
[b2,ind2]=sort(tvall419pm);
numsorted419pm=numsall2(ind2);
tvsorted419pm=tvall419pm(ind2);

tvall419apv=[tv419apvA tv419apvB];
numsall2=[zeros(1,length([tv419apvA])) ones(1,length([tv419apvB]))];
[b2,ind2]=sort(tvall419apv);
numsorted419apv=numsall2(ind2);
tvsorted419apv=tvall419apv(ind2);


tvall420am=[tv420amA tv420amB];
numsall2=[zeros(1,length([tv420amA])) ones(1,length([tv420amB]))];
[b2,ind2]=sort(tvall420am);
numsorted420am=numsall2(ind2);
tvsorted420am=tvall420am(ind2);

tvall420pm=[tv420pmA tv420pmB];
numsall2=[zeros(1,length([tv420pmA])) ones(1,length([tv420pmB]))];
[b2,ind2]=sort(tvall420pm);
numsorted420pm=numsall2(ind2);
tvsorted420pm=tvall420pm(ind2);

tvall420apv=[tv420apvA tv420apvB];
numsall2=[zeros(1,length([tv420apvA])) ones(1,length([tv420apvB]))];
[b2,ind2]=sort(tvall420apv);
numsorted420apv=numsall2(ind2);
tvsorted420apv=tvall420apv(ind2);

tvall421pm=[tv421pmA tv421pmB];
numsall2=[zeros(1,length([tv421pmA])) ones(1,length([tv421pmB]))];
[b2,ind2]=sort(tvall421pm);
numsorted421pm=numsall2(ind2);
tvsorted421pm=tvall421pm(ind2);

tvall421apv=[tv421apvA tv421apvB];
numsall2=[zeros(1,length([tv421apvA])) ones(1,length([tv421apvB]))];
[b2,ind2]=sort(tvall421apv);
numsorted421apv=numsall2(ind2);
tvsorted421apv=tvall421apv(ind2);

tvall422pre=[tv422preA tv422preB];
numsall2=[zeros(1,length([tv422preA])) ones(1,length([tv422preB]))];
[b2,ind2]=sort(tvall422pre);
numsorted422pre=numsall2(ind2);
tvsorted422pre=tvall422pre(ind2);

tvall422reverse=[tv422reverseA tv422reverseB];
numsall2=[zeros(1,length([tv422reverseA])) ones(1,length([tv422reverseB]))];
[b2,ind2]=sort(tvall422reverse);
numsorted422reverse=numsall2(ind2);
tvsorted422reverse=tvall422reverse(ind2);

tvall422revrev=[tv422revrevA tv422revrevB];
numsall2=[zeros(1,length([tv422revrevA])) ones(1,length([tv422revrevB]))];
[b2,ind2]=sort(tvall422revrev);
numsorted422revrev=numsall2(ind2);
tvsorted422revrev=tvall422revrev(ind2);

tvall423=[tv423A tv423B];
numsall2=[zeros(1,length([tv423A])) ones(1,length([tv423B]))];
[b2,ind2]=sort(tvall423);
numsorted423=numsall2(ind2);
tvsorted423=tvall423(ind2);

tvall423wnoff=[tv423wnoffA tv423wnoffB];
numsall2=[zeros(1,length([tv423wnoffA])) ones(1,length([tv423wnoffB]))];
[b2,ind2]=sort(tvall423wnoff);
numsorted423wnoff=numsall2(ind2);
tvsorted423wnoff=tvall423wnoff(ind2);

% probability of B as a function of time
figure;hold on;
plot(runningaverage(tvsorted413am,50),runningaverage(numsorted413am,50),'*','Color','k')
plot(runningaverage(tvsorted413,50),runningaverage(numsorted413,50),'*','Color','k')
plot(runningaverage(tvsorted414,50),runningaverage(numsorted414,50),'*','Color','k')
plot(runningaverage(tvsorted416wn,50),runningaverage(numsorted416wn,50),'*','Color','b')
plot(runningaverage(tvsorted415wn4,50),runningaverage(numsorted415wn4,50),'*','Color','b')
plot(runningaverage(tvsorted416apv,50),runningaverage(numsorted416apv,50),'*','Color','r')
plot(runningaverage(tvsorted417apv,50),runningaverage(numsorted417apv,50),'*','Color','r')
plot(runningaverage(tvsorted418apv,50),runningaverage(numsorted418apv,50),'*','Color','r')
plot(runningaverage(tvsorted417,50),runningaverage(numsorted417,50),'*','Color','b')
% plot(runningaverage(tvsorted417pm,50),runningaverage(numsorted417pm,50),'*','Color','b')
plot(runningaverage(tvsorted418am,50),runningaverage(numsorted418am,50),'*','Color','b')
plot(runningaverage(tvsorted418pm,50),runningaverage(numsorted418pm,50),'*','Color','b')
plot(runningaverage(tvsorted417,50),runningaverage(numsorted417,50),'*','Color','b')
plot(runningaverage(tvsorted419am,50),runningaverage(numsorted419am,50),'*','Color','b')
plot(runningaverage(tvsorted419pm,50),runningaverage(numsorted419pm,50),'*','Color','b')
plot(runningaverage(tvsorted419apv,50),runningaverage(numsorted419apv,50),'*','Color','r')
plot(runningaverage(tvsorted420am,50),runningaverage(numsorted420am,50),'*','Color','b')
plot(runningaverage(tvsorted420pm,50),runningaverage(numsorted420pm,50),'*','Color','b')
plot(runningaverage(tvsorted420apv,30),runningaverage(numsorted420apv,30),'*','Color','r')
plot(runningaverage(tvsorted421pm,50),runningaverage(numsorted421pm,50),'*','Color','b')
plot(runningaverage(tvsorted421apv,30),runningaverage(numsorted421apv,30),'*','Color','r')
plot(runningaverage(tvsorted422pre,30),runningaverage(numsorted422pre,30),'*','Color','b')
plot(runningaverage(tvsorted422reverse,50),runningaverage(numsorted422reverse,50),'*','Color','r')
plot(runningaverage(tvsorted422revrev,50),runningaverage(numsorted422revrev,50),'*','Color','k')
plot(runningaverage(tvsorted423,20),runningaverage(numsorted423,20),'*','Color','b')
plot(runningaverage(tvsorted423wnoff,30),runningaverage(numsorted423wnoff,30),'*','Color','k')


figure;hold on;
plot(tv414B,pitch414B(300,:),'*')
plot(tv414pmB,pitch414pmB(300,:),'*')
plot(tv415wn4B,pitch415wn4B(300,:),'*')
plot(tv416amB,pitch416amB(300,:),'*')
plot(tv416apvB,pitch416apvB(300,:),'*','Color','r')
plot(tv417B,pitch417B(300,:),'*')
plot(tv418apvB,pitch418apvB(300,:),'*','Color','r')
plot(tv419apvB,pitch419apvB(300,:),'*','Color','r')
plot(tv422revrevB,pitch422revrevB(300,:),'*','Color','r')
plot(tv423B,pitch423B(300,:),'*')
plot(tv423wnoffB,pitch423wnoffB(300,:),'*')




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

% probes in
tvall326=[tv326A tv326B];
numsall2=[zeros(1,length([tv326A])) ones(1,length([tv326B]))];
[b2,ind2]=sort(tvall326);
numsorted326=numsall2(ind2);
tvsorted326=tvall326(ind2);

tvall326apv=[tv326apvA tv326apvB];
numsall2=[zeros(1,length([tv326apvA])) ones(1,length([tv326apvB]))];
[b2,ind2]=sort(tvall326apv);
numsorted326apv=numsall2(ind2);
tvsorted326apv=tvall326apv(ind2);

tvallWN=[tv328wnA tv328wnB];
numsall2=[zeros(1,length([tv328wnA])) ones(1,length([tv328wnB]))];
[b2,ind2]=sort(tvallWN);
numsortedWN=numsall2(ind2);
tvsortedWN=tvallWN(ind2);



% probability of B as a function of time
figure;hold on;
plot(runningaverage(tvsorted326apv,30),runningaverage(numsorted326apv,30),'*','Color','r')
plot(runningaverage(tvsorted326,30),runningaverage(numsorted326,30),'*','Color','k')
plot(runningaverage(tvsortedWN,50),runningaverage(numsortedWN,50),'*','Color','r')



