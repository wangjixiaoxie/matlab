                                                                                                                                                                              % pu3bk24
% syntax bird - AP5 in RA

% Bilateral RA implantation surgery on Friday December 4, 2009
% Started singing December 6, 2009
% Branch point is a sweep + very short stack (A)
    % Transition to sweep (B) + intro note (D) - ~85%
    % Transition to sweep (B) + high stack note (C) - ~15%
% segment at 1000

% 12.7.09 at 2:45pm - launchpad - template to hit intro note A-->BD
% 12.8.09 at 2pm - probes in
% 12.9.09 - started singing in the morning
%   ~52% D and 48% C
%       3:04pm - 4mMapv on at 1.0uL/min 
%              44% D
%              56% C
%       7:16pm - ACSF on at 1.5uL/min
%       7:30pm - 0.6uL/min
% 12.10.09 
%      1:20pm - 4mMapv on at 1.0uL/min
%      5:38pm - acsf on at 1.5uL/min
%      6:10pm - 0.6uL/min

% 12.13.09 morning - ACSF after WN - 77% C and 23% D
% 12.13.09 -AP5 on - 68% C and 32% D
    % looks like a reversion
% clog, probes out, didn't sing for a day or so
% 12.15.09 - 89% C and 11% D

% 12.16.09 - did control pitch shift experiment


% 3.02.10 - began screening again
% 3.03.10 - overcame recording issues, recording song
%       7:30pm - 59% C and 41% D (Is it stable - wait a day before WN)



tv1209pre=[tv1209preC tv1209preD];
numsall2=[zeros(1,length(tv1209preC)) ones(1,length(tv1209preD))];
[b2,ind2]=sort(tv1209pre);
numsorted1209pre=numsall2(ind2);
tvsorted1209pre=tv1209pre(ind2);

tv1209apv=[tv1209apvC tv1209apvD];
numsall2=[zeros(1,length(tv1209apvC)) ones(1,length(tv1209apvD))];
[b2,ind2]=sort(tv1209apv);
numsorted1209apv=numsall2(ind2);
tvsorted1209apv=tv1209apv(ind2);

tv1212=[tv1212C tv1212D];
numsall2=[zeros(1,length(tv1212C)) ones(1,length(tv1212D))];
[b2,ind2]=sort(tv1212);
numsorted1212=numsall2(ind2);
tvsorted1212=tv1212(ind2);

tv1213apv=[tv1213apvC tv1213apvD];
numsall2=[zeros(1,length(tv1213apvC)) ones(1,length(tv1213apvD))];
[b2,ind2]=sort(tv1213apv);
numsorted1213apv=numsall2(ind2);
tvsorted1213apv=tv1213apv(ind2);


figure;hold on;
% Pre - 1208_probein
    plot(runningaverage(tvsorted1209pre(1:120),50),runningaverage(numsorted1209pre(1:120),50),'*','Color','k')
    plot(runningaverage(tvsorted1209pre(121:end),50),runningaverage(numsorted1209pre(121:end),50),'*','Color','k')
% APV - 1209_4mMapv
    plot(runningaverage(tvsorted1209apv,50),runningaverage(numsorted1209apv,50),'*','Color','r')
% WN on - second half (1211wnon_newprobe)- TW has the rest of the data   
    plot(runningaverage(tvsorted1212-60,50),runningaverage(numsorted1212,50),'*','Color','k')
% APV - 1213_4mMAPV
    plot(runningaverage(tvsorted1213apv-60,50),runningaverage(numsorted1213apv,50),'*','Color','r')

% July 2010 summary analysis
figure;hold on;
subplot(211);
imagesc(t,f,log(avZ));syn;ylim([0,1e4]);xlim([-1.3 0.2])
%
tv1213apvC=timing3(fv1213apvC);
tv1213apvD=timing3(fv1213apvD);
tv1213apv=[tv1213apvC tv1213apvD];
numsall2=[zeros(1,length(tv1213apvC)) ones(1,length(tv1213apvD))];
[b2,ind2]=sort(tv1213apv);
numsorted1213apv=numsall2(ind2);
tvsorted1213apv=tv1213apv(ind2);

subplot(212);hold on;
% Pre
plot(runningaverage(tvsorted1209pre,60),runningaverage(1-numsorted1209pre,60),'*','Color','k')
plot(runningaverage(tvsorted1210pre,20),runningaverage(1-numsorted1210pre,20),'*','Color','k')
plot(runningaverage(tvsorted1210Bacsf,20),runningaverage(1-numsorted1210Bacsf,20),'*','Color','k')
plot(runningaverage(tvsorted1210apv,20),runningaverage(1-numsorted1210apv,20),'*','Color','g')
% WN on
plot(runningaverage(tvsorted1211acsf,15),runningaverage(1-numsorted1211acsf,15),'*','Color','r')
plot(runningaverage(tvsorted1213apv,15),runningaverage(1-numsorted1213apv,15),'*','Color','g')
ylim([0 1])
% FF variability (using non-catch trials too because non-targeted; note
% that apv file ignores first 29 labeled notes to allow time for apv to
% diffuse into brain)
mean(std(pitch1213acsfCallforFF(220:280,:)'))/mean(mean(pitch1213acsfCallforFF(220:280,:)'))
mean(std(pitch1213apvCallforFF(220:280,:)'))/mean(mean(pitch1213apvCallforFF(220:280,:)'))
mean(std(pitchpreC(220:280,:)'))/mean(mean(pitchpreC(220:280,:)'))
mean(std(pitch1210apvC(220:280,:)'))/mean(mean(pitch1210apvC(220:280,:)'))