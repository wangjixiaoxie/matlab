bk67bk44
% Syntax lesion candidate
% Sept 2010 - TW did syntax learning screening - looks great

% 9.15.10 - pitch learning screening
    % LOOKS LIKE CIRCADIAN
    % 1:36pm - WN on - hit above 2485Hz - 80th prctile for the same morning
        % accident - meant to hit below
    % 4:18pm (4:54pm computer time) - hit below 2500Hz - 75th percentile...
    
    
% TW - lesions, etc.

% 10.04.10 at 8:25pm - downward syntax shift
    % Reduce probability of transition to low stack note (B)
    % Starting at 62% B for 10.03
    % 59% for 10.04
    % Current template looks good (except around 4:30pm on 10.03 - 2 misses)
    % 10.06 morning - 45% B
    % 10.06 evening - 45% B
 % Plan to keep WN on   
 % 10.07.10 - looks like even more learning
    % 10.07 afternoon/evening - 43.5% B -- wn off
    
    
    
 % 10.15.10 - repeat experiment
 % /bk67bk44/1011_wnoff
% label notes
[timevals,distnt] = jcrepeatdist('a','batchnotes');
% Tested 10.17 - looks good -minor problems with ambiguous first note in
    % the series of repeats
    % baseline median is 6, baseline mean is 6.6
        % MinRepeat = 5 (consider lowering this to 4)
        % MaxRepeat = 20
        % RepReset = 0.2 
        % RepRefrac = 0.08 
        % TrigRefrac = 0.08
% 10.18.10 - dawn - wn on
    % looks good (still worth considering lowering min to 4)
    % lowered MIN to 4 at 11am
% 10.19.10 - 9:30pm - wn off

% 10.22.10 - PLAN REPEAT EXPANSION
% 2:30pm - wn on for 30min to test strategy
% baseline mean is 7.1
% pvsly hit all 4 or more
% now hit all 4 or less (5 or less?)
% learned to increase to 8.7
[timevalsPre1022,distntPre1022] = jcrepeatdist('a','batchnotes2');  

figure;hold on;
plot(runningaverage(timevalspre1022,50),runningaverage(distntpre1022,50),'*')
plot(runningaverage(timevalswn1023,50),runningaverage(distntwn1023,50),'*','Color','r')
plot(runningaverage(timevalspost1024,20),runningaverage(distntpost1024,20),'*','Color','b')


[timevalsPost1020,distntPost1020] = jcrepeatdist('a','batchnotes');  

 [a,b]=hist(distntPRE1017,0:1:20);
[c,d]=hist(distntWN1018(190:end),0:1:20);
%[e,f]=hist(distntPost1020,0:1:20);
a=a/length(distntPRE1017);
c=c/length(distntWN1018(190:end));
%e=e/length(distntPost1019);





%%%%%% Make figures

clear all
load /cardinal3/SyntaxBirds/bk67bk44/RepeatData1018.mat
figure;hold on;
subplot(322);hold on;
stairs(b,a,'Linewidth',3)
stairs(d,c,'color','r','linewidth',3)
%stairs(f,e,'color','k','linewidth',3)
xlim([0 20])
subplot(324);hold on;
stairs(b,c-a,'color','r','Linewidth',3)
plot([0 20],[0 0],'k')
winwidth=50;
subplot(326);plot(runningaverage(timevalsPRE1017-min(timevalsPRE1017),winwidth),runningaverage(distntPRE1017,winwidth),'*') % 1015_evtafreptest
hold on;plot(runningaverage(timevalsWN1018(1:189)-min(timevalsPRE1017),winwidth),runningaverage(distntWN1018(1:189),winwidth),'*','Color','r') % 1017_wnonrepeats
hold on;plot(runningaverage(timevalsWN1018(190:end)-min(timevalsPRE1017),winwidth),runningaverage(distntWN1018(190:end),winwidth),'*','Color','r') % 1017_wnonrepeats
hold on;plot(runningaverage(timevalsPost1020-min(timevalsPRE1017),winwidth),runningaverage(distntPost1020,winwidth),'*','Color','k') % 1017_wnonrepeats   
plot([0 100],[mean(distntPRE1017) mean(distntPRE1017)])
%%%%

%%% load bk1bk29 data
clear all
load /cardinal3/SyntaxBirds/bk1bk29/Repeats1017.mat
subplot(321);hold on;
stairs(b,a,'Linewidth',3)
stairs(d,c,'color','r','linewidth',3)
%stairs(f,e,'color','k','linewidth',3)
xlim([0 20])
subplot(323);hold on;
stairs(b,c-a,'color','r','Linewidth',3)
plot([0 20],[0 0],'k')
winwidth=50;
subplot(325);plot(runningaverage(timevalsPRE1016-min(timevalsPRE1016),winwidth),runningaverage(distntPRE1016,winwidth),'*') % 1015_evtafreptest
hold on;plot(runningaverage(timevalsWN1017(1:206)-min(timevalsPRE1016),winwidth),runningaverage(distntWN1017(1:206),winwidth),'*','Color','r') % 1017_wnonrepeats
hold on;plot(runningaverage(timevalsWN1017(207:end)-min(timevalsPRE1016),winwidth),runningaverage(distntWN1017(207:end),winwidth),'*','Color','r') % 1017_wnonrepeats
hold on;plot(runningaverage(timevalsPost1019-min(timevalsPRE1016),winwidth),runningaverage(distntPost1019,winwidth),'*','Color','k') % 1017_wnonrepeats
plot([0 100],[mean(distntPRE1016) mean(distntPRE1016)])