bk1bk29
% Syntax lesion 
% Sept 2010 - TW did syntax learning screening - looks great

% 10.05.10 at dawn - planned to do syntax shift but p(B) is too high
    % Increase probability of transition to low stack note (B)
    % Starting at 81% B for 10.03
    % Current template looks bad
    
    
% 10.15.10 - REPEAT REDUCTION EXPERIMENT
        % /bk1bk29/1011_wnoff
        % label notes
        [timevals,distnt] = jcrepeatdist('a','batchnotes');
        % Tested 10.16 - looks good -minor problems with ambiguous first note in
        % the series of repeats

    % WN on 10.17 at dawn
    % baseline median is 4, baseline mean is 4.9
        % MinRepeat = 4 (consider lowering this to 3)
        % MaxRepeat = 15
        % RepReset = 0.2 
        % RepRefrac = 0.06 

        % TrigRefrac = 0.06
    % WN off 10.18 at 8:30pm

% REPEAT EXPANSION EXPERIMENT

    % WN on 10.22 at dawn
        % baseline mean is 5.1
        % MinRepeat = 1 
        % MaxRepeat = 4 (consider lowering this to 3)
        % RepReset = 0.2 
        % RepRefrac = 0.06 
        % TrigRefrac = 0.06

    % WN off 10.23 at 9pm



%[timevalsPost1019,distntPost1019] = jcrepeatdist('a','batchnotes');   

[a,b]=hist(distntPRE1016,0:1:20);
[c,d]=hist(distntWN1017(207:end),0:1:20);
[e,f]=hist(distntPost1019,0:1:20);
a=a/length(distntPRE1016);
c=c/length(distntWN1017(207:end));
e=e/length(distntPost1019);
figure;hold on;
subplot(211);hold on;
stairs(b,a,'Linewidth',3)
stairs(d,c,'color','r','linewidth',3)
%stairs(f,e,'color','k','linewidth',3)
xlim([0 20])
winwidth=50;
subplot(212);plot(runningaverage(timevalsPRE1016-min(timevalsPRE1016),winwidth),runningaverage(distntPRE1016,winwidth),'*') % 1015_evtafreptest
hold on;plot(runningaverage(timevalsWN1017(1:206)-min(timevalsPRE1016),winwidth),runningaverage(distntWN1017(1:206),winwidth),'*','Color','r') % 1017_wnonrepeats
hold on;plot(runningaverage(timevalsWN1017(207:end)-min(timevalsPRE1016),winwidth),runningaverage(distntWN1017(207:end),winwidth),'*','Color','r') % 1017_wnonrepeats
hold on;plot(runningaverage(timevalsPost1019-min(timevalsPRE1016),winwidth),runningaverage(distntPost1019,winwidth),'*','Color','k') % 1017_wnonrepeats