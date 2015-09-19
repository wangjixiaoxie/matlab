bk91w60
% NatNeuro control candidate
% Stim candidate

% To Do:
%   - transfer EvTaf4

% 9.02.10 - screen
%   1:20pm - made template
%   How does he respond to WN? 
% Not singing much - wait til I get back to start WN

% 9.07.10
% singing
% problem with EvTaf4 - move to bigbird?
% checking template

% 9.08.10
    % check wn notch - looks PERFECT
    % hit below 2328Hz (median) at 8:35pm
% 9.09.10
    % hit below 2350Hz (median) at 6:18pm
    % hit below 2360Hz at 8:12pm
% 9.10.10
    % hit below 2340Hz at 2:40pm - learning has subsided?
            % slightly increased wn volume at 2:48pm
% 9.12.10 - wn off at some point - turned off by TW


% Look at song until end of 9.11.10
    % not much learning
% 9.14.10 - wn on at noon - hit below 2340Hz - louder volume - notched
    % 1:51pm - hit below 2360Hz
    % 2:42pm - hit below 2340Hz
    
    % 8:10pm - looks like a big decrease in hit rate - check tomorrow
% 9.15.10
    % 10:30am - hit below 2360Hz
% 9.16.10
    % 12:01pm - hit below 2380Hz
% 9.17.10
    % 1:50pm - WN off

%%%%%%%%% STIM
% 9.24.10
    % Surgery - implanted LMAN bilateral stim array - see notebook for details
% 9.27.10 - began singing again in morning

% 9.29.10 - plugged in
% 10.04.10 - song at 9:45am
% 10.06.10 - stim test
    % 12:30pm - ??R, ??L, 30ms delay, 90ms stim, 70uA on both stimulators
        % R and L were not going where the should be for the double-implant
        % connector - so I re-soldered new connections
        
    % 5pm - 23R, 14L, 30ms delay, 90ms stim, 70uA on both stimulators
    % HUGE deflection upwards starting slightly after syllable onset
% 10.07.10 - 10ms delay, 50uA -- to get effect earlier and less prominent
    % worked well
    % 12:50pm - changed config to 23R, 23L to test
    % 1:58pm - new tmp for WN

% 10.11.10 dawn - wn on hit below 2430Hz
    % 10ms delay, 90ms duration, 23R, 23L, 50uA
    % 10:30am - changed config to 24R, 24L
    % 11:10am - new stimulator (although old one probably worked fine) 
    % 11:10am - changed delay to 0ms
    % 5:31pm - changed template to have more stim hits
% 10.12.10
    % MISTAKE - I had stim at 80% - no learning, reduced singing
    % 3:45pm - wn off - stim at 20%
    % wn back on at dark - hit below 2430Hz
% 10.13.10 - WN day 1
    % looks like learning...
    % 6pm - increased stim amplitude from 50uA to 70uA to ensure effect
        % BIG stim effect
    % 6:45pm - improved WN template to hit more reliably
    % 6:55pm - catch trial fraction to 0.5 - amplitude to 50uA
    % 8:37pm - threshold to hit below 2470Hz
% 10.14.10 - WN day 2
    % 8:50am - catch trial fraction back to 0.8
    % 10:30am - threshold to hit below 2420Hz (it was hitting everything b/c new template is quite good)
    % 1pm - 2400Hz
    % 3:20pm - 2410Hz
    % 3:40pm - 2430Hz
    % 5:30pm - 2450Hz - catch trial frac to 0.6
    % 10pm - catch trial frac to 0.8
% 10.15.10 - WN day 3
    % 10:10am - catch trial frac to 0.7 (need to improve tmp)
    % 10:40am - 2480Hz
    % 11am - new stim template - should do better - frac to 0.8
    % 6:45pm - 2490Hz (conservative)
    
% 10.18.10 
    % wn off at 1:10pm
    
 % PROBLEM - stim effect on left side appears to be dying
% 10.19.10
    % 4:45pm - switched stimulators - thus stim effect on right side is now dying
    

%%%%%%%%
%%%%%%%%
% # July 2011
% 7.6.11 - started screening
% 7.7.11 - stim 100prct (accident), 20msdel, 30msdur, R34,L34, 40uA
    % 20ms del is bad (too late), 40uA may be too little
% 7.8.11 - stim 20prct, same params except 0ms del 
    % works better but not well
% 7.8.11 - try R14, L14 (0ms del, 30ms dur, 40uA)
% 7.8.11 - try R23, L23 (0ms del, 30ms dur, 40uA)
% 7.9.11 - Return to R34, L34 (best), keep 0ms del, 30ms dur, increase to 80uA
% Conclusion - no stim effect
    
    
    
    
% Baseline is 9.14 morning
% Learning is 9.14 evening and 9.15 morning - batchHitMiss
% fvBaseline=findwnoteJC('batchnotes','b','','',0,[6000 8100],8500,1,'obs0',0);
% for i=1:length(fvBaseline)
%     shiftedBaseline(i,:)=fvBaseline(i).datt;
% end
% pitchBaseline=jc_pitchmat1024(shiftedBaseline,1024,1020,1,6000,8100,[1],'obs0',1);
fvWNon=findwnoteNN('batchnotes','b','','',0,[6000 8100],8500,1,'obs0',0);
        recognized=[];
        hitrange=[];
        catchtrial=[];
        for i=1:length(fvWNon)
            shiftedWNon(i,:)=fvWNon(i).datt;
            recognized(i)=fvWNon(i).CATCH>-1;
            hitrange(i)=fvWNon(i).TEMPLATE==1;
            catchtrial(i)=fvWNon(i).CATCH==1;
        end
        indEscapeAbove=recognized & hitrange & catchtrial;
        indHitAbove=recognized & hitrange & ~catchtrial;
pitchWNon=jc_pitchmat1024(shiftedWNon,1024,1020,1,6000,8100,[1],'obs0',1);    
figure;plot(median(pitchWNon(:,indEscapeAbove)'))
hold on;plot(median(pitchWNon(:,indHitAbove)'),'r')
% PRE-PROCESSING
    % How much does pitch at t depend on pitch of the most recent 100
    % performances.
        % experiment2 --- [900 950]
        % experiment1 --- 
    ptwindow=[260 310];
    histwindow=100; 
    residWNon=[];
    for i=1:size(pitchWNon,2)
        residWNon(:,i)=pitchWNon(:,i)-mean(pitchBaseline')';
    end
    residWNon=residWNon(ptwindow,:); % only look at reasonable part

% REGRESSION
        b=zeros(size(pitchWNon,2)-histwindow,histwindow);
        clear bint
    % Rows of b are the syllable renditions being predicted ("predictees")
    % Columns of b are the syllable renditions doing the predicting ("predictors")
        % Final column is the most recent syllable rendition
        for i=histwindow+1:size(pitchWNon,2)
            [b(i-histwindow,1:histwindow),bint(i).data]=regress(residWNon(:,i),residWNon(:,i-histwindow:i-1));
            i
        end    
        regcoefs=b;
% POST-PROCESSING
    clear recentcoefs recentescapesabove recenthitsabove mrecentescapes mrecenthits mrecentabove
    for i=1:size(regcoefs,1)
        recentcoefs(i,:)=regcoefs(i,:);
        recentescapesabove(i,:)=indEscapeAbove(i:i+histwindow-1);
        recenthitsabove(i,:)=indHitAbove(i:i+histwindow-1);
    end


    for i=1:histwindow
        aind=find(recentescapesabove(:,i));
        if ~isempty(aind)
            holder=(recentcoefs(aind,i));
            mrecentescapes(i)=mean(holder(abs(holder)<1)); % Only look at the regression coefficients that are not outliers
        end
        bind=find(recenthitsabove(:,i));
        if ~isempty(bind)
            holder=(recentcoefs(bind,i));
            mrecenthits(i)=mean(holder(abs(holder)<1)); % Only look at the regression coefficients that are not outliers
        end
    %     abind=find(recentescapesabove(:,i)+recenthitsabove(:,i));
    %     if ~isempty(aind) | ~isempty(bind)
    %         mrecentabove(i)=median(recentcoefs(abind,i));
    %     end
    end
    mrecent=mean(recentcoefs);
  % Flip so that first point is most recent syllable.
  fmrecentescapes=fliplr(mrecentescapes);
  fmrecenthits=fliplr(mrecenthits);
  fmrecent=fliplr(mrecent);
% PLOT
figure;hold on;
plot(runningaverage(fmrecentescapes,10),'b')
plot(runningaverage(fmrecenthits,10),'r')
figure;hold on;
plot(fmrecentescapes,'b')
plot(fmrecenthits,'r')


% load /cardinal3/bk80w28/Analyzed910.mat
    % Beautiful!!!!!
        figure;plot(runningaverage(t(indEscapeAbove),50),runningaverage(mb(indEscapeAbove),50),'r')
        hold on;plot(runningaverage(t(indHitAbove),50),runningaverage(mb(indHitAbove),50),'b')


%%% 09.13.10
clear recentcoefs recentescapesabove recenthitsabove mrecentescapes mrecenthits mrecentabove
k=100;
        indEscapeAbove=recognized & hitrange & catchtrial;
        indHitAbove=recognized & hitrange & ~catchtrial;
pitchWNon=jc_pitchmat1024(shiftedWNon,1024,1020,1,6000,8100,[1],'obs0',1);    
figure;plot(median(pitchWNon(:,indEscapeAbove)'))
hold on;plot(median(pitchWNon(:,indHitAbove)'),'r')


% PRE-PROCESSING
    % How much does pitch at t depend on pitch of the most recent 100
    % performances.
    ptwindow=[850:1100];
    histwindow=100; 
    residWNon=[];
    for i=1:size(pitchWNon,2)
        residWNon(:,i)=pitchWNon(:,i)-mean(pitchBaseline')';
    end
    residWNon=residWNon(ptwindow,:); % only look at reasonable part

% REGRESSION
        b=zeros(size(pitchWNon,2)-histwindow,histwindow);
        clear bint
    % Rows of b are the syllable renditions being predicted ("predictees")
    % Columns of b are the syllable renditions doing the predicting ("predictors")
        % Final column is the most recent syllable rendition
        for i=histwindow+1:size(pitchWNon,2)
            [b(i-histwindow,1:histwindow),bint(i).data]=regress(residWNon(:,i),residWNon(:,i-histwindow:i-1));
            i
        end    
        regcoefs=b;
% POST-PROCESSING
    clear recentcoefs recentescapesabove recenthitsabove mrecentescapes mrecenthits mrecentabove
    for i=1:size(regcoefs,1)
        recentcoefs(i,:)=regcoefs(i,:);
        recentescapesabove(i,:)=indEscapeAbove(i:i+histwindow-1);
        recenthitsabove(i,:)=indHitAbove(i:i+histwindow-1);
    end


    for i=1:histwindow
        aind=find(recentescapesabove(:,i));
        if ~isempty(aind)
            holder=(recentcoefs(aind,i));
            mrecentescapes(i)=mean(holder(abs(holder)<1)); % Only look at the regression coefficients that are not outliers
        end
        bind=find(recenthitsabove(:,i));
        if ~isempty(bind)
            holder=(recentcoefs(bind,i));
            mrecenthits(i)=mean(holder(abs(holder)<1)); % Only look at the regression coefficients that are not outliers
        end
    %     abind=find(recentescapesabove(:,i)+recenthitsabove(:,i));
    %     if ~isempty(aind) | ~isempty(bind)
    %         mrecentabove(i)=median(recentcoefs(abind,i));
    %     end
    end
    mrecent=mean(recentcoefs);
  % Flip so that first point is most recent syllable.
  fmrecentescapes=fliplr(mrecentescapes);
  fmrecenthits=fliplr(mrecenthits);
  fmrecent=fliplr(mrecent);
% PLOT
figure;hold on;
plot(runningaverage(fmrecentescapes,10),'b')
plot(runningaverage(fmrecenthits,10),'r')
figure;hold on;
plot(fmrecentescapes,'b')
plot(fmrecenthits,'r')


% load /cardinal3/bk80w28/Analyzed910.mat
    % Beautiful!!!!!
        figure;plot(runningaverage(t(indEscapeAbove),50),runningaverage(mb(indEscapeAbove),50),'r')
        hold on;plot(runningaverage(t(indHitAbove),50),runningaverage(mb(indHitAbove),50),'b')


% load /cardinal3/bk80w28/Analyzed910.mat
    % Beautiful!!!!!
        figure;plot(runningaverage(t(indEscapeAbove),50),runningaverage(mb(indEscapeAbove),50),'r')
        hold on;plot(runningaverage(t(indHitAbove),50),runningaverage(mb(indHitAbove),50),'b')



    
    