% bk80w28
% 09.07 - wn test for notch - too loud
%   template test - first threshold crossing
%   hit below 2285Hz at 8pm
% 09.08 
%   hit below 2300Hz at 10:06am
%   hit below 2320Hz at 12:55pm
%   hit below 2335Hz at 3:25pm
%   hit below 2350Hz at 4:20pm
%   hit below 2365Hz at 5pm
%   hit below 2380Hz at 8:15pm
%  looks like he learned a lot in the morning and then it has leveled off
% 09.09
%   hit below 2420Hz
% 09.10
%   WN off at 12:25pm - allow recovery and then repeat the same experiment

% 09.12 - median for entire day is 2300Hz
%   10:30pm - wn on, hit below 2300Hz
% 09.13
%   10:38am - hit below 2330Hz
%   10:46am - hit below 2380Hz
%   11:05am - adjusted template slightly to recognize BTAF better
%   12:23pm - hit below 2400Hz
%   1:45pm - hit below 2420Hz
% Labeled through 4:20pm
%   6:43pm - wn off

% Plan -tues- 09.14 - decay
%      -weds- 09.15 - decay
%      -thrs- 09.16 - WN on
%      -fri-  09.17 - WN on - asymptote
%      -sat-  09.18 - WN on - asymptote
%     ------decay or asymptote at Asilomar



% Look at averages on both exps:
figure;plot((median(pitchWNon(600:end,end-100:end)')-median(pitchBaseline(600:end,:)'))/max((median(pitchWNon(850:1100,end-100:end)')-median(pitchBaseline(850:1100,:)'))),'k')
hold on;plot((median(pitchWNon(600:end,indEscapeAbove)')-median(pitchBaseline(600:end,:)'))/max(median(pitchWNon(850:1100,indEscapeAbove)')-median(pitchBaseline(850:1100,:)')),'r')


% Experiment 1
    % baseline: /cardinal3/bk80w28/0903_baseline/batchnotes
        % a random sampling of songs from 9.07.10
    % wn: /cardinal3/bk80w28/0907_wnon/batchnotes
        % every song in the folder - starts at dark on 9.07.10 and ends at 10am on 9.09.10
% Experiment 2
    % baseline: /cardinal3/bk80w28/0910_wnoff/batchnotes
        % a random sampling of songs from 9.12.10
    % wn: /cardinal3/bk80w28/0913_wnon/batchnotes
        % every song in the folder - starts at dark on 9.12.10 and ends at 4:12pm on 9.13.10
% Experiment 3
    % baseline: /cardinal3/bk80w28/0915_tmptest/???
    % wn /cardinal3/bk80w28/0916_wnon
        % starts at dark on 9.15.10 and ends at 10am on 9.17.10
        
% Experiment 1
        ptwindow=[350:400];
        fvBaseline=findwnoteJC('batchnotes','b','','',0,[6000 8100],8500,1,'obs0',0);
        for i=1:length(fvBaseline)
            shiftedBaseline(i,:)=fvBaseline(i).datt;
        end
        pitchBaseline=jc_pitchmat1024(shiftedBaseline,1024,1020,1,6000,8100,[1],'obs0',1);    
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

% Experiment 2 - adjusted because hits between 'a' and 'b' make it
% difficult to accurately estimate the onset of 'b'.
        ptwindow=[840:920]; % or 850-1100 works
        fvBaseline=findwnoteJC('batchnotes','a','-','b',0,[6000 8100],8500,1,'obs0',0);
        for i=1:length(fvBaseline)
            shiftedBaseline(i,:)=fvBaseline(i).datt;
        end
        pitchBaseline=jc_pitchmat1024(shiftedBaseline,1024,1020,1,6000,8100,[1],'obs0',1);
        fvWNon=findwnoteNN('batchnotes','-','','b',0,[6000 8100],8500,1,'obs0',0);
        fvWNonB=findwnoteNN('batchnotes','b','','',0,[6000 8100],8500,1,'obs0',0);
                recognized=[];
                hitrange=[];
                catchtrial=[];
                for i=1:length(fvWNon)
                    shiftedWNon(i,:)=fvWNon(i).datt;
                    recognized(i)=fvWNonB(i).CATCH>-1;
                    hitrange(i)=fvWNonB(i).TEMPLATE==1;
                    catchtrial(i)=fvWNonB(i).CATCH==1;
                end
                indEscapeAbove=recognized & hitrange & catchtrial;
                indHitAbove=recognized & hitrange & ~catchtrial;
        pitchWNon=jc_pitchmat1024(shiftedWNon,1024,1020,1,6000,8100,[1],'obs0',1);    
        figure;plot(median(pitchWNon(:,indEscapeAbove)'))
        hold on;plot(median(pitchWNon(:,indHitAbove)'),'r')


% PRE-PROCESSING
    % How much does pitch at t depend on pitch of the most recent 100
    % performances.
    histwindow=100; 
    residWNon=[];
    for i=1:size(pitchWNon,2)
        residWNon(:,i)=pitchWNon(:,i)-mean(pitchBaseline')';
    end
    residWNon=residWNon(ptwindow,:); % only look at reasonable part

% REGRESSION
%         b=zeros(size(residWNon,2)-histwindow,histwindow);
%         clear bint
%     % Rows of b are the syllable renditions being predicted ("predictees")
%     % Columns of b are the syllable renditions doing the predicting ("predictors")
%         % Final column is the most recent syllable rendition
%         for i=histwindow+1:size(residWNon,2)
%             [b(i-histwindow,1:histwindow),bint(i).data]=regress(residWNon(:,i),residWNon(:,i-histwindow:i-1));
%             i
%         end    
%         regcoefs=b;
 b=zeros(size(residWNon,2)-histwindow,histwindow);
        for i=histwindow+1:size(residWNon,2)
            for j=1:100
                x=mean((residWNon(:,i)./residWNon(:,i-j)));
                b(i-histwindow,j)=x;
            end
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
            mrecentescapes(i)=median(holder); % Only look at the regression coefficients that are not outliers
        end
        bind=find(recenthitsabove(:,i));
        if ~isempty(bind)
            holder=(recentcoefs(bind,i));
            mrecenthits(i)=median(holder); % Only look at the regression coefficients that are not outliers
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


%%% 09.13.10
clear recentcoefs recentescapesabove recenthitsabove mrecentescapes mrecenthits mrecentabove
k=100;
for i=k+1:200%size(regcoefs,1)
    recentcoefs(i-k,:)=regcoefs(i,i-k:i-1);
    recentescapesabove(i-k,:)=indEscapeAbove(i-k:i-1);
    recenthitsabove(i-k,:)=indHitAbove(i-k:i-1);
end
for i=1:k
    aind=find(recentescapesabove(:,i));
    if ~isempty(aind)
        mrecentescapes(i)=mean(recentcoefs(aind,i));
    end
    bind=find(recenthitsabove(:,i));
    if ~isempty(bind)
        mrecenthits(i)=mean(recentcoefs(bind,i));
    end
    abind=find(recentescapesabove(:,i)+recenthitsabove(:,i));
    if ~isempty(aind) | ~isempty(bind)
        mrecentabove(i)=mean(recentcoefs(abind,i));
    end
end
mrecent=mean(recentcoefs);






fvals80w7Ahit=findwnoteJC('batch807Anotes','a','','',0,[2000 2700],8500,1,'obs0');
fvals807Amiss=findwnoteJC('batch807Anotes','b','','',0,[2000 2700],8500,1,'obs0');
for i=1:length(fvals807Ahit)
shifted807Ahit(i,:)=fvals807Ahit(i).datt;
end
for i=1:length(fvals807Amiss)
shifted807Amiss(i,:)=fvals807Amiss(i).datt;
end
pitch807Ahit=jc_pitchmat1024(shifted807Ahit,1024,1020,1,2000,2700,[1],'obs0',1);
pitch807Amiss=jc_pitchmat1024(shifted807Amiss,1024,1020,1,2000,2700,[1],'obs0',1);
pitch807A=[pitch807Ahit,pitch807Amiss];
figure;plot(mean(pitchBase'))
hold on;plot(mean(pitch807A'),'r')
hold on;plot(mean(pitch807Amiss'),'g')
hold on;plot(mean(pitch807Ahit'),'k')

for i=k+1:200%size(regcoefs,1)
    recentcoefs(i-k,:)=regcoefs(i,i-k:i-1);
    recentescapesabove(i-k,:)=indEscapeAbove(i-k:i-1);
    recenthitsabove(i-k,:)=indHitAbove(i-k:i-1);
end
for i=1:k
    aind=find(recentescapesabove(:,i));
    if ~isempty(aind)
        mrecentescapes(i)=mean(recentcoefs(aind,i));
    end
    bind=find(recenthitsabove(:,i));
    if ~isempty(bind)
        mrecenthits(i)=mean(recentcoefs(bind,i));
    end
    abind=find(recentescapesabove(:,i)+recenthitsabove(:,i));
    if ~isempty(aind) | ~isempty(bind)
        mrecentabove(i)=mean(recentcoefs(abind,i));
    end
end
mrecent=mean(recentcoefs);






fvals80w7Ahit=findwnoteJC('batch807Anotes','a','','',0,[2000 2700],8500,1,'obs0');
fvals807Amiss=findwnoteJC('batch807Anotes','b','','',0,[2000 2700],8500,1,'obs0');
for i=1:length(fvals807Ahit)
shifted807Ahit(i,:)=fvals807Ahit(i).datt;
end
for i=1:length(fvals807Amiss)
shifted807Amiss(i,:)=fvals807Amiss(i).datt;
end
pitch807Ahit=jc_pitchmat1024(shifted807Ahit,1024,1020,1,2000,2700,[1],'obs0',1);
pitch807Amiss=jc_pitchmat1024(shifted807Amiss,1024,1020,1,2000,2700,[1],'obs0',1);
pitch807A=[pitch807Ahit,pitch807Amiss];
figure;plot(mean(pitchBase'))
hold on;plot(mean(pitch807A'),'r')
hold on;plot(mean(pitch807Amiss'),'g')
hold on;plot(mean(pitch807Ahit'),'k')
