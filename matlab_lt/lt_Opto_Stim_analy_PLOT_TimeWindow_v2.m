function [StatsStruct, Params]=lt_Opto_Stim_analy_PLOT_TimeWindow(StatsStruct,Params,RunStats,KeepOutliers)
%% LT 4/7/15 - v2 - overwrites StatsStruct to save space.  Naming convention changed to assist automaticity. finds pitch outliers based on variant of tukey's method, within userdefined windows.
% = 0, saves outlier inds, but does not remove from data.
% RunStats=1 lt_Opto_Stim_analy_PLOT_TimeWindow_Statistics_v2 to plot and run stats



%% LT 1/15/15 - takes PC data from lt_Opto_Stim_analy_PLOT_Compare2 and allows you to look at certain time windows (in PC);
% INPUTS:
% StatsStruct
% Params
% (both outputs of lt_Opto_Stim_analy_PLOT_Compare2)

%% AUTO PARAMS

NumFields=length(Params.FieldsToCheck); % which trial classes? (e.g. stim vs. no stim);
PreDur=1000*Params.PreDur;
StimDur=Params.StimDur;

        if ~exist('KeepOutliers','var');
            KeepOutliers=1; % will include outliers in data
        end

%% FIRST, PLOT PCs to allow user to determine ideal tiem windows


if ~isfield(Params,'TimeWindowList') || ~isfield(Params,'TimeField'); % only run if timefields not already specified in params
    tPC=Params.tf_bins.tPC;
    tSP=Params.tf_bins.tSP;
    fSP=Params.tf_bins.fSP;
    
    for iii=1:NumFields;
        
        fieldname=Params.FieldsToCheck{iii};
        NumSyls=size(StatsStruct.(fieldname).PC,2); % num trials
        
        figure; hold on;
        
        % 1) PC
        PC=StatsStruct.(fieldname).PC;
        
        
        h1(1)=subplot(4,1,1); hold on;
        plot(tPC*1000,PC,'LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade
        plot(tPC*1000,mean(PC'),'Linewidth',2)
        ylabel('Frequency (Hz)')
        title([fieldname ': Pitch contours, Specgram, W.Entropy, and Amplitude. n=' num2str(NumSyls) '.']);
        xlim([tPC(1)*1000 tPC(end)*1000])
        
        % 2) MEAN SPEC and AMPLITUDE
        sp_mean=StatsStruct.(fieldname).sp_mean;
        sm_log_mean=StatsStruct.(fieldname).sm_log_mean;
        
        sm=StatsStruct.(fieldname).sm;
        sm_log=log(sm);
        
        
        % First, take log of sp
        % first, convert any sp values of 0 to non-zero(to the lowest value present);
        % solves problem of taking log of 0
        pp=find(sp_mean>0);
        mntmp = min(min(sp_mean(pp)));
        pp=find(sp_mean==0);
        sp_mean(pp) = mntmp;
        
        % second, take log
        sptemp=log(sp_mean);
        sptemp = sptemp - min(min(sptemp));
        sptemp = uint8(((2^8) - 1)*(sptemp./max(max(sptemp)))); % SAVE SOME MEMORY 8X less than 64 bit double
        
        h1(2)=subplot(4,1,2); hold on;
        imagesc(tSP*1000, fSP, sptemp);
        ylabel('Frequency (hz)');
        axis([tSP(1) tSP(end) fSP(1) fSP(end)]);
        
        
        
        % 3) Weiner Entropy
        h1(3)=subplot(4,1,3); hold on;
        WE=StatsStruct.(fieldname).WEntropyTimecourse;
        
        plot(tSP*1000,WE,'LineStyle','--','Color',[0.6 0.6 0.6]);
        plot(tSP*1000,mean(WE,2),'Linewidth',2,'Color','r'); % mean
        
        ylabel('(log) Weiner Entropy (-inf to 0)');
        
        
        
        % 4) Amplitude
        tSM=linspace(tSP(1),tSP(end),length(sm_log_mean));
        h1(4)=subplot(4,1,4); hold on;
        
        % individual contours
        plot(tSM*1000,sm_log,'LineStyle','--','Color',[0.6 0.6 0.6]); % individual contours
        plot(tSM*1000,sm_log_mean','Linewidth',2); % mean
        ylabel('Smoothed Amplitude (log scale)')
        xlim([tSP(1) tSP(end)])
        xlabel(['Time (ms); aligned at ' num2str(PreDur) 'ms'])
        
        
        lt_Opto_Stim_analy_PLOT_Compare2_PLOTSTIMS
        linkaxes(h1,'x');
        
    end
end


%% PLOT STATS
% Plan, input temporal windows, and this function will spit out stats
% on pitch, ampl, and WE in those windows

tPC=Params.tf_bins.tPC;

% QUERY USER for time windows, unless that info is already in params
if ~isfield(Params,'TimeWindowList') || ~isfield(Params,'TimeField'); % continue if not already fields
    
    % Query user to find epochs to analyze further.
    NumTimeWind=input('How many time epochs to analyze? (0 to skip) ');
    
    % get time windows
    for i=1:NumTimeWind;
        Params.TimeWindowList{i}=input(['Enter time window (in ms) for epoch # ' num2str(i) ' (e.g. [250 270])']);
    end
    % convert to strings - field names
    for i=1:length(Params.TimeWindowList);
        TimeWind = Params.TimeWindowList{i};
        if ceil(Params.TimeWindowList{i})~=Params.TimeWindowList{i}; % means that there is fraction
            Params.TimeField{i}=['Wind' num2str(floor(TimeWind(1))) 'fr_' num2str(floor(TimeWind(2))) 'frms']; % put "fr" at end to signify fraction
        else
            Params.TimeField{i}=['Wind' num2str(TimeWind(1)) '_' num2str(TimeWind(2)) 'ms'];
        end
    end
end



%% REPLOT PLOTS, BUT NOW WITH TIME WINDOW LINES
close all;

tPC=Params.tf_bins.tPC;
tSP=Params.tf_bins.tSP;
fSP=Params.tf_bins.fSP;


for iii=1:NumFields;
    
    fieldname=Params.FieldsToCheck{iii};
    NumSyls=size(StatsStruct.(fieldname).PC,2); % num trials
    
    figure; hold on;
    
    % 1) PC
    PC=StatsStruct.(fieldname).PC;
    
    
    h1(1)=subplot(4,1,1); hold on;
    plot(tPC*1000,PC,'LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade
    plot(tPC*1000,mean(PC'),'Linewidth',2)
    ylabel('Frequency (Hz)')
    title([fieldname ': Pitch contours, Specgram, W.Entropy, and Amplitude. n=' num2str(NumSyls) '.']);
    xlim([tPC(1)*1000 tPC(end)*1000])
    
    
    % add time window lines
    for k=1:length(Params.TimeWindowList);
        timewind=Params.TimeWindowList{k};
        
        line([timewind(1) timewind(1)], ylim,'Color','r');
        line([timewind(2) timewind(2)], ylim,'Color','r');
    end
    
    
    % 2) MEAN SPEC and AMPLITUDE
    sp_mean=StatsStruct.(fieldname).sp_mean;
    sm_log_mean=StatsStruct.(fieldname).sm_log_mean;
    
        sm=StatsStruct.(fieldname).sm;
        sm_log=log(sm);
    
    
    % First, take log of sp
    % first, convert any sp values of 0 to non-zero(to the lowest value present);
    % solves problem of taking log of 0
    pp=find(sp_mean>0);
    mntmp = min(min(sp_mean(pp)));
    pp=find(sp_mean==0);
    sp_mean(pp) = mntmp;
    
    % second, take log
    sptemp=log(sp_mean);
    sptemp = sptemp - min(min(sptemp));
    sptemp = uint8(((2^8) - 1)*(sptemp./max(max(sptemp)))); % SAVE SOME MEMORY 8X less than 64 bit double
    
    h1(2)=subplot(4,1,2); hold on;
    imagesc(tSP*1000, fSP, sptemp);
    ylabel('Frequency (hz)');
    axis([tSP(1) tSP(end) fSP(1) fSP(end)]);
    
    
    
    % 3) Weiner Entropy
    h1(3)=subplot(4,1,3); hold on;
    WE=StatsStruct.(fieldname).WEntropyTimecourse;
    
    plot(tSP*1000,WE,'LineStyle','--','Color',[0.6 0.6 0.6]);
    plot(tSP*1000,mean(WE,2),'Linewidth',2,'Color','r'); % mean
    
    ylabel('(log) Weiner Entropy (-inf to 0)');
    
    
    
    % 4) Amplitude
    tSM=linspace(tSP(1),tSP(end),length(sm_log_mean));
    h1(4)=subplot(4,1,4); hold on;
    
    % individual contours
    plot(tSM*1000,sm_log,'LineStyle','--','Color',[0.6 0.6 0.6]); % individual contours
    plot(tSM*1000,sm_log_mean','Linewidth',2); % mean
    ylabel('Smoothed Amplitude (log scale)')
    xlim([tSP(1) tSP(end)])
    xlabel(['Time (ms); aligned at ' num2str(PreDur) 'ms'])
    
    
    lt_Opto_Stim_analy_PLOT_Compare2_PLOTSTIMS
    linkaxes(h1,'x');
    
    % add time window lines
    for k=1:length(Params.TimeWindowList);
        timewind=Params.TimeWindowList{k};
        
        line([timewind(1) timewind(1)], ylim,'Color','r');
        line([timewind(2) timewind(2)], ylim,'Color','r');
    end
    
end






%% for each time window, extract data
for i=1:length(Params.TimeWindowList);
    TimeField=Params.TimeField{i};
    TimeWind = Params.TimeWindowList{i};
    
    for ii=1:length(Params.FieldsToCheck); % for each window, analyze each trial type
        fieldname = Params.FieldsToCheck{ii};
        NumSyls=length(StatsStruct.(fieldname).TimeSinceLastTrig); % num trials;
        
        % EXTRACT INFO
        % Sample size
        StatsStruct.(fieldname).WINDOWED.(TimeField).n=NumSyls;
        
        
        % Time of day info
        for ii=1:NumSyls; % renditions
            StatsStruct.(fieldname).WINDOWED.(TimeField).Time.datenum(ii)=...
                StatsStruct.(fieldname).datenum(ii); % datenum
            
            [A B]=lt_convert_datenum_to_hour(StatsStruct.(fieldname).WINDOWED.(TimeField).Time.datenum(ii));
            StatsStruct.(fieldname).WINDOWED.(TimeField).Time.hours(ii)=B.hours; % units of hours
            StatsStruct.(fieldname).WINDOWED.(TimeField).Time.days(ii)=B.days; % units of days
            StatsStruct.(fieldname).WINDOWED.(TimeField).Time.date{ii}=A.ddmmmyyyy; % the date
        end
        
        
        % PITCH CONTOUR
        t1=find(tPC-TimeWind(1)/1000>0,1,'first'); % start index
        t2=find(tPC-TimeWind(2)/1000<0,1,'last'); % last index
        
        
        StatsStruct.(fieldname).WINDOWED.(TimeField).Pitch.vals=mean(StatsStruct.(fieldname).PC(t1:t2,:)); % take time mean of PC in that window
        StatsStruct.(fieldname).WINDOWED.(TimeField).Pitch.mean=mean(StatsStruct.(fieldname).WINDOWED.(TimeField).Pitch.vals);
        StatsStruct.(fieldname).WINDOWED.(TimeField).Pitch.median=median(StatsStruct.(fieldname).WINDOWED.(TimeField).Pitch.vals);
        StatsStruct.(fieldname).WINDOWED.(TimeField).Pitch.SD=std(StatsStruct.(fieldname).WINDOWED.(TimeField).Pitch.vals);
        StatsStruct.(fieldname).WINDOWED.(TimeField).Pitch.SEM=StatsStruct.(fieldname).WINDOWED.(TimeField).Pitch.SD/...
            sqrt(StatsStruct.(fieldname).WINDOWED.(TimeField).n-1);
        
        % AMPLITUDE
        t1=find(tSM-TimeWind(1)/1000>0,1,'first'); % start index
        t2=find(tSM-TimeWind(2)/1000<0,1,'last'); % last index
        
        StatsStruct.(fieldname).WINDOWED.(TimeField).Ampl.vals_log=mean(log(StatsStruct.(fieldname).sm(t1:t2,:)));
        StatsStruct.(fieldname).WINDOWED.(TimeField).Ampl.mean=mean(StatsStruct.(fieldname).WINDOWED.(TimeField).Ampl.vals_log);
        StatsStruct.(fieldname).WINDOWED.(TimeField).Ampl.median=median(StatsStruct.(fieldname).WINDOWED.(TimeField).Ampl.vals_log);
        StatsStruct.(fieldname).WINDOWED.(TimeField).Ampl.SD=std(StatsStruct.(fieldname).WINDOWED.(TimeField).Ampl.vals_log);
        StatsStruct.(fieldname).WINDOWED.(TimeField).Ampl.SEM=StatsStruct.(fieldname).WINDOWED.(TimeField).Ampl.SD...
            /sqrt(StatsStruct.(fieldname).WINDOWED.(TimeField).n-1);
        
        % ENTROPY
        t1=find(tSP-TimeWind(1)/1000>0,1,'first'); % start index
        t2=find(tSP-TimeWind(2)/1000<0,1,'last'); % last index
        
        % if window so small that there is no time bin inside, then take
        % bin immediately larger then desired window
        if t1>t2; % this can happen because i am trying to get inner values.
            tmp=t2;
            t2=t1;
            t1=tmp; % this gets smallest outer bins
            disp('NOTE: in getting entropy time bins, time window too small, so got smallest window enclosing desired times');
        end
        
        
        if size(StatsStruct.(fieldname).WEntropyTimecourse(t1:t2,:),1)>1; % i.e. if more than one time bin, then take temporal mean, otherwise take that bin's value.
            StatsStruct.(fieldname).WINDOWED.(TimeField).WEntropy.vals=mean(StatsStruct.(fieldname).WEntropyTimecourse(t1:t2,:));
        else
            StatsStruct.(fieldname).WINDOWED.(TimeField).WEntropy.vals=StatsStruct.(fieldname).WEntropyTimecourse(t1:t2,:);
        end
        
        StatsStruct.(fieldname).WINDOWED.(TimeField).WEntropy.mean=mean(StatsStruct.(fieldname).WINDOWED.(TimeField).WEntropy.vals);
        StatsStruct.(fieldname).WINDOWED.(TimeField).WEntropy.median=median(StatsStruct.(fieldname).WINDOWED.(TimeField).WEntropy.vals);
        
        StatsStruct.(fieldname).WINDOWED.(TimeField).WEntropy.SD=std(StatsStruct.(fieldname).WINDOWED.(TimeField).WEntropy.vals);
        StatsStruct.(fieldname).WINDOWED.(TimeField).WEntropy.SEM=StatsStruct.(fieldname).WINDOWED.(TimeField).WEntropy.SD/...
            sqrt(StatsStruct.(fieldname).WINDOWED.(TimeField).n-1);
        
    end
end


%% FIND OUTLIERS FOR EACH TIME WINDOW, BASED ON PITCH CONTOUR
% will simply save the inds of the outliers
disp('Finding outliers based on PC... will save Inds but not remove data');

Outliers={};
OutlInds={};

for i=1:length(Params.TimeWindowList);
    TimeField=Params.TimeField{i};
    TimeWind = Params.TimeWindowList{i};
    
    for ii=1:length(Params.FieldsToCheck); % for each window, analyze each trial type
        fieldname = Params.FieldsToCheck{ii};
        NumSyls=length(StatsStruct.(fieldname).TimeSinceLastTrig); % num trials;
        
        
        % determine if trials are outliers based on their position in
        % spread of all data
        PC_mat=StatsStruct.(fieldname).PC; % cols are trials, rows are timepoints
        
        % which timepoints to check?
        t1=find(tPC-TimeWind(1)/1000>0,1,'first'); % start index
        t2=find(tPC-TimeWind(2)/1000<0,1,'last'); % last index
        
        rows_to_check=[t1 t2];
        num_rows=t2-t1+1;
        
        Outliers{i}{ii}=[];
        
        for iii=1:num_rows; % for each column, perform check.
            % empirically determined 2 to be a good iqr value to use.
            % (based on pu37wh20)
            [~, B, C]=lt_db_tukey_outlier_90tile(PC_mat',rows_to_check(1)+iii-1,3);
            Outliers{i}{ii}=[Outliers{i}{ii} B' C']; % accumulate index of outliers
        end
        OutlInds{i}{ii}=unique(Outliers{i}{ii});
        
        
        % save outliers
        StatsStruct.(fieldname).WINDOWED.(TimeField).OutlierInds=OutlInds{i}{ii};
        
    end
end


% count how many outliers total
outliers={};
for ii=1:length(Params.FieldsToCheck); % for each window, analyze each trial type
    fieldname=Params.FieldsToCheck{ii};
    
    outliers{ii}=[];
    for i=1:length(Params.TimeWindowList);
        
        outliers{ii}=[outliers{ii} OutlInds{i}{ii}];
        
        
    end
    unique_outliers{ii}=unique(outliers{ii});
    
    % save outliers
    StatsStruct.(fieldname).OutlierInds_AllWindsCombn=unique_outliers{ii};
    
    
    % get title for plots
    tmp=100*length(unique_outliers{ii})/(length(StatsStruct.(Params.FieldsToCheck{ii}).TimeSinceLastTrig)); % fraction that are outliers
    
    outlier_title{ii}=([Params.FieldsToCheck{ii} ' has ' num2str(length(unique_outliers{ii})) ' outliers out of ' ...
        num2str(length(StatsStruct.(Params.FieldsToCheck{ii}).TimeSinceLastTrig)) ' total trials (' num2str(tmp,2) '%)']);
end


% Plot outliers
figure; hold on;
% plotcolors=lt_make_plot_colors(length(Params.TimeWindowList),0,0); % one for each time field.

for iii=1:NumFields;
    
    fieldname=Params.FieldsToCheck{iii};
    %     NumSyls=size(StatsStruct.(fieldname).PC,2); % num trials
    
    subplot(NumFields,1,iii); hold on;
    title([outlier_title{iii}]);
    
    % first plot all trials, not just outliers
    PC=StatsStruct.(fieldname).PC;
    
    %     plot(tPC*1000,PC,'LineStyle','--','Color',[0.8 0.8 0.8]) % plot all pitch contours in light shade
    plot(PC,'LineStyle','--','Color',[0.7 0.7 0.7]) % plot all pitch contours
    
    % then plot outliers
    for j=1:length(Params.TimeWindowList);
        TimeField=Params.TimeField{j};
        TimeWind = Params.TimeWindowList{j};
        
        outinds=OutlInds{j}{iii};
        
        if length(outinds)>0;
            
            % 1) PC
            PC=StatsStruct.(fieldname).PC;
            
            
            %             plot(tPC*1000,PC(:,outinds),'-','Color',[rand rand rand]) % plot all pitch contours in light shade
            plot(PC(:,outinds),'-','Color',[0.4 0.8 0.4]) % plot all pitch contours in light shade
            
            ylabel('Frequency (Hz)')
            %             xlim([tPC(1)*1000 tPC(end)*1000])
            
            
            
            
        end
    end
    % add time window lines
    for k=1:length(Params.TimeWindowList);
        timewind=Params.TimeWindowList{k};
        
        t1=find(tPC-timewind(1)/1000>0,1,'first'); % start index
        t2=find(tPC-timewind(2)/1000<0,1,'last'); % last index
        
        %         line([timewind(1) timewind(1)], ylim,'Color','r');
        %         line([timewind(2) timewind(2)], ylim,'Color','r');
        
        line([t1 t1], ylim,'Color','k', 'LineWidth', 2);
        line([t2 t2], ylim,'Color','k', 'LineWidth',2);
        
    end
    
end

lt_subtitle('Outliers, based on pitch contours');






%% PLOT TIMEWINDOW STATS
% X=[];
% Xse=[];
%
% PlotColors=lt_make_plot_colors(NumFields,0);
% TimeWindFields=fieldnames(StatsStruct.(fieldname).WINDOWED);
% xx=1:length(Params.TimeWindowList); % x coordinates
%
% % get means and stats
% for i=1:length(Params.FieldsToCheck);
%     fieldname = Params.FieldsToCheck{i};
%     NumSyls=length(StatsStruct.(fieldname).TimeSinceLastTrig);
%
%     for ii=1:length(Params.TimeWindowList);
%         TimeWind = Params.TimeWindowList{ii};
%         TimeField=['Wind' num2str(TimeWind(1)) '_' num2str(TimeWind(2)) 'ms'];
%
%         X.pitch(ii)=StatsStruct.(fieldname).WINDOWED.(TimeField).Pitch.mean;
%         X.ampl(ii)=StatsStruct.(fieldname).WINDOWED.(TimeField).Ampl.mean;
%         X.entr(ii)=StatsStruct.(fieldname).WINDOWED.(TimeField).WEntropy.mean;
%
%         Xse.pitch(ii)=StatsStruct.(fieldname).WINDOWED.(TimeField).Pitch.SEM;
%         Xse.ampl(ii)=StatsStruct.(fieldname).WINDOWED.(TimeField).Ampl.SEM;
%         Xse.entr(ii)=StatsStruct.(fieldname).WINDOWED.(TimeField).WEntropy.SEM;
%
%     end
% end
%
%     % PITCH
%     subplot(3,1,1); hold on;
%     errorbar(xx+0.02*(i-2),X.pitch,Xse.pitch,'o','Color',PlotColors{i},'MarkerSize',8);
%     set(gca,'XTick',1:length(Params.TimeWindowList),'XTickLabel',TimeWindFields);
%     title('Mean pitch +/- SEM');
%     ylabel('Freq (hz)');
%
%
%     % AMPLIT
%     subplot(3,1,2); hold on;
%     errorbar(xx+0.02*(i-2),X.ampl,Xse.ampl,'o','Color',PlotColors{i},'MarkerSize',8);
%     title('Mean Amplitude +/- SEM');
%     ylabel('Amplitude (arbitrary log units)');
%     set(gca,'XTick',1:length(Params.TimeWindowList),'XTickLabel',TimeWindFields);
%
%     % ENTROPY
%     subplot(3,1,3); hold on;
%     hfig5(i)=errorbar(xx+0.02*(i-2),X.entr,Xse.entr,'o','Color', PlotColors{i},'MarkerSize',8);
%     title('Mean W. Entropy +/- SEM');
%     ylabel('W. Entropy (range -inf to 0)');
%     set(gca,'XTick',1:length(Params.TimeWindowList),'XTickLabel',TimeWindFields);
%
% end
%
% lt_subtitle('Statistics in temporal windows');
%
% legend(hfig5,Params.FieldsToCheck);





%% PLOT AND PERFORM STATISTICS
if RunStats==1;
    if length(Params.FieldsToCheck)==2; % only do comaprison stats if have 2 trial types.
        [StatsStruct,Params]=lt_Opto_Stim_analy_PLOT_TimeWindow_Statistics_v2(StatsStruct,Params,KeepOutliers);
    end
end


%% PLOT TIMEWINDOW STATS
if (0)
PlotColors=lt_make_plot_colors(NumFields,0);
TimeWindFields=Params.TimeField;
NumTimeWinds=length(TimeWindFields);

% PITCH
figure; hold on;
for i=1:NumTimeWinds;
    timewind=TimeWindFields{i};
    subplot(1,NumTimeWinds,i); hold on;
    title(timewind);
    for ii=1:NumFields; % for each trail type
        fieldname=Params.FieldsToCheck{ii};
        
        % plot individual vals
        Y=StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.vals;
        X=0.5*ii-0.4+0.3*rand(length(Y),1); % random values to disperse vals
        plot(X,Y,'o','Color',PlotColors{ii});
        ylabel('Hz');
        
        % plot mean vals
        hfig5(ii)=errorbar(0.5*ii-0.25,StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.mean,...
            StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.SEM,'s','Color','k','MarkerFaceColor',PlotColors{ii}, 'MarkerSize',8);
    end
end

legend(hfig5,Params.FieldsToCheck);
lt_subtitle('Pitch (MEAN), all trials, sorted by time window and trial');

% PITCH - using median
figure; hold on;
for i=1:NumTimeWinds;
    timewind=TimeWindFields{i};
    subplot(1,NumTimeWinds,i); hold on;
    title(timewind);
    for ii=1:NumFields; % for each trail type
        fieldname=Params.FieldsToCheck{ii};
        
        % plot individual vals
        Y=StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.vals;
        X=0.5*ii-0.4+0.3*rand(length(Y),1); % random values to disperse vals
        plot(X,Y,'o','Color',PlotColors{ii});
        ylabel('Hz');
        
        % plot mean vals
        hfig5(ii)=errorbar(0.5*ii-0.25,StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.median,...
            StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.SEM,'s','Color','k','MarkerFaceColor',PlotColors{ii}, 'MarkerSize',8);
    end
end

legend(hfig5,Params.FieldsToCheck);
lt_subtitle('Pitch (MEDIAN), all trials, sorted by time window and trial');



% AMPLITUDE
figure; hold on;
for i=1:NumTimeWinds;
    timewind=TimeWindFields{i};
    subplot(1,NumTimeWinds,i); hold on;
    for ii=1:NumFields; % for each trail type
        fieldname=Params.FieldsToCheck{ii};
        
        % plot individual vals
        Y=StatsStruct.(fieldname).WINDOWED.(timewind).Ampl.vals_log;
        X=0.5*ii-0.4+0.3*rand(length(Y),1); % random values to disperse vals
        plot(X,Y,'o','Color',PlotColors{ii});
        ylabel('Amplitude (log)');
        % plot mean vals
        hfig6(ii)=errorbar(0.5*ii-0.25,StatsStruct.(fieldname).WINDOWED.(timewind).Ampl.mean,...
            StatsStruct.(fieldname).WINDOWED.(timewind).Ampl.SEM,'s','Color','k','MarkerFaceColor',PlotColors{ii}, 'MarkerSize',8);
    end
end

legend(hfig6,Params.FieldsToCheck);
lt_subtitle('Amplitude, all trials, sorted by time window and trial');



% ENTROPY
figure; hold on;
for i=1:NumTimeWinds;
    timewind=TimeWindFields{i};
    subplot(1,NumTimeWinds,i); hold on;
    for ii=1:NumFields; % for each trail type
        fieldname=Params.FieldsToCheck{ii};
        
        % plot individual vals
        Y=StatsStruct.(fieldname).WINDOWED.(timewind).WEntropy.vals;
        X=0.5*ii-0.4+0.3*rand(length(Y),1); % random values to disperse vals
        plot(X,Y,'o','Color',PlotColors{ii});
        ylabel('W. Entropy (log)');
        % plot mean vals
        hfig7(ii)=errorbar(0.5*ii-0.25,StatsStruct.(fieldname).WINDOWED.(timewind).WEntropy.mean,...
            StatsStruct.(fieldname).WINDOWED.(timewind).WEntropy.SEM,'s','Color','k','MarkerFaceColor',PlotColors{ii}, 'MarkerSize',8);
    end
end

legend(hfig7,Params.FieldsToCheck);
lt_subtitle('Entropy, all trials, sorted by time window and trial');



%% PLOT AS ABOVE, BUT WITHOUT SCATTER OF DATA - just means and SEM

PlotColors=lt_make_plot_colors(NumFields,0);
TimeWindFields=Params.TimeField;
NumTimeWinds=length(TimeWindFields);

% PITCH
figure; hold on;
for i=1:NumTimeWinds;
    timewind=TimeWindFields{i};
    subplot(1,NumTimeWinds,i); hold on;
    title(timewind); ylabel('Hz');
    
    for ii=1:NumFields; % for each trail type
        fieldname=Params.FieldsToCheck{ii};
        
        % plot mean vals
        hfig5(ii)=errorbar(0.5*ii-0.25,StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.mean,...
            StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.SEM,'s','Color','k','MarkerFaceColor',PlotColors{ii}, 'MarkerSize',8);
    end
end

legend(hfig5,Params.FieldsToCheck);
lt_subtitle('Pitch (MEAN/SEM), sorted by time window and trial');


% PITCH - median
figure; hold on;
for i=1:NumTimeWinds;
    timewind=TimeWindFields{i};
    subplot(1,NumTimeWinds,i); hold on;
    title(timewind); ylabel('Hz');
    
    for ii=1:NumFields; % for each trail type
        fieldname=Params.FieldsToCheck{ii};
        
        % plot mean vals
        hfig5(ii)=errorbar(0.5*ii-0.25,StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.median,...
            StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.SEM,'s','Color','k','MarkerFaceColor',PlotColors{ii}, 'MarkerSize',8);
    end
end

legend(hfig5,Params.FieldsToCheck);
lt_subtitle('Pitch (MEDIAN/SEM), sorted by time window and trial');



% AMPLITUDE
figure; hold on;
for i=1:NumTimeWinds;
    timewind=TimeWindFields{i};
    subplot(1,NumTimeWinds,i); hold on;
    ylabel('Amplitude (log)');
    
    for ii=1:NumFields; % for each trail type
        fieldname=Params.FieldsToCheck{ii};
        
        % plot mean vals
        hfig6(ii)=errorbar(0.5*ii-0.25,StatsStruct.(fieldname).WINDOWED.(timewind).Ampl.mean,...
            StatsStruct.(fieldname).WINDOWED.(timewind).Ampl.SEM,'s','Color','k','MarkerFaceColor',PlotColors{ii}, 'MarkerSize',8);
    end
end

legend(hfig6,Params.FieldsToCheck);
lt_subtitle('Amplitude, sorted by time window and trial');



% ENTROPY
figure; hold on;
for i=1:NumTimeWinds;
    timewind=TimeWindFields{i};
    subplot(1,NumTimeWinds,i); hold on;
    ylabel('W. Entropy (log)');
    for ii=1:NumFields; % for each trail type
        fieldname=Params.FieldsToCheck{ii};
        
        % plot mean vals
        hfig7(ii)=errorbar(0.5*ii-0.25,StatsStruct.(fieldname).WINDOWED.(timewind).WEntropy.mean,...
            StatsStruct.(fieldname).WINDOWED.(timewind).WEntropy.SEM,'s','Color','k','MarkerFaceColor',PlotColors{ii}, 'MarkerSize',8);
    end
end

legend(hfig7,Params.FieldsToCheck);
lt_subtitle('Entropy, sorted by time window and trial');




%% PLOT HISTOGRAM
if (0); % I usually ignore, so just not plot for now.
    numbins=20;
    
    % determine subplot sizes
    [~, row_plots, col_plots]=lt_get_subplot_size(NumTimeWinds,NumTimeWinds);
    
    
    % PITCH
    figure; hold on;
    for i=1:length(TimeWindFields);
        timefield=TimeWindFields{i};
        subplot(row_plots,col_plots,i); hold on;
        title(timefield);
        
        % determine what edges to use for hist (based on max and min of all trial
        % types)
        XX=[];
        for ii=1:NumFields; % all trial types
            fieldname=Params.FieldsToCheck{ii};
            XX=[XX StatsStruct.(fieldname).WINDOWED.(timefield).Pitch.vals];
        end
        xmin=min(XX);
        xmax=max(XX);
        centers=linspace(xmin-20,xmax+20,numbins);
        
        
        vals=[];
        % get histogram and plot
        for ii=1:NumFields; % all trial types
            fieldname=Params.FieldsToCheck{ii};
            
            % PITCH
            vals=StatsStruct.(fieldname).WINDOWED.(timefield).Pitch.vals;
            [N,~]= hist(vals,centers);
            
            % get probability density
            Npdf=N/sum(N);
            
            % plot
            hpdf1(ii)=plot(centers,Npdf,'-','Color',PlotColors{ii});
            
            % Mark mean and sem
            Xmean=StatsStruct.(fieldname).WINDOWED.(timefield).Pitch.mean;
            Xsem=StatsStruct.(fieldname).WINDOWED.(timefield).Pitch.SEM;
            
            line([Xmean Xmean],[0.0075*(ii-1) 0.0075*ii],'Color',PlotColors{ii},'LineWidth',2);
            line([Xmean-Xsem Xmean+Xsem],[0.0075*ii 0.0075*ii]-0.0037,'Color',PlotColors{ii},'LineWidth',2);
            
            
        end
        legend(hpdf1,Params.FieldsToCheck);
    end
    lt_subtitle('Pitch');
    
    
    % AMPLITUDE
    figure; hold on;
    for i=1:length(TimeWindFields);
        timefield=TimeWindFields{i};
        subplot(row_plots,col_plots,i); hold on;
        title(timefield);
        
        % determine what edges to use for hist (based on max and min of all trial
        % types)
        XX=[];
        for ii=1:NumFields; % all trial types
            fieldname=Params.FieldsToCheck{ii};
            XX=[XX StatsStruct.(fieldname).WINDOWED.(timefield).Ampl.vals_log];
        end
        xmin=min(XX);
        xmax=max(XX);
        centers=linspace(xmin-1,xmax+1,numbins);
        
        
        vals=[];
        % get histogram and plot
        for ii=1:NumFields; % all trial types
            fieldname=Params.FieldsToCheck{ii};
            
            % PITCH
            vals=StatsStruct.(fieldname).WINDOWED.(timefield).Ampl.vals_log;
            [N,~]= hist(vals,centers);
            
            % get probability density
            Npdf=N/sum(N);
            
            % plot
            hpdf2(ii)=plot(centers,Npdf,'-','Color',PlotColors{ii});
            
            % Mark mean and sem
            Xmean=StatsStruct.(fieldname).WINDOWED.(timefield).Ampl.mean;
            Xsem=StatsStruct.(fieldname).WINDOWED.(timefield).Ampl.SEM;
            
            line([Xmean Xmean],[0.0075*(ii-1) 0.0075*ii],'Color',PlotColors{ii},'LineWidth',2);
            line([Xmean-Xsem Xmean+Xsem],[0.0075*ii 0.0075*ii]-0.0037,'Color',PlotColors{ii},'LineWidth',2);
            
            
        end
        legend(hpdf2,Params.FieldsToCheck);
    end
    lt_subtitle('Amplitude');
    
    
    % ENTROPY
    figure; hold on;
    for i=1:length(TimeWindFields);
        timefield=TimeWindFields{i};
        subplot(row_plots,col_plots,i); hold on;
        title(timefield);
        
        % determine what edges to use for hist (based on max and min of all trial
        % types)
        XX=[];
        for ii=1:NumFields; % all trial types
            fieldname=Params.FieldsToCheck{ii};
            XX=[XX StatsStruct.(fieldname).WINDOWED.(timefield).WEntropy.vals];
        end
        xmin=min(XX);
        xmax=max(XX);
        centers=linspace(xmin-1,xmax+1,numbins);
        
        
        vals=[];
        % get histogram and plot
        for ii=1:NumFields; % all trial types
            fieldname=Params.FieldsToCheck{ii};
            
            % PITCH
            vals=StatsStruct.(fieldname).WINDOWED.(timefield).WEntropy.vals;
            [N,~]= hist(vals,centers);
            
            % get probability density
            Npdf=N/sum(N);
            
            % plot
            hpdf3(ii)=plot(centers,Npdf,'-','Color',PlotColors{ii});
            
            % Mark mean and sem
            Xmean=StatsStruct.(fieldname).WINDOWED.(timefield).WEntropy.mean;
            Xsem=StatsStruct.(fieldname).WINDOWED.(timefield).WEntropy.SEM;
            
            line([Xmean Xmean],[0.0075*(ii-1) 0.0075*ii],'Color',PlotColors{ii},'LineWidth',2);
            line([Xmean-Xsem Xmean+Xsem],[0.0075*ii 0.0075*ii]-0.0037,'Color',PlotColors{ii},'LineWidth',2);
            
            
        end
        legend(hpdf3,Params.FieldsToCheck);
    end
    lt_subtitle('Entropy');
    
end

end
%% SAVE
disp('Saving...')

cd(Params.savefolder);
tstamp=lt_get_timestamp(0);

% save
save('StatsStruct','StatsStruct');
save('Params','Params');

% note completion
DoneNote=['DONE_TimeWindow_' tstamp '.txt'];

fid1=fopen(DoneNote,'w');
fclose(fid1);

% save figs
% try cd('FIGURES/TimeWindow');
% catch err
%     mkdir('FIGURES/TimeWindow');
%     cd('FIGURES/TimeWindow');
% end
% 
% lt_save_all_figs;
% 
% cd('../../');

disp('Done!')

end


