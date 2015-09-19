function [Params, AllDays_PlotLearning]= lt_seq_dep_pitch_PlotLearning_Musc(Params, AllDays_PlotLearning, saveON)

%     Params.PlotLearning.Lag_time=1.7; % time from switch to musc
%
%     Params.PlotLearning.PBS_window=[-1.5 0]; % time before switch for PBS
%       Params.PlotLearning.Dur_of_PBS_dat_that_counts_as_MUSC=0.5; half
%       hour from musc to pbs, still call musc data.
%     saveON=1;

%% PARAMS1

% == DEFAULTS
if ~isfield(Params.PlotLearning, 'Lag_time');
    Params.PlotLearning.Lag_time=1.7; % time from switch to musc
end

if ~isfield(Params.PlotLearning, 'PBS_window');
    Params.PlotLearning.PBS_window=[-1.5 0]; % time before switch for PBS
end

if ~isfield(Params.PlotLearning, 'Dur_of_PBS_dat_that_counts_as_MUSC');
    Params.PlotLearning.Dur_of_PBS_dat_that_counts_as_MUSC=0;
end

if ~isfield('saveON', 'var');
    saveON=1;
end

if ~isfield(Params.PlotLearning, 'timeline')
    Params.PlotLearning.timeline={};
end




% ==
Lag_time=Params.PlotLearning.Lag_time; % how many hours post MUSC start to collect data?
PBS_window=Params.PlotLearning.PBS_window; % [-2 0] means take PBS data in window -2hrs to 0.

ExptCondition_codes={'PBS','MUSC'}; % i.e. if filename says "PBS", then is code 1 (musc=2 and so forth)
Params.PlotLearning.ExptCondition_codes=ExptCondition_codes;

%     % convert to codes
%     PBS_code=find(strcmp(ExptCondition_codes, 'PBS'));
%     MUSC_code=find(strcmp(ExptCondition_codes, 'PBS'));

plotcols=lt_make_plot_colors(length(ExptCondition_codes), 0, 0);

SylFields_Unique=Params.PlotLearning.SylFields_Unique; % all seq and syls (assume 1st day has all possible syls)

NumDays=length(AllDays_PlotLearning.DataMatrix.(SylFields_Unique{1}).FFvals);
FirstDay=Params.SeqFilter.FirstDay;
plotWNdays=Params.PlotLearning.plotWNdays;
LastDay=Params.SeqFilter.LastDay;



% globals
global WNTimeOnInd
WNTimeOnInd=Params.PlotLearning.WNTimeOnInd;

global WNTimeOffInd
WNTimeOffInd=Params.PlotLearning.WNTimeOffInd;

global DaysToMarkInds
DaysToMarkInds=Params.PlotLearning.DaysToMarkInds;

plotWNdays=Params.PlotLearning.plotWNdays;

% %% THROW OUT DAYS WITH NO MUSCIMOL EFFECT
%
% for i=Params.PlotLearning.MuscimolDaysToThrowout_Inds; % days
%
%     for ii=1:length(SylFields_Unique);
%         syl=SylFields_Unique{ii};
%
%         [AllDays_PlotLearning.DataMatrix_MUSC.(syl)(:)]
%



% ================= CHUNKING DAYS (early consolid, late consolid, bidir)
% ====== Figure out the days that are in those three caregories
% First check that timeline info exists
if isfield(Params.PlotLearning, 'timeline');
    try
        eventtimes=lt_convert_EventTimes_to_RelTimes(FirstDay, {Params.PlotLearning.timeline.consolid_start});
        Params.PlotLearning.timeline.consolid_start_Inds=eventtimes.FinalValue;
    catch err
    end
    
    try
        eventtimes=lt_convert_EventTimes_to_RelTimes(FirstDay, {Params.PlotLearning.timeline.consolid_end});
        Params.PlotLearning.timeline.consolid_end_Inds=eventtimes.FinalValue;
    catch err
    end
    
    try
        eventtimes=lt_convert_EventTimes_to_RelTimes(FirstDay, {Params.PlotLearning.timeline.bidir_start});
        Params.PlotLearning.timeline.bidir_start_Inds=eventtimes.FinalValue;
    catch err
    end
    
    try
        eventtimes=lt_convert_EventTimes_to_RelTimes(FirstDay, {Params.PlotLearning.timeline.bidir_end});
        Params.PlotLearning.timeline.bidir_end_Inds=eventtimes.FinalValue;
    catch err
    end
    
%     % ====== CONSOLIDATION PERIOD
%     if isfield(Params.PlotLearning.timeline, 'consolid_start_Inds')
%         consolid_start=Params.PlotLearning.timeline.consolid_start_Inds; % start
%         if isfield(Params.PlotLearning.timeline, 'consolid_end_Inds')
%             consolid_end=Params.PlotLearning.timeline.consolid_end_Inds; % end (if defined)
%         
%         
%         
%         else
%             consolid_end=Params.PlotLearning.WNTimeOffInd; % if not defined, then use the last day of WN)
%         end
%         
%         duration_consolidation=consolid_end-consolid_start+1;
%         midday_consolid=consolid_start + ceil(duration_consolidation/2)-1;
%         
%         Params.PlotLearning.timeline.days_consolid_early=consolid_start:midday_consolid; % take first half of data
%         Params.PlotLearning.timeline.days_consolid_late=midday_consolid+1:consolid_end; % take second half of data
%         
%         
%     end
    
    % ====== CONSOLIDATION PERIOD
    if isfield(Params.PlotLearning.timeline, 'consolid_start_Inds')
        consolid_start=Params.PlotLearning.timeline.consolid_start_Inds; % start
        if isfield(Params.PlotLearning.timeline, 'consolid_end_Inds')
            consolid_end=Params.PlotLearning.timeline.consolid_end_Inds; % end (if defined)

        duration_consolidation=consolid_end-consolid_start+1;
        midday_consolid=consolid_start + ceil(duration_consolidation/2)-1;
        
        Params.PlotLearning.timeline.days_consolid_early=consolid_start:midday_consolid; % take first half of data
        Params.PlotLearning.timeline.days_consolid_late=midday_consolid+1:consolid_end; % take second half of data
            
            
        else
            consolid_end=Params.PlotLearning.WNTimeOffInd; % if not defined, then use the last day of WN)
            
         duration_consolidation=consolid_end-consolid_start+1;

         Params.PlotLearning.timeline.days_consolid_early=consolid_start:consolid_end; % take all data btw consol start and last WN day
        end
        
    end
%     
%     % ====== BIDIR LEARNING PERIOD
%     if isfield(Params.PlotLearning.timeline, 'bidir_start_Inds')
%         
%         bidir_start=Params.PlotLearning.timeline.bidir_start_Inds; % start
%         if isfield(Params.PlotLearning.timeline, 'bidir_end_Inds')
%             bidir_end=Params.PlotLearning.timeline.bidir_end_Inds; % end (if defined)
%             
%             
%         else
%             bidir_end=Params.PlotLearning.WNTimeOffInd; % if not defined, then use the last day of WN)
%         end
%         
%         duration_bidir=bidir_end-bidir_start+1;
%         midday_bidir=bidir_start + ceil(duration_bidir/2)-1;
%
%         Params.PlotLearning.timeline.days_bidir_early=bidir_start:midday_bidir; % take first half of data
%         Params.PlotLearning.timeline.days_bidir_late=midday_bidir+1:bidir_end; % take second half of data
%
%     end

% ====== BIDIR LEARNING PERIOD
if isfield(Params.PlotLearning.timeline, 'bidir_start_Inds')
    
    bidir_start=Params.PlotLearning.timeline.bidir_start_Inds; % start
    if isfield(Params.PlotLearning.timeline, 'bidir_end_Inds')
        bidir_end=Params.PlotLearning.timeline.bidir_end_Inds; % end (if defined)
        duration_bidir=bidir_end-bidir_start+1;
        midday_bidir=bidir_start + ceil(duration_bidir/2)-1;
        
        Params.PlotLearning.timeline.days_bidir_early=bidir_start:midday_bidir; % take first half of data
        Params.PlotLearning.timeline.days_bidir_late=midday_bidir+1:bidir_end; % take second half of data
        
        
    else
        bidir_end=Params.PlotLearning.WNTimeOffInd; % if not defined, then use the last day of WN)
        
        duration_bidir=bidir_end-bidir_start+1;
        
        Params.PlotLearning.timeline.days_bidir_early=bidir_start:bidir_end; % take first half of data
        
    end
end

else
    disp('Skipping - no "timeline" in Params');
end


%% FILTER OUT RELEVANT MUSC AND PBS DATA (BASED ON TEMPORAL WINDOWS)

% == 1) extract relevant data (i.e. lags from muscimol switch)

for i=1:length(SylFields_Unique);
    syl=SylFields_Unique{i};
    
    AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow=cell([1, NumDays]);
    AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow=cell([1, NumDays]);
    
    
    for ii=1:NumDays;
        
        % check if day has data
        if isempty(AllDays_PlotLearning.DataMatrix.(syl).FFvals{ii});
            continue
        end
        
        % === EXTRACT DATA, for both PBS and MUSC
        FFvals_PBS=cell2mat(AllDays_PlotLearning.DataMatrix.(syl).FFvals{ii});
        Tvals_PBS=cell2mat(AllDays_PlotLearning.DataMatrix.(syl).Tvals{ii});
        
        FFvals_MUSC=cell2mat(AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals{ii});
        Tvals_MUSC=cell2mat(AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals{ii});
        
        
        % == for each condition (e.g. PBS or MUSC), filter out relevant data (musc
        % is lag after switch, PBS is in window before switch.
        for k=1:length(ExptCondition_codes);
            
            % -- For muscimol, get data after lag after start, to account
            % for lag
            if any(strcmp(ExptCondition_codes{k}, 'MUSC'));
                if ~isempty(FFvals_MUSC)
                    
                    % THROW OUT THIS DAY IF IS BAD MUSCIMOL
                    if any(Params.PlotLearning.MuscimolDaysToThrowout_Inds==ii);
                        continue
                    end
                    
                    % When did the switch happen?
                    MUSC_start=Params.PlotLearning.MuscimolSchedule_ByDayInds{ii}.start;
                    
                    
                    % sort data by time
                    [Tvals_MUSC, inds] = sort(Tvals_MUSC);
                    FFvals_MUSC=FFvals_MUSC(inds);
                    
                    % get the data coming after the lag
                    % convert time to hours
                    [~, DataTimes] = lt_convert_datenum_to_hour(Tvals_MUSC);
                    
                    IndsToKeep=find(DataTimes.hours>MUSC_start+Lag_time);
                    
                    % == extract the inds that pass lag time
                    Tvals_WithinTimeWindow=Tvals_MUSC(IndsToKeep);
                    FFvals_WithinTimeWindow=FFvals_MUSC(IndsToKeep);
                    
                    % ==== add some of PBS data coming right after switch,
                    % if that is desired
                    if Params.PlotLearning.Dur_of_PBS_dat_that_counts_as_MUSC>0;
                        [~, DataTimes] = lt_convert_datenum_to_hour(Tvals_PBS);
                        
                        pbs_hours=DataTimes.hours;
                        start_time=Params.PlotLearning.MuscimolSchedule_ByDayInds{ii}.end-0.2; % in case songs come right at switch
                        end_time=Params.PlotLearning.MuscimolSchedule_ByDayInds{ii}.end+...
                            Params.PlotLearning.Dur_of_PBS_dat_that_counts_as_MUSC;
                        
                        Inds_of_pbs_to_take=pbs_hours>start_time & pbs_hours<end_time;
                        
                        Tvals_WithinTimeWindow=[Tvals_WithinTimeWindow Tvals_PBS(Inds_of_pbs_to_take)];
                        FFvals_WithinTimeWindow=[FFvals_WithinTimeWindow FFvals_PBS(Inds_of_pbs_to_take)];
                    end
                    
                    
                    % == save those values
                    AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{ii}=FFvals_WithinTimeWindow;
                    AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{ii}=Tvals_WithinTimeWindow;
                    
                end
            end
            
            
            
            % -- FOR PBS, TAKE DATA IN WINDOW PRECEDING SWITCH
            if any(strcmp(ExptCondition_codes{k}, 'PBS'));
                
                % When did the switch happen?
                if length(Params.PlotLearning.MuscimolSchedule_ByDayInds)>=ii && ~isempty(Params.PlotLearning.MuscimolSchedule_ByDayInds{ii}) && ~any(Params.PlotLearning.MuscimolDaysToThrowout_Inds==ii); % today has switch, and is not day to thrwo out.
                    % if list of muscimol schedule contains today, and if today
                    % is not empty, and if I did not want to throw today out,
                    % then find the MUSC switch time.
                    MUSC_start=Params.PlotLearning.MuscimolSchedule_ByDayInds{ii}.start;
                else
                    % if no switch today, then use the median
                    % switch time
                    MUSC_start=Params.PlotLearning.MuscimolSchedule_MedianStartTime;
                end
                
                % sort data by time
                [Tvals_PBS, inds] = sort(Tvals_PBS);
                FFvals_PBS=FFvals_PBS(inds);
                
                % get the data coming in window 2 hrs before switch
                % convert time to hours
                [~, DataTimes] = lt_convert_datenum_to_hour(Tvals_PBS);
                
                IndsToKeep=find(DataTimes.hours>MUSC_start+PBS_window(1) & DataTimes.hours<MUSC_start+PBS_window(2));
                
                % == extract the inds that pass lag time
                Tvals_WithinTimeWindow=Tvals_PBS(IndsToKeep);
                FFvals_WithinTimeWindow=FFvals_PBS(IndsToKeep);
                
                % == save those values
                AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{ii}=FFvals_WithinTimeWindow;
                AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{ii}=Tvals_WithinTimeWindow;
            end
        end
    end
end

%% PREPROCESSING - what days actually have muscimol?

% what days have muscimol?
DaysWithMuscimol=[];
for i=1:NumDays;
    if length(AllDays_PlotLearning.DataMatrix_MUSC.(SylFields_Unique{1}).FFvals_WithinTimeWindow)>=i; % if last day doesn't have muscimol then it won't even be in the cell.
        
        if ~isempty(AllDays_PlotLearning.DataMatrix_MUSC.(SylFields_Unique{1}).FFvals_WithinTimeWindow{i});
            DaysWithMuscimol=[DaysWithMuscimol i];
        end
    end
end



%% PREPROCESSING - GET BASELINE (and deviation from baseline) stats for timte windowed muscimol and PBS

for i=1:length(SylFields_Unique);
    syl=SylFields_Unique{i};
    
    Yvals_musc=[];
    Yvals_pbs=[];
    for ii=Params.SeqFilter.BaselineDays;
        % == PBS
        try
            Yvals_pbs=[Yvals_pbs AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{ii}];
        catch err
        end
        
        % == musc
        try
            Yvals_musc=[Yvals_musc AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{ii}];
        catch err
        end
    end
    
    % === GET STATS
    % = PBS
    AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_WithinTimeWindow=Yvals_pbs;
    AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow=mean(Yvals_pbs);
    
    % == MUSC
    AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).rawFF_WithinTimeWindow=Yvals_musc;
    AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF_WithinTimeWindow=mean(Yvals_musc);
    
    % ===== FOR ALL DAYS, GET DEVIATION FROM BASELINE
    AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow=cell([1, NumDays]);
    AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow=cell([1, NumDays]);
    
    for ii=1:NumDays;
        
        % == PBS
        %             if ~isempty(AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{ii});
        AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{ii}=AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{ii} ...
            - AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
        
        AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow(ii)=mean(AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{ii});
        AllDays_PlotLearning.DataMatrix.(syl).semFF_DevFromBase_WithinTimeWindow(ii)=lt_sem(AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{ii});
        %             end
        
        % == MUSC
        %                     if
        %                     ~isempty(AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{ii});
        %                     % DONT HAVE EMPTY CLAUSE, NOW WILL HAVE NAN FOR DAYS
        %                     WITH NO DATA.
        %
        AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{ii}=AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{ii} ...
            - AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF_WithinTimeWindow;
        
        AllDays_PlotLearning.DataMatrix_MUSC.(syl).meanFF_DevFromBase_WithinTimeWindow(ii)=mean(AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{ii});
        AllDays_PlotLearning.DataMatrix_MUSC.(syl).semFF_DevFromBase_WithinTimeWindow(ii)=lt_sem(AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{ii});
        
        %                     end
    end
    
end





%% PLOT LEARNING, 1) raw data, day means, and CVs - one figure for each syllable

for i=1:length(SylFields_Unique);
    syl=SylFields_Unique{i};
    
    % ====== All data each syl a plot --> mult days (for both conditions)
    lt_figure; hold on; grid on;
    %         title(syl);
    %         xlabel('days'); ylabel('FF (hz)');
    
    
    % == 1) all renditions
    lt_subplot(2,2,1); hold on; grid on;
    title(syl);
    
    for ii=1:NumDays;
        
        if isempty(AllDays_PlotLearning.DataMatrix.(syl).FFvals{ii})
            % if day empty, skip
            continue
        end
        
        % ===================================== PBS
        Y=cell2mat(AllDays_PlotLearning.DataMatrix.(syl).FFvals{ii});
        X=cell2mat(AllDays_PlotLearning.DataMatrix.(syl).Tvals{ii});
        
        % convert time to days
        tmp=lt_convert_EventTimes_to_RelTimes(FirstDay, X);
        X=tmp.FinalValue;
        
        plot(X, Y, 'ob');
        
        % == ANNOTATE CATCH SONG STATUS (if data are available)
        if isfield(AllDays_PlotLearning.DataMatrix.(syl), 'CatchSongStatus');
            if ~isempty(AllDays_PlotLearning.DataMatrix.(syl).CatchSongStatus{ii});
                CatchSongvals=AllDays_PlotLearning.DataMatrix.(syl).CatchSongStatus{ii};
                
                % --- Plot (replot catch song vals in different shape
                X_catchsong=X(logical(CatchSongvals));
                Y_catchsong=Y(logical(CatchSongvals));
                
                plot(X_catchsong, Y_catchsong, 'ok', 'MarkerFaceColor','k','MarkerSize',8);
            end
        end
        
        
        % === OVERLAY with 1) data within windows for analysis, and 2)
        % lines indication swiches
        Y_windowed=AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{ii};
        X_windowed=AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{ii};
        
        % convert time to days
        if ~isempty(X_windowed)
            tmp=lt_convert_EventTimes_to_RelTimes(FirstDay, X_windowed);
            X_windowed=tmp.FinalValue;
            
            plot(X_windowed, Y_windowed, 'ob', 'MarkerFaceColor', 'b');
            
            % mean
            Ymean_windowed=mean(Y_windowed);
            Ysem_windowed=lt_sem(Y_windowed);
            errorbar(ceil(max(X))-rand/5, Ymean_windowed, Ysem_windowed, 'sb', 'MarkerFaceColor', 'b', 'MarkerSize', 9);
        end
        
        
        %     % plot mean
        %     Ymean=mean(Y);
        %     Ysem=lt_sem(Y);
        %     errorbar(ceil(max(X))-rand/5, Ymean, Ysem, 's', 'color', plotcols{k}, 'MarkerSize', 9);
        
        % ============================================== MUSC
        if ~isempty(AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals{ii})
            Y=cell2mat(AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals{ii});
            X=cell2mat(AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals{ii});
            
            % convert time to days
            tmp=lt_convert_EventTimes_to_RelTimes(FirstDay, X);
            X=tmp.FinalValue;
            
            plot(X, Y, 'or');
            
            
            % == ANNOTATE CATCH SONG STATUS (if data are available)
            if isfield(AllDays_PlotLearning.DataMatrix_MUSC.(syl), 'CatchSongStatus');
                if ~isempty(AllDays_PlotLearning.DataMatrix_MUSC.(syl).CatchSongStatus{ii});
                    CatchSongvals=AllDays_PlotLearning.DataMatrix_MUSC.(syl).CatchSongStatus{ii};
                    
                    % --- Plot (replot catch song vals in different
                    % shape)
                    X_catchsong=X(logical(CatchSongvals));
                    Y_catchsong=Y(logical(CatchSongvals));
                    
                    plot(X_catchsong, Y_catchsong, 'ok', 'MarkerFaceColor','k','MarkerSize',8);
                end
            end
            
            % ============= OVERLAY with 1) data within windows for analysis, and 2)
            % lines indication swiches
            Y_windowed=AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{ii};
            X_windowed=AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{ii};
            
            % convert time to days
            if ~isempty(X_windowed)
                tmp=lt_convert_EventTimes_to_RelTimes(FirstDay, X_windowed);
                X_windowed=tmp.FinalValue;
                
                plot(X_windowed, Y_windowed, 'or', 'MarkerFaceColor', 'r');
                
                % mean
                Ymean_windowed=mean(Y_windowed);
                Ysem_windowed=lt_sem(Y_windowed);
                errorbar(ceil(max(X))-rand/5, Ymean_windowed, Ysem_windowed, 'sr', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
                
            end
        end
        
        % ======================== Put lines indicating experiment switches
        if ~isempty(Params.PlotLearning.MuscimolSchedule_ByDayInds{ii})
            Xline=Params.PlotLearning.MuscimolSchedule_ByDayInds{ii}.start;
            Xline=ii+Xline/24; % convert to days;
            line([Xline Xline], ylim);
            
            Xline=Params.PlotLearning.MuscimolSchedule_ByDayInds{ii}.end;
            Xline=ii+Xline/24; % convert to days;
            line([Xline Xline], ylim);
        end
        
    end
    
    % ============================= 2) JUST MEANS, SEMS
    hsplot(1)=subplot(2,2,2); hold on; grid on;
    xlabel('days'); ylabel('FF (hz)');
    
    % --------------------- PBS
    Xall=[];
    Yall.mean=[];
    Yall.sem=[];
    Yall.cv=[];
    Yall.N=[];
    
    % - collect data
    for ii=1:NumDays;
        if isempty(AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{ii});
            continue
        end
        
        Yvals=AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{ii};
        
        Yall.mean=[Yall.mean mean(Yvals)];
        Yall.sem=[Yall.sem lt_sem(Yvals)];
        Yall.N=[Yall.N length(Yvals)];
        Ystd=std(Yvals);
        Yall.cv=[Yall.cv Ystd/mean(Yvals)];
        
        Xall=[Xall ii];
    end
    
    % == Plot all days
    shadedErrorBar(Xall, Yall.mean, Yall.sem, {'Color', 'b'},1);
    lt_plot(Xall, Yall.mean, {'Color', 'b'});
    ylabel('FF (hz) (SEM)');
    
    % == Plot N and CV
    hsplot(2)=lt_subplot(4,1,3); hold on;
    lt_plot(Xall, Yall.cv, {'Color', 'b'});
    ylabel('CV');
    
    hsplot(3)=lt_subplot(4,1,4); hold on;
    lt_plot(Xall, Yall.N, {'Color', 'b'});
    ylabel('N');
    
    % -------------------- MUSC
    subplot(2,2,2); hold on; grid on;
    
    Xall=[];
    Yall.mean=[];
    Yall.sem=[];
    Yall.cv=[];
    Yall.N=[];
    
    for ii=1:NumDays;
        
        if isempty(AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{ii});
            continue
        end
        
        Yvals=AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{ii};
        
        Yall.mean=[Yall.mean mean(Yvals)];
        Yall.sem=[Yall.sem lt_sem(Yvals)];
        Yall.N=[Yall.N length(Yvals)];
        Ystd=std(Yvals);
        Yall.cv=[Yall.cv Ystd/mean(Yvals)];
        
        Xall=[Xall ii];
    end
    
    % == Plot all days
    if length(Xall)>1;
        shadedErrorBar(Xall, Yall.mean, Yall.sem, {'Color', 'r'},1);
    end
    lt_plot(Xall, Yall.mean, {'Color', 'r'});
    ylabel('FF (hz) (SEM)');
    
    
    % == Plot N and CV
    lt_subplot(4,1,3); hold on;
    lt_plot(Xall, Yall.cv, {'Color', 'r'});
    ylabel('CV');
    
    lt_subplot(4,1,4); hold on;
    lt_plot(Xall, Yall.N, {'Color', 'r'});
    ylabel('N');
    
    % -- annotate for this syl
    linkaxes(hsplot,'x');
    %     lt_subtitle(syl);
end



%% === PLOT RELATIVE TO BASELINE, AND TAKING INTO ACCOUNT BASELINE - one
% plot for all syls.
% NOTE: TO CHANGE - recalculate minus baseline, using tempoarl windows.


for j=1:length(Params.SeqFilter.SylLists.FieldsInOrder)
    syllist=Params.SeqFilter.SylLists.FieldsInOrder{j};
    
    lt_figure; hold on;
    title('day means, PBS (lines), and MUSC (dots)');
    xlabel('days');
    ylabel('pitch, deviation from baseline (hz)');
    plotcols=lt_make_plot_colors(length(syllist), 0, 0);
    
    hplot=[];
    
    % -- DEVIATION FROM BASELINE
    for i=1:length(syllist);
        
        syl=syllist{i};
        
        % == PBS
        
        Y=AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow;
        Ysem=AllDays_PlotLearning.DataMatrix.(syl).semFF_DevFromBase_WithinTimeWindow;
        X=1:length(Y);
        
        %     shadedErrorBar(X, Y, Ysem, {'Color', plotcols{i}},1);
        hplot(i)=plot(X, Y,'-', 'Color', plotcols{i},'LineWidth',2);
        
        % == MUSC
        Y=AllDays_PlotLearning.DataMatrix_MUSC.(syl).meanFF_DevFromBase_WithinTimeWindow;
        Ysem=AllDays_PlotLearning.DataMatrix_MUSC.(syl).semFF_DevFromBase_WithinTimeWindow;
        X=1:length(Y);
        
        X=X(~isnan(Y));
        Y=Y(~isnan(Y));
        
        %     shadedErrorBar(X, Y, Ysem, {'Color', plotcols{i}},1);
        %     lt_plot(X, Y, {'Color', plotcols{i}});
        plot(X, Y,'--', 'Color', plotcols{i},'LineWidth',2);
        
        
        Fn_AnnotateWNLines(plotWNdays,ylim);
        lt_plot_zeroline;
    end
    
    legend(hplot, syllist);
end

% === REVERSION (i.e. MUSCIMOL data as percent of PBS data, first taking
% into account baseline effect of MUSCIMOL);

% -- Collect data
for i=1:length(SylFields_Unique);
    syl=SylFields_Unique{i};
    
    % - PBS data
    Y_pbs=AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow;
    
    % - muscimol data
    Y_musc=AllDays_PlotLearning.DataMatrix_MUSC.(syl).meanFF_DevFromBase_WithinTimeWindow;
    
    % == muscimol pitch as difference from PBS pitch (using baseline
    % subtraacted pitch)
    AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_MuscMinusPBS_MinusBase=Y_musc-Y_pbs;
    
    % == muscimol pitch as percent of PBS pitch (using baseline subtracted
    % data)
    AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_MuscDividePBS_MinusBase=Y_musc./Y_pbs;
    
    % == "AFP bias" - defined as PBS minus MUSC
    AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_PBSMinusMusc_MinusBase=Y_pbs-Y_musc;
    
    
end

% -- PLOT

plotcols=lt_make_plot_colors(length(SylFields_Unique),0,0);

lt_figure; hold on;
hplot=[];
for i=1:length(SylFields_Unique);
    syl=SylFields_Unique{i};
    
    % == PLOT (difference)
    lt_subplot(2,1,1); hold on;
    title('Pitch - muscimol minus PBS (taking into accnt baseline)');
    ylabel('Hz');
    xlabel('days');
    
    X=find(~isnan(AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_MuscMinusPBS_MinusBase));
    lt_plot(X, AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_MuscMinusPBS_MinusBase(X), {'LineStyle','-','Color',plotcols{i}});
    
    lt_plot_zeroline
    Fn_AnnotateWNLines(plotWNdays,ylim)
    
    % == PLOT (reversion) - important things are switching direction around
    % y=1; should mark things with low learning - reversion less
    % informative
    lt_subplot(2,1,2); hold on;
    title('Pitch - muscimol divided by PBS (taking into accnt baseline) (fill: z>0.75) (square: z>2)');
    ylabel('Hz');
    xlabel('days');
    
    X=find(~isnan(AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_MuscDividePBS_MinusBase));
    % plot all as open circles
    hplot(i)=lt_plot(X, AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_MuscDividePBS_MinusBase(X), {'LineStyle','-', 'MarkerFaceColor','none','Color',plotcols{i}, 'MarkerSize', 3});
    
    % plot only those that are showing actual learning (e.g. z-score >0.75
    % from baseline)
    % -- z>0.75
    X=find(abs(AllDays_PlotLearning.DataMatrix.(syl).meanFF_zscore)>0.75);
    if ~isempty(X);
        hplot(i)=lt_plot(X, AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_MuscDividePBS_MinusBase(X), {'LineStyle','-', 'MarkerFaceColor',plotcols{i},'Color',plotcols{i},'MarkerSize',7});
    end
    
    % -- z>2
    X=find(abs(AllDays_PlotLearning.DataMatrix.(syl).meanFF_zscore)>2);
    if ~isempty(X);
        hplot(i)=lt_plot(X, AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_MuscDividePBS_MinusBase(X), {'LineStyle','-', 'Marker','s', 'MarkerFaceColor',plotcols{i},'Color','k','MarkerSize',7});
    end
    
    ylim([-1 3]);
    
    line(xlim, [1 1],'LineStyle','--','Color','k');
    Fn_AnnotateWNLines(plotWNdays,ylim)
    
end

legend(hplot, SylFields_Unique);


%% == PLOT MUSCIMOL AND PBS ON SEPARATE PLOTS - IN PROGRESS

% ====== PBS ONLY

datamatfield='DataMatrix';

% ==================================== PBS
% PLOT pitch deviation and zscore - each with own axis.
lt_figure; hold on;
for j=1:length(Params.SeqFilter.SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
    FieldsList=Params.SeqFilter.SylLists.FieldsToPlot{j};
    
    plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
    
    % PLOT MEAN PITCH SUBTRACTING BASELINE
    lt_subplot(length(Params.SeqFilter.SylLists.FieldsToPlot),2,-1+j*2); hold on;
    h=[];
    for jj=1:length(FieldsList); % how many fields within this set?
        syl=FieldsList{jj}; % actual syl name (e.g. 'a')
        
        % compile pitch means and SEM for each day, for this syl.
        shadedErrorBar(1:length(AllDays_PlotLearning.(datamatfield).(syl).meanFF_DevFromBase_WithinTimeWindow), AllDays_PlotLearning.(datamatfield).(syl).meanFF_DevFromBase_WithinTimeWindow, ...
            AllDays_PlotLearning.(datamatfield).(syl).semFF,{'Color',plot_colors{jj},'LineWidth',1.5},1);
        
        h(jj)=plot(AllDays_PlotLearning.(datamatfield).(syl).meanFF_DevFromBase_WithinTimeWindow,'o','Color',plot_colors{jj},'MarkerFaceColor',plot_colors{jj},'MarkerSize',6);
    end
    
    legend(h,FieldsList)
    title(['Day mean Pitch (minus baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    ylabel('FF (hz) (rel to baseline)','FontSize',12,'FontWeight','bold')
    xlabel('days','FontSize',12,'FontWeight','bold')
    
    lt_plot_zeroline
    Fn_AnnotateWNLines(plotWNdays,ylim)
    
    
    %     % PLOT Z-score
    %     lt_subplot(length(Params.SeqFilter.SylLists.FieldsToPlot),2,j*2); hold on;
    %     h=[];
    %     for jj=1:length(FieldsList); % how many fields within this set?
    %         syl=FieldsList{jj}; % actual syl name (e.g. 'a')
    %
    %         % Plot
    %         plot(AllDays_PlotLearning.(datamatfield).(syl).meanFF_zscore,'-','Color',plot_colors{jj},'LineWidth',1.5)
    %         h(jj)=plot(AllDays_PlotLearning.(datamatfield).(syl).meanFF_zscore,'o','Color',plot_colors{jj},'MarkerFaceColor',plot_colors{jj},'MarkerSize',6);
    %     end
    %
    % annotate
    %     legend(h,FieldsList)
    %     title(['Z-scored mean pitch (relative to baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    %     ylabel('FF (Z-scored to baseline)','FontSize',12,'FontWeight','bold')
    %     xlabel('days','FontSize',12,'FontWeight','bold')
    %
    %     lt_plot_zeroline;
    %     Fn_AnnotateWNLines(plotWNdays,ylim)
    
end

lt_subtitle('PBS');


% =========================================== MUSC ONLY
datamatfield='DataMatrix_MUSC';

% PLOT pitch deviation and zscore - each with own axis.
lt_figure; hold on;
for j=1:length(Params.SeqFilter.SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
    FieldsList=Params.SeqFilter.SylLists.FieldsToPlot{j};
    
    plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
    
    % PLOT MEAN PITCH SUBTRACTING BASELINE
    lt_subplot(length(Params.SeqFilter.SylLists.FieldsToPlot),2,-1+j*2); hold on;
    h=[];
    for jj=1:length(FieldsList); % how many fields within this set?
        syl=FieldsList{jj}; % actual syl name (e.g. 'a')
        
        % compile pitch means and SEM for each day, for this syl.
        shadedErrorBar(1:length(AllDays_PlotLearning.(datamatfield).(syl).meanFF_DevFromBase_WithinTimeWindow), AllDays_PlotLearning.(datamatfield).(syl).meanFF_DevFromBase_WithinTimeWindow, ...
            real(AllDays_PlotLearning.(datamatfield).(syl).semFF_DevFromBase_WithinTimeWindow),{'Color',plot_colors{jj},'LineWidth',1.5},1);
        
        h(jj)=plot(AllDays_PlotLearning.(datamatfield).(syl).meanFF_DevFromBase_WithinTimeWindow,'o','Color',plot_colors{jj},'MarkerFaceColor',plot_colors{jj},'MarkerSize',6);
    end
    
    legend(h,FieldsList)
    title(['Day mean Pitch (minus baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
    ylabel('FF (hz) (rel to baseline)','FontSize',12,'FontWeight','bold')
    xlabel('days','FontSize',12,'FontWeight','bold')
    
    lt_plot_zeroline
    Fn_AnnotateWNLines(plotWNdays,ylim)
    
end

lt_subtitle('MUSC');

%% ONE FIGURE FOR EACH SYL, PLOT EACH DAY AS ONE POINT, AFP bias VS. GENERALIZATION/LEARNING

syllist=SylFields_Unique;
ThrowOutBaseline=1;

if (0); % ==== DON'T PLOT ALL INDIVIDUALLY, JUST PLOT ONE FIGURE FOR ALL.
    for i=1:length(syllist);
        syl=syllist{i};
        
        lt_figure; hold on;
        
        % === PLOT each day as one day
        
        Learning_array=AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow; % all days, FF minus baseline
        AFP_Bias_array=AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_PBSMinusMusc_MinusBase;
        
        % -- THROW OUT BASELINE?
        if ThrowOutBaseline==1;
            Learning_array(Params.SeqFilter.BaselineDays)=[];
            AFP_Bias_array(Params.SeqFilter.BaselineDays)=[];
        end
        
        % only plot days with both musc and pbs
        Learning_array=Learning_array(~isnan(AFP_Bias_array));
        AFP_Bias_array=AFP_Bias_array(~isnan(AFP_Bias_array));
        
        % -- PLOT
        lt_regress(AFP_Bias_array, Learning_array, 1, 0);
        
        title(syl);
        xlabel('Learning (hz minus baseline, PBS)');
        ylabel('AFP Bias (Pitch, PBS minus MUSC)');
        
        % lines
        line(xlim, [0 0], 'Color','b');
        line([0 0], ylim, 'Color','b');
        
    end
end

% ==================== REPLOT, ALL IN ONE FIGURE
lt_figure; hold on;
hplot=[];

for i=1:length(syllist);
    syl=syllist{i};
    
    hplot(i)=lt_subplot(ceil(length(syllist)/3), 3, i); hold on;
    
    % === PLOT each day as one day
    Learning_array=AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow; % all days, FF minus baseline
    AFP_Bias_array=AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_PBSMinusMusc_MinusBase;
    
    % -- THROW OUT BASELINE?
    if ThrowOutBaseline==1;
        Learning_array(Params.SeqFilter.BaselineDays)=[];
        AFP_Bias_array(Params.SeqFilter.BaselineDays)=[];
    end
    
    % get day inds to overlay on plot
    DayInds_array=find(~isnan(AFP_Bias_array));
    
    
    % only plot days with both musc and pbs
    Learning_array=Learning_array(~isnan(AFP_Bias_array));
    AFP_Bias_array=AFP_Bias_array(~isnan(AFP_Bias_array));
    
    % -- PLOT
    try
        lt_regress(AFP_Bias_array, Learning_array, 1, 0);
    catch err
    end
    
    % -- plot day inds
    for ii=1:length(DayInds_array);
        text(Learning_array(ii)+10, AFP_Bias_array(ii), num2str(DayInds_array(ii)),'FontSize',8, 'Color','r');
    end
    
    title(syl);
    %     xlabel('Learning (hz minus baseline, PBS)');
    %     ylabel('Reversion (Pitch, MUSC minus PBS)');
    
    % lines
    line(xlim, [0 0], 'Color','b');
    line([0 0], ylim, 'Color','b');
    
    % limits
    xlim([-250 150]);
    ylim([-100 100]);
end
lt_subtitle('AFP Bias (Pitch, PBS minus MUSC) vs Learning (hz minus baseline, PBS)');

linkaxes(hplot,'xy');



% ========== REPLOT, ALL IN ONE FIG, BUT WITH AFP BIAS AS PERCENT OF CURRENT
% LEARNING
% == REPLOT, ALL IN ONE FIGURE
lt_figure; hold on;
hplot=[];

for i=1:length(syllist);
    syl=syllist{i};
    
    hplot(i)=lt_subplot(ceil(length(syllist)/3), 3, i); hold on;
    
    % === PLOT each day as one day
    Learning_array=AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow; % all days, FF minus baseline
    AFP_Bias_array=AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_PBSMinusMusc_MinusBase;
    
    
    % -- THROW OUT BASELINE?
    if ThrowOutBaseline==1;
        Learning_array(Params.SeqFilter.BaselineDays)=[];
        AFP_Bias_array(Params.SeqFilter.BaselineDays)=[];
    end
    
    % get day inds to overlay on plot
    DayInds_array=find(~isnan(AFP_Bias_array));
    
    % only plot days with both musc and pbs
    Learning_array=Learning_array(~isnan(AFP_Bias_array));
    AFP_Bias_array=AFP_Bias_array(~isnan(AFP_Bias_array));
    
    % == normalize AFP bias by current learning (i.e. divide)
    AFP_Bias_array=AFP_Bias_array./Learning_array;
    
    % -- PLOT
    try
        lt_regress(AFP_Bias_array, Learning_array, 1, 0);
    catch err
    end
    
    % -- plot day inds
    for ii=1:length(DayInds_array);
        text(Learning_array(ii)+10, AFP_Bias_array(ii), num2str(DayInds_array(ii)),'FontSize',8, 'Color','r');
    end
    
    title(syl);
    %     xlabel('Learning (hz minus baseline, PBS)');
    %     ylabel('Reversion (Pitch, MUSC minus PBS)');
    
    % lines
    line(xlim, [0 0], 'Color','b');
    line([0 0], ylim, 'Color','b');
    
    % limits
    xlim([-250 150]);
    ylim([-1 1]);
end
lt_subtitle('AFP Bias (divided by current learning) vs Learning');

linkaxes(hplot,'xy');



%% ONE FIGURE FOR ONE DAY (OR DAY BIN), ALL SYLS, EACH AS A DOT, AFP VS. GENERALIZATION

lt_figure; hold on;
hsplot=[];
for i=1:length(DaysWithMuscimol);
    
    day=DaysWithMuscimol(i);
    
    hsplot(i)=lt_subplot(ceil(length(DaysWithMuscimol)/3),3,i); hold on;
    title(['day: ' num2str(day)]);
    
    % === Similar syls
    for i=1:length(Params.SeqFilter.SylLists.SylsSame);
        syl = Params.SeqFilter.SylLists.SylsSame{i};
        
        % collect stats
        learning=AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow(day); % just for today
        AFP_bias=AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_PBSMinusMusc_MinusBase(day);
        
        % plot one dot
        lt_plot(learning, AFP_bias, {'Color', 'b'});
        
        % plot text with syl ID
        text(learning+1, AFP_bias, syl,'FontSize',8, 'Color','b');
    end
    
    
    % === Different syls
    for i=1:length(Params.SeqFilter.SylLists.SylsDifferent);
        syl = Params.SeqFilter.SylLists.SylsDifferent{i};
        
        % collect stats
        learning=AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow(day); % just for today
        AFP_bias=AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_PBSMinusMusc_MinusBase(day);
        
        % plot one dot
        lt_plot(learning, AFP_bias, {'Color', 'r'});
        
        % plot text with syl ID
        text(learning+1, AFP_bias, syl,'FontSize',8, 'Color','r');
    end
    
    
    
    % === Target syl
    syl = Params.SeqFilter.SylLists.TargetSyls{1};
    
    % collect stats
    learning=AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow(day); % just for today
    AFP_bias=AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_PBSMinusMusc_MinusBase(day);
    
    % plot one dot
    lt_plot(learning, AFP_bias, {'Color', 'k'});
    
    % plot text with syl ID
    text(learning+1, AFP_bias, syl,'FontSize',8, 'Color','k');
    
    
    % === annotate this day
    % lines
    line(xlim, [0 0], 'Color','g');
    line([0 0], ylim, 'Color','g');
    
    xlim([-250 250]);
    ylim([-150 150]);
    
end
lt_subtitle('AFP bias vs. Learning - one plot each day');
linkaxes(hsplot,'xy');



% ============= REPLOT, but normalizing AFP bias to current learning

lt_figure; hold on;
hsplot=[];
for i=1:length(DaysWithMuscimol);
    
    day=DaysWithMuscimol(i);
    
    hsplot(i)=lt_subplot(ceil(length(DaysWithMuscimol)/3),3,i); hold on;
    title(['day: ' num2str(day)]);
    
    % === Similar syls
    for i=1:length(Params.SeqFilter.SylLists.SylsSame);
        syl = Params.SeqFilter.SylLists.SylsSame{i};
        
        % collect stats
        learning=AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow(day); % just for today
        AFP_bias=AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_PBSMinusMusc_MinusBase(day);
        
        AFP_bias=AFP_bias/learning;
        
        % plot one dot
        lt_plot(learning, AFP_bias, {'Color', 'b'});
        
        % plot text with syl ID
        text(learning+1, AFP_bias, syl,'FontSize',8, 'Color','b');
    end
    
    
    % === Different syls
    for i=1:length(Params.SeqFilter.SylLists.SylsDifferent);
        syl = Params.SeqFilter.SylLists.SylsDifferent{i};
        
        % collect stats
        learning=AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow(day); % just for today
        AFP_bias=AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_PBSMinusMusc_MinusBase(day);
        AFP_bias=AFP_bias/learning;
        
        % plot one dot
        lt_plot(learning, AFP_bias, {'Color', 'r'});
        
        % plot text with syl ID
        text(learning+1, AFP_bias, syl,'FontSize',8, 'Color','r');
    end
    
    
    
    % === Target syl
    syl = Params.SeqFilter.SylLists.TargetSyls{1};
    
    % collect stats
    learning=AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow(day); % just for today
    AFP_bias=AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_PBSMinusMusc_MinusBase(day);
    AFP_bias=AFP_bias/learning;
    
    % plot one dot
    lt_plot(learning, AFP_bias, {'Color', 'k'});
    
    % plot text with syl ID
    text(learning+1, AFP_bias, syl,'FontSize',8, 'Color','k');
    
    
    % === annotate this day
    % lines
    line(xlim, [0 0], 'Color','g');
    line([0 0], ylim, 'Color','g');
    
    xlim([-250 250]);
    ylim([-1 1]);
    
end
lt_subtitle('AFP bias (norm to laerning) vs. Learning - one plot each day');
linkaxes(hsplot,'xy');


%% === DOES THE GENERALIZATION TO SOME NON-TARGETS REFLECT CONSOLIDATION?
% 1) Supports if AFP bias does not seem to drive learning (as it does
% with target. Above scatter plots. Also, plot the correlation between
% AFP bias and current recent learning (or next day learning) - would
% expect target to show strong correlation, while non-targets wont.

% 2) Also, ask whether generalization correlates with consolidation -
% i.e. consolidation = Musc/PBS at target.  Can ask about MP-sepcific
% generalization (i.e. MUSC pitch) and PBS generalization.

% === MORE TO DO:
% 1) make sure all days have good MUSC time window
% 2) plot learning for MUSC and PBS separately.


%% BAR PLOT, ONE DAY, FROM MOST AFP INFLUENCE (E.G. TARGET) TO LEAST.  ARGUE FOR INCREASING SPECIFICITY ACROSS LEARNING?
% === TO DO: PLOT bar for 1) PBS, 2) MUSC together.

% === PLOT IN ORDER: [target, sim, diff]
syllist=[Params.SeqFilter.SylLists.TargetSyls{1} Params.SeqFilter.SylLists.SylsSame Params.SeqFilter.SylLists.SylsDifferent]; % order I want to plot syls in bar plot

lt_figure; hold on;
h=[];
% one plot for each day
for i=1:length(DaysWithMuscimol);
    day=DaysWithMuscimol(i);
    
    AFP_bias_array=[];
    
    for ii=1:length(syllist);
        
        syl=syllist{ii};
        
        AFP_bias_array(ii)=AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_PBSMinusMusc_MinusBase(day); % AFP bias
        
    end
    
    
    % == PLOT FOR THIS DAY
    h(i)=lt_subplot(ceil(length(DaysWithMuscimol)/3), 3, i); hold on;
    title(['day: ' num2str(day)]);
    
    bar(AFP_bias_array);
    
    set(h,'XTick',1:length(syllist));
    set(h, 'XTickLabel', syllist);
    
end
lt_subtitle('AFP bias for syls in order: [target, similar, diff]');
linkaxes(h, 'y')


% === PLOT IN ORDER: [motifs]
syllist=[];
for i=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
    syllist=[syllist Params.SeqFilter.SylLists.FieldsInOrder{i}];
end

lt_figure; hold on;
h=[];
% one plot for each day
for i=1:length(DaysWithMuscimol);
    day=DaysWithMuscimol(i);
    
    AFP_bias_array=[];
    
    for ii=1:length(syllist);
        
        syl=syllist{ii};
        
        AFP_bias_array(ii)=AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_PBSMinusMusc_MinusBase(day); % AFP bias
        
    end
    
    
    % == PLOT FOR THIS DAY
    h(i)=lt_subplot(ceil(length(DaysWithMuscimol)/3), 3, i); hold on;
    title(['day: ' num2str(day)]);
    
    bar(AFP_bias_array);
    
    set(h,'XTick',1:length(syllist));
    set(h, 'XTickLabel', syllist);
    
end
lt_subtitle('AFP bias for syls by motif: [motif1 motif2 ...]');
linkaxes(h, 'y')

% ========== NOTE
% the above indicates that syls move together (in terms of AFP bias).

%% === SCATTER OF AFP BIAS FOR NON-TARGETS VS. TARGET.  ONE PLOT FOR EACH DAY

% lt_figure; hold on;
%
%
%
% for i=DaysWithMuscimol;
%
%     lt_subplot(length(DaysWithMuscimol),3,i); hold on;
%     title(['day: ' num2str(i)]);
%
%     % ===
%
% end
%




%% PLOT DOTS comparing 1) PBS, 2) MUSC, for all syls.  Divide up by day
% Goal: Look at spread within days, and change across days

% ====================== BY DAY, ACTUAL FF, ALL PLOTS SAME AXIS.
plots_per_fig=9;
[num_figures, ~, ~]=lt_get_subplot_size(length(DaysWithMuscimol), plots_per_fig); % max 9 figures per plot

% Initiate corrent num of figures;
hfig=[];
for k=1:num_figures;
    hfig(k) = lt_figure;
end

% Get other things
plotcols=lt_make_plot_colors(length(SylFields_Unique), 0, 0);
count=1;

% Plot
for i=1:length(DaysWithMuscimol);
    day=DaysWithMuscimol(i);
    
    % ------------------------------ Decide which figure number to plot on
    fignum=ceil(count/plots_per_fig);
    figure(hfig(fignum)); % go to that figure;
    
    % decide which subplot to plot on
    subplotnum=mod(count, plots_per_fig);
    if subplotnum==0;
        subplotnum=plots_per_fig;
    end
    hsplot(count)=lt_subplot(ceil(plots_per_fig/3), 3, subplotnum); hold on;
    title(['day: ' num2str(day)]);
    
    count=count+1;
    % --------------------------------------------------------------------
    
    % ==== PLOT TARGET
    targsyl=Params.SeqFilter.SylLists.TargetSyls{1};
    
    FF_pbs_TARG=AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase_WithinTimeWindow(day);
    FF_musc_TARG=AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).meanFF_DevFromBase_WithinTimeWindow(day);
    
    lt_plot([1 2], [FF_pbs_TARG FF_musc_TARG], {'LineStyle', '-'});
    
    % ==== PLOT ALL SYLS
    for j=1:length(SylFields_Unique);
        syl=SylFields_Unique{j};
        
        % skip if this is the target
        if strcmp(syl, targsyl);
            continue
        end
        
        FF_pbs=AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow(day);
        FF_musc=AllDays_PlotLearning.DataMatrix_MUSC.(syl).meanFF_DevFromBase_WithinTimeWindow(day);
        
        % plot differently for similar vs. different syls
        if any(strcmp(syl, Params.SeqFilter.SylLists.SylsSame))
            % similar
            plot([1 2], [FF_pbs FF_musc], '-o','Color',plotcols{j});
        elseif any(strcmp(syl, Params.SeqFilter.SylLists.SylsDifferent))
            % different
            plot([1 2], [FF_pbs FF_musc], '-s','Color',plotcols{j});
        else
            plot([1 2], [FF_pbs FF_musc], '-.','Color',plotcols{j});
        end
    end
    
    % --- format plot
    xlim([-0.5 3.5]);
    lt_plot_zeroline;
    
end

% Give all figures same title
for i=1:num_figures;
    figure(hfig(i));
    lt_subtitle('PBS vs. MUSC, FF minus baseline');
end

linkaxes(hsplot, 'y');

disp(['DAYS TO MARK: ' num2str(cell2mat(Params.PlotLearning.DaysToMarkInds))])


% ======================= ALL DAYS, Normalized to maximum learning (of all
% days).
plots_per_fig=9;
[num_figures, ~, ~]=lt_get_subplot_size(length(DaysWithMuscimol), plots_per_fig); % max 9 figures per plot

% Initiate corrent num of figures;
hfig=[];
for k=1:num_figures;
    hfig(k) = lt_figure;
end

% Get other things
plotcols=lt_make_plot_colors(length(SylFields_Unique), 0, 0);
count=1;

% === Find maximum learning by target (over all days)
FF_pbs_targ_max=[];
for i=1:length(DaysWithMuscimol);
    day=DaysWithMuscimol(i);
    
    FF_pbs_TARG=AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase_WithinTimeWindow(day);
    FF_pbs_targ_max=[FF_pbs_targ_max FF_pbs_TARG];
end
[~, ind]=max(abs(FF_pbs_targ_max));
FF_pbs_targ_max=FF_pbs_targ_max(ind);


% Plot
for i=1:length(DaysWithMuscimol);
    day=DaysWithMuscimol(i);
    
    % ------------------------------ Decide which figure number to plot on
    fignum=ceil(count/plots_per_fig);
    figure(hfig(fignum)); % go to that figure;
    
    % decide which subplot to plot on
    subplotnum=mod(count, plots_per_fig);
    if subplotnum==0;
        subplotnum=plots_per_fig;
    end
    hsplot(count)=lt_subplot(ceil(plots_per_fig/3), 3, subplotnum); hold on;
    title(['day: ' num2str(day)]);
    
    count=count+1;
    % --------------------------------------------------------------------
    
    % ==== PLOT TARGET
    targsyl=Params.SeqFilter.SylLists.TargetSyls{1};
    
    FF_pbs_TARG=AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase_WithinTimeWindow(day);
    FF_musc_TARG=AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).meanFF_DevFromBase_WithinTimeWindow(day);
    
    FF_pbs_TARG=FF_pbs_TARG/FF_pbs_targ_max;
    FF_musc_TARG=FF_musc_TARG/FF_pbs_targ_max;
    
    lt_plot([1 2], [FF_pbs_TARG FF_musc_TARG], {'LineStyle', '-'});
    
    % ==== PLOT ALL SYLS
    for j=1:length(SylFields_Unique);
        syl=SylFields_Unique{j};
        
        % skip if this is the target
        if strcmp(syl, targsyl);
            continue
        end
        
        FF_pbs=AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow(day);
        FF_musc=AllDays_PlotLearning.DataMatrix_MUSC.(syl).meanFF_DevFromBase_WithinTimeWindow(day);
        
        FF_pbs=FF_pbs/FF_pbs_targ_max;
        FF_musc=FF_musc/FF_pbs_targ_max;
        
        % plot differently for similar vs. different syls
        if any(strcmp(syl, Params.SeqFilter.SylLists.SylsSame))
            % similar
            plot([1 2], [FF_pbs FF_musc], '-o','Color',plotcols{j});
        elseif any(strcmp(syl, Params.SeqFilter.SylLists.SylsDifferent))
            % different
            plot([1 2], [FF_pbs FF_musc], '-s','Color',plotcols{j});
        else
            plot([1 2], [FF_pbs FF_musc], '-.','Color',plotcols{j});
        end
    end
    
    % --- format plot
    xlim([-0.5 3.5]);
    ylim([-1 1]);
    lt_plot_zeroline;
    
end

% Give all figures same title
for i=1:num_figures;
    figure(hfig(i));
    lt_subtitle('PBS vs. MUSC, Fraction of max learning by target');
end

linkaxes(hsplot, 'y');

disp(['DAYS TO MARK: ' num2str(cell2mat(Params.PlotLearning.DaysToMarkInds))])



%% =================== CHUNKING DAYS (consolid and bidir, early and late)
% REQUIRES TIMELINE INPUTED IN PARAMS
if isfield(Params.PlotLearning, 'timeline');
    
    % -- WHICH epochs to plot?
    EpochField_list={};
    if isfield(Params.PlotLearning.timeline, 'days_consolid_early');
        if isfield(Params.PlotLearning.timeline, 'days_consolid_late');
            EpochField_list={'days_consolid_early', 'days_consolid_late'};
        else
            EpochField_list= {'days_consolid_early'};
        end
    end
    
    if isfield(Params.PlotLearning.timeline, 'days_bidir_early');
        EpochField_list=[EpochField_list {'days_bidir_early'}];
        
    end
    
    if isfield(Params.PlotLearning.timeline, 'days_bidir_late');
        EpochField_list=[EpochField_list {'days_bidir_late'}];
        
    end
    
    % --------------
    
    for k=1:length(EpochField_list);
        epoch_field=EpochField_list{k};
        
        % Run
        days_list=Params.PlotLearning.timeline.(epoch_field);
        days_list=intersect(DaysWithMuscimol, days_list); % get only days that have muscimol
        
        AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).days_list=days_list;
        
        % ==== COLLECT DATA ACROSS DAYS
        for i=1:length(SylFields_Unique);
            syl=SylFields_Unique{i};
            
            FFvals_pbs=[];
            Tvals_pbs=[];
            FFvals_musc=[];
            Tvals_musc=[];
            for ii=1:length(days_list);
                day=days_list(ii);
                
                % ======== USING DAY MEANS
                FF_pbs=AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow(day);
                FF_musc=AllDays_PlotLearning.DataMatrix_MUSC.(syl).meanFF_DevFromBase_WithinTimeWindow(day);
                
                AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).FF_pbs(ii)=FF_pbs;
                AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).FF_musc(ii)=FF_musc;
                
                % ========== USING ALL RENDS
                FFvals_pbs=[FFvals_pbs AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{day}];
                Tvals_pbs=[Tvals_pbs AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{day}];
                FFvals_musc=[FFvals_musc AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{day}];
                Tvals_musc=[Tvals_musc AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day}];
            end
            
            % ==== STORE ALL RENDS DATA
            AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).ALLRENDS.FFvals_pbs=FFvals_pbs;
            AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).ALLRENDS.Tvals_pbs=Tvals_pbs;
            AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).ALLRENDS.FFvals_musc=FFvals_musc;
            AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).ALLRENDS.Tvals_musc=Tvals_musc;
            
            
            %             % --- mean across days (using day means)
            %             AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).meanFF_pbs=nanmean(AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).FF_pbs);
            %             AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).semFF_pbs=lt_sem(AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).FF_pbs);
            %
            %             AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).meanFF_musc=nanmean(AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).FF_musc);
            %             AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).semFF_musc=lt_sem(AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).FF_musc);
            
            % --- mean across days (using all rends)
            AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).meanFF_pbs=nanmean(AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).ALLRENDS.FFvals_pbs);
            AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).semFF_pbs=lt_sem(AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).ALLRENDS.FFvals_pbs);
            
            AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).meanFF_musc=nanmean(AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).ALLRENDS.FFvals_musc);
            AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).semFF_musc=lt_sem(AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).ALLRENDS.FFvals_musc);
            
        end
    end
    
%     lt_figure; hold on;
%     hsplot=[];
%     % ==== PLOT (USING DAY MEANS)
%     for k=1:length(EpochField_list);
%         epoch_field=EpochField_list{k};
%         
%         % --- TARGET
%         hsplot(k)=lt_subplot(1, length(EpochField_list), k); hold on;
%         title(epoch_field);
%         
%         targsyl=Params.SeqFilter.SylLists.TargetSyls{1};
%         
%         meanFF_pbs=nanmean(AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(targsyl).FF_pbs);
%         meanFF_musc=nanmean(AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(targsyl).FF_musc);
%         
%         lt_plot([1 2], [meanFF_pbs meanFF_musc], {'LineStyle', '-', 'MarkerSize', 7});
%         
%         hplot=[];
%         SylsPlotted={};
%         % --- ALL SYLS
%         for j=1:length(SylFields_Unique);
%             syl=SylFields_Unique{j};
%             
%             % skip if this is the target
%             if strcmp(syl, targsyl);
%                 continue
%             end
%             
%             SylsPlotted=[SylsPlotted {syl}];
%             
%             meanFF_pbs=nanmean(AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).FF_pbs);
%             meanFF_musc=nanmean(AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).FF_musc);
%             
%             % plot differently for similar vs. different syls
%             if any(strcmp(syl, Params.SeqFilter.SylLists.SylsSame))
%                 % similar
%                 hplot(j)=plot([1 2], [meanFF_pbs meanFF_musc], '-o','Color',plotcols{j}, 'MarkerSize', 7);
%             elseif any(strcmp(syl, Params.SeqFilter.SylLists.SylsDifferent))
%                 % different
%                 hplot(j)=plot([1 2], [meanFF_pbs meanFF_musc], '-s','Color',plotcols{j}, 'MarkerSize', 7);
%             else
%                 hplot(j)=plot([1 2], [meanFF_pbs meanFF_musc], '-d','Color',plotcols{j}, 'MarkerSize', 7);
%             end
%         end
%         
%         % --- format plot
%         xlim([-0.5 3.5]);
%         lt_plot_zeroline;
%         
%     end
%     
%     lt_subtitle('mean of day means');
%     
%     linkaxes(hsplot, 'y');
%     try
%     legend(hplot, SylFields_Unique)
%     catch err
%     end
    
    % ==== PLOT (USING MEAN OF RENDS)
    lt_figure; hold on;
    hsplot=[];
    
    for k=1:length(EpochField_list);
        epoch_field=EpochField_list{k};
        
        % --- TARGET
        hsplot(k)=lt_subplot(1, length(EpochField_list), k); hold on;
        title(epoch_field);
        
        targsyl=Params.SeqFilter.SylLists.TargetSyls{1};
        
        ffvals_pbs=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(targsyl).ALLRENDS.FFvals_pbs;
        meanFF_pbs=mean(ffvals_pbs);
        semff_pbs=lt_sem(ffvals_pbs);
        
        ffvals_musc=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(targsyl).ALLRENDS.FFvals_musc;
        meanFF_musc=mean(ffvals_musc);
        semFF_musc=lt_sem(ffvals_musc);
        
        errorbar([1 2], [meanFF_pbs meanFF_musc], [semff_pbs semFF_musc], 'o-', 'Color','k','MarkerFaceColor','k','MarkerSize', 8);
        
        hplot=[];
        SylsPlotted={};
        
        % --- ALL SYLS
        for j=1:length(SylFields_Unique);
            syl=SylFields_Unique{j};
            
            % skip if this is the target
            if strcmp(syl, targsyl);
                continue
            end
            
            SylsPlotted=[SylsPlotted {syl}];
            
            ffvals_pbs=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).ALLRENDS.FFvals_pbs;
            meanFF_pbs=mean(ffvals_pbs);
            semff_pbs=lt_sem(ffvals_pbs);
            
            ffvals_musc=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).ALLRENDS.FFvals_musc;
            meanFF_musc=mean(ffvals_musc);
            semFF_musc=lt_sem(ffvals_musc);
            
            % plot differently for similar vs. different syls
            if any(strcmp(syl, Params.SeqFilter.SylLists.SylsSame))
                % similar
                hplot(j)=errorbar([1 2], [meanFF_pbs meanFF_musc], [semff_pbs semFF_musc], 'o-', 'Color',plotcols{j}, 'MarkerFaceColor', plotcols{j}, 'MarkerSize', 7);
            elseif any(strcmp(syl, Params.SeqFilter.SylLists.SylsDifferent))
                % different
                hplot(j)=errorbar([1 2], [meanFF_pbs meanFF_musc], [semff_pbs semFF_musc], 's-', 'Color',plotcols{j},'MarkerSize', 7);
            else
                hplot(j)=errorbar([1 2], [meanFF_pbs meanFF_musc], [semff_pbs semFF_musc], 'd-', 'Color',plotcols{j},'MarkerSize', 7);
            end
        end
        
        % --- format plot
        xlim([-0.5 3.5]);
        lt_plot_zeroline;
        
    end
    
    lt_subtitle('mean of all rends');
    try
    linkaxes(hsplot, 'y');
    
    legend(hplot, SylFields_Unique)
    catch err
    end
    
end

%% PLOT IN ORDER OF MOTIFS (bar plot, one for each day)

count=1;
SubplotsPerFig=6;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

for i=1:length(DaysWithMuscimol);
    dayind=DaysWithMuscimol(i);
    
    [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
    
    % =================== PBS
    FFmean_All=[];
    FFsem_All=[];
    Syllables={};
    for j=1:length(SylFields_Unique);
        syl=SylFields_Unique{j};
        
        ffvals=AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{dayind};
        
        ffmean=mean(ffvals);
        ffsem=lt_sem(ffvals);
        
        % --- collect
        FFmean_All(j)=ffmean;
        FFsem_All(j)=ffsem;
        
        if length(syl)==1;
            Syllables=[Syllables syl];
        else
            syl_lower=regexp(syl, '[A-Z]', 'match');
            Syllables=[Syllables lower(syl_lower)];
        end
    end
    
    % --- PLOT
    [hBar hErrorbar]=barwitherr(FFsem_All, [1:length(FFmean_All)]-0.2,FFmean_All);
    set(hBar, 'FaceColor', 'b');
    set(hBar, 'BarWidth', 0.5);
    
    hold on;
    
    % ==================== MUSC
    FFmean_All=[];
    FFsem_All=[];
    for j=1:length(SylFields_Unique);
        syl=SylFields_Unique{j};
        
        ffvals=AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{dayind};
        
        ffmean=mean(ffvals);
        ffsem=lt_sem(ffvals);
        
        FFmean_All(j)=ffmean;
        FFsem_All(j)=ffsem;
        
    end
    
    % --- PLOT
    [hBar hErrorbar]=barwitherr(FFsem_All, [1:length(FFmean_All)]+0.2, FFmean_All);
    set(hBar, 'FaceColor', 'r');
    set(hBar, 'BarWidth', 0.5);
    
    
    % ==================== GLOBAL
    set(gca, 'XTick', 1:length(Syllables));
    set(gca, 'XTickLabel', Syllables);
    title(['day ' num2str(dayind)]);
    
end
lt_subtitle('PBS (blue) vs. MUSC (red)')



%% PLOT BARS, but for consolidation

if isfield(Params.PlotLearning.timeline, 'consolid_start');
    lt_figure; hold on;
    
    % ===== CONSOLIDATION EARLY
    datafield='days_consolid_early';
    % ----------------- PBS
    lt_subplot(1,2,1); hold on;
    Syllables={};
    FFmean_All=[];
    FFsem_All=[]; % across syllables
    for j=1:length(SylFields_Unique);
        syl=SylFields_Unique{j};
        
        ffmean=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(datafield).(syl).meanFF_pbs;
        ffsem=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(datafield).(syl).semFF_pbs;
        
        % --colllect across syls
        FFmean_All(j)=ffmean;
        FFsem_All(j)=ffsem;
        
        % get syllable ID
        if length(syl)==1;
            Syllables=[Syllables syl];
        else
            syl_lower=regexp(syl, '[A-Z]', 'match');
            Syllables=[Syllables lower(syl_lower)];
        end
        
        % --- PLOT
        [hBar hErrorbar]=barwitherr(FFsem_All, [1:length(FFmean_All)]-0.2, FFmean_All);
        set(hBar, 'FaceColor', 'b');
        set(hBar, 'BarWidth', 0.5);
        title('early consolidation');
        
    end
    
    % ------------------ MUSC
    Syllables={};
    FFmean_All=[];
    FFsem_All=[]; % across syllables
    for j=1:length(SylFields_Unique);
        syl=SylFields_Unique{j};
        
        ffmean=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(datafield).(syl).meanFF_musc;
        ffsem=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(datafield).(syl).semFF_musc;
        
        % --colllect across syls
        FFmean_All(j)=ffmean;
        FFsem_All(j)=ffsem;
        
        % get syllable ID
        if length(syl)==1;
            Syllables=[Syllables syl];
        else
            syl_lower=regexp(syl, '[A-Z]', 'match');
            Syllables=[Syllables lower(syl_lower)];
        end
        
        % --- PLOT
        hold on;
        [hBar hErrorbar]=barwitherr(FFsem_All, [1:length(FFmean_All)]+0.2, FFmean_All);
        set(hBar, 'FaceColor', 'r');
        set(hBar, 'BarWidth', 0.5);
        hold on;
    end
    hold on;
    set(gca, 'XTick', 1:length(Syllables));
    set(gca, 'XTickLabel', Syllables);
    
    % ===== CONSOLIDATION LATE
    datafield='days_consolid_late';
    if isfield(Params.PlotLearning.timeline, datafield)
        % ----------------- PBS
        lt_subplot(1,2,2); hold on;
        Syllables={};
        FFmean_All=[];
        FFsem_All=[]; % across syllables
        for j=1:length(SylFields_Unique);
            syl=SylFields_Unique{j};
            
            ffmean=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(datafield).(syl).meanFF_pbs;
            ffsem=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(datafield).(syl).semFF_pbs;
            
            % --colllect across syls
            FFmean_All(j)=ffmean;
            FFsem_All(j)=ffsem;
            
            % get syllable ID
            if length(syl)==1;
                Syllables=[Syllables syl];
            else
                syl_lower=regexp(syl, '[A-Z]', 'match');
                Syllables=[Syllables lower(syl_lower)];
            end
            
            % --- PLOT
            [hBar hErrorbar]=barwitherr(FFsem_All, [1:length(FFmean_All)]-0.2, FFmean_All);
            set(hBar, 'FaceColor', 'b');
            set(hBar, 'BarWidth', 0.5);
            
        end
        
        % ------------------ MUSC
        Syllables={};
        FFmean_All=[];
        FFsem_All=[]; % across syllables
        for j=1:length(SylFields_Unique);
            syl=SylFields_Unique{j};
            
            ffmean=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(datafield).(syl).meanFF_musc;
            ffsem=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(datafield).(syl).semFF_musc;
            
            % --colllect across syls
            FFmean_All(j)=ffmean;
            FFsem_All(j)=ffsem;
            
            % get syllable ID
            if length(syl)==1;
                Syllables=[Syllables syl];
            else
                syl_lower=regexp(syl, '[A-Z]', 'match');
                Syllables=[Syllables lower(syl_lower)];
            end
            
            % --- PLOT
            hold on;
            [hBar hErrorbar]=barwitherr(FFsem_All, [1:length(FFmean_All)]+0.2, FFmean_All);
            set(hBar, 'FaceColor', 'r');
            set(hBar, 'BarWidth', 0.5);
            hold on;title('late consolidation');
            
        end
        hold on;
        set(gca, 'XTick', 1:length(Syllables));
        set(gca, 'XTickLabel', Syllables);
    end
    
    % % ==================== GLOBAL
    lt_subtitle('rends mean, across days')
    
end

%% PLOT BARS, but for bidir

if isfield(Params.PlotLearning.timeline, 'days_bidir_early');
    lt_figure; hold on;
    
    % ===== BIDIR EARLY
    datafield='days_bidir_early';
    % ----------------- PBS
    lt_subplot(1,2,1); hold on;
    Syllables={};
    FFmean_All=[];
    FFsem_All=[]; % across syllables
    for j=1:length(SylFields_Unique);
        syl=SylFields_Unique{j};
        
        ffmean=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(datafield).(syl).meanFF_pbs;
        ffsem=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(datafield).(syl).semFF_pbs;
        
        % --colllect across syls
        FFmean_All(j)=ffmean;
        FFsem_All(j)=ffsem;
        
        % get syllable ID
        if length(syl)==1;
            Syllables=[Syllables syl];
        else
            syl_lower=regexp(syl, '[A-Z]', 'match');
            Syllables=[Syllables lower(syl_lower)];
        end
        
        % --- PLOT
        [hBar hErrorbar]=barwitherr(FFsem_All, [1:length(FFmean_All)]-0.2, FFmean_All);
        set(hBar, 'FaceColor', 'b');
        set(hBar, 'BarWidth', 0.5);
        
    end
    
    % ------------------ MUSC
    Syllables={};
    FFmean_All=[];
    FFsem_All=[]; % across syllables
    for j=1:length(SylFields_Unique);
        syl=SylFields_Unique{j};
        
        ffmean=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(datafield).(syl).meanFF_musc;
        ffsem=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(datafield).(syl).semFF_musc;
        
        % --colllect across syls
        FFmean_All(j)=ffmean;
        FFsem_All(j)=ffsem;
        
        % get syllable ID
        if length(syl)==1;
            Syllables=[Syllables syl];
        else
            syl_lower=regexp(syl, '[A-Z]', 'match');
            Syllables=[Syllables lower(syl_lower)];
        end
        
        % --- PLOT
        hold on;
        [hBar hErrorbar]=barwitherr(FFsem_All, [1:length(FFmean_All)]+0.2, FFmean_All);
        set(hBar, 'FaceColor', 'r');
        set(hBar, 'BarWidth', 0.5);
        hold on;
    end
    hold on;
    title('early bidir');
    set(gca, 'XTick', 1:length(Syllables));
    set(gca, 'XTickLabel', Syllables);
end

% ===== BIDIR LATE
if isfield(Params.PlotLearning.timeline, 'days_bidir_late');
    
    datafield='days_bidir_late';
    % ----------------- PBS
    lt_subplot(1,2,2); hold on;
    Syllables={};
    FFmean_All=[];
    FFsem_All=[]; % across syllables
    for j=1:length(SylFields_Unique);
        syl=SylFields_Unique{j};
        
        ffmean=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(datafield).(syl).meanFF_pbs;
        ffsem=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(datafield).(syl).semFF_pbs;
        
        % --colllect across syls
        FFmean_All(j)=ffmean;
        FFsem_All(j)=ffsem;
        
        % get syllable ID
        if length(syl)==1;
            Syllables=[Syllables syl];
        else
            syl_lower=regexp(syl, '[A-Z]', 'match');
            Syllables=[Syllables lower(syl_lower)];
        end
        
        % --- PLOT
        [hBar hErrorbar]=barwitherr(FFsem_All, [1:length(FFmean_All)]-0.2, FFmean_All);
        set(hBar, 'FaceColor', 'b');
        set(hBar, 'BarWidth', 0.5);
        
    end
    
    % ------------------ MUSC
    Syllables={};
    FFmean_All=[];
    FFsem_All=[]; % across syllables
    for j=1:length(SylFields_Unique);
        syl=SylFields_Unique{j};
        
        ffmean=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(datafield).(syl).meanFF_musc;
        ffsem=AllDays_PlotLearning.EpochData.PBS_and_MUSC.(datafield).(syl).semFF_musc;
        
        % --colllect across syls
        FFmean_All(j)=ffmean;
        FFsem_All(j)=ffsem;
        
        % get syllable ID
        if length(syl)==1;
            Syllables=[Syllables syl];
        else
            syl_lower=regexp(syl, '[A-Z]', 'match');
            Syllables=[Syllables lower(syl_lower)];
        end
        
        % --- PLOT
        hold on;
        [hBar hErrorbar]=barwitherr(FFsem_All, [1:length(FFmean_All)]+0.2, FFmean_All);
        set(hBar, 'FaceColor', 'r');
        set(hBar, 'BarWidth', 0.5);
        hold on;
        title('late bidir');
        
    end
    hold on;
    set(gca, 'XTick', 1:length(Syllables));
    set(gca, 'XTickLabel', Syllables);
    lt_subtitle('rends mean, across days')
    
end
%% SAVE
if saveON==1;
    timestampSv=lt_get_timestamp(0);
    cd(Params.PlotLearning.savedir);
    
    save('Params','Params');
    save('AllDays_PlotLearning','AllDays_PlotLearning');
    
    % write a text file that tells you when files were made
    fid1=fopen(['DONE_PlotLearning_Musc_' timestampSv '.txt'],'w');
    fclose(fid1);
    
    try
        cd FIGURES/PlotLearning_Musc
    catch err
        mkdir FIGURES/PlotLearning_Musc;
        cd FIGURES/PlotLearning_Musc;
    end
    
    try
    lt_save_all_figs
    catch err
    end
    
    cd ../../
    
end



end

function Fn_AnnotateWNLines(plotWNdays,ylim)

global WNTimeOnInd
global WNTimeOffInd
global DaysToMarkInds

if plotWNdays==1;
    line([WNTimeOnInd-0.5 WNTimeOnInd-0.5],ylim,'LineStyle','--','Color','r'); % minus 0.5 since the datapoint is on the day, so want the line to be before datapojnt.
    line([WNTimeOffInd+0.5 WNTimeOffInd+0.5],ylim,'LineStyle','--','Color','r')
    
    for i=1:length(DaysToMarkInds);
        line([DaysToMarkInds{i} DaysToMarkInds{i}],ylim,'LineStyle','--','Color','k');
    end
end
end
