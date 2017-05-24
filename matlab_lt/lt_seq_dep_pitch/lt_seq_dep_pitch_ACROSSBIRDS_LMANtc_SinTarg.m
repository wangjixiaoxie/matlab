function [OUTPUT DATSTRUCT_FirstWNDayIndEqualsOne] = lt_seq_dep_pitch_ACROSSBIRDS_LMANTimeCourse(SeqDepPitch_AcrossBirds, Params, OnlyConsolPeriod_bidir, BinRelConsolDay1, TakeAverageWithinExpt, NumDaysInBin, GetAutoConsolWindows)
%% LT 3/4/16 - modified from LMANtimecourse, now for only the one targ phase

NumBirds=length(SeqDepPitch_AcrossBirds.birds);
%     NumDaysInBin=2;
BinToAnalyze=2;

%% == plot mean day for start of matinained shift period


tmp=length(Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart)
inds=3:3:tmp;
ConsolWindows=Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart(inds);

lt_figure; hold on;

tmp=cell2mat(ConsolWindows);

onsetDays=tmp(1:2:end);
offsetDays=tmp(2:2:end);
durations=offsetDays-onsetDays+1

lt_plot_annotation(1, ['onset (std) = ' num2str(mean(onsetDays)) '(' num2str(std(onsetDays)) ')']);
lt_plot_annotation(2, ['offset (std) = ' num2str(mean(offsetDays)) '(' num2str(std(offsetDays)) ')']);
lt_plot_annotation(3, ['durations (std) = ' num2str(mean(durations)) '(' num2str(std(durations)) ')']);


%% DISPLAY STATS FOR WHEN THERE WERE INACTIVATION DAYS FOR ALL BIRDS
disp(' ');
disp('WN day 1  --  Days with musc (WNday1=1)  --  Bidir day 1(WNday1=1)');

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        days_with_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
        WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        days_with_musc_relWNday1=days_with_musc-WNday1+1;
        
        bidir_day1=nan;
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart_Ind');
            bidir_day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind+1;
            bidir_day1=bidir_day1-WNday1+1;
        end
        
        %
        %         if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_start_Inds');
        %             bidir_day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_start_Inds(1);
        %             bidir_day1=bidir_day1-WNday1+1;
        %         end
        %
        
        disp([birdname '-' exptname ':  ' num2str(WNday1) ' -- ' mat2str(days_with_musc_relWNday1) ' -- ' num2str(bidir_day1)]);
        
    end
end


%% PLOT TIMELINE FOR ALL EXPERIMENTS [IN PROGRESS, ADDING REVERSION AT TARGET, SEE LINE WITH: FFminusBase_MUSC]
% EACH ROW IS ONE EXPERIMENT, ALL LOCKED TO WN ON. SHOW EACH INACTIVATION
% DAY. SHOW BIDIR. SHOW CV DROP AND REVERSION AT TARGET, AND PERCENT MAX
% LEARNING BY TARGET

% WN day 1 =  ind1

lt_figure; hold on;
BirdnamesExptname_all={};
title('Timeline (muscimol in bold) (samedir in blue)');
xlabel('day (WN day 1 = 1)');



CountExpt=1;
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        days_with_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
        WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        days_with_musc_relWNday1=days_with_musc-WNday1+1;
        
        
               
        % --- same dir start and end
        samedir_day1=nan;
        samedir_lastday=nan;
        
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'SameDir_OneDayBeforeStart');
            samedir_day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.SameDir_OneDayBeforeStart_Ind+1;
            samedir_day1=samedir_day1-WNday1+1;
            
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'SameDir_LastDay_Ind');
                samedir_lastday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.SameDir_LastDay_Ind;
            else
                samedir_lastday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            end
            samedir_lastday=samedir_lastday-WNday1+1;
            
        end
        
        %         if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_start_Inds');
        %             bidir_day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_start_Inds(1);
        %             bidir_day1=bidir_day1-WNday1+1;
        %
        %             if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_end_Inds');
        %                 bidir_lastday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_end_Inds(end);
        %             else
        %                 bidir_lastday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
        %             end
        %             bidir_lastday=bidir_lastday-WNday1+1;
        %         end
        
        % -- last day
            WN_last=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            WN_last=WN_last-WNday1+1;
        
        % ---- first day (relative to WN on day)
        firstday=-(WNday1-2);
        
        
        % ======= REVERSION AT TARGET
        FFminusBase_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).meanFF_DevFromBase_WithinTimeWindow;
        FFminusBase_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase_WithinTimeWindow;
        
        
        % === plot
        % - line for all days with any data
        X=firstday:WN_last;
        Y=CountExpt*ones(1,length(X));
        
        plot(X, Y, '-ok', 'LineWidth', 1);
        
        
        
        % --- line for end of WN
        line([WN_last+0.5 WN_last+0.5], [CountExpt-0.4 CountExpt+0.4], 'Color','k');

        
        % -- line for start and end of bidir
%         line([bidir_day1-0.5 bidir_day1-0.5], [CountExpt-0.4 CountExpt+0.4]);
%         line([bidir_lastday+0.5 bidir_lastday+0.5], [CountExpt-0.4 CountExpt+0.4]);
        
        line([samedir_day1-0.5 samedir_day1-0.5], [CountExpt-0.4 CountExpt+0.4]);
        line([samedir_lastday+0.5 samedir_lastday+0.5], [CountExpt-0.4 CountExpt+0.4]);
        
        
        
        % -- fill in musc days
        X=days_with_musc_relWNday1;
        Y=CountExpt*ones(1, length(X));
        
        plot(X, Y, 'ok', 'MarkerFaceColor','k')
        
        
        
        % -- collect bird name and expt
        BirdnamesExptname_all=[BirdnamesExptname_all [birdname(1:4) '-' exptname(end-4:end)]];
        
        CountExpt=CountExpt+1;
    end
end

set(gca, 'YTick', 1:CountExpt-1);
set(gca, 'YTickLabel', BirdnamesExptname_all);


% ---- WN on line
line([0.5 0.5], ylim, 'Color','r');



%%  for all days plot shift + reversion at target [IN PROGRESS - DOES NOTHING]

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        days_with_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
        WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        days_with_musc_relWNday1=days_with_musc-WNday1+1;
        
        % --- same dir start and end
        samedir_day1=nan;
        samedir_lastday=nan;
        
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'SameDir_OneDayBeforeStart');
            samedir_day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.SameDir_OneDayBeforeStart_Ind+1;
            samedir_day1=samedir_day1-WNday1+1;
            
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'SameDir_LastDay_Ind');
                samedir_lastday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.SameDir_LastDay_Ind;
            else
                samedir_lastday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            end
            samedir_lastday=samedir_lastday-WNday1+1;
            
        end
        
        
        % -- last day
            WN_last=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            WN_last=WN_last-WNday1+1;
        
        % ---- first day (relative to WN on day)
        firstday=-(WNday1-2);
        
        
        % -- collect bird name and expt
        BirdnamesExptname_all=[BirdnamesExptname_all [birdname(1:4) '-' exptname(end-4:end)]];
        
        CountExpt=CountExpt+1;
    end
end



%% [EXTRACT[ ===== DURING CONSOLIDATION PERIOD OF SAME-DIR, SEE CONSOLIDATION?

disp(' ----- ')

ExptCount=1;

DATSTRUCT_FirstWNDayIndEqualsOne=struct;

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        % ==== GET ALL DAYS, MAKE THEM RELATIVE TO BIDIR DAY 1 [getting
        % dates]
        days_with_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
        WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        WNlastday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
        
        
        days_with_musc_RelToWNDay1=days_with_musc-WNday1+1;
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION, 'MultiDirSyls');
            othertarg=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls{2};
        elseif isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION, 'SameDirSyls');
            othertarg= SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SameDirSyls{2};
        else
            othertarg=[];
        end
        
        % -- only keep othertarg if it is same type;
        if ~isempty(othertarg)
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(othertarg).similar_to_targ==0
                othertarg=[];
            end
        end
           
        % --------- FIGURE OUT WHAT ARE SAME-TYPES
        SameTypeSyls = [SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STDS ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STSS];
        
       
    
    
        % ==== collect data
        
        DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(ExptCount).FFminusBase_Mean_PBS=nan(1, WNlastday-WNday1+1);
        DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(ExptCount).FFminusBase_Mean_MUSC=nan(1, WNlastday-WNday1+1);
        
        DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(ExptCount).HitRate_PBS=nan(1, WNlastday-WNday1+1);
        DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(ExptCount).HitRate_MUSC=nan(1, WNlastday-WNday1+1);

        DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(ExptCount).ffSTD_PBS=nan(1, WNlastday-WNday1+1);
        DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(ExptCount).ffSTD_MUSC=nan(1, WNlastday-WNday1+1);

        
        DATSTRUCT_FirstWNDayIndEqualsOne.secondtarg(ExptCount).FFminusBase_Mean_PBS=nan(1, WNlastday-WNday1+1);
        DATSTRUCT_FirstWNDayIndEqualsOne.secondtarg(ExptCount).FFminusBase_Mean_MUSC=nan(1, WNlastday-WNday1+1);

        DATSTRUCT_FirstWNDayIndEqualsOne.meanOfSameType(ExptCount).FFminusBase_Mean_PBS=nan(1, WNlastday-WNday1+1);
        DATSTRUCT_FirstWNDayIndEqualsOne.meanOfSameType(ExptCount).FFminusBase_Mean_MUSC=nan(1, WNlastday-WNday1+1);
       
        
        
        % ====== PBS - COLLECT ALL DAYS THAT ARE POSSIBLE
        AllWNNDays_ActualInds=WNday1:WNlastday;
        for j=1:length(AllWNNDays_ActualInds);
            day=AllWNNDays_ActualInds(j);
            postWNDay=day-WNday1+1;
            
            if length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase_WithinTimeWindow)<day
                continue
            end
            
            
            if isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase_WithinTimeWindow(day))
                continue;
            end
            
            % --- TARG
            ffmean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase_WithinTimeWindow(day);
            
            DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(ExptCount).FFminusBase_Mean_PBS(postWNDay)=ffmean;
            
            hits=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).HitRateStatus_WithinTimeWindow{day});
            totalrends=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).HitRateStatus_WithinTimeWindow{day});
            HitRate=hits/totalrends;
            
            DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(ExptCount).HitRate_PBS(postWNDay)=HitRate;
            
            % Daily S.D.
            ffstd=std(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals_DevFromBase_WithinTimeWindow{day});
            DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(ExptCount).ffSTD_PBS(postWNDay)=ffstd;

            
            
            % ---- SECOND TARG
            if ~isempty(othertarg)
            ffmean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(othertarg).meanFF_DevFromBase_WithinTimeWindow(day);
            
            DATSTRUCT_FirstWNDayIndEqualsOne.secondtarg(ExptCount).FFminusBase_Mean_PBS(postWNDay)=ffmean;
            end
            
            
            % ----- MEAN OF SAME TYPES (take mean of means)
            if ~isempty(SameTypeSyls)
                ffmeanall = [];
                for k=1:length(SameTypeSyls)
                   syltmp = SameTypeSyls{k};
                   ffmean = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syltmp).meanFF_DevFromBase_WithinTimeWindow(day);
                   
                   ffmeanall = [ffmeanall ffmean];
                end
                
                
                DATSTRUCT_FirstWNDayIndEqualsOne.meanOfSameType(ExptCount).FFminusBase_Mean_PBS(postWNDay)=...
                    mean(ffmeanall);
            end
            
        end
        
        
        
        % ===== MUSC - COLLECT ONLY DAYS WITH MUSC INACTIVATION
        inds=days_with_musc_RelToWNDay1>0 & days_with_musc_RelToWNDay1 <= (WNlastday-WNday1+1); % these are days during bidir
        
        ActualDayInds=days_with_musc(inds);
        PostWNInds=days_with_musc_RelToWNDay1(inds);
        
        for j=1:length(ActualDayInds)
            day=ActualDayInds(j);
            postWNDay=PostWNInds(j);
            
            % ---- collect data for this day -TARG
            ffmean_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).meanFF_DevFromBase_WithinTimeWindow(day);
            DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(ExptCount).FFminusBase_Mean_MUSC(postWNDay)=ffmean_MUSC;
            
                        % hit rate
            hits=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).HitRateStatus_WithinTimeWindow{day});
            totalrends=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).HitRateStatus_WithinTimeWindow{day});
            HitRate=hits/totalrends;
            
            DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(ExptCount).HitRate_MUSC(postWNDay)=HitRate;

            
            % Daily S.D.
            ffstd=std(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).FFvals_DevFromBase_WithinTimeWindow{day});
            DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(ExptCount).ffSTD_MUSC(postWNDay)=ffstd;

            
                        % ---- SECOND TARG
            if ~isempty(othertarg)
            ffmean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(othertarg).meanFF_DevFromBase_WithinTimeWindow(day);
            
            DATSTRUCT_FirstWNDayIndEqualsOne.secondtarg(ExptCount).FFminusBase_Mean_MUSC(postWNDay)=ffmean;
            end
            
            
            % ----- MEAN OF SAME TYPES (take mean of means)
            if ~isempty(SameTypeSyls)
                ffmeanall = [];
                for k=1:length(SameTypeSyls)
                   syltmp = SameTypeSyls{k};
                   ffmean = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syltmp).meanFF_DevFromBase_WithinTimeWindow(day);
                   
                   ffmeanall = [ffmeanall ffmean];
                end
                
                
                DATSTRUCT_FirstWNDayIndEqualsOne.meanOfSameType(ExptCount).FFminusBase_Mean_MUSC(postWNDay)=...
                    mean(ffmeanall);
            end

        end
        
        
        % ==== collect dates to put lines
        ImportantInds=[];
        try
            multidirStart=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind+1;
            multidirEnd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_LastDay_Ind;
            ImportantInds=[ImportantInds multidirStart multidirEnd];
        catch err
        end
        try
            samedirStart= SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.SameDir_OneDayBeforeStart_Ind+1;
            samedirEnd= SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.SameDir_LastDay_Ind;
            ImportantInds=[ImportantInds samedirStart samedirEnd];
        end
        ImportantInds=ImportantInds-WNday1+1;
        
        % ====== SAVE INFORMATION ABOUT THIS EXPERIMENT
        DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(ExptCount).birdname=birdname;
        DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(ExptCount).exptname=exptname;
        DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(ExptCount).targ_learn_dir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(ExptCount).ImportantInds_relWNday1=ImportantInds;
        DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(ExptCount).SameTypeSyls=SameTypeSyls;
        
        ExptCount=ExptCount+1;
        
    end
end


%% =================== RECALCALATE MAINTAINED SHIFT WINDOW?
% if so, then does following.
% uses a 6-day window, slides it starting from day 3. one gets window where
% mean of each day is within all-day mean +/- 0.75 SD(mean sd of renditions
% across days), then stops. then keeps adding days to end until breaks that
% rule.
% if 6 day doesn't work, then uses 5 day, and so on, smallest window to try
% is 4
targSylTypeList={'firsttarget'};

DayWindowSizes=[6 5 4];

for lll=1:length(targSylTypeList)
    targsylType=targSylTypeList{lll};

for i=1:length(DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget);
    ChosenConsolDays=[];
    
    for k=1:length(DayWindowSizes);
        
        numDaysInWind=DayWindowSizes(k);
        
        for kk=3:1000;
            % starting from day 3, keep sliding window
            
            DaysWind=kk:kk+numDaysInWind-1;
            
            % stop if the last day in window is past the last day with data
            if DaysWind(end)>length(DATSTRUCT_FirstWNDayIndEqualsOne.(targsylType)(i).FFminusBase_Mean_PBS);
                break % break out of loop. this window size did not get anything useful. try next smallest window size
            end
            
            % CHECK THIS WINDOW
            FFvals=DATSTRUCT_FirstWNDayIndEqualsOne.(targsylType)(i).FFminusBase_Mean_PBS(DaysWind);
            meanFF_overall=mean(FFvals);
            meanSTD=mean(DATSTRUCT_FirstWNDayIndEqualsOne.(targsylType)(i).ffSTD_PBS(DaysWind));
            numZ=0.75;
            
            FFlims=[meanFF_overall-numZ*meanSTD meanFF_overall+numZ*meanSTD];
            
            % if all ffvals are within the limits, then keep this starting day.
            % keep adding days until fails the requirement of being within
            % window.
            
            func=@(x)(x>FFlims(1) && x<FFlims(2));
            withinLims=arrayfun(func, FFvals);
            
            
            while all(withinLims)==1
                ChosenConsolDays=[DaysWind(1) DaysWind(end)];
                DaysWind=[DaysWind DaysWind(end)+1];
                
                if DaysWind(end)>length(DATSTRUCT_FirstWNDayIndEqualsOne.(targsylType)(i).FFminusBase_Mean_PBS);
                    break % break out of loop. DONE
                end
                
                % CHECK THIS WINDOW
                FFvals=DATSTRUCT_FirstWNDayIndEqualsOne.(targsylType)(i).FFminusBase_Mean_PBS(DaysWind);
                meanFF_overall=mean(FFvals);
                meanSTD=mean(DATSTRUCT_FirstWNDayIndEqualsOne.(targsylType)(i).ffSTD_PBS(DaysWind));
                numZ=0.75;
                
                FFlims=[meanFF_overall-numZ*meanSTD meanFF_overall+numZ*meanSTD];
                
                % if all ffvals are within the limits, then keep this starting day.
                % keep adding days until fails the requirement of being within
                % window.
                
                func=@(x)(x>FFlims(1) && x<FFlims(2));
                withinLims=arrayfun(func, FFvals);
            end
            
            if ~isempty(ChosenConsolDays)
                break
            end
            
            
        end
        
        if ~isempty(ChosenConsolDays)
            break
        end
        
    end
    
    Params.LMANTimeCourse.ConsolPer_IndsFromStart(i).(targsylType)=ChosenConsolDays;
    
end
end


%% ==== replace the consol periods that were hand entered
if GetAutoConsolWindows==1;
    
for i=1:length(DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget);

    birdname=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(i).birdname;
    exptname=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(i).exptname;
    
    % ====
    
    birdInd=find(strcmp(Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart, birdname));
    exptInd=find(strcmp(Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart, exptname));
    
    ind=intersect(birdInd+1, exptInd);
    
    if ~isempty(ind);
        consolPeriod_OLD=Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart{ind+1};
        
        consolPeriod_NEW_ft=Params.LMANTimeCourse.ConsolPer_IndsFromStart(i).firsttarget;

        if isempty(consolPeriod_NEW_ft)
        disp([birdname '-' exptname ': auto consol period EMPTY!!' ])
            continue
        end
        
        consolPeriod_NEW_olap=consolPeriod_NEW_ft;
        disp([birdname '-' exptname ': ' num2str(consolPeriod_OLD) ' (OLD) ===== ' num2str(consolPeriod_NEW_olap) ' (NEW, OLAP!!)'])
        
        % REPLACE THOSE VALUES
        Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart{ind+1}=consolPeriod_NEW_olap;
        
    end

end

end

%% ====== PLOT [EACH EXPT ONE FIG --> FF, HIT RATE, CONSOL]

for i=1:length(DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget);
    
    birdname=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(i).birdname;
    exptname=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(i).exptname;
    
    lt_figure; hold on;
    
    % +++++ 1) FF of each targ
    lt_subplot(3,1,1); hold on;
    title([birdname '-' exptname]);
    xlabel('day (day1 = bidir day 1)');
    ylabel('shift from baseline (hz)');
    
    % == first targ
    ffPBS=DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(i).FFminusBase_Mean_PBS;
    ffMUSC=DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(i).FFminusBase_Mean_MUSC;
    inds=find(~isnan(ffPBS));
    
    lt_plot(inds, ffPBS(inds), {'Color', 'b', 'LineStyle', '-'});
    lt_plot(inds, ffMUSC(inds), {'Color', 'r', 'LineStyle', '-'});
    
    
    lt_plot_zeroline;
    
    % line for consol start and end
    birdInd=find(strcmp(Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart, birdname));
    exptInd=find(strcmp(Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart, exptname));
    
    ind=intersect(birdInd+1, exptInd);
    
    if ~isempty(ind);
        consolPeriod=Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart{ind+1};
        
        line([consolPeriod(1)-0.4 consolPeriod(1)-0.4], ylim, 'Color','m');
        line([consolPeriod(2)+0.4 consolPeriod(2)+0.4], ylim, 'Color','m');
    end
    
    % +++++++ 2) CONSOLIDATION
    lt_subplot(3,1,2); hold on;
    ylabel('consolidation');
    
    % == extract separation
    % -- PBS
    ffPBS=DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(i).FFminusBase_Mean_PBS;
    
    
    % -- MUSC
    ffMUSC=DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(i).FFminusBase_Mean_MUSC;
    
    
    % - plot
    inds=find(~isnan(ffMUSC));
    
    Consolidation = ffMUSC(inds)./ffPBS(inds);
    
    lt_plot(inds, Consolidation, {'Color', 'k', 'LineStyle', '-'});
    
    line(xlim, [1 1], 'LineStyle','--')
    
    

    xlim([0 50]); ylim([-0.2 1.2]); lt_plot_zeroline;
    
    
    
    % +++++++++++++++++ 3) HIT RATE (First target);
    lt_subplot(3,1,3); hold on;
    title('first targ');
    ylabel('hit rate');
    
    hrPBS=DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(i).HitRate_PBS;
    hrMUSC=DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(i).HitRate_MUSC;
    
    inds=find(~isnan(hrPBS));
    lt_plot(inds, hrPBS(inds), {'Color', 'b', 'LineStyle', '-'});
    lt_plot(inds, hrMUSC(inds), {'Color', 'r'});
    
    lt_plot_zeroline;
    line(xlim, [1 1]);
    
    

    
    
end


%% ====== PLOT

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

for i=1:length(DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget);
    
    birdname=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(i).birdname;
    exptname=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(i).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([birdname '-' exptname]);
    xlabel('day (day1 = WN day 1)');
    ylabel('shift from baseline (hz)');
    
    
    % == first targ
    ffPBS=DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(i).FFminusBase_Mean_PBS;
    ffMUSC=DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(i).FFminusBase_Mean_MUSC;
    inds=find(~isnan(ffPBS));
    
    lt_plot(inds, ffPBS(inds), {'Color', 'b', 'LineStyle', '-'});
    lt_plot(inds, ffMUSC(inds), {'Color', 'r', 'LineStyle', '-'});
    
    
    lt_plot_zeroline;
    
    % line for important inds
    importantInds=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(i).ImportantInds_relWNday1;
    if ~isempty(importantInds);
        for j=1:length(importantInds);
            if mod(j,2)==1; % then if odd, start of epoch
                            line([importantInds(j)-0.5 importantInds(j)-0.5], ylim, 'Color','k');

            else
                %then is even, so is end of epoch
            line([importantInds(j)+0.5 importantInds(j)+0.5], ylim, 'Color','k');
            end
        end
    end
    
    
    % line for consol start and end
    birdInd=find(strcmp(Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart, birdname));
    exptInd=find(strcmp(Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart, exptname));
    
    ind=intersect(birdInd+1, exptInd);
    
    if ~isempty(ind);
        consolPeriod=Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart{ind+1};
        
        line([consolPeriod(1)-0.4 consolPeriod(1)-0.4], ylim, 'Color','m');
        line([consolPeriod(2)+0.4 consolPeriod(2)+0.4], ylim, 'Color','m');
        
                % ================= PUT MEAN PICH +/- 0.75 STD DURING MAINTAINED
        % SHIFT PERIOD
        numZ=0.75;
        
        % - first targ
        meanFF=mean(DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(i).FFminusBase_Mean_PBS(consolPeriod(1):consolPeriod(2)));
        meanStdFF=mean(DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(i).ffSTD_PBS(consolPeriod(1):consolPeriod(2)));
        
        Xlim=xlim;
               
        shadedErrorBar(Xlim, [meanFF meanFF], [numZ*meanStdFF numZ*meanStdFF],{'Color','b'},1);
        

    end
    
    
    
end

lt_subtitle('red=WN end;, magen=consol period');




%% === extract consolidation days


for j=1:length(DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget); % num expts
    birdname=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).birdname;
    exptname=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).exptname;
    
    
    birdInd=find(strcmp(Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart, birdname));
    exptInd=find(strcmp(Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart, exptname));
    
    ind=intersect(birdInd+1, exptInd);
    
    if ~isempty(ind);
        consolPeriod=Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart{ind+1};
        
        DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).consolPeriod=consolPeriod;
        
    else
        disp(['WNARNING - no consol period entered for ' birdname '-' exptname ' - using all days']);
        
    end
end


%% ==== PLOT MEAN ACROSS ALL EXPERIMENTS
% ==== DON'T CARE ABOUT CONSOLIDATION - SIMPLY BIN IN DAYS RELATIVE TO
% LEARNING START

% === to add:
% 1) stop when reach end of bidir
% 3) only run if bird has data for all bins

if BinRelConsolDay1==0;
    
    
    maxday=length([DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget.FFminusBase_Mean_PBS]);
    DayBins_StartEdges=1:NumDaysInBin:maxday;
    
    OUTPUT=struct; % one cell for each bin
    Exptcount=1;
    
    for i=1:length(DayBins_StartEdges);
        
        OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC=[];
        OUTPUT.firsttarget.FFrelBaseInBin(i).PBS=[];
        OUTPUT.firsttarget.FFrelBaseInBin(i).MPbias=[];
%         
        OUTPUT.firsttarget.otherstats(i).hitratePBS=[];
        OUTPUT.firsttarget.otherstats(i).hitrateMUSC=[];

        
        OUTPUT.INFORMATION(i).experimentNum=[];
        
        OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC=[];
        OUTPUT.secondtarg.FFrelBaseInBin(i).PBS=[];

        OUTPUT.meanOfSameType.FFrelBaseInBin(i).MUSC=[];
        OUTPUT.meanOfSameType.FFrelBaseInBin(i).PBS=[];
        
        DaysInBin=DayBins_StartEdges(i):DayBins_StartEdges(i)+NumDaysInBin-1;
        
        for day=DaysInBin;
            
            
            
            
            % == for this day, collect all MUSC. for days with musc, also
            % collect PBS
            
            for j=1:length(DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget); % num expts
                
                % ++++++++++++++++++++++++++++++++++++++
                % === only continue if vector is long enough to today
                if length(DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(j).FFminusBase_Mean_MUSC)<day;
                    continue;
                end
                
                % === only continue if first target has MUSC data for today
                if isnan(DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(j).FFminusBase_Mean_MUSC(day));
                    continue;
                end
                
                % ==== make sure this day is in consolidation period (if
                % desired)
                if OnlyConsolPeriod_bidir==1;
                    birdname=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).birdname;
                    exptname=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).exptname;
                    
                    
                    birdInd=find(strcmp(Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart, birdname));
                    exptInd=find(strcmp(Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart, exptname));
                    
                    ind=intersect(birdInd+1, exptInd);
                    
                    if ~isempty(ind);
                        consolPeriod=Params.LMANTimeCourse.SamedirConsolPeriod_IndsFromSamedirStart{ind+1};
                        
                        if day<consolPeriod(1) | day>consolPeriod(2)
                            continue
                        end
                        
                    else
                        disp(['WNARNING - no consol period entered for ' birdname '-' exptname ' - using all days']);
                        
                    end
                end
                
                
                
                % ++++++++++++++++++++++++++
                
                
                
                % ============== first target
                ffMUSC=DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(j).FFminusBase_Mean_MUSC(day); % MUSC
                ffPBS=DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(j).FFminusBase_Mean_PBS(day); % PBS
                
                % --- flip sign if learning was neg
                ffMUSC=ffMUSC*DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                ffPBS=ffPBS*DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                
                % ---- REVERSION (MP BIAS)
                MPbias=ffMUSC/ffPBS;

                % --- hr
                hrPBS=DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(j).HitRate_PBS(day);
                hrMUSC=DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(j).HitRate_MUSC(day);

                % ---- OUTPUT
                OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC=[OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC ffMUSC];
                OUTPUT.firsttarget.FFrelBaseInBin(i).PBS=[OUTPUT.firsttarget.FFrelBaseInBin(i).PBS ffPBS];
                
                OUTPUT.firsttarget.FFrelBaseInBin(i).MPbias=[OUTPUT.firsttarget.FFrelBaseInBin(i).MPbias MPbias];
                
                OUTPUT.firsttarget.otherstats(i).hitratePBS=[OUTPUT.firsttarget.otherstats(i).hitratePBS hrPBS];
                OUTPUT.firsttarget.otherstats(i).hitrateMUSC=[OUTPUT.firsttarget.otherstats(i).hitrateMUSC hrMUSC];


                 % =============== second targ
                ffMUSC=DATSTRUCT_FirstWNDayIndEqualsOne.secondtarg(j).FFminusBase_Mean_MUSC(day); % MUSC
                ffPBS=DATSTRUCT_FirstWNDayIndEqualsOne.secondtarg(j).FFminusBase_Mean_PBS(day); % PBS
                
                % --- flip sign if learning was neg
                ffMUSC=ffMUSC*DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                ffPBS=ffPBS*DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                                                  
                OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC=[OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC ffMUSC];
                OUTPUT.secondtarg.FFrelBaseInBin(i).PBS=[OUTPUT.secondtarg.FFrelBaseInBin(i).PBS ffPBS];
               
                
                % ======================== Mean of same type
                ffMUSC=DATSTRUCT_FirstWNDayIndEqualsOne.meanOfSameType(j).FFminusBase_Mean_MUSC(day); % MUSC
                ffPBS=DATSTRUCT_FirstWNDayIndEqualsOne.meanOfSameType(j).FFminusBase_Mean_PBS(day); % PBS
                
                % --- flip sign if learning was neg
                ffMUSC=ffMUSC*DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                ffPBS=ffPBS*DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                                                  
                OUTPUT.meanOfSameType.FFrelBaseInBin(i).MUSC=[OUTPUT.meanOfSameType.FFrelBaseInBin(i).MUSC ffMUSC];
                OUTPUT.meanOfSameType.FFrelBaseInBin(i).PBS=[OUTPUT.meanOfSameType.FFrelBaseInBin(i).PBS ffPBS];
                
                
                % ===== INFORMATION
                
                OUTPUT.INFORMATION(i).experimentNum=[OUTPUT.INFORMATION(i).experimentNum j];
                
                Exptcount=Exptcount+1;
            end
            
        end
    end
    
elseif BinRelConsolDay1==1;
    %% ==== PLOT MEAN ACROSS ALL EXPERIMENTS FOR BIDIR  [consol day 1 = bin 1]
    
    % === to add:
    % 1) stop when reach end of bidir
    % 3) only run if bird has data for all bins
    
    
    maxday=length([DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget.FFminusBase_Mean_PBS]);
    DayBins_StartEdges=1:NumDaysInBin:maxday;
    
    OUTPUT=struct; % one cell for each bin
    Exptcount=1;
    
    for i=1:length(DayBins_StartEdges);
        
        OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC=[];
        OUTPUT.firsttarget.FFrelBaseInBin(i).PBS=[];
        OUTPUT.firsttarget.FFrelBaseInBin(i).MPbias=[];
        
        OUTPUT.firsttarget.otherstats(i).hitratePBS=[];
        OUTPUT.firsttarget.otherstats(i).hitrateMUSC=[];
        
        
        OUTPUT.INFORMATION(i).experimentNum=[];
        
        OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC=[];
        OUTPUT.secondtarg.FFrelBaseInBin(i).PBS=[];

        OUTPUT.meanOfSameType.FFrelBaseInBin(i).MUSC=[];
        OUTPUT.meanOfSameType.FFrelBaseInBin(i).PBS=[];
        
        DaysInBin=DayBins_StartEdges(i):DayBins_StartEdges(i)+NumDaysInBin-1;
        
        for day=DaysInBin;
            
            
            
            
            % == for this day, collect all MUSC. for days with musc, also
            % collect PBS
            
            for j=1:length(DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget); % num expts
                
                birdname=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).birdname;
                exptname=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).exptname;
                
                
                % === CONVERT TO DAY RELATIVE TO START OF CONSOLIDATION (CONSOL DAY 1 =
                % IND 1)
                if isempty(DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).consolPeriod)
                    
                    disp(['[skipping] No consol period for ' birdname '-' exptname])
                    continue
                end
                
                
                consolDay1_relBidirDay1=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).consolPeriod(1);
                day_rel_consol=day-consolDay1_relBidirDay1+1;
                
                % which day bin would this be in?
                daybin=ceil(day_rel_consol/NumDaysInBin);
                
                
                % ++++++++++++++++++++++++++++++++++++++ FILTERS
                % == make sure day rel consol is not <1
                if day_rel_consol<1
                    continue
                end
                
                % ==== make sure this day is in consolidation period (if
                % desired)
                if OnlyConsolPeriod_bidir==1;
                    consolPeriod=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).consolPeriod;
                    
                    if day>consolPeriod(2)
                        continue
                    end
                end
                
                % === only continue if vector is long enough to today
                if length(DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(j).FFminusBase_Mean_MUSC)<day;
                    continue;
                end
                
                % === only continue if first target has MUSC data for today
                if isnan(DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(j).FFminusBase_Mean_MUSC(day));
                    continue;
                end
                
                
                
                % ++++++++++++++++++++++++++
                % ============== first target
                ffMUSC=DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(j).FFminusBase_Mean_MUSC(day); % MUSC
                ffPBS=DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(j).FFminusBase_Mean_PBS(day); % PBS
                
                % --- flip sign if learning was neg
                ffMUSC=ffMUSC*DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                ffPBS=ffPBS*DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                
                % ---- REVERSION (MP BIAS)
                MPbias=ffMUSC/ffPBS;
                
                hrPBS=DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(j).HitRate_PBS(day);
                hrMUSC=DATSTRUCT_FirstWNDayIndEqualsOne.firsttarget(j).HitRate_MUSC(day);
                
                % ---- OUTPUT
                OUTPUT.firsttarget.FFrelBaseInBin(daybin).MUSC=[OUTPUT.firsttarget.FFrelBaseInBin(daybin).MUSC ffMUSC];
                OUTPUT.firsttarget.FFrelBaseInBin(daybin).PBS=[OUTPUT.firsttarget.FFrelBaseInBin(daybin).PBS ffPBS];
                
                OUTPUT.firsttarget.FFrelBaseInBin(daybin).MPbias=[OUTPUT.firsttarget.FFrelBaseInBin(daybin).MPbias MPbias];
                                
                OUTPUT.firsttarget.otherstats(daybin).hitratePBS=[OUTPUT.firsttarget.otherstats(daybin).hitratePBS hrPBS];
                OUTPUT.firsttarget.otherstats(daybin).hitrateMUSC=[OUTPUT.firsttarget.otherstats(daybin).hitrateMUSC hrMUSC];
                
                
                % =============== second targ
                ffMUSC=DATSTRUCT_FirstWNDayIndEqualsOne.secondtarg(j).FFminusBase_Mean_MUSC(day); % MUSC
                ffPBS=DATSTRUCT_FirstWNDayIndEqualsOne.secondtarg(j).FFminusBase_Mean_PBS(day); % PBS
                
                % --- flip sign if learning was neg
                ffMUSC=ffMUSC*DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                ffPBS=ffPBS*DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                                                  
                OUTPUT.secondtarg.FFrelBaseInBin(daybin).MUSC=[OUTPUT.secondtarg.FFrelBaseInBin(daybin).MUSC ffMUSC];
                OUTPUT.secondtarg.FFrelBaseInBin(daybin).PBS=[OUTPUT.secondtarg.FFrelBaseInBin(daybin).PBS ffPBS];
               
                
                
                
                % =============== MEAN OF SAME TYPE
                ffMUSC=DATSTRUCT_FirstWNDayIndEqualsOne.meanOfSameType(j).FFminusBase_Mean_MUSC(day); % MUSC
                ffPBS=DATSTRUCT_FirstWNDayIndEqualsOne.meanOfSameType(j).FFminusBase_Mean_PBS(day); % PBS
                
                % --- flip sign if learning was neg
                ffMUSC=ffMUSC*DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                ffPBS=ffPBS*DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                                                  
                OUTPUT.meanOfSameType.FFrelBaseInBin(daybin).MUSC=[OUTPUT.meanOfSameType.FFrelBaseInBin(daybin).MUSC ffMUSC];
                OUTPUT.meanOfSameType.FFrelBaseInBin(daybin).PBS=[OUTPUT.meanOfSameType.FFrelBaseInBin(daybin).PBS ffPBS];
                
                
                
                
                % ===== INFORMATION
                
                OUTPUT.INFORMATION(daybin).experimentNum=[OUTPUT.INFORMATION(daybin).experimentNum j];
                
                
                Exptcount=Exptcount+1;
            end
            
        end
    end
end


%% ===== PLOT (DAY BINS REL BIDIR START DAY 1)
%
% lt_figure; hold on;
% title('mean of MP biases (within each extp) (blue=second targ), (bk=first targ)');
% xlabel([ 'day bin (rel bidir day 1) (binsize: ' num2str(NumDaysInBin)]);
% ylabel('MP bias (mean of value for each expt');
%
% X=1:length(DayBins_StartEdges);
% Yfirsttarg=[];
% Ysecondtarg=[];
% Yfirsttarg_sem=[];
% Ysecondtarg_sem=[];
%
%
% for i=1:length(X);
%
%     if isempty(OUTPUT.firsttarget.FFrelBaseInBin(i).MPbias)
%     continue
%     end
%       %=== first targ
%       % -- plot each point
%       y=OUTPUT.firsttarget.FFrelBaseInBin(i).MPbias;
%
%       plot(i, y, 'ok');
%
%       % -- collect mean
%         Yfirsttarg(i)=mean(y);
%         Yfirsttarg_sem(i)=lt_sem(y);
%
%         lt_plot_bar(i, mean(y), {'Errors', lt_sem(y), 'Color','k'});
%
%         % === second targ
%        % -- plot each point
%       y=OUTPUT.secondtarg.FFrelBaseInBin(i).MPbias;
%
%       plot(i, y, 'ob');
%
%       % -- collect mean
%         Ysecondtarg(i)=mean(y);
%         Ysecondtarg_sem(i)=lt_sem(y);
%
%                 lt_plot_bar(i, mean(y), {'Errors', lt_sem(y), 'Color','b'});
%
% end



%% === if an experiment has >1 sample in a bin, take average of those samples

NumBins=length(OUTPUT.INFORMATION);

OUTPUT_WithinExptAvg=struct;

for i=1:NumBins
    
    % ==== extract data vectors
%     MPbias_ft=OUTPUT.firsttarget.FFrelBaseInBin(i).MPbias;
    MUSC_ft=OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC;
    PBS_ft=OUTPUT.firsttarget.FFrelBaseInBin(i).PBS;
    
    hr_first_PBS=OUTPUT.firsttarget.otherstats(i).hitratePBS;
    hr_first_MUSC=OUTPUT.firsttarget.otherstats(i).hitrateMUSC;

    exptNums=OUTPUT.INFORMATION(i).experimentNum;
    
        
    MUSC_st=OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC;
    PBS_st=OUTPUT.secondtarg.FFrelBaseInBin(i).PBS;
    
    MUSC_meansametype = OUTPUT.meanOfSameType.FFrelBaseInBin(i).MUSC;
    PBS_meansametype = OUTPUT.meanOfSameType.FFrelBaseInBin(i).PBS;
    
    
    % === prepare output
    Exptnum_out=[];
    MUSC_ft_out=[];
    PBS_ft_out=[];
    
    MUSC_st_out=[];
    PBS_st_out=[];
    
    MUSC_meansametype_out = [];
    PBS_meansametype_out = [];
    
    hr_first_PBS_out=[];
    hr_first_MUSC_out=[];

    % go thru each exptnum, if there are multiple cases, then avearge them
    MaxExptNum=max(exptNums);
    
    for k=1:MaxExptNum;
        
        inds=find(exptNums==k);
        
        if ~isempty(inds)
                % average them
                musc_ft=mean(MUSC_ft(inds));
                pbs_ft=mean(PBS_ft(inds));
                musc_st=mean(MUSC_st(inds));
                pbs_st=mean(PBS_st(inds));
                musc_meanst=mean(MUSC_meansametype(inds));
                pbs_meanst=mean(PBS_meansametype(inds));
                
                
            hr_ft_pbs=mean(hr_first_PBS(inds));
            hr_ft_musc=mean(hr_first_MUSC(inds));
            
            % ===== add to new vectors
            %             MPbias_WithinExptAveraged=[MPbias_WithinExptAveraged mpbias];
            Exptnum_out=[Exptnum_out k];
            MUSC_ft_out=[MUSC_ft_out musc_ft];
            PBS_ft_out=[PBS_ft_out pbs_ft];
            MUSC_st_out=[MUSC_st_out musc_st];
            PBS_st_out=[PBS_st_out pbs_st];
            MUSC_meansametype_out = [MUSC_meansametype_out musc_meanst];
            PBS_meansametype_out = [PBS_meansametype_out pbs_meanst];
            
            hr_first_PBS_out=[hr_first_PBS_out hr_ft_pbs];
            hr_first_MUSC_out=[hr_first_MUSC_out hr_ft_musc];
            
            
        end
        
    end
    
    
    % ====== REPLACE OLD VECTORS
    %     OUTPUT_WithinExptAvg.firsttarget.FFrelBaseInBin(i).MPbias=MPbias_WithinExptAveraged;
    OUTPUT_WithinExptAvg.INFORMATION(i).experimentNum=Exptnum_out;
    
    OUTPUT_WithinExptAvg.firsttarget.FFrelBaseInBin(i).MUSC=MUSC_ft_out;
    OUTPUT_WithinExptAvg.firsttarget.FFrelBaseInBin(i).PBS=PBS_ft_out;
    
    OUTPUT_WithinExptAvg.firsttarget.otherstats(i).hitratePBS=hr_first_PBS_out;
    OUTPUT_WithinExptAvg.firsttarget.otherstats(i).hitrateMUSC=hr_first_MUSC_out;

    OUTPUT_WithinExptAvg.secondtarg.FFrelBaseInBin(i).MUSC=MUSC_st_out;
    OUTPUT_WithinExptAvg.secondtarg.FFrelBaseInBin(i).PBS=PBS_st_out;

    OUTPUT_WithinExptAvg.meanOfSameType.FFrelBaseInBin(i).MUSC=MUSC_meansametype_out;
    OUTPUT_WithinExptAvg.meanOfSameType.FFrelBaseInBin(i).PBS=PBS_meansametype_out;
end




%% ==== ONLY KEEP ONE VALUE PER EXPERIMENT PER BIN

if TakeAverageWithinExpt==1
    OUTPUT=OUTPUT_WithinExptAvg;   
end
    clear OUTPUT_WithinExptAvg;   
    
    



%% PLOT (RAW PBS AND MUSC, for targ and new targ);
lt_figure; hold on;

lt_subplot(2,1,1); hold on;
title('Raw PBS and MUSC, from bidir day 1, all experiments');
xlabel([ 'day bin  (binsize: ' num2str(NumDaysInBin) ')']);
ylabel('FF rel baseline (in targ dir)');

X=1:length(DayBins_StartEdges);
subdivisions=0.1; % for plotting diff things in one ind

for i=1:length(X);
    
    if isempty(OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC)
        continue
    end
    
    
    %=== first targ
    x=[i-subdivisions*2 i-subdivisions];
    % -- plot each point
    yPBS=OUTPUT.firsttarget.FFrelBaseInBin(i).PBS;
    yMUSC=OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC;
    
    plot(x, [yPBS' yMUSC']', '-ok');
    
    % -- means
    if size([yPBS' yMUSC'],1)>1
        %       lt_plot_bar(x-0.05, mean([yPBS' yMUSC']), {'Color', 'k', 'Errors', lt_sem([yPBS' yMUSC'])});
        lt_plot_bar(x-subdivisions/4, mean([yPBS' yMUSC']), {'Color', 'k'});
    end
    
    % stats
%     [~, p]=ttest(yPBS, yMUSC);
    p=signrank(yPBS, yMUSC);
    if p<0.1;
        lt_plot_text(x(1), max([yPBS yMUSC]), num2str(p, '%3.2g'));
    end
    
%     % === second targ
%     x=[i+subdivisions i+2*subdivisions];
%     
%     % -- plot each point
%     yPBS=OUTPUT.secondtarg.FFrelBaseInBin(i).PBS;
%     yMUSC=OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC;
%     
%     plot(x, [yPBS' yMUSC']', '-ob');
%     
%     % -- means
%     if size([yPBS' yMUSC'],1)>1
%         %       lt_plot_bar(x+0.05, mean([yPBS' yMUSC']), {'Color', 'b', 'Errors', lt_sem([yPBS' yMUSC'])});
%         lt_plot_bar(x+subdivisions/4, mean([yPBS' yMUSC']), {'Color', 'b'});
%     end
%     
%     % stats
%     [~, p]=ttest(yPBS, yMUSC);
%     if p<0.1;
%         lt_plot_text(x(1), max([yPBS yMUSC]), num2str(p, '%3.2g'));
%     end
    
    
    
end


% ========= [ only mean, each mean norm to PBS]
lt_subplot(2,1,2); hold on;
title('Raw PBS and MUSC, from bidir day 1, all experiments [mean, norm so PBS is 1]');
xlabel([ 'day bin (binsize: ' num2str(NumDaysInBin) ')']);
ylabel('FF rel baseline (in targ dir)');

for i=1:length(X);
    
    if isempty(OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC)
        continue
    end
    
    %=== first targ
    x=[i-subdivisions*2 i-subdivisions];
    % -- plot each point
    yPBS=OUTPUT.firsttarget.FFrelBaseInBin(i).PBS;
    yMUSC=OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC;
    
    
    % -- means
    if size([yPBS' yMUSC'],1)>1
        lt_plot_bar(x-0.05, mean([yPBS' yMUSC'])./mean([yPBS' yPBS']), {'Color', 'k', 'Errors', lt_sem([yPBS' yMUSC'])./mean([yPBS' yPBS'])});
    end
    
    
%     % === second targ
%     x=[i+subdivisions i+2*subdivisions];
%     
%     % -- plot each point
%     yPBS=OUTPUT.secondtarg.FFrelBaseInBin(i).PBS;
%     yMUSC=OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC;
%     
%     % -- means
%     if size([yPBS' yMUSC'],1)>1
%         lt_plot_bar(x-0.05, -mean([yPBS' yMUSC'])./mean([yPBS' yPBS']), {'Color', 'b', 'Errors', lt_sem([yPBS' yMUSC'])./mean([yPBS' yPBS'])});
%     end
    
end

%% PLOT (consolidation score, diff from warren et al.?)
lt_figure; hold on;

title('consolidation');
xlabel([ 'day bin  (binsize: ' num2str(NumDaysInBin) ')']);
ylabel('consolidation');


ConsolidationCell={};

for i=1:length(X);
    
    if isempty(OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC)
        continue
    end
    
    %=== first targ
    x=[i];
    % -- plot each point
    yPBS=OUTPUT.firsttarget.FFrelBaseInBin(i).PBS;
    yMUSC=OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC;
    
    consolidationAll=yMUSC./yPBS;
    
    plot(x+0.1,consolidationAll, 'ok');

    
    % -- means
    if length(consolidationAll)>1
        lt_plot_bar(x, mean(consolidationAll), {'Color', 'k', 'Errors', lt_sem(consolidationAll)});
    end
    
    % --- plot text
    if length(consolidationAll)>1
        N=length(consolidationAll);
        ymean=mean(consolidationAll);
        ysem=lt_sem(consolidationAll);

    lt_plot_text(x, 1.2*mean(consolidationAll), ['n=' num2str(N) '; mean(sem)=' num2str(ymean) '(' num2str(ysem) ')'], 'm');
    end
    
    
    % stats
%     [~, p]=ttest(consolidationAll, 0.849);
    p=signrank(consolidationAll, 0.849);
    if p<0.1;
        lt_plot_text(x(1), max(consolidationAll), num2str(p, '%3.2g'));
    end
    
    ConsolidationCell{i}=consolidationAll;
    
    
end

% ==== COMPARE diff bins using 1) rank sum for bin 3 vs bin 1 and 2) anova
% for diff in means for all bins (either just 1 -3 as in warren or all
% bins)
p=ranksum(ConsolidationCell{3}, ConsolidationCell{1});

consolidationVals=[];
binVals=[];
for i=1:length(ConsolidationCell);
    consolidationVals=[consolidationVals ConsolidationCell{i}];
    binVals=[binVals i*ones(1, length(ConsolidationCell{i}))];
end

consolidationVals=[];
binVals=[];
for i=1:3;
    consolidationVals=[consolidationVals ConsolidationCell{i}];
    binVals=[binVals i*ones(1, length(ConsolidationCell{i}))];
end
anovan(consolidationVals, {binVals})


% line for consolidation from warren et al.,
line(xlim, [0.85 0.85])

%% ==================================== AFP bias at second target predicts consolidation at first target?
lt_figure; hold on;

% ---- 1) 
exptnums=OUTPUT.INFORMATION(BinToAnalyze).experimentNum;

firsttarg_consol=OUTPUT.firsttarget.FFrelBaseInBin(BinToAnalyze).MUSC./OUTPUT.firsttarget.FFrelBaseInBin(BinToAnalyze).PBS;

secondtarg_AFPbias=OUTPUT.secondtarg.FFrelBaseInBin(BinToAnalyze).PBS-OUTPUT.secondtarg.FFrelBaseInBin(BinToAnalyze).MUSC;


inds=~isnan(secondtarg_AFPbias);

lt_subplot(3,2,1); hold on;
xlabel('AFP bias, 2nd targ');
ylabel('consol, 1st targ');

lt_regress(firsttarg_consol(inds), secondtarg_AFPbias(inds), 1, 0, 1, 1, 'k');

for i=1:length(firsttarg_consol);
    if inds(i)==1
        
        birdname=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(exptnums(i)).birdname;
         exptname=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(exptnums(i)).exptname;
       
         lt_plot_text(secondtarg_AFPbias(i), firsttarg_consol(i), [birdname '-' exptname(end-3:end)], 'b');
    end
end

line(xlim, [0 0]);
line([0 0], ylim);
line(xlim, [1 1]);
        


% ---- 1) 
exptnums=OUTPUT.INFORMATION(BinToAnalyze).experimentNum;

firsttarg_consol=OUTPUT.firsttarget.FFrelBaseInBin(BinToAnalyze).MUSC./OUTPUT.firsttarget.FFrelBaseInBin(BinToAnalyze).PBS;

secondtarg_consol=OUTPUT.secondtarg.FFrelBaseInBin(BinToAnalyze).MUSC./OUTPUT.secondtarg.FFrelBaseInBin(BinToAnalyze).PBS;


inds=~isnan(secondtarg_consol);

lt_subplot(3,2,2); hold on;
xlabel('consol, 2nd targ');
ylabel('consol, 1st targ');

lt_regress(firsttarg_consol(inds), secondtarg_consol(inds), 1, 0, 1, 1, 'k');

for i=1:length(firsttarg_consol);
    if inds(i)==1
        
        birdname=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(exptnums(i)).birdname;
         exptname=DATSTRUCT_FirstWNDayIndEqualsOne.INFORMATION(exptnums(i)).exptname;
       
         lt_plot_text(secondtarg_consol(i), firsttarg_consol(i), [birdname '-' exptname(end-3:end)], 'b');
    end
end

line(xlim, [0 0]);
line([0 0], ylim);
line(xlim, [1 1]);
        





