function [OUTPUT, DATSTRUCT_FirstBidirDayIndEqualsOne]=lt_seq_dep_pitch_ACROSSBIRDS_LMANTimeCourse(SeqDepPitch_AcrossBirds, Params, OnlyConsolPeriod_bidir, BinRelConsolDay1, TakeAverageWithinExpt, NumDaysInBin, GetAutoConsolWindows)
%% LT 12/22/15

NumBirds=length(SeqDepPitch_AcrossBirds.birds);
%     NumDaysInBin=2;


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
title('Timeline (muscimol in bold) (bidir in blue)');
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
        
        % --- bidir start and end
        bidir_day1=nan;
        bidir_lastday=nan;
        
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart_Ind');
            bidir_day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind+1;
            bidir_day1=bidir_day1-WNday1+1;
            
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_LastDay_Ind');
                bidir_lastday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_LastDay_Ind;
            else
                bidir_lastday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            end
            bidir_lastday=bidir_lastday-WNday1+1;
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
        if ~isnan(bidir_lastday);
            WN_last=bidir_lastday;
        else
            WN_last=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            WN_last=WN_last-WNday1+1;
        end
        
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
        
        
        
        % -- line for start and end of bidir
        line([bidir_day1-0.5 bidir_day1-0.5], [CountExpt-0.4 CountExpt+0.4]);
        line([bidir_lastday+0.5 bidir_lastday+0.5], [CountExpt-0.4 CountExpt+0.4]);
        
        
        % -- fill in musc days
        X=days_with_musc_relWNday1;
        Y=CountExpt*ones(1, 1:length(X));
        
        plot(X, Y, 'ok', 'MarkerFaceColor','k')
        
        % --- line for end of WN
        line([WN_last+0.5 WN_last+0.5], [CountExpt-0.4 CountExpt+0.4], 'Color','r');
        
        
        % -- collect bird name and expt
        BirdnamesExptname_all=[BirdnamesExptname_all [birdname(1:4) '-' exptname(end-4:end)]];
        
        CountExpt=CountExpt+1;
    end
end

set(gca, 'YTick', 1:CountExpt-1);
set(gca, 'YTickLabel', BirdnamesExptname_all);


% ---- WN on line
line([0.5 0.5], ylim, 'Color','r');



%%  for all days plot shift + reversion at target [THINKING OF PUTTING IN CODE ABOVE INSTEAD. SEE ABOVE, STILL IN PROGRESS]


for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        days_with_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
        WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        days_with_musc_relWNday1=days_with_musc-WNday1+1;
        
        % --- bidir start and end
        bidir_day1=nan;
        bidir_lastday=nan;
        
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart_Ind');
            bidir_day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind+1;
            bidir_day1=bidir_day1-WNday1+1;
            
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_LastDay_Ind');
                bidir_lastday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_LastDay_Ind;
            else
                bidir_lastday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            end
            bidir_lastday=bidir_lastday-WNday1+1;
        end
        
        
        % -- last day
        if ~isnan(bidir_lastday);
            WN_last=bidir_lastday;
        else
            WN_last=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            WN_last=WN_last-WNday1+1;
        end
        
        % ---- first day (relative to WN on day)
        firstday=-(WNday1-2);
        
        
        
        
        
        % -- collect bird name and expt
        BirdnamesExptname_all=[BirdnamesExptname_all [birdname(1:4) '-' exptname(end-4:end)]];
        
        CountExpt=CountExpt+1;
    end
end



%% ===== [EXTRACT DATSTRUCT] DURING CONSOLIDATION PERIOD OF 2-DIR, SEE CONSOLIDATION?

disp(' ----- ')
disp(' BIDIR CONSOLIDATION ');

ExptCount=1;


DATSTRUCT_FirstBidirDayIndEqualsOne=struct;

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        % ==== GET ALL DAYS, MAKE THEM RELATIVE TO BIDIR DAY 1 [getting
        % dates]
        days_with_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
        WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart_Ind');
            bidir_day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind+1;
        else
            bidir_day1=nan;
        end
        
        %         if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_start_Inds');
        %             bidir_day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_start_Inds(1);
        %         else
        %             bidir_day1=nan;
        %         end
        
        if isnan(bidir_day1)
            disp([birdname '-' exptname '- SKIPPED bidir consolidation (no bidir data)']);
            continue
        end
        
        days_with_musc_RelToBidirDay1=days_with_musc-bidir_day1+1;
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_LastDay_Ind');
            bidir_last_day=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_LastDay_Ind;
        else
            bidir_last_day=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
        end
        %         bidir_last_day = length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase_WithinTimeWindow); % last day with data, could be after WN end
        
        
        bidir_last_day_actual = bidir_last_day;
        %         if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_end_Inds');
        %         bidir_last_day_actual=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_end_Inds;
        % %         disp([ num2str(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_end_Inds) ' - ' num2str(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd)]);
        %
        %         else
        %         bidir_last_day_actual=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
        %         end
        
        
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.num_multidir_syls~=2;
            disp([birdname '-' exptname '- SKIPPED num multidir syls is not 2']);
            continue
        end
        
        % ========================== DO THIS EXPERIMENT
        % ==== figure out syls
        syl1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls{1};
        syl2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls{2};
        
        if strcmp(syl1, targsyl)
            othersyl=syl2;
        elseif strcmp(syl2,targsyl)
            othersyl=syl1;
        else
            disp('PROBLEM!!!');
            dafascawr3a;
        end
        
        
        % ==== collect data
        DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(ExptCount).FFminusBase_Mean_PBS=nan(1, bidir_last_day-bidir_day1+1);
        DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(ExptCount).FFminusBase_Mean_MUSC=nan(1, bidir_last_day-bidir_day1+1);
        DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(ExptCount).FFminusBase_Mean_PBS=nan(1, bidir_last_day-bidir_day1+1);
        DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(ExptCount).FFminusBase_Mean_MUSC=nan(1, bidir_last_day-bidir_day1+1);
        
        DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(ExptCount).ffSTD_PBS=nan(1, bidir_last_day-bidir_day1+1);
        DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(ExptCount).ffSTD_MUSC=nan(1, bidir_last_day-bidir_day1+1);
        DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(ExptCount).ffSTD_PBS=nan(1, bidir_last_day-bidir_day1+1);
        DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(ExptCount).ffSTD_MUSC=nan(1, bidir_last_day-bidir_day1+1);
        
        
        DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(ExptCount).HitRate_PBS=nan(1, bidir_last_day-bidir_day1+1);
        DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(ExptCount).HitRate_MUSC=nan(1, bidir_last_day-bidir_day1+1);
        DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(ExptCount).HitRate_PBS=nan(1, bidir_last_day-bidir_day1+1);
        DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(ExptCount).HitRate_MUSC=nan(1, bidir_last_day-bidir_day1+1);
        
        
        DATSTRUCT_FirstBidirDayIndEqualsOne.bothTargs(ExptCount).HitRate_Combined_PBS=nan(1, bidir_last_day-bidir_day1+1);
        DATSTRUCT_FirstBidirDayIndEqualsOne.bothTargs(ExptCount).HitRate_Combined_MUSC=nan(1, bidir_last_day-bidir_day1+1);
        
        
        
        % ====== PBS - COLLECT ALL DAYS THAT ARE POSSIBLE
        AllBidirDays_ActualInds=bidir_day1:bidir_last_day;
        for j=1:length(AllBidirDays_ActualInds);
            day=AllBidirDays_ActualInds(j);
            postbidirDay=day-bidir_day1+1;
            
            if isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase_WithinTimeWindow(day))
                continue;
            end
            
            % --- TARG
            ffmean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase_WithinTimeWindow(day);
            
            DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(ExptCount).FFminusBase_Mean_PBS(postbidirDay)=ffmean;
            
            
            % hit rate
            hits=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).HitRateStatus_WithinTimeWindow{day});
            totalrends=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).HitRateStatus_WithinTimeWindow{day});
            HitRate=hits/totalrends;
            
            %             DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(ExptCount).HitRate_hits(postbidirDay)=hits;
            %             DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(ExptCount).HitRate_totalrends(postbidirDay)=totalrends;
            DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(ExptCount).HitRate_PBS(postbidirDay)=HitRate;
            
            totalrends1=totalrends; % to collect for global hit rate
            hits1=hits;
            
            % Daily S.D.
            ffstd=std(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals_DevFromBase_WithinTimeWindow{day});
            DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(ExptCount).ffSTD_PBS(postbidirDay)=ffstd;
            
            
            
            % --- OTHER TARG
            ffmean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(othersyl).meanFF_DevFromBase_WithinTimeWindow(day);
            DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(ExptCount).FFminusBase_Mean_PBS(postbidirDay)=ffmean;
            
            % hit rate
            hits=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(othersyl).HitRateStatus_WithinTimeWindow{day});
            totalrends=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(othersyl).HitRateStatus_WithinTimeWindow{day});
            HitRate=hits/totalrends;
            
            DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(ExptCount).HitRate_PBS(postbidirDay)=HitRate;
            
            totalrends2=totalrends; % to collect for global hit rate
            hits2=hits;
            
            % Daily S.D.
            ffstd=std(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(othersyl).FFvals_DevFromBase_WithinTimeWindow{day});
            DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(ExptCount).ffSTD_PBS(postbidirDay)=ffstd;

            
            
            % ----- COMBINED BOTH TARGETS
            HitRate_Combined_PBS=(hits1 + hits2)/(totalrends1 + totalrends2);
            DATSTRUCT_FirstBidirDayIndEqualsOne.bothTargs(ExptCount).HitRate_Combined_PBS(postbidirDay)=HitRate_Combined_PBS;
            
        end
        
        
        
        % ===== MUSC - COLLECT ONLY DAYS WITH MUSC INACTIVATION
        inds=days_with_musc_RelToBidirDay1>0 & days_with_musc_RelToBidirDay1 <= (bidir_last_day-bidir_day1+1); % these are days during bidir
        
        ActualDayInds=days_with_musc(inds);
        PostBidirInds=days_with_musc_RelToBidirDay1(inds);
        
        for j=1:length(ActualDayInds)
            day=ActualDayInds(j);
            postbidirDay=PostBidirInds(j);
            
            % ---- collect data for this day -TARG
            ffmean_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).meanFF_DevFromBase_WithinTimeWindow(day);
            DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(ExptCount).FFminusBase_Mean_MUSC(postbidirDay)=ffmean_MUSC;
            
            % hit rate
            hits=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).HitRateStatus_WithinTimeWindow{day});
            totalrends=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).HitRateStatus_WithinTimeWindow{day});
            HitRate=hits/totalrends;
            
            DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(ExptCount).HitRate_MUSC(postbidirDay)=HitRate;
            
            totalrends1=totalrends; % to collect for global hit rate
            hits1=hits;
            
            % Daily S.D.
            ffstd=std(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).FFvals_DevFromBase_WithinTimeWindow{day});
            DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(ExptCount).ffSTD_MUSC(postbidirDay)=ffstd;
            
            
            
            % ---- collect data for this day - OTHERTARG
            ffmean_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(othersyl).meanFF_DevFromBase_WithinTimeWindow(day);
            DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(ExptCount).FFminusBase_Mean_MUSC(postbidirDay)=ffmean_MUSC;
            
            % hit rate
            hits=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(othersyl).HitRateStatus_WithinTimeWindow{day});
            totalrends=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(othersyl).HitRateStatus_WithinTimeWindow{day});
            HitRate=hits/totalrends;
            
            DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(ExptCount).HitRate_MUSC(postbidirDay)=HitRate;
            
            totalrends2=totalrends; % to collect for global hit rate
            hits2=hits;
            
            % Daily S.D.
            ffstd=std(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(othersyl).FFvals_DevFromBase_WithinTimeWindow{day});
            DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(ExptCount).ffSTD_MUSC(postbidirDay)=ffstd;

            
            % ----- COMBINED BOTH TARGETS
            HitRate_Combined_MUSC=(hits1 + hits2)/(totalrends1 + totalrends2);
            DATSTRUCT_FirstBidirDayIndEqualsOne.bothTargs(ExptCount).HitRate_Combined_MUSC(postbidirDay)=HitRate_Combined_MUSC;
            
        end
        
        
        
        
        
        % ====== SAVE INFORMATION ABOUT THIS EXPERIMENT
        if (0)
%         generalization=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(othersyl).LEARNING.learning_metric.mean...
%             /SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
%         if isnan(generalization)
%             disp('PROBLEM - generalization is nan [using consolid days instead]!!!');
%             generalization=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(othersyl).LEARNING.consolid_start_meanFF_minbase...
%                 /SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.consolid_start_meanFF_minbase;
%         end
        end
        generalization=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(othersyl).meanFF_pbs/...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(targsyl).meanFF_pbs;
        
      
        
        % --- MP generalization
        generalizationMP=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(othersyl).meanFF_musc/...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(targsyl).meanFF_musc;
        disp([birdname '-' exptname ':'  num2str(generalization) '(gener) -- ' num2str(generalizationMP) '(generMP)'])

        %         disp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(othersyl).LEARNING.learning_metric.mean...
        %             /SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean);
        %         disp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(othersyl).LEARNING.consolid_start_meanFF_minbase...
        %                 /SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.consolid_start_meanFF_minbase);
        
        acoustic_distance=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(othersyl).eucldist_from_targ_zscore;
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).eucldist_from_targ_zscore~=0;
            % then targsyl is not actually targsyl!!
            disp('PROBLEM !! targsyl is not actually targsyl?');
        end
        
        correlation = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(othersyl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl);
        
        
        DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(ExptCount).day1_FromStartWN=bidir_day1-WNday1+1;
        DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(ExptCount).birdname=birdname;
        DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(ExptCount).exptname=exptname;
        DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(ExptCount).bidir_last_day_actual_RelInds=bidir_last_day_actual-bidir_day1+1;
        DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(ExptCount).targ_learn_dir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        
        % --- other stats
        DATSTRUCT_FirstBidirDayIndEqualsOne.OtherStats(ExptCount).generalization_OtherOverTarg=generalization;
        DATSTRUCT_FirstBidirDayIndEqualsOne.OtherStats(ExptCount).generalization_MP_OtherOverTarg=generalizationMP;
        DATSTRUCT_FirstBidirDayIndEqualsOne.OtherStats(ExptCount).AcDist_OtherToTarg=acoustic_distance;
        DATSTRUCT_FirstBidirDayIndEqualsOne.OtherStats(ExptCount).Corr_Song=correlation;
        
        
        
        
        
        % ========= CALCULATE SEPARATION
        % --- PBS
        ffPBS_targ1=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(ExptCount).FFminusBase_Mean_PBS;
        ffPBS_targ2=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(ExptCount).FFminusBase_Mean_PBS;
        
        % subtract, then multiply by sign of learning by targ
        signOfLearn=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(ExptCount).targ_learn_dir;
        
        Separation_PBS=(ffPBS_targ1-ffPBS_targ2).*signOfLearn;
        DATSTRUCT_FirstBidirDayIndEqualsOne.separation(ExptCount).PBS_hzInLearnDir=Separation_PBS;
        
        
        % --- MUSC
        ffMUSC_targ1=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(ExptCount).FFminusBase_Mean_MUSC;
        ffMUSC_targ2=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(ExptCount).FFminusBase_Mean_MUSC;
        
        Separation_MUSC=(ffMUSC_targ1-ffMUSC_targ2).*signOfLearn;
        DATSTRUCT_FirstBidirDayIndEqualsOne.separation(ExptCount).MUSC_hzInLearnDir=Separation_MUSC;
        
        
        
        
        
        
        % ====== TO DO: PLOT EACH EXPERIEMNT - VERIFY CORRECT
        % ====== GET AVERAGE BY USING DAY BINS
        
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

targSylTypeList={'firsttarget','secondtarget'};

DayWindowSizes=[6 5 4];
if (0) % use version below - here does each syl separately.
for lll=1:length(targSylTypeList)
    targsylType=targSylTypeList{lll};

for i=1:length(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget);
    ChosenConsolDays=[];
    
    for k=1:length(DayWindowSizes);
        
        numDaysInWind=DayWindowSizes(k);
        
        for kk=3:1000;
            % starting from day 3, keep sliding window
            
            DaysWind=kk:kk+numDaysInWind-1;
            
            % stop if the last day in window is past the last day with data
            if DaysWind(end)>length(DATSTRUCT_FirstBidirDayIndEqualsOne.(targsylType)(i).FFminusBase_Mean_PBS);
                break % break out of loop. this window size did not get anything useful. try next smallest window size
            end
            
            % CHECK THIS WINDOW
            FFvals=DATSTRUCT_FirstBidirDayIndEqualsOne.(targsylType)(i).FFminusBase_Mean_PBS(DaysWind);
            meanFF_overall=mean(FFvals);
            meanSTD=mean(DATSTRUCT_FirstBidirDayIndEqualsOne.(targsylType)(i).ffSTD_PBS(DaysWind));
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
                
                if DaysWind(end)>length(DATSTRUCT_FirstBidirDayIndEqualsOne.(targsylType)(i).FFminusBase_Mean_PBS);
                    break % break out of loop. DONE
                end
                
                % CHECK THIS WINDOW
                FFvals=DATSTRUCT_FirstBidirDayIndEqualsOne.(targsylType)(i).FFminusBase_Mean_PBS(DaysWind);
                meanFF_overall=mean(FFvals);
                meanSTD=mean(DATSTRUCT_FirstBidirDayIndEqualsOne.(targsylType)(i).ffSTD_PBS(DaysWind));
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
end

%% ================ VERSION 2 - FIRST FIND START OF WINDOW GOOD FOR BOTH
% SYLS, THEN GET EXTENDED WINDOW
DayWindowSizes=[6 5 4];
targsylType='firsttarget';
for i=1:length(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget);
    ChosenConsolDays=[];
    
    for k=1:length(DayWindowSizes);
        
        numDaysInWind=DayWindowSizes(k);
        
        for kk=3:1000;
            % starting from day 3, keep sliding window
            
            DaysWind=kk:kk+numDaysInWind-1;
            
            % stop if the last day in window is past the last day with data
            if DaysWind(end)>length(DATSTRUCT_FirstBidirDayIndEqualsOne.(targsylType)(i).FFminusBase_Mean_PBS);
                break % break out of loop. this window size did not get anything useful. try next smallest window size
            end
            
            % CHECK THIS WINDOW
            % - first target
            targsylType='firsttarget';
            FFvals=DATSTRUCT_FirstBidirDayIndEqualsOne.(targsylType)(i).FFminusBase_Mean_PBS(DaysWind);
            meanFF_overall=mean(FFvals);
            meanSTD=mean(DATSTRUCT_FirstBidirDayIndEqualsOne.(targsylType)(i).ffSTD_PBS(DaysWind));
            numZ=0.75;
            
            FFlims=[meanFF_overall-numZ*meanSTD meanFF_overall+numZ*meanSTD];
            
            % if all ffvals are within the limits, then keep this starting day.
            % keep adding days until fails the requirement of being within
            % window.
            
            func=@(x)(x>FFlims(1) && x<FFlims(2));
            withinLims1=arrayfun(func, FFvals);
            
            % - second targ
            targsylType='secondtarget';
            FFvals=DATSTRUCT_FirstBidirDayIndEqualsOne.(targsylType)(i).FFminusBase_Mean_PBS(DaysWind);
            meanFF_overall=mean(FFvals);
            meanSTD=mean(DATSTRUCT_FirstBidirDayIndEqualsOne.(targsylType)(i).ffSTD_PBS(DaysWind));
            numZ=0.75;
            
            FFlims=[meanFF_overall-numZ*meanSTD meanFF_overall+numZ*meanSTD];
            
            % if all ffvals are within the limits, then keep this starting day.
            % keep adding days until fails the requirement of being within
            % window.
            
            func=@(x)(x>FFlims(1) && x<FFlims(2));
            withinLims2=arrayfun(func, FFvals);
            
            while all(withinLims1)==1 & all(withinLims2)==1
                ChosenConsolDays=[DaysWind(1) DaysWind(end)];
                DaysWind=[DaysWind DaysWind(end)+1];
                
                if DaysWind(end)>length(DATSTRUCT_FirstBidirDayIndEqualsOne.(targsylType)(i).FFminusBase_Mean_PBS);
                    break % break out of loop. DONE
                end
                
                % CHECK THIS WINDOW
                % - first target
                targsylType='firsttarget';
                FFvals=DATSTRUCT_FirstBidirDayIndEqualsOne.(targsylType)(i).FFminusBase_Mean_PBS(DaysWind);
                meanFF_overall=mean(FFvals);
                meanSTD=mean(DATSTRUCT_FirstBidirDayIndEqualsOne.(targsylType)(i).ffSTD_PBS(DaysWind));
                numZ=0.75;
                
                FFlims=[meanFF_overall-numZ*meanSTD meanFF_overall+numZ*meanSTD];
                
                % if all ffvals are within the limits, then keep this starting day.
                % keep adding days until fails the requirement of being within
                % window.
                
                func=@(x)(x>FFlims(1) && x<FFlims(2));
                withinLims1=arrayfun(func, FFvals);
                
                % - second targ
                targsylType='secondtarget';
                FFvals=DATSTRUCT_FirstBidirDayIndEqualsOne.(targsylType)(i).FFminusBase_Mean_PBS(DaysWind);
                meanFF_overall=mean(FFvals);
                meanSTD=mean(DATSTRUCT_FirstBidirDayIndEqualsOne.(targsylType)(i).ffSTD_PBS(DaysWind));
                numZ=0.75;
                
                FFlims=[meanFF_overall-numZ*meanSTD meanFF_overall+numZ*meanSTD];
                
                % if all ffvals are within the limits, then keep this starting day.
                % keep adding days until fails the requirement of being within
                % window.
                
                func=@(x)(x>FFlims(1) && x<FFlims(2));
                withinLims2=arrayfun(func, FFvals);
            end
            
            if ~isempty(ChosenConsolDays)
                break
            end
            
            
        end
        
        if ~isempty(ChosenConsolDays)
            break
        end
        
    end
    
    Params.LMANTimeCourse.ConsolPer_IndsFromStart(i).firsttarget=ChosenConsolDays;
    Params.LMANTimeCourse.ConsolPer_IndsFromStart(i).secondtarget=ChosenConsolDays;
    
end

%% ==== replace the consol periods that were hand entered
if GetAutoConsolWindows==1;
    
for i=1:length(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget);

    birdname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).birdname;
    exptname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).exptname;
    
    % ====
    
    birdInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, birdname));
    exptInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, exptname));
    
    ind=intersect(birdInd+1, exptInd);
    
    if ~isempty(ind);
        consolPeriod_OLD=Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1};
        
        consolPeriod_NEW_ft=Params.LMANTimeCourse.ConsolPer_IndsFromStart(i).firsttarget;
        consolPeriod_NEW_st=Params.LMANTimeCourse.ConsolPer_IndsFromStart(i).secondtarget;
        
        if isempty(consolPeriod_NEW_ft) | isempty(consolPeriod_NEW_st)
        disp([birdname '-' exptname ': auto consol period EMPTY!!' ])
            continue
        end
        
        consolPeriod_NEW_olap=intersect([consolPeriod_NEW_ft(1):consolPeriod_NEW_ft(2)], [consolPeriod_NEW_st(1): consolPeriod_NEW_st(2)]);
        consolPeriod_NEW_olap=[min(consolPeriod_NEW_olap) max(consolPeriod_NEW_olap)];
        disp([birdname '-' exptname ': ' num2str(consolPeriod_OLD) ' (OLD) - ' num2str(consolPeriod_NEW_ft) ' (NEW,ft) - ' num2str(consolPeriod_NEW_st) ' (NEW,st ===== ' num2str(consolPeriod_NEW_olap) ' (NEW, OLAP!!)'])
        
        % REPLACE THOSE VALUES
        Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1}=consolPeriod_NEW_olap;
        
    end

end

end


%% ====== PLOT [EACH EXPT ONE FIG --> FF, HIT RATE, CONSOL]

for i=1:length(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget);
    
    birdname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).birdname;
    exptname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).exptname;
    
    lt_figure; hold on;
    
    % +++++ 1) FF of each targ
    lt_subplot(5,1,1); hold on;
    title([birdname '-' exptname]);
    xlabel('day (day1 = bidir day 1)');
    ylabel('shift from baseline (hz)');
    
    % == first targ
    ffPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(i).FFminusBase_Mean_PBS;
    ffMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(i).FFminusBase_Mean_MUSC;
    inds=find(~isnan(ffPBS));
    
    lt_plot(inds, ffPBS(inds), {'Color', 'b', 'LineStyle', '-'});
    lt_plot(inds, ffMUSC(inds), {'Color', 'r', 'LineStyle', '-'});
    
    % == second targ
    ffPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(i).FFminusBase_Mean_PBS;
    ffMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(i).FFminusBase_Mean_MUSC;
    inds=find(~isnan(ffPBS));
    
    lt_plot(inds, ffPBS(inds), {'Color', 'b', 'MarkerFaceColor','none','LineStyle', '-'});
    lt_plot(inds, ffMUSC(inds), {'Color', 'r', 'MarkerFaceColor','none', 'LineStyle', '-'});
    
    
    lt_plot_zeroline;
    
    % line for end of bidir
    bidir_last_day_actual=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).bidir_last_day_actual_RelInds;
    line([bidir_last_day_actual+0.5 bidir_last_day_actual+0.5], ylim, 'Color','r');
    
    % line for consol start and end
    
    birdInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, birdname));
    exptInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, exptname));
    
    ind=intersect(birdInd+1, exptInd);
    
    if ~isempty(ind);
        consolPeriod=Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1};
        
        line([consolPeriod(1)-0.4 consolPeriod(1)-0.4], ylim, 'Color','m');
        line([consolPeriod(2)+0.4 consolPeriod(2)+0.4], ylim, 'Color','m');
    end
    
    % +++++++ 2) CONSOLIDATION
    lt_subplot(5,1,2); hold on;
    ylabel('consolidation');
    
    % == extract separation
    % -- PBS
    Separation_PBS=DATSTRUCT_FirstBidirDayIndEqualsOne.separation(i).PBS_hzInLearnDir;
    
    %     inds=find(~isnan(Separation_PBS));
    %     lt_plot(inds, Separation_PBS(inds), {'Color', 'b', 'LineStyle', '-'});
    
    % -- MUSC
    Separation_MUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.separation(i).MUSC_hzInLearnDir;
    
    
    % - plot
    inds=find(~isnan(Separation_MUSC));
    
    Consolidation = Separation_MUSC(inds)./Separation_PBS(inds);
    
    lt_plot(inds, Consolidation, {'Color', 'k', 'LineStyle', '-'});
    
    line(xlim, [1 1], 'LineStyle','--')
    
    
    % line for end of bidir
    bidir_last_day_actual=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).bidir_last_day_actual_RelInds;
    line([bidir_last_day_actual+0.5 bidir_last_day_actual+0.5], ylim, 'Color','r');
    
    % line for consol start and end
    
    birdInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, birdname));
    exptInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, exptname));
    
    ind=intersect(birdInd+1, exptInd);
    
    if ~isempty(ind);
        consolPeriod=Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1};
        
        line([consolPeriod(1)-0.4 consolPeriod(1)-0.4], ylim, 'Color','m');
        line([consolPeriod(2)+0.4 consolPeriod(2)+0.4], ylim, 'Color','m');
    end
    
    xlim([0 50]); ylim([-0.2 1.2]); lt_plot_zeroline;
    
    
    
    % +++++++++++++++++ 3) HIT RATE (First target);
    lt_subplot(5,1,3); hold on;
    title('first targ');
    ylabel('hit rate');
    
    hrPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(i).HitRate_PBS;
    hrMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(i).HitRate_MUSC;
    
    inds=find(~isnan(hrPBS));
    lt_plot(inds, hrPBS(inds), {'Color', 'b', 'LineStyle', '-'});
    lt_plot(inds, hrMUSC(inds), {'Color', 'r'});
    
    lt_plot_zeroline;
    line(xlim, [1 1]);
    
    
    % +++++++++++++++++ 4) HIT RATE (second target);
    lt_subplot(5,1,4); hold on;
    title('second targ');
    ylabel('hit rate');
    
    hrPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(i).HitRate_PBS;
    hrMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(i).HitRate_MUSC;
    
    inds=find(~isnan(hrPBS));
    
    lt_plot(inds, hrPBS(inds), {'Color', 'b', 'MarkerFaceColor','none','LineStyle', '-'});
    lt_plot(inds, hrMUSC(inds), {'Color', 'r', 'MarkerFaceColor','none'});
    
    
    lt_plot_zeroline;
    line(xlim, [1 1]);
    
    
    % +++++++++++++++ 5) COMBINED HIT RATE (average)
    lt_subplot(5,1,5); hold on;
    title('both combined (actual rends, not average)');
    ylabel('hit rate');
    
    hrPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.bothTargs(i).HitRate_Combined_PBS;
    hrMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.bothTargs(i).HitRate_Combined_MUSC;
    
    inds=find(~isnan(hrPBS));
    
    lt_plot(inds, hrPBS(inds), {'Color', 'b', 'MarkerFaceColor','none','LineStyle', '-'});
    lt_plot(inds, hrMUSC(inds), {'Color', 'r', 'MarkerFaceColor','none'});
    
    
    lt_plot_zeroline;
    line(xlim, [1 1]);
    
    
end




%% ====== PLOT [EACH TARG, + CONSOLIDATION METRIC]

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

for i=1:length(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget);
    
    birdname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).birdname;
    exptname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([birdname '-' exptname]);
    xlabel('day (day1 = bidir day 1)');
    ylabel('shift from baseline (hz)');
    xlim([0 50])
    
    % == first targ
    ffPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(i).FFminusBase_Mean_PBS;
    ffMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(i).FFminusBase_Mean_MUSC;
    inds=find(~isnan(ffPBS));
    
    lt_plot(inds, ffPBS(inds), {'Color', 'b', 'LineStyle', '-'});
    lt_plot(inds, ffMUSC(inds), {'Color', 'r', 'LineStyle', '-'});
    
    
    % == second targ
    ffPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(i).FFminusBase_Mean_PBS;
    ffMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(i).FFminusBase_Mean_MUSC;
    inds=find(~isnan(ffPBS));
    
    lt_plot(inds, ffPBS(inds), {'Color', 'b', 'MarkerFaceColor','none','LineStyle', '-'});
    lt_plot(inds, ffMUSC(inds), {'Color', 'r', 'MarkerFaceColor','none', 'LineStyle', '-'});
    
    
    lt_plot_zeroline;
    
    % line for end of bidir
    bidir_last_day_actual=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).bidir_last_day_actual_RelInds;
    line([bidir_last_day_actual+0.5 bidir_last_day_actual+0.5], ylim, 'Color','r');
    
    % line for consol start and end
    
    birdInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, birdname));
    exptInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, exptname));
    
    ind=intersect(birdInd+1, exptInd);
    
    if ~isempty(ind);
        consolPeriod=Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1};
        
        line([consolPeriod(1)-0.4 consolPeriod(1)-0.4], ylim, 'Color','m');
        line([consolPeriod(2)+0.4 consolPeriod(2)+0.4], ylim, 'Color','m');
        
        % ================= PUT MEAN PICH +/- 0.75 STD DURING MAINTAINED
        % SHIFT PERIOD
        numZ=0.75;
        
        % - first targ
        meanFF=mean(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(i).FFminusBase_Mean_PBS(consolPeriod(1):consolPeriod(2)));
        meanStdFF=mean(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(i).ffSTD_PBS(consolPeriod(1):consolPeriod(2)));
        
        Xlim=xlim;
               
        shadedErrorBar(Xlim, [meanFF meanFF], [numZ*meanStdFF numZ*meanStdFF],{'Color','b'},1);
        
        % - second targ
        meanFF=mean(DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(i).FFminusBase_Mean_PBS(consolPeriod(1):consolPeriod(2)));
        meanStdFF=mean(DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(i).ffSTD_PBS(consolPeriod(1):consolPeriod(2)));
        
        Xlim=xlim;
               
        shadedErrorBar(Xlim, [meanFF meanFF], [numZ*meanStdFF numZ*meanStdFF],{'Color','b'},1);
        
    end
    
    
    
    
    % ============ CONSOLIDATION
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([birdname '-' exptname]);
    xlabel('day (day1 = bidir day 1)');
    ylabel('consolidation');
    
    % == extract separation
    % -- PBS
    Separation_PBS=DATSTRUCT_FirstBidirDayIndEqualsOne.separation(i).PBS_hzInLearnDir;
    
    %     inds=find(~isnan(Separation_PBS));
    %     lt_plot(inds, Separation_PBS(inds), {'Color', 'b', 'LineStyle', '-'});
    
    % -- MUSC
    Separation_MUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.separation(i).MUSC_hzInLearnDir;
    
    
    % - plot
    inds=find(~isnan(Separation_MUSC));
    
    Consolidation = Separation_MUSC(inds)./Separation_PBS(inds);
    
    lt_plot(inds, Consolidation, {'Color', 'k', 'LineStyle', '-'});
    
    line(xlim, [1 1], 'LineStyle','--')
    
    
    % line for end of bidir
    bidir_last_day_actual=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).bidir_last_day_actual_RelInds;
    line([bidir_last_day_actual+0.5 bidir_last_day_actual+0.5], ylim, 'Color','r');
    
    % line for consol start and end
    
    birdInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, birdname));
    exptInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, exptname));
    
    ind=intersect(birdInd+1, exptInd);
    
    if ~isempty(ind);
        consolPeriod=Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1};
        
        line([consolPeriod(1)-0.4 consolPeriod(1)-0.4], ylim, 'Color','m');
        line([consolPeriod(2)+0.4 consolPeriod(2)+0.4], ylim, 'Color','m');
    end
    
    xlim([0 50]); ylim([-0.2 1.2]); lt_plot_zeroline;
    
    
end

lt_subtitle('red=WN end;, magen=consol period');


%% ==== PLOT (IN TERMS OF SEPARATION (I.E COMBINING INFO FOR THE TWO TARGETS)

figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

for i=1:length(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget);
    
    birdname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).birdname;
    exptname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([birdname '-' exptname]);
    xlabel('day (day1 = bidir day 1)');
    ylabel('Separation (hz)');
    
    
    
    % == extract separation
    % -- PBS
    Separation_PBS=DATSTRUCT_FirstBidirDayIndEqualsOne.separation(i).PBS_hzInLearnDir;
    
    inds=find(~isnan(Separation_PBS));
    lt_plot(inds, Separation_PBS(inds), {'Color', 'b', 'LineStyle', '-'});
    
    % -- MUSC
    Separation_MUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.separation(i).MUSC_hzInLearnDir;
    
    
    
    % - plot
    inds=find(~isnan(Separation_MUSC));
    lt_plot(inds, Separation_MUSC(inds), {'Color', 'r', 'LineStyle', '-'});
    
    
    
    lt_plot_zeroline;
    
    % line for end of bidir
    bidir_last_day_actual=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).bidir_last_day_actual_RelInds;
    line([bidir_last_day_actual+0.5 bidir_last_day_actual+0.5], ylim, 'Color','r');
    
    % line for consol start and end
    
    birdInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, birdname));
    exptInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, exptname));
    
    ind=intersect(birdInd+1, exptInd);
    
    if ~isempty(ind);
        consolPeriod=Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1};
        
        line([consolPeriod(1)-0.4 consolPeriod(1)-0.4], ylim, 'Color','m');
        line([consolPeriod(2)+0.4 consolPeriod(2)+0.4], ylim, 'Color','m');
    end
    
    
end

lt_subtitle('red=WN end;, magen=consol period');


%% ==== PLOT (IN TERMS OF CONSOLIDATION (I.E COMBINING INFO FOR THE TWO TARGETS)
if (0);
    figcount=1;
    subplotrows=3;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    
    for i=1:length(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget);
        
        birdname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).birdname;
        exptname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).exptname;
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname '-' exptname]);
        xlabel('day (day1 = bidir day 1)');
        ylabel('consolidation');
        
        
        
        % == extract separation
        % -- PBS
        Separation_PBS=DATSTRUCT_FirstBidirDayIndEqualsOne.separation(i).PBS_hzInLearnDir;
        
        %     inds=find(~isnan(Separation_PBS));
        %     lt_plot(inds, Separation_PBS(inds), {'Color', 'b', 'LineStyle', '-'});
        
        % -- MUSC
        Separation_MUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.separation(i).MUSC_hzInLearnDir;
        
        
        % - plot
        inds=find(~isnan(Separation_MUSC));
        
        Consolidation = Separation_MUSC(inds)./Separation_PBS(inds);
        
        lt_plot(inds, Consolidation, {'Color', 'k', 'LineStyle', '-'});
        
        line(xlim, [1 1], 'LineStyle','--')
        
        
        % line for end of bidir
        bidir_last_day_actual=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).bidir_last_day_actual_RelInds;
        line([bidir_last_day_actual+0.5 bidir_last_day_actual+0.5], ylim, 'Color','r');
        
        % line for consol start and end
        
        birdInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, birdname));
        exptInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, exptname));
        
        ind=intersect(birdInd+1, exptInd);
        
        if ~isempty(ind);
            consolPeriod=Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1};
            
            line([consolPeriod(1)-0.4 consolPeriod(1)-0.4], ylim, 'Color','m');
            line([consolPeriod(2)+0.4 consolPeriod(2)+0.4], ylim, 'Color','m');
        end
        
        xlim([0 50]); ylim([-0.2 1.2]); lt_plot_zeroline;
        
    end
    
    lt_subtitle('red=WN end;, magen=consol period');
    
end

%% ==== PLOT (IN TERMS OF SEPARATION (I.E COMBINING INFO FOR THE TWO TARGETS) [plot values]
if (0)
    figcount=1;
    subplotrows=3;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    
    for i=1:length(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget);
        
        birdname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).birdname;
        exptname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).exptname;
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname '-' exptname]);
        xlabel('day (day1 = bidir day 1)');
        ylabel('Separation (hz)');
        
        
        
        % == extract separation
        % -- PBS
        Separation_PBS=DATSTRUCT_FirstBidirDayIndEqualsOne.separation(i).PBS_hzInLearnDir;
        
        inds=find(~isnan(Separation_PBS));
        lt_plot(inds, Separation_PBS(inds), {'Color', 'b', 'LineStyle', '-'});
        
        for k=inds
            lt_plot_text(k, Separation_PBS(k), num2str(Separation_PBS(k)), 'b');
        end
        
        % -- MUSC
        Separation_MUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.separation(i).MUSC_hzInLearnDir;
        
        for k=inds
            lt_plot_text(k, Separation_MUSC(k), num2str(Separation_MUSC(k)), 'r');
        end
        
        
        % - plot
        inds=find(~isnan(Separation_MUSC));
        lt_plot(inds, Separation_MUSC(inds), {'Color', 'r', 'LineStyle', '-'});
        
        
        
        lt_plot_zeroline;
        
        % line for end of bidir
        bidir_last_day_actual=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).bidir_last_day_actual_RelInds;
        line([bidir_last_day_actual+0.5 bidir_last_day_actual+0.5], ylim, 'Color','r');
        
        % line for consol start and end
        
        birdInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, birdname));
        exptInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, exptname));
        
        ind=intersect(birdInd+1, exptInd);
        
        if ~isempty(ind);
            consolPeriod=Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1};
            
            line([consolPeriod(1)-0.4 consolPeriod(1)-0.4], ylim, 'Color','m');
            line([consolPeriod(2)+0.4 consolPeriod(2)+0.4], ylim, 'Color','m');
        end
        
        
    end
    
    lt_subtitle('red=WN end;, magen=consol period');
    
    
    
end

%% ====== PLOT [no WN lines]
%
% figcount=1;
% subplotrows=3;
% subplotcols=2;
% fignums_alreadyused=[];
% hfigs=[];
%
% for i=1:length(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget);
%
% birdname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).birdname;
% exptname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(i).exptname;
%
%         [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% title([birdname '-' exptname]);
% xlabel('day (day1 = bidir day 1)');
% ylabel('shift from baseline (hz)');
%
%
%     % == first targ
%     ffPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(i).FFminusBase_Mean_PBS;
%     ffMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(i).FFminusBase_Mean_MUSC;
%     inds=find(~isnan(ffPBS));
%
%     lt_plot(inds, ffPBS(inds), {'Color', 'b', 'LineStyle', '-'});
%     lt_plot(inds, ffMUSC(inds), {'Color', 'r', 'LineStyle', '-'});
%
%     % == second targ
%     ffPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(i).FFminusBase_Mean_PBS;
%     ffMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(i).FFminusBase_Mean_MUSC;
%     inds=find(~isnan(ffPBS));
%
%     lt_plot(inds, ffPBS(inds), {'Color', 'b', 'MarkerFaceColor','none','LineStyle', '-'});
%     lt_plot(inds, ffMUSC(inds), {'Color', 'r', 'MarkerFaceColor','none', 'LineStyle', '-'});
%
%     lt_plot_zeroline;
%
%     % line for consol start and end
%
%     birdInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, birdname));
%     exptInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, exptname));
%
%     ind=intersect(birdInd+1, exptInd);
%
%     if ~isempty(ind);
%         consolPeriod=Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1};
%
%         line([consolPeriod(1)-0.4 consolPeriod(1)-0.4], ylim, 'Color','m');
%         line([consolPeriod(2)+0.4 consolPeriod(2)+0.4], ylim, 'Color','m');
%     end
%
%
% end
%
% lt_subtitle('red=WN end;, magen=consol period');
%

%% === extract consolidation days


for j=1:length(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget); % num expts
    birdname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).birdname;
    exptname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).exptname;
    
    
    birdInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, birdname));
    exptInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, exptname));
    
    ind=intersect(birdInd+1, exptInd);
    
    if ~isempty(ind);
        consolPeriod=Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1};
        
        DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).consolPeriod=consolPeriod;
        
    else
        disp(['WNARNING - no consol period entered for ' birdname '-' exptname ' - using all days']);
        
    end
end


%% ==== PLOT MEAN ACROSS ALL EXPERIMENTS FOR BIDIR
% ==== DON'T CARE ABOUT CONSOLIDATION - SIMPLY BIN IN DAYS RELATIVE TO
% LEARNING START

% === to add:
% 1) stop when reach end of bidir
% 3) only run if bird has data for all bins

if BinRelConsolDay1==0;
    
    
    maxday=length([DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget.FFminusBase_Mean_PBS]);
    DayBins_StartEdges=1:NumDaysInBin:maxday;
    
    OUTPUT=struct; % one cell for each bin
    Exptcount=1;
    
    for i=1:length(DayBins_StartEdges);
        
        OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC=[];
        OUTPUT.firsttarget.FFrelBaseInBin(i).PBS=[];
        OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC=[];
        OUTPUT.secondtarg.FFrelBaseInBin(i).PBS=[];
        OUTPUT.firsttarget.FFrelBaseInBin(i).MPbias=[];
        OUTPUT.secondtarg.FFrelBaseInBin(i).MPbias=[];
        
        OUTPUT.firsttarget.otherstats(i).hitratePBS=[];
        OUTPUT.firsttarget.otherstats(i).hitrateMUSC=[];
        OUTPUT.secondtarg.otherstats(i).hitratePBS=[];
        OUTPUT.secondtarg.otherstats(i).hitrateMUSC=[];
        
        OUTPUT.combined.otherstats(i).hitratePBS=[];
        OUTPUT.combined.otherstats(i).hitrateMUSC=[];
       
        
        
        OUTPUT.separation(i).PBS=[];
        OUTPUT.separation(i).MUSC=[];
        
        OUTPUT.INFORMATION(i).experimentNum=[];
        OUTPUT.INFORMATION(i).birdname={};
        OUTPUT.INFORMATION(i).exptname={};
        
        DaysInBin=DayBins_StartEdges(i):DayBins_StartEdges(i)+NumDaysInBin-1;
        
        for day=DaysInBin;
            
            
            
            
            % == for this day, collect all MUSC. for days with musc, also
            % collect PBS
            
            for j=1:length(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget); % num expts
                
                % ++++++++++++++++++++++++++++++++++++++
                % === only continue if vector is long enough to today
                if length(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(j).FFminusBase_Mean_MUSC)<day;
                    continue;
                end
                
                % === only continue if first target has MUSC data for today
                if isnan(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(j).FFminusBase_Mean_MUSC(day));
                    continue;
                end
                
                % ==== make sure this day is in consolidation period (if
                % desired)
                if OnlyConsolPeriod_bidir==1;
                    birdname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).birdname;
                    exptname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).exptname;
                    
                    
                    birdInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, birdname));
                    exptInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, exptname));
                    
                    ind=intersect(birdInd+1, exptInd);
                    
                    if ~isempty(ind);
                        consolPeriod=Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1};
                        
                        if day<consolPeriod(1) | day>consolPeriod(2)
                            continue
                        end
                        
                    else
                        disp(['WNARNING - no consol period entered for ' birdname '-' exptname ' - using all days']);
                        
                    end
                end
                
                
                
                % ++++++++++++++++++++++++++
                
                % ==================== BOTH TARGS
                hrCombinedPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.bothTargs(j).HitRate_Combined_PBS(day);
                hrCombinedMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.bothTargs(j).HitRate_Combined_MUSC(day);

                OUTPUT.combined.otherstats(i).hitratePBS=[OUTPUT.combined.otherstats(i).hitratePBS hrCombinedPBS];
                OUTPUT.combined.otherstats(i).hitrateMUSC=[OUTPUT.combined.otherstats(i).hitrateMUSC hrCombinedMUSC];
                
                
                % ============== first target
                ffMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(j).FFminusBase_Mean_MUSC(day); % MUSC
                ffPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(j).FFminusBase_Mean_PBS(day); % PBS
                
                % --- flip sign if learning was neg
                ffMUSC=ffMUSC*DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                ffPBS=ffPBS*DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                
                % ---- REVERSION (MP BIAS)
                MPbias=ffMUSC/ffPBS;
                                
                % --- HIT RATE
                hrPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(j).HitRate_PBS(day);
                hrMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(j).HitRate_MUSC(day);
                         
                % ---- OUTPUT
                OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC=[OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC ffMUSC];
                OUTPUT.firsttarget.FFrelBaseInBin(i).PBS=[OUTPUT.firsttarget.FFrelBaseInBin(i).PBS ffPBS];
                
                OUTPUT.firsttarget.FFrelBaseInBin(i).MPbias=[OUTPUT.firsttarget.FFrelBaseInBin(i).MPbias MPbias];
                
                OUTPUT.firsttarget.otherstats(i).hitratePBS=[OUTPUT.firsttarget.otherstats(i).hitratePBS hrPBS];
                OUTPUT.firsttarget.otherstats(i).hitrateMUSC=[OUTPUT.firsttarget.otherstats(i).hitrateMUSC hrMUSC];
                    

                
                
                
                
                % ============= second target
                ffMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(j).FFminusBase_Mean_MUSC(day); % MUSC
                ffPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(j).FFminusBase_Mean_PBS(day); % PBS
                
                % --- flip sign if learning was neg
                ffMUSC=ffMUSC*DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                ffPBS=ffPBS*DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                
                % ---- REVERSION (MP BIAS)
                MPbias=ffMUSC/ffPBS;
                
                                % --- HIT RATE
                hrPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(j).HitRate_PBS(day);
                hrMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(j).HitRate_MUSC(day);

                
                % ---- OUTPUT
                OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC=[OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC ffMUSC];
                OUTPUT.secondtarg.FFrelBaseInBin(i).PBS=[OUTPUT.secondtarg.FFrelBaseInBin(i).PBS ffPBS];
                
                OUTPUT.secondtarg.FFrelBaseInBin(i).MPbias=[OUTPUT.secondtarg.FFrelBaseInBin(i).MPbias MPbias];
                                
                OUTPUT.secondtarg.otherstats(i).hitratePBS=[OUTPUT.secondtarg.otherstats(i).hitratePBS hrPBS];
                OUTPUT.secondtarg.otherstats(i).hitrateMUSC=[OUTPUT.secondtarg.otherstats(i).hitrateMUSC hrMUSC];

                
                % ======== SEPARATION
                separ_PBS=DATSTRUCT_FirstBidirDayIndEqualsOne.separation(j).PBS_hzInLearnDir(day);
                OUTPUT.separation(i).PBS = [OUTPUT.separation(i).PBS separ_PBS];
                
                separ_MUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.separation(j).MUSC_hzInLearnDir(day);
                OUTPUT.separation(i).MUSC = [OUTPUT.separation(i).MUSC separ_MUSC];
                
                
                % ===== INFORMATION
                
                OUTPUT.INFORMATION(i).experimentNum=[OUTPUT.INFORMATION(i).experimentNum j];
                OUTPUT.INFORMATION(i).birdname=[OUTPUT.INFORMATION(i).birdname birdname];
                OUTPUT.INFORMATION(i).exptname=[OUTPUT.INFORMATION(i).exptname exptname];
                Exptcount=Exptcount+1;
            end
            
        end
    end
    
elseif BinRelConsolDay1==1;
    %% ==== PLOT MEAN ACROSS ALL EXPERIMENTS FOR BIDIR  [consol day 1 = bin 1]
    % === to add:
    % 1) stop when reach end of bidir
    % 3) only run if bird has data for all bins
    
    
    maxday=length([DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget.FFminusBase_Mean_PBS]);
    DayBins_StartEdges=1:NumDaysInBin:maxday;
    
    OUTPUT=struct; % one cell for each bin
    Exptcount=1;
    
    for i=1:length(DayBins_StartEdges);
        
        OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC=[];
        OUTPUT.firsttarget.FFrelBaseInBin(i).PBS=[];
        OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC=[];
        OUTPUT.secondtarg.FFrelBaseInBin(i).PBS=[];
        OUTPUT.firsttarget.FFrelBaseInBin(i).MPbias=[];
        OUTPUT.secondtarg.FFrelBaseInBin(i).MPbias=[];
        
        OUTPUT.firsttarget.otherstats(i).hitratePBS=[];
        OUTPUT.firsttarget.otherstats(i).hitrateMUSC=[];
        OUTPUT.secondtarg.otherstats(i).hitratePBS=[];
        OUTPUT.secondtarg.otherstats(i).hitrateMUSC=[];
        
        OUTPUT.combined.otherstats(i).hitratePBS=[];
        OUTPUT.combined.otherstats(i).hitrateMUSC=[];
        
        
        OUTPUT.separation(i).PBS=[];
        OUTPUT.separation(i).MUSC=[];
        
        OUTPUT.INFORMATION(i).experimentNum=[];
        OUTPUT.INFORMATION(i).birdname={};
        OUTPUT.INFORMATION(i).exptname={};
        
        DaysInBin=DayBins_StartEdges(i):DayBins_StartEdges(i)+NumDaysInBin-1;
        
        for day=DaysInBin;
            
            
            
            
            % == for this day, collect all MUSC. for days with musc, also
            % collect PBS
            
            for j=1:length(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget); % num expts
                
                birdname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).birdname;
                exptname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).exptname;
                
                
                % === CONVERT TO DAY RELATIVE TO START OF CONSOLIDATION (CONSOL DAY 1 =
                % IND 1)
                if isempty(DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).consolPeriod)
                    
                    disp(['[skipping] No consol period for ' birdname '-' exptname])
                    continue
                end
                
                
                consolDay1_relBidirDay1=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).consolPeriod(1);
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
                    consolPeriod=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).consolPeriod;
                    
                    if day>consolPeriod(2)
                        continue
                    end
                end
                
                % === only continue if vector is long enough to today
                if length(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(j).FFminusBase_Mean_MUSC)<day;
                    continue;
                end
                
                % === only continue if first target has MUSC data for today
                if isnan(DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(j).FFminusBase_Mean_MUSC(day));
                    continue;
                end
                
                
                % ==================== BOTH TARGS
                hrCombinedPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.bothTargs(j).HitRate_Combined_PBS(day);
                hrCombinedMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.bothTargs(j).HitRate_Combined_MUSC(day);
                
                OUTPUT.combined.otherstats(daybin).hitratePBS=[OUTPUT.combined.otherstats(daybin).hitratePBS hrCombinedPBS];
                OUTPUT.combined.otherstats(daybin).hitrateMUSC=[OUTPUT.combined.otherstats(daybin).hitrateMUSC hrCombinedMUSC];
                
                
                % ++++++++++++++++++++++++++
                % ============== first target
                ffMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(j).FFminusBase_Mean_MUSC(day); % MUSC
                ffPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(j).FFminusBase_Mean_PBS(day); % PBS
                
                % --- flip sign if learning was neg
                ffMUSC=ffMUSC*DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                ffPBS=ffPBS*DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                
                % ---- REVERSION (MP BIAS)
                MPbias=ffMUSC/ffPBS;
                
                % --- HIT RATE
                hrPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(j).HitRate_PBS(day);
                hrMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.firsttarget(j).HitRate_MUSC(day);
               
                % ---- OUTPUT
                OUTPUT.firsttarget.FFrelBaseInBin(daybin).MUSC=[OUTPUT.firsttarget.FFrelBaseInBin(daybin).MUSC ffMUSC];
                OUTPUT.firsttarget.FFrelBaseInBin(daybin).PBS=[OUTPUT.firsttarget.FFrelBaseInBin(daybin).PBS ffPBS];
                OUTPUT.firsttarget.FFrelBaseInBin(daybin).MPbias=[OUTPUT.firsttarget.FFrelBaseInBin(daybin).MPbias MPbias];
                
                OUTPUT.firsttarget.otherstats(daybin).hitratePBS=[OUTPUT.firsttarget.otherstats(daybin).hitratePBS hrPBS];
                OUTPUT.firsttarget.otherstats(daybin).hitrateMUSC=[OUTPUT.firsttarget.otherstats(daybin).hitrateMUSC hrMUSC];
                
                
                
                % ============= second target
                ffMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(j).FFminusBase_Mean_MUSC(day); % MUSC
                ffPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(j).FFminusBase_Mean_PBS(day); % PBS
                
                % --- flip sign if learning was neg
                ffMUSC=ffMUSC*DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                ffPBS=ffPBS*DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(j).targ_learn_dir;
                
                % ---- REVERSION (MP BIAS)
                MPbias=ffMUSC/ffPBS;
                
                % --- HIT RATE
                hrPBS=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(j).HitRate_PBS(day);
                hrMUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.secondtarget(j).HitRate_MUSC(day);

                % ---- OUTPUT
                OUTPUT.secondtarg.FFrelBaseInBin(daybin).MUSC=[OUTPUT.secondtarg.FFrelBaseInBin(daybin).MUSC ffMUSC];
                OUTPUT.secondtarg.FFrelBaseInBin(daybin).PBS=[OUTPUT.secondtarg.FFrelBaseInBin(daybin).PBS ffPBS];
                
                OUTPUT.secondtarg.FFrelBaseInBin(daybin).MPbias=[OUTPUT.secondtarg.FFrelBaseInBin(daybin).MPbias MPbias];
                
                OUTPUT.secondtarg.otherstats(daybin).hitratePBS=[OUTPUT.secondtarg.otherstats(daybin).hitratePBS hrPBS];
                OUTPUT.secondtarg.otherstats(daybin).hitrateMUSC=[OUTPUT.secondtarg.otherstats(daybin).hitrateMUSC hrMUSC];
                
                
                % ======== SEPARATION
                separ_PBS=DATSTRUCT_FirstBidirDayIndEqualsOne.separation(j).PBS_hzInLearnDir(day);
                OUTPUT.separation(daybin).PBS = [OUTPUT.separation(daybin).PBS separ_PBS];
                
                separ_MUSC=DATSTRUCT_FirstBidirDayIndEqualsOne.separation(j).MUSC_hzInLearnDir(day);
                OUTPUT.separation(daybin).MUSC = [OUTPUT.separation(daybin).MUSC separ_MUSC];
                
                
                % ===== INFORMATION
                
                OUTPUT.INFORMATION(daybin).experimentNum=[OUTPUT.INFORMATION(daybin).experimentNum j];
                OUTPUT.INFORMATION(daybin).birdname=[OUTPUT.INFORMATION(daybin).birdname birdname];
                OUTPUT.INFORMATION(daybin).exptname=[OUTPUT.INFORMATION(daybin).exptname exptname];
                
                
                
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
    MUSC_ft=OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC;
    PBS_ft=OUTPUT.firsttarget.FFrelBaseInBin(i).PBS;
    
    MUSC_st=OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC;
    PBS_st=OUTPUT.secondtarg.FFrelBaseInBin(i).PBS;
    
    SeparationMUSC=OUTPUT.separation(i).MUSC;
    SeparationPBS=OUTPUT.separation(i).PBS;
    
    
    hr_first_PBS=OUTPUT.firsttarget.otherstats(i).hitratePBS;
    hr_first_MUSC=OUTPUT.firsttarget.otherstats(i).hitrateMUSC;
    hr_second_PBS=OUTPUT.secondtarg.otherstats(i).hitratePBS;
    hr_second_MUSC=OUTPUT.secondtarg.otherstats(i).hitrateMUSC;
    
    hr_combined_PBS=OUTPUT.combined.otherstats(i).hitratePBS;
    hr_combined_MUSC=OUTPUT.combined.otherstats(i).hitrateMUSC;
    
    
    exptNums=OUTPUT.INFORMATION(i).experimentNum;
    
    
    % === prepare output
    %     MPbias_WithinExptAveraged=[];
    Exptnum_out=[];
    MUSC_ft_out=[];
    PBS_ft_out=[];
    MUSC_st_out=[];
    PBS_st_out=[];
    
    SeparationMUSC_out=[];
    SeparationPBS_out=[];
    
    hr_first_PBS_out=[];
    hr_first_MUSC_out=[];
    hr_second_PBS_out=[];
    hr_second_MUSC_out=[];
    
    hr_comb_PBS_out=[];
    hr_comb_MUSC_out=[];
    
    
    % go thru each exptnum, if there are multiple cases, then avearge them
    MaxExptNum=max(exptNums);
    
    for k=1:MaxExptNum;
        
        inds=find(exptNums==k);
        
        if ~isempty(inds)
%             if length(inds)>1
                % average them
                %                 mpbias=mean(MPbias_ft(inds));
                musc_ft=mean(MUSC_ft(inds));
                pbs_ft=mean(PBS_ft(inds));
                
                musc_st=mean(MUSC_st(inds));
                pbs_st=mean(PBS_st(inds));
                
                separation_musc=mean(SeparationMUSC(inds));
                separation_pbs=mean(SeparationPBS(inds));
                
                hr_ft_pbs=mean(hr_first_PBS(inds));
                hr_ft_musc=mean(hr_first_MUSC(inds));
                hr_st_pbs=mean(hr_second_PBS(inds));
                hr_st_musc=mean(hr_second_MUSC(inds));
                
                hr_comb_pbs=mean(hr_combined_PBS(inds));
                hr_comb_musc=mean(hr_combined_MUSC(inds));
                
            % ===== add to new vectors
            %             MPbias_WithinExptAveraged=[MPbias_WithinExptAveraged mpbias];
            Exptnum_out=[Exptnum_out k];
            MUSC_ft_out=[MUSC_ft_out musc_ft];
            PBS_ft_out=[PBS_ft_out pbs_ft];
            MUSC_st_out=[MUSC_st_out musc_st];
            PBS_st_out=[PBS_st_out pbs_st];
            
            SeparationMUSC_out=[SeparationMUSC_out separation_musc];
            SeparationPBS_out=[SeparationPBS_out separation_pbs];
            
            hr_first_PBS_out=[hr_first_PBS_out hr_ft_pbs];
            hr_first_MUSC_out=[hr_first_MUSC_out hr_ft_musc];
            hr_second_PBS_out=[hr_second_PBS_out hr_st_pbs];
            hr_second_MUSC_out=[hr_second_MUSC_out hr_st_musc];
            
            hr_comb_PBS_out=[hr_comb_PBS_out hr_comb_pbs];
            hr_comb_MUSC_out=[hr_comb_MUSC_out hr_comb_musc];
            
        end
        
    end
    
    
    % ====== REPLACE OLD VECTORS
    %     OUTPUT_WithinExptAvg.firsttarget.FFrelBaseInBin(i).MPbias=MPbias_WithinExptAveraged;
    OUTPUT_WithinExptAvg.INFORMATION(i).experimentNum=Exptnum_out;
    
    OUTPUT_WithinExptAvg.firsttarget.FFrelBaseInBin(i).MUSC=MUSC_ft_out;
    OUTPUT_WithinExptAvg.firsttarget.FFrelBaseInBin(i).PBS=PBS_ft_out;
    
    OUTPUT_WithinExptAvg.secondtarg.FFrelBaseInBin(i).MUSC=MUSC_st_out;
    OUTPUT_WithinExptAvg.secondtarg.FFrelBaseInBin(i).PBS=PBS_st_out;
    
    OUTPUT_WithinExptAvg.separation(i).MUSC=SeparationMUSC_out;
    OUTPUT_WithinExptAvg.separation(i).PBS=SeparationPBS_out;
    
    OUTPUT_WithinExptAvg.firsttarget.otherstats(i).hitratePBS=hr_first_PBS_out;
    OUTPUT_WithinExptAvg.firsttarget.otherstats(i).hitrateMUSC=hr_first_MUSC_out;
    OUTPUT_WithinExptAvg.secondtarg.otherstats(i).hitratePBS=hr_second_PBS_out;
    OUTPUT_WithinExptAvg.secondtarg.otherstats(i).hitrateMUSC=hr_second_MUSC_out;

    OUTPUT_WithinExptAvg.combined.otherstats(i).hitratePBS=hr_comb_PBS_out;
    OUTPUT_WithinExptAvg.combined.otherstats(i).hitrateMUSC=hr_comb_MUSC_out;

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
    [~, p]=ttest(yPBS, yMUSC);
    [p]=signrank(yPBS, yMUSC);
    if p<0.1;
        lt_plot_text(x(1), max([yPBS yMUSC]), num2str(p, '%3.2g'));
    end
    
            % -- plot text (reversion)
    reverVals=yMUSC-yPBS;
    lt_plot_text(x(1), 1.3*max([yPBS yMUSC]), ['PBSmean(sem)=' num2str(mean(yPBS)) '(' num2str(lt_sem(yPBS)) ')'] ,'m')
    lt_plot_text(x(1), 1.2*max([yPBS yMUSC]), ['MUSCmean(sem)=' num2str(mean(yMUSC)) '(' num2str(lt_sem(yMUSC)) ')'], 'm')

    % === second targ
    x=[i+subdivisions i+2*subdivisions];
    
    % -- plot each point
    yPBS=OUTPUT.secondtarg.FFrelBaseInBin(i).PBS;
    yMUSC=OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC;
    
    plot(x, [yPBS' yMUSC']', '-ob');
    
    % -- means
    if size([yPBS' yMUSC'],1)>1
        %       lt_plot_bar(x+0.05, mean([yPBS' yMUSC']), {'Color', 'b', 'Errors', lt_sem([yPBS' yMUSC'])});
        lt_plot_bar(x+subdivisions/4, mean([yPBS' yMUSC']), {'Color', 'b'});
    end
    
    % stats
%     [~, p]=ttest(yPBS, yMUSC);
    [p]=signrank(yPBS, yMUSC);
    if p<0.1;
        lt_plot_text(x(1), max([yPBS yMUSC]), num2str(p, '%3.2g'));
    end
    
        
            % -- plot text (reversion)
    reverVals=yMUSC-yPBS;
    lt_plot_text(x(1), 1.3*max([yPBS yMUSC]), ['PBSmean(sem)=' num2str(mean(yPBS)) '(' num2str(lt_sem(yPBS)) ')'] ,'r')
    lt_plot_text(x(1), 1.2*max([yPBS yMUSC]), ['MUSCmean(sem)=' num2str(mean(yMUSC)) '(' num2str(lt_sem(yMUSC)) ')'], 'r')

    
    
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
    
    
    % === second targ
    x=[i+subdivisions i+2*subdivisions];
    
    % -- plot each point
    yPBS=OUTPUT.secondtarg.FFrelBaseInBin(i).PBS;
    yMUSC=OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC;
    
    % -- means
    if size([yPBS' yMUSC'],1)>1
        lt_plot_bar(x-0.05, -mean([yPBS' yMUSC'])./mean([yPBS' yPBS']), {'Color', 'b', 'Errors', lt_sem([yPBS' yMUSC'])./mean([yPBS' yPBS'])});
    end
    
end



%% PLOT (SEPARATION)
lt_figure; hold on;

lt_subplot(3,1,1); hold on;
title('Raw SEPARATION, from bidir day 1, all experiments');
xlabel([ 'day bin (binsize: ' num2str(NumDaysInBin) ')']);
ylabel('SEPARATION (hz, in targ dir)');

X=1:length(DayBins_StartEdges);
subdivisions=0.1; % for plotting diff things in one ind

for i=1:length(X);
    
    if isempty(OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC)
        continue
    end
    
    % ==== SEPARATION
    x=[i-subdivisions*2 i-subdivisions];
    
    
    Y=[OUTPUT.separation(i).PBS' OUTPUT.separation(i).MUSC'];
    
    plot(x, Y', '-ok');
    
    % -- means
    if size(Y,1)>1
        lt_plot_bar(x-subdivisions/4, mean(Y), {'Color', 'k'});
    end
    
    
    
    % stats
    [p]=signrank(Y(:,1), Y(:,2));
    if p<0.1;
        lt_plot_text(x(1), max(max(Y)), num2str(p, '%3.2g'));
    end
end

lt_subplot(3,1,2); hold on;
title('Raw SEPARATION, from bidir day 1, all experiments');
xlabel([ 'day bin (binsize: ' num2str(NumDaysInBin) ')']);
ylabel('SEPARATION (hz, in targ dir)');

X=1:length(DayBins_StartEdges);
subdivisions=0.1; % for plotting diff things in one ind

for i=1:length(X);
    
    if isempty(OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC)
        continue
    end
    
    % ==== SEPARATION
    x=[i-subdivisions*2 i-subdivisions];
    
    
    Y=[OUTPUT.separation(i).PBS' OUTPUT.separation(i).MUSC'];
    
    Ymeans=mean(Y);
    Ysems=lt_sem(Y);
    
    % -- means
    if size(Y,1)>1
        lt_plot_bar(x-subdivisions/4, mean(Y), {'Color', 'k', 'Errors', Ysems});
    end
    
    % stats
    [p]=signrank(Y(:,1), Y(:,2));
    if p<0.1;
        lt_plot_text(x(1), max(max(Y)), num2str(p, '%3.2g'));
    end
end


lt_subplot(3,1,3); hold on;
title('Raw SEPARATION, from bidir day 1, all experiments');
xlabel([ 'day bin (binsize: ' num2str(NumDaysInBin) ')']);
ylabel('SEPARATION (hz, in targ dir)');

X=1:length(DayBins_StartEdges);
subdivisions=0.1; % for plotting diff things in one ind

for i=1:length(X);
    
    if isempty(OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC)
        continue
    end
    
    % ==== SEPARATION
    x=[i-subdivisions*2 i-subdivisions];
    
    
    Y=[OUTPUT.separation(i).PBS' OUTPUT.separation(i).MUSC'];
    Norm=Y(:,1);
    
    Y=Y./repmat(Norm, 1, size(Y,2));
    
    
    Ymeans=mean(Y);
    Ysems=lt_sem(Y);
    
    
    % -- means
    if size(Y,1)>1
        lt_plot_bar(x-subdivisions/4, mean(Y), {'Color', 'k', 'Errors', Ysems});
    end
    
    
    
    % stats
    [p]=signrank(Y(:,1), Y(:,2));
    if p<0.1;
        lt_plot_text(x(1), max(max(Y)), num2str(p, '%3.2g'));
    end
end



%% DOES MAGNITUDE OF CONSOLIDATION NEGATIVELY CORRELATE WITH MEASURES OF SIMILARITY?

if NumDaysInBin==3;
    DayBinForTwoDir=2; %
elseif NumDaysInBin==2;
    DayBinForTwoDir=3;
end
lt_figure; hold on;

% === COLLECT MEAN CONSOLIDATION
Separation_PBS=OUTPUT.firsttarget.FFrelBaseInBin(DayBinForTwoDir).PBS-OUTPUT.secondtarg.FFrelBaseInBin(DayBinForTwoDir).PBS;
Separation_MP=OUTPUT.firsttarget.FFrelBaseInBin(DayBinForTwoDir).MUSC-OUTPUT.secondtarg.FFrelBaseInBin(DayBinForTwoDir).MUSC;

Consolid_mean=Separation_MP./Separation_PBS;

% === COLLECT GENERALIZATION SCORE
exptnums=OUTPUT.INFORMATION(DayBinForTwoDir).experimentNum;
Generalization=[];
for i=1:length(exptnums)
    exptNum=exptnums(i);
    disp([DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(exptNum).birdname '-' DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(exptNum).exptname]);
    
    Generalization=[Generalization DATSTRUCT_FirstBidirDayIndEqualsOne.OtherStats(exptNum).generalization_OtherOverTarg];
    
end

% === COLLECT GENERALIZATION(MP) SCORE
exptnums=OUTPUT.INFORMATION(DayBinForTwoDir).experimentNum;
GeneralizationMP=[];
for i=1:length(exptnums)
    exptNum=exptnums(i);
    disp([DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(exptNum).birdname '-' DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(exptNum).exptname]);
    
    GeneralizationMP=[GeneralizationMP DATSTRUCT_FirstBidirDayIndEqualsOne.OtherStats(exptNum).generalization_MP_OtherOverTarg];
    
end


% ======= COLLECT first targ consolidation score
Consol_firsttarg=OUTPUT.firsttarget.FFrelBaseInBin(DayBinForTwoDir).MUSC./OUTPUT.firsttarget.FFrelBaseInBin(DayBinForTwoDir).PBS;


% === COLLECT ACOUSTIC DISTANCE SCORE
exptnums=OUTPUT.INFORMATION(DayBinForTwoDir).experimentNum;
Acoustic=[];
for i=1:length(exptnums)
    exptNum=exptnums(i);
    %     disp([DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(exptNum).birdname '-' DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(exptNum).exptname]);
    
    Acoustic=[Acoustic DATSTRUCT_FirstBidirDayIndEqualsOne.OtherStats(exptNum).AcDist_OtherToTarg];
    
end

% === COLLECT ACOUSTIC DISTANCE SCORE
exptnums=OUTPUT.INFORMATION(DayBinForTwoDir).experimentNum;
CorrSong=[];
for i=1:length(exptnums)
    exptNum=exptnums(i);
    %     disp([DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(exptNum).birdname '-' DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(exptNum).exptname]);
    
    CorrSong=[CorrSong DATSTRUCT_FirstBidirDayIndEqualsOne.OtherStats(exptNum).Corr_Song];
    
end


% === collect bird names
exptnums=OUTPUT.INFORMATION(DayBinForTwoDir).experimentNum;
Birdnames={};

for j=1:length(exptnums);
    exptnum=exptnums(j);
    
    birdname=DATSTRUCT_FirstBidirDayIndEqualsOne.INFORMATION(exptnum).birdname;
    Birdnames=[Birdnames birdname];
    
end

% ---- PLOT
lt_subplot(3,2,1); hold on;
xlabel('generalization');
ylabel('consolidation');

lt_regress(Consolid_mean, Generalization, 1, 0, 1, 1)

for j=1:length(Birdnames);
    lt_plot_text(Generalization(j), Consolid_mean(j), Birdnames{j}, 'k');
end

lt_subplot(3,2,2); hold on;
xlabel('Acoustic');
ylabel('consolidation');

lt_regress(Consolid_mean, Acoustic, 1, 0, 1, 1)

lt_subplot(3,2,3); hold on;
xlabel('Corr(song)');
ylabel('consolidation');

lt_regress(Consolid_mean, CorrSong, 1, 0, 1, 1)


lt_subplot(3,2,4); hold on;
xlabel('generalization (MP))');
ylabel('consolidation');

lt_regress(Consolid_mean, GeneralizationMP, 1, 0, 1, 1)

% --
lt_subplot(3,2,5); hold on;
xlabel('generalization');
ylabel('consolidation (first)');

lt_regress(Consol_firsttarg, Generalization, 1, 0, 1, 1)
for j=1:length(Birdnames);
    lt_plot_text(Generalization(j), Consol_firsttarg(j), Birdnames{j}, 'k');
end

% --
lt_subplot(3,2,6); hold on;
xlabel('generalization (MP)');
ylabel('consolidation (first)');

lt_regress(Consol_firsttarg, GeneralizationMP, 1, 0, 1, 1)
for j=1:length(Birdnames);
    lt_plot_text(GeneralizationMP(j), Consol_firsttarg(j), Birdnames{j}, 'k');
end




%% OUTPUT STATS for use in comparison with samedir analysis







