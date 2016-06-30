function lt_seq_dep_pitch_ACROSSBIRDS_TimeCourse(SeqDepPitch_AcrossBirds, Params, DayToLockAllTo)
%% LT 12/22/15

NumBirds=length(SeqDepPitch_AcrossBirds.birds);


%% === PARAMS
NumDaysInBin=Params.TimeCourse.NumDaysInBin;
AccountForEmptyDays_WNon=Params.TimeCourse.AccountForEmptyDays_WNon;
DaysToPlot=Params.TimeCourse.DaysToPlot;
PlotOnlyGaplessExpt=Params.TimeCourse.PlotOnlyGaplessExpt;
ThrowOutPoorLearner=Params.TimeCourse.ThrowOutPoorLearner;

%%
if ThrowOutPoorLearner==1;
    
filter='learning_metric';
[SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);

    
    
end

%% PLOT TIMELINE FOR ALL EXPERIMENTS
% EACH ROW IS ONE EXPERIMENT, ALL LOCKED TO WN ON.
% NOTE: SKIPPED - COLLECTING BIDIR INFORMATION - CHANGE THAT.

% lt_figure; hold on;
% BirdnamesExptname_all={};
% title('Timeline (bidir in blue)');
% xlabel('day (WN day 1 = 1)');
%
%
% CountExpt=1;
% for i=1:NumBirds;
%     birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
%     numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%
%
%     for ii=1:numexpts;
%
%         exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
%
%         targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
%
%         WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
%
%         % --- bidir start and end
%         bidir_day1=nan;
%         bidir_lastday=nan;
%         if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart_Ind'); % STOPPED HERE!!!!!!!!!!!!!!!!!!!!!
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
%
%         % -- last day
%         if ~isnan(bidir_lastday);
%             WN_last=bidir_lastday;
%         else
%             WN_last=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
%             WN_last=WN_last-WNday1+1;
%         end
%
%         % ---- first day (relative to WN on day)
%         firstday=-(WNday1-2);
%
%
%         % === plot
%         % - line for all days with any data
%         X=firstday:WN_last;
%         Y=CountExpt*ones(1,length(X));
%
%         plot(X, Y, '-ok', 'LineWidth', 1);
%
%
%
%         % -- line for start and end of bidir
%         line([bidir_day1-0.5 bidir_day1-0.5], [CountExpt-0.4 CountExpt+0.4]);
%         line([bidir_lastday+0.5 bidir_lastday+0.5], [CountExpt-0.4 CountExpt+0.4]);
%
%
% %         % -- fill in musc days
% %         X=days_with_musc_relWNday1;
% %         Y=CountExpt*ones(1, 1:length(X));
% %
% %         plot(X, Y, 'ok', 'MarkerFaceColor','k')
%
%         % --- line for end of WN
%         line([WN_last+0.5 WN_last+0.5], [CountExpt-0.4 CountExpt+0.4], 'Color','r');
%
%
%         % -- collect bird name and expt
%         BirdnamesExptname_all=[BirdnamesExptname_all [birdname(1:4) '-' exptname(end-4:end)]];
%
%         CountExpt=CountExpt+1;
%     end
% end
%
% set(gca, 'YTick', 1:CountExpt-1);
% set(gca, 'YTickLabel', BirdnamesExptname_all);
%
%
% % ---- WN on line
% line([0.5 0.5], ylim, 'Color','r');



%% ===== GET ALL DATA LOCKED TO SPECIFIC DAY

disp(' ----- ')
% disp(' BIDIR CONSOLIDATION ');
%
ExptCount=1;


DATSTRUCT=struct;

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        % === GET DAY INDS RELATIVE TO A LOCK DAY (DAY1)
        
        if strcmp(DayToLockAllTo, 'WNday1');
            
            FirstDay_AllElseLockedToThis=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            
            % get consolidation days if exist
            birdInd=find(strcmp(Params.TimeCourse.ConsolidationList, birdname));
            exptInd=find(strcmp(Params.TimeCourse.ConsolidationList, exptname));
            
            ind=intersect(birdInd+1, exptInd);
            
            if ~isempty(ind);
                consolPeriod=Params.TimeCourse.ConsolidationList{ind+1};
            else
                consolPeriod=[];
            end
            
            if AccountForEmptyDays_WNon==1;
                % then push back if no songs first few WN days
                FirstDay_AllElseLockedToThis=FirstDay_AllElseLockedToThis+SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode;
                consolPeriod=consolPeriod-SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode;
            end
            
            
        elseif strcmp(DayToLockAllTo, 'ConsolidDay1');
            
            % get consolidation days if exist
            birdInd=find(strcmp(Params.TimeCourse.ConsolidationList, birdname));
            exptInd=find(strcmp(Params.TimeCourse.ConsolidationList, exptname));
            
            ind=intersect(birdInd+1, exptInd);
            
            if ~isempty(ind);
                consolPeriod=Params.TimeCourse.ConsolidationList{ind+1};
            else
                disp(['skipped ' birdname '-' exptname ' - NO CONSOLID PERIOD ENTERED']);
                continue
            end
            
            WNOnDay=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            FirstDay_AllElseLockedToThis=WNOnDay+consolPeriod(1)-1;
            ConsolEnd_RelInds=consolPeriod(2)-consolPeriod(1)+1;
        end
        
        
        
        % ---- GET OTHER STUFF
        LastDayWithData_ActualInds=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase);
        LastDayWithData_RelInds=LastDayWithData_ActualInds-FirstDay_AllElseLockedToThis+1;
        
        targLearnDir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        
        % ++++++++++++++++++++++++++++++++++++++++++++++++++ COLLECT DATA
        % is this LMAN expt? if so, use within time window data
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1;
            datafieldFF='meanFF_DevFromBase_WithinTimeWindow';
        else
            datafieldFF='meanFF_DevFromBase';
        end
        
        % ===== FIRST TARGET
        DATSTRUCT.firsttarget(ExptCount).FFminusBase_Mean=nan(1, LastDayWithData_ActualInds-FirstDay_AllElseLockedToThis+1);
        
        % ----
        ffmeans_actualInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).(datafieldFF);
        ffmeans_RelInds=ffmeans_actualInds(FirstDay_AllElseLockedToThis:LastDayWithData_ActualInds);
        
        tmp=size(ffmeans_RelInds);
        if tmp(1)>1;
            ffmeans_RelInds=ffmeans_RelInds';
        end
        ffmeans_RelInds=ffmeans_RelInds.*targLearnDir;
        
        if strcmp(DayToLockAllTo, 'ConsolidDay1');
            ffmeans_RelInds=ffmeans_RelInds(1:ConsolEnd_RelInds);
        end
        DATSTRUCT.firsttarget(ExptCount).FFminusBase_Mean=ffmeans_RelInds;
        
        
        % ====== NONTARGETS (PUT INTO CLASSES) ============================
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        SAMETYPE=[];
        DIFFTYPE=[];
        ALL=[];
        
        DATSTRUCT.nontargets_Same(ExptCount).FFminusBase_Mean=nan(1, LastDayWithData_ActualInds-FirstDay_AllElseLockedToThis+1);
        DATSTRUCT.nontargets_Diff(ExptCount).FFminusBase_Mean=nan(1, LastDayWithData_ActualInds-FirstDay_AllElseLockedToThis+1);
        DATSTRUCT.nontargets_All(ExptCount).FFminusBase_Mean=nan(1, LastDayWithData_ActualInds-FirstDay_AllElseLockedToThis+1);
        
        for k=1:length(SylsUnique);
            syl=SylsUnique{k};
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
                continue;
            end
            
            % ====== extract ff information
            ffmeans_actualInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafieldFF);
            ffmeans_RelInds=ffmeans_actualInds(FirstDay_AllElseLockedToThis:LastDayWithData_ActualInds);
            tmp=size(ffmeans_RelInds);
            if tmp(1)>1;
                ffmeans_RelInds=ffmeans_RelInds';
            end
            ffmeans_RelInds=ffmeans_RelInds.*targLearnDir;
            if strcmp(DayToLockAllTo, 'ConsolidDay1');
                ffmeans_RelInds=ffmeans_RelInds(1:ConsolEnd_RelInds);
            end
            
            
            % --- SAME TYPE
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
                SAMETYPE=[SAMETYPE; ffmeans_RelInds];
                % --- DIFF TYPE
            elseif SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ==0;
                DIFFTYPE=[DIFFTYPE; ffmeans_RelInds];
            end
            % --- ALL
            ALL=[ALL; ffmeans_RelInds];
            
            
        end
        
        % ==== OUTPUT - nontargets =================================
        % get mean and sem
        DATSTRUCT.nontargets_Same(ExptCount).FFminusBase_Mean=nanmean(SAMETYPE, 1);
        DATSTRUCT.nontargets_Diff(ExptCount).FFminusBase_Mean=nanmean(DIFFTYPE, 1);
        DATSTRUCT.nontargets_All(ExptCount).FFminusBase_Mean=nanmean(ALL, 1);
        
        DATSTRUCT.nontargets_Same(ExptCount).FFminusBase_Mean_AllSyls=SAMETYPE;
        DATSTRUCT.nontargets_Diff(ExptCount).FFminusBase_Mean_AllSyls=DIFFTYPE;
        DATSTRUCT.nontargets_All(ExptCount).FFminusBase_Mean_AllSyls=ALL;
        
        
        % ====== SAVE INFORMATION ABOUT THIS EXPERIMENT
        
        DATSTRUCT.INFORMATION(ExptCount).birdname=birdname;
        DATSTRUCT.INFORMATION(ExptCount).exptname=exptname;
        
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart_Ind');
            BidirStartActual=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind+1;
            BidirStart_Rel=BidirStartActual-FirstDay_AllElseLockedToThis+1;
            DATSTRUCT.INFORMATION(ExptCount).BidirDay1_RelInd=BidirStart_Rel;
        end
        
        
        DATSTRUCT.INFORMATION(ExptCount).LastWNDay_RelInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd-FirstDay_AllElseLockedToThis+1;
        
        DATSTRUCT.INFORMATION(ExptCount).LastDayWithData_RelInds=LastDayWithData_RelInds;
        DATSTRUCT.INFORMATION(ExptCount).targ_learn_dir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        DATSTRUCT.INFORMATION(ExptCount).consolPeriod_RelInds=consolPeriod;
        
        if strcmp(DayToLockAllTo, 'ConsolidDay1');
            DATSTRUCT.INFORMATION(ExptCount).ConsolLastDay_RelInds=ConsolEnd_RelInds;
        end
        
        ExptCount=ExptCount+1;
        
    end
end


%% ====== PLOT [targ and other targs]

figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

for i=1:length(DATSTRUCT.firsttarget);
    
    birdname=DATSTRUCT.INFORMATION(i).birdname;
    exptname=DATSTRUCT.INFORMATION(i).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([birdname '-' exptname]);
    xlabel('day (rel WN day 1)');
    ylabel('shift from baseline (hz)');
    
    
    % == first targ
    ffMeans=DATSTRUCT.firsttarget(i).FFminusBase_Mean;
    inds=find(~isnan(ffMeans));
    
    lt_plot(inds, ffMeans(inds), {'Color', 'k', 'LineStyle', '-'});
    
    
    % ==== other nontargs
    % === SAME
    ffMeans=DATSTRUCT.nontargets_Same(i).FFminusBase_Mean;
    inds=find(~isnan(ffMeans));
    
    plot(inds, ffMeans(inds), '-ob');
    
    % === DIFF
    ffMeans=DATSTRUCT.nontargets_Diff(i).FFminusBase_Mean;
    inds=find(~isnan(ffMeans));
    
    plot(inds, ffMeans(inds), '-or');
    
    % === all
    ffMeans=DATSTRUCT.nontargets_All(i).FFminusBase_Mean;
    inds=find(~isnan(ffMeans));
    
    plot(inds, ffMeans(inds), '--ok');
    
    
    
    lt_plot_zeroline;
    
    % line for end of bidir
    if strcmp(DayToLockAllTo, 'WNday1');
        if isfield(DATSTRUCT.INFORMATION, 'BidirDay1_RelInd')
            bidir_start=DATSTRUCT.INFORMATION(i).BidirDay1_RelInd;
            if ~isempty(bidir_start)
                line([bidir_start-0.5 bidir_start-0.5], ylim, 'Color','r');
            end
        end
        %         bidir_last_day_actual=DATSTRUCT.INFORMATION(i).LastDay_RelInds_BidirIfExist_OtherwiseWN;
        %         line([bidir_last_day_actual+0.5 bidir_last_day_actual+0.5], ylim, 'Color','r');
    end
    
    % line for consol start and end
    if strcmp(DayToLockAllTo, 'ConsolidDay1');
        consolEnd=DATSTRUCT.INFORMATION(i).ConsolLastDay_RelInds;
        line([consolEnd+0.4 consolEnd+0.4], ylim, 'Color', 'm');
    else
        consolPeriod=DATSTRUCT.INFORMATION(i).consolPeriod_RelInds;
        if ~isempty(consolPeriod);
            line([consolPeriod(1)-0.4 consolPeriod(1)-0.4], ylim, 'Color', 'm');
            line([consolPeriod(2)+0.4 consolPeriod(2)+0.4], ylim, 'Color', 'm');
        end
    end
    
    if strcmp(DayToLockAllTo, 'SingleDirConsolid');
        disp('WARNING - NOT YET CODE DONE!!!');
        FXdsaVCADFEWACEW
        birdname=DATSTRUCT.INFORMATION(i).birdname;
        exptname=DATSTRUCT.INFORMATION(i).exptname;
        
        % --- need to enter consolid dates in params
        birdInd=find(strcmp(Params.LMANTimeCourse_v2.SingleDirConsolid, birdname));
        exptInd=find(strcmp(Params.LMANTimeCourse_v2.SingleDirConsolid, exptname));
        
        ind=intersect(birdInd+1, exptInd);
        
        consolPeriod=Params.LMANTimeCourse_v2.SingleDirConsolid{ind+1};
        
        % convert it to relative to consol start
        consolEnd=consolPeriod(2)-consolPeriod(1)+1;
        
        line([consolEnd+0.5 consolEnd+0.5], ylim, 'Color','m');
    end
    
    % line for WN off
    WNoffRelInd=DATSTRUCT.INFORMATION(i).LastWNDay_RelInds;
    line([WNoffRelInd+0.5 WNoffRelInd+0.5], ylim, 'Color','k');
    
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
    
    
end

lt_subtitle('black: targ; blue: sametype; red: difftype');



%% ==== PLOT ACROSS DAYS,  LOCKED TO START [ RAW] [expt means]
lt_figure; hold on;

numexpts=length(DATSTRUCT.firsttarget);

maxdays=[];
for i=1:numexpts;
    maxdays=max([maxdays length(DATSTRUCT.firsttarget(i).FFminusBase_Mean)]);
end



% === TARG
color='k';
lt_subplot(3,1,1); hold on;
title('TARGS');

Yall=nan(numexpts, maxdays); % [expt x days];
% 1) collect and plot indivs
for i=1:length(DATSTRUCT.firsttarget);
    
    % for each experiment
    Yff=DATSTRUCT.firsttarget(i).FFminusBase_Mean;
    
    % ==== does this experiment have continuous data (if desired) ====
    if ~isempty(DaysToPlot)
        if length(Yff)>=DaysToPlot
            Yff=Yff(1:DaysToPlot);
        end
    end
    
    if PlotOnlyGaplessExpt==1;
        if any(isnan(Yff)) || length(Yff)<DaysToPlot;
            continue
        end
    end
    % ==========================================================
    
    plot(1:length(Yff), Yff, '-o', 'Color',color);
    
    Yall(i,1:length(Yff))=Yff;
end
% 2) plot mean
Ymean=nanmean(Yall,1);
Ysem=lt_sem(Yall);

lt_plot_bar(1:length(Ymean), Ymean, {'Errors', Ysem, 'Color', color});


% === SAME TYPE
color='b';
lt_subplot(3,1,2); hold on;
title('SAME');

Yall=nan(numexpts, maxdays); % [expt x days];
% 1) collect and plot indivs
for i=1:length(DATSTRUCT.firsttarget);
    % for each experiment
    Yff=DATSTRUCT.nontargets_Same(i).FFminusBase_Mean;
    % ==== does this experiment have continuous data (if desired) ====
    if ~isempty(DaysToPlot)
        if length(Yff)>=DaysToPlot
            Yff=Yff(1:DaysToPlot);
        end
    end
    
    if PlotOnlyGaplessExpt==1;
        if any(isnan(Yff)) || length(Yff)<DaysToPlot;
            continue
        end
    end
    % ==========================================================
    
    plot(1:length(Yff), Yff, '-o', 'Color',color);
    
    Yall(i,1:length(Yff))=Yff;
end
% 2) plot mean
Ymean=nanmean(Yall,1);
Ysem=lt_sem(Yall);

lt_plot_bar(1:length(Ymean), Ymean, {'Errors', Ysem, 'Color', color});



% === DIFF TYPE
color='r';
lt_subplot(3,1,3); hold on;
title('DIFF');

Yall=nan(numexpts, maxdays); % [expt x days];
% 1) collect and plot indivs
for i=1:length(DATSTRUCT.firsttarget);
    % for each experiment
    Yff=DATSTRUCT.nontargets_Diff(i).FFminusBase_Mean;
    
    % ==== does this experiment have continuous data (if desired) ====
    if ~isempty(DaysToPlot)
        if length(Yff)>=DaysToPlot
            Yff=Yff(1:DaysToPlot);
        end
    end
    
    if PlotOnlyGaplessExpt==1;
        if any(isnan(Yff)) || length(Yff)<DaysToPlot;
            continue
        end
    end
    % ==========================================================
    
    plot(1:length(Yff), Yff, '-o', 'Color',color);
    
    Yall(i,1:length(Yff))=Yff;
end
% 2) plot mean
Ymean=nanmean(Yall,1);
Ysem=lt_sem(Yall);

lt_plot_bar(1:length(Ymean), Ymean, {'Errors', Ysem, 'Color', color});



%% ==== PLOT ACROSS DAYS,  LOCKED TO START [ RAW] [all rends]
if PlotOnlyGaplessExpt==1;
    hplots=[];

lt_figure; hold on;

numexpts=length(DATSTRUCT.firsttarget);

maxdays=[];
for i=1:numexpts;
    maxdays=max([maxdays length(DATSTRUCT.firsttarget(i).FFminusBase_Mean)]);
end



% === TARG
color='k';
hplot=lt_subplot(3,1,1); hold on;
hplots=[hplots hplot];
title('TARGS');

Yall=nan(numexpts, maxdays); % [expt x days];
% 1) collect and plot indivs
for i=1:length(DATSTRUCT.firsttarget);
    
    % for each experiment
    Yff=DATSTRUCT.firsttarget(i).FFminusBase_Mean;
    
    % ==== does this experiment have continuous data (if desired) ====
    if ~isempty(DaysToPlot)
        if length(Yff)>=DaysToPlot
            Yff=Yff(1:DaysToPlot);
        end
    end
    
    if PlotOnlyGaplessExpt==1;
        if any(isnan(Yff)) || length(Yff)<DaysToPlot;
            continue
        end
    end
    % ==========================================================
    
    plot(1:length(Yff), Yff, '-o', 'Color',color);
    
    Yall(i,1:length(Yff))=Yff;
end
% 2) plot mean
Ymean=nanmean(Yall,1);
Ysem=lt_sem(Yall);

lt_plot_bar(1:length(Ymean), Ymean, {'Errors', Ysem, 'Color', color});


% === SAME TYPE
color='b';
hplot=lt_subplot(3,1,2); hold on;
hplots=[hplots hplot];
title('SAME');

Yall=[]; % [expt x days];
% 1) collect and plot indivs
for i=1:length(DATSTRUCT.firsttarget);
    % for each experiment
    Yff=DATSTRUCT.nontargets_Same(i).FFminusBase_Mean_AllSyls;
    % ==== does this experiment have continuous data (if desired) ====
        if size(Yff,2)>=DaysToPlot
            Yff=Yff(:, 1:DaysToPlot);
        end
    
    if PlotOnlyGaplessExpt==1;
        if any(any(isnan(Yff))) || size(Yff,2)<DaysToPlot;
            continue
        end
    end
    % ==========================================================
    
    plot(1:size(Yff,2), Yff, '-o', 'Color',color);
    
    Yall=[Yall; Yff];
end
% 2) plot mean
Ymean=nanmean(Yall,1);
Ysem=lt_sem(Yall);

lt_plot_bar(1:length(Ymean), Ymean, {'Errors', Ysem, 'Color', color});



% === DIFF TYPE
color='r';
hplot=lt_subplot(3,1,3); hold on;
hplots=[hplots hplot];
title('DIFF');

Yall=[]; % [expt x days];
% 1) collect and plot indivs
for i=1:length(DATSTRUCT.firsttarget);
    % for each experiment
    Yff=DATSTRUCT.nontargets_Diff(i).FFminusBase_Mean_AllSyls;
    % ==== does this experiment have continuous data (if desired) ====
        if size(Yff,2)>=DaysToPlot
            Yff=Yff(:, 1:DaysToPlot);
        end
    
    if PlotOnlyGaplessExpt==1;
        if any(any(isnan(Yff))) || size(Yff,2)<DaysToPlot;
            continue
        end
    end
    % ==========================================================
    
    plot(1:size(Yff,2), Yff, '-o', 'Color',color);
    
    Yall=[Yall; Yff];
end
% 2) plot mean
Ymean=nanmean(Yall,1);
Ysem=lt_sem(Yall);

lt_plot_bar(1:length(Ymean), Ymean, {'Errors', Ysem, 'Color', color});


linkaxes(hplots, 'x');
end


%% ==== PLOT ACROSS DAYS,  LOCKED TO START [MEANS] [all rends]
if PlotOnlyGaplessExpt==1;
    hplots=[];

lt_figure; hold on;

numexpts=length(DATSTRUCT.firsttarget);

maxdays=[];
for i=1:numexpts;
    maxdays=max([maxdays length(DATSTRUCT.firsttarget(i).FFminusBase_Mean)]);
end



% === TARG
color='k';
hplot=lt_subplot(3,1,1); hold on;
hplots=[hplots hplot];
title('TARGS');
ylabel('Pitch shift (hz)');
xlabel('day (during WN)');

Yall=nan(numexpts, maxdays); % [expt x days];
% 1) collect and plot indivs
for i=1:length(DATSTRUCT.firsttarget);
    
    % for each experiment
    Yff=DATSTRUCT.firsttarget(i).FFminusBase_Mean;
    
    % ==== does this experiment have continuous data (if desired) ====
    if ~isempty(DaysToPlot)
        if length(Yff)>=DaysToPlot
            Yff=Yff(1:DaysToPlot);
        end
    end
    
    if PlotOnlyGaplessExpt==1;
        if any(isnan(Yff)) || length(Yff)<DaysToPlot;
            continue
        end
    end
    % ==========================================================
    
%     plot(1:length(Yff), Yff, '-o', 'Color',color);
%     
    Yall(i,1:length(Yff))=Yff;
end
% 2) plot mean
Ymean=nanmean(Yall,1);
Ysem=lt_sem(Yall);

lt_plot_bar(1:length(Ymean), Ymean, {'Errors', Ysem, 'Color', color});


% === SAME TYPE
color='b';
hplot=lt_subplot(3,1,2); hold on;
hplots=[hplots hplot];
title('SAME');
ylabel('Pitch shift (hz)');
xlabel('day (during WN)');

Yall=[]; % [expt x days];
% 1) collect and plot indivs
for i=1:length(DATSTRUCT.firsttarget);
    % for each experiment
    Yff=DATSTRUCT.nontargets_Same(i).FFminusBase_Mean_AllSyls;
    % ==== does this experiment have continuous data (if desired) ====
        if size(Yff,2)>=DaysToPlot
            Yff=Yff(:, 1:DaysToPlot);
        end
    
    if PlotOnlyGaplessExpt==1;
        if any(any(isnan(Yff))) || size(Yff,2)<DaysToPlot;
            continue
        end
    end
    % ==========================================================
    
%     plot(1:size(Yff,2), Yff, '-o', 'Color',color);
%     
    Yall=[Yall; Yff];
end
% 2) plot mean
Ymean=nanmean(Yall,1);
Ysem=lt_sem(Yall);

lt_plot_bar(1:length(Ymean), Ymean, {'Errors', Ysem, 'Color', color});



% === DIFF TYPE
color='r';
hplot=lt_subplot(3,1,3); hold on;
hplots=[hplots hplot];
title('DIFF');
ylabel('Pitch shift (hz)');
xlabel('day (during WN)');

Yall=[]; % [expt x days];
% 1) collect and plot indivs
for i=1:length(DATSTRUCT.firsttarget);
    % for each experiment
    Yff=DATSTRUCT.nontargets_Diff(i).FFminusBase_Mean_AllSyls;
    % ==== does this experiment have continuous data (if desired) ====
        if size(Yff,2)>=DaysToPlot
            Yff=Yff(:, 1:DaysToPlot);
        end
    
    if PlotOnlyGaplessExpt==1;
        if any(any(isnan(Yff))) || size(Yff,2)<DaysToPlot;
            continue
        end
    end
    % ==========================================================
    
%     plot(1:size(Yff,2), Yff, '-o', 'Color',color);
%     
    Yall=[Yall; Yff];
end
% 2) plot mean
Ymean=nanmean(Yall,1);
Ysem=lt_sem(Yall);

lt_plot_bar(1:length(Ymean), Ymean, {'Errors', Ysem, 'Color', color});
% ---- FOR EACH DAY, IS IT DIFF FROM 0?
numdays=size(Yall, 2);
for j=1:numdays
    [h, p]=ttest(Yall(:,j), 0, 'Tail', 'left');
%     p=signrank(Yall(:,j));
    disp(p);
    if p<0.05 
        lt_plot_text(j, max([0, nanmean(Yall(:,j))]), '*', 'm');
    end
end





linkaxes(hplots, 'x');
end




%% +++++++++++++++++++++++++++++++++++++++++++++++++++++=
%% ==== PLOT MEAN ACROSS ALL EXPERIMENTS [BINNED]
% ==== DON'T CARE ABOUT CONSOLIDATION - SIMPLY BIN IN DAYS RELATIVE TO
% LEARNING START

% === to add:
% 1) stop when reach end of bidir
% 3) only run if bird has data for all bins

maxday=[];
for i=1:length(DATSTRUCT.firsttarget);
    maxday=max([maxday length(DATSTRUCT.firsttarget(i).FFminusBase_Mean)]);
end
DayBins_StartEdges=1:NumDaysInBin:maxday;

OUTPUT=struct; % one cell for each bin

for i=1:length(DayBins_StartEdges);
    
    OUTPUT.firsttarget(i).FFrelBaseInBin=[];
    OUTPUT.nontargets_Same(i).FFrelBaseInBin=[];
    OUTPUT.nontargets_Diff(i).FFrelBaseInBin=[];
    OUTPUT.nontargets_All(i).FFrelBaseInBin=[];
    
    DaysInBin=DayBins_StartEdges(i):DayBins_StartEdges(i)+NumDaysInBin-1;
    
    for day=DaysInBin;
        
        % == for this day, collect all MUSC. for days with musc, also
        % collect PBS
        
        for j=1:length(DATSTRUCT.firsttarget); % num expts
            
            
            
            % ++++++++++++++++++++++++++++++++++++++
            % === only continue if vector is long enough to today
            if length(DATSTRUCT.firsttarget(j).FFminusBase_Mean)<day;
                continue;
            end
            
            if isnan(DATSTRUCT.firsttarget(j).FFminusBase_Mean(day));
                continue;
            end
            
            % ==== make sure this day is in consolidation period (if
            % desired)
            %             if strcmp(DayToLockAllTo, 'SingleDirConsolid');
            %             birdname=DATSTRUCT.INFORMATION(j).birdname;
            %             exptname=DATSTRUCT.INFORMATION(j).exptname;
            %
            %                 % --- need to enter consolid dates in params
            %             birdInd=find(strcmp(Params.LMANTimeCourse_v2.SingleDirConsolid, birdname));
            %             exptInd=find(strcmp(Params.LMANTimeCourse_v2.SingleDirConsolid, exptname));
            %
            %             ind=intersect(birdInd+1, exptInd);
            %
            %                 consolPeriod=Params.LMANTimeCourse_v2.SingleDirConsolid{ind+1};
            %
            %                 % convert it to relative to consol start
            %                 consolEnd=consolPeriod(2)-consolPeriod(1)+1;
            %                     if day>consolEnd
            %                         continue
            %                     end
            %             end
            
            
            
            % ============== first target
            ffMeans=DATSTRUCT.firsttarget(j).FFminusBase_Mean(day);
            % ---- OUTPUT
            OUTPUT.firsttarget(i).FFrelBaseInBin=[OUTPUT.firsttarget(i).FFrelBaseInBin ffMeans];
            
            
            % =========== OTHER TARGS -------------------------------------
            if ~isempty(DATSTRUCT.nontargets_Same(j).FFminusBase_Mean)
                % === same type
                ffMeans=DATSTRUCT.nontargets_Same(j).FFminusBase_Mean(day);
                % ---- OUTPUT
                OUTPUT.nontargets_Same(i).FFrelBaseInBin=[OUTPUT.nontargets_Same(i).FFrelBaseInBin ffMeans];
            end
            
            % === diff type
            if ~isempty(DATSTRUCT.nontargets_Diff(j).FFminusBase_Mean)
                ffMeans=DATSTRUCT.nontargets_Diff(j).FFminusBase_Mean(day);
                % ---- OUTPUT
                OUTPUT.nontargets_Diff(i).FFrelBaseInBin=[OUTPUT.nontargets_Diff(i).FFrelBaseInBin ffMeans];
            end
            
            % === all type
            if ~isempty(DATSTRUCT.nontargets_All(j).FFminusBase_Mean)
            ffMeans=DATSTRUCT.nontargets_All(j).FFminusBase_Mean(day);
            % ---- OUTPUT
            OUTPUT.nontargets_All(i).FFrelBaseInBin=[OUTPUT.nontargets_All(i).FFrelBaseInBin ffMeans];
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
%     if isempty(OUTPUT.firsttarget(i).FFrelBaseInBin.MPbias)
%     continue
%     end
%       %=== first targ
%       % -- plot each point
%       y=OUTPUT.firsttarget(i).FFrelBaseInBin.MPbias;
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
%
%
%
%% PLOT (RAW PBS AND MUSC, for targ and new targ);
lt_figure; hold on;

lt_subplot(2,1,1); hold on;
title('Raw,  all experiments');
xlabel([ 'day bin (binsize: ' num2str(NumDaysInBin) ')']);
ylabel('FF rel baseline (in targ dir)');

X=1:length(DayBins_StartEdges);
subdivisions=0.1; % for plotting diff things in one ind

for i=1:length(X);
    
    %=== first targ
    x=[i-subdivisions*3];
    % -- plot each point
    YmeansInBin=OUTPUT.firsttarget(i).FFrelBaseInBin;
    
    if isempty(YmeansInBin);
        continue
    end
    
    plot(x, YmeansInBin, 'ok');
    
    % -- means
    if length(YmeansInBin)>1
        lt_plot_bar(x-subdivisions/4, mean(YmeansInBin), {'Color', 'k', 'BarWidth', 0.1});
    end
    
    
    % ====== NON TARGS --------------------
    % == same type
    x=[i-subdivisions];
    % -- plot each point
    YmeansInBin=OUTPUT.nontargets_Same(i).FFrelBaseInBin;
    if isempty(YmeansInBin);
        continue
    end
    plot(x, YmeansInBin, 'ob');
    
    % -- means
    if length(YmeansInBin)>1
        lt_plot_bar(x-subdivisions/4, mean(YmeansInBin), {'Color', 'b', 'BarWidth', 0.1});
    end
    
    
    %       % == diff type
    x=[i];
    % -- plot each point
    YmeansInBin=OUTPUT.nontargets_Diff(i).FFrelBaseInBin;
     if isempty(YmeansInBin);
        continue
    end
   
    plot(x, YmeansInBin, 'or');
    
    % -- means
    if length(YmeansInBin)>1
        lt_plot_bar(x-subdivisions/4, mean(YmeansInBin), {'Color', 'r', 'BarWidth', 0.1});
    end
    
    
    %         % === all
    x=i+subdivisions;
    % -- plot each point
    YmeansInBin=OUTPUT.nontargets_All(i).FFrelBaseInBin;
     if isempty(YmeansInBin);
        continue
    end
   
    plot(x, YmeansInBin, 'ok');
    
    % -- means
    if length(YmeansInBin)>1
        lt_plot_bar(x-subdivisions/4, mean(YmeansInBin), {'Color', 'k', 'BarWidth', 0.1});
    end
    
    
    
    lt_plot_annotation('', {'bk=targ','bl=same','red=diff','mag=all'},'k');
    
    
    %       % === second targ
    %       x=[i+subdivisions i+2*subdivisions];
    %
    %       % -- plot each point
    %       YmeansInBin=OUTPUT.secondtarg.FFrelBaseInBin(i).PBS;
    %       yMUSC=OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC;
    %
    %       plot(x, [YmeansInBin' yMUSC'], '-ob');
    %
    %       % -- means
    %       if size([YmeansInBin' yMUSC'],1)>1
    % %       lt_plot_bar(x+0.05, mean([YmeansInBin' yMUSC']), {'Color', 'b', 'Errors', lt_sem([YmeansInBin' yMUSC'])});
    %       lt_plot_bar(x+subdivisions/4, mean([YmeansInBin' yMUSC']), {'Color', 'b'});
    %       end
    
end


% ========= [ only mean, each mean norm to PBS]
% lt_subplot(2,1,2); hold on;
% title('Raw, all experiments [mean, norm]');
% xlabel([ 'day bin (binsize: ' num2str(NumDaysInBin) ')']);
% ylabel('FF rel baseline (in targ dir)');
%
% for i=1:length(X);
%
%     if isempty(OUTPUT.firsttarget(i).FFrelBaseInBin)
%         continue
%     end
%
%       %=== first targ
%       x=[i-subdivisions*2];
%       % -- plot each point
%       YmeansInBin=OUTPUT.firsttarget(i).FFrelBaseInBin;
%
%       % -- means
%       if length(YmeansInBin)>1
%       lt_plot_bar(x-0.05, mean(YmeansInBin)./mean([YmeansInBin' YmeansInBin']), {'Color', 'k', 'Errors', lt_sem([YmeansInBin' yMUSC'])./mean([YmeansInBin' YmeansInBin'])});
%       end
%
%
% %       % === second targ
% %       x=[i+subdivisions i+2*subdivisions];
% %
% %       % -- plot each point
% %       YmeansInBin=OUTPUT.secondtarg.FFrelBaseInBin(i).PBS;
% %       yMUSC=OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC;
% %
% %       % -- means
% %       if size([YmeansInBin' yMUSC'],1)>1
% %       lt_plot_bar(x-0.05, -mean([YmeansInBin' yMUSC'])./mean([YmeansInBin' YmeansInBin']), {'Color', 'b', 'Errors', lt_sem([YmeansInBin' yMUSC'])./mean([YmeansInBin' YmeansInBin'])});
% %       end
%
% end





%% ==================== COLLECT, BUT HAVE FIRST WN DAY BE LOCKED TO DAY 1 (INCLUDE BIDIR, BUT ANNOTATE AS SUCH)
%


