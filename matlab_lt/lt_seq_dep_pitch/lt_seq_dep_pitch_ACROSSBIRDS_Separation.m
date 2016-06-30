function [SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_Separation(SeqDepPitch_AcrossBirds, PARAMS, OnlyAnalyzeIfHasAllDaysData, DaysToPlot, ThrowOutIfDriveMore,ThrowOut_pu11_exptwithbadday7)
%% LT 9/21/15 - LOOK AT NON-TARGETS AND WHETHER GENERALIZATION IS GREATER EARLY ON VS. LATER ON.

% E.g.:
% OnlyAnalyzeIfHasAllDaysData=1; % if 1, then throws out expts with nans or not enough days.
% DaysToPlot=1:15; % i.e. day 1 is WN day 1
% ThrowOutIfDriveMore=1; % throws out expts where drive more pitch change
% during the days
% ThrowOut_pu11_exptwithbadday7 =1; throws out pu11wh87','SeqDepPitchLMAN2'
    % ONLY WORKS if ThrowOutIfDriveMore=1;
    
%% Params

NumBirds=length(SeqDepPitch_AcrossBirds.birds);
NumExperiments = length(SeqDepPitch_AcrossBirds.AllExpt);

if ThrowOut_pu11_exptwithbadday7==1
DriveMoreList={'pu64bk13', 'SeqDepPitchShift2', 12, ...
    'gr41gr90', 'SeqDepPitchShift', 13, ...
    'gr41gr90', 'SeqDepPitchShift2', 17, ...
    'pu11wh87','SeqDepPitchShift2', 11,...
    'pu11wh87','SeqDepPitchLMAN2', 1, ...
    'pu53wh88', 'SeqDepPitchShift3',13}; % i.e. value triplets {birdname, exptname, LastDayBeforeDriveMore_WNdayInds}
else
 DriveMoreList={'pu64bk13', 'SeqDepPitchShift2', 12, ...
    'gr41gr90', 'SeqDepPitchShift', 13, ...
    'gr41gr90', 'SeqDepPitchShift2', 17, ...
    'pu11wh87','SeqDepPitchShift2', 11,...
    'pu53wh88', 'SeqDepPitchShift3',13}; % i.e. value triplets {birdname, exptname, LastDayBeforeDriveMore_WNdayInds}
end

%% ==== EACH EXPERIMENT ONE DATAPOINT (NEW METHOD USING PREEXTRACTED EXPERIMENT DATA)
count=1;
SubplotsPerFig=8;
subplotrows=2;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];

% to collect for mean
clear STRUCT_OUT;
STRUCT_OUT.Targets=[];
STRUCT_OUT.Nontargets=[];
STRUCT_OUT.Similar=[];
STRUCT_OUT.Diff=[];

% === PLOT AND COLLECT DATA
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('each datapoint is one experiment');

for i=1:NumExperiments;
    birdnum=SeqDepPitch_AcrossBirds.AllExpt(i).birdnum;
    exptnum=SeqDepPitch_AcrossBirds.AllExpt(i).exptnum;
    
    birdname=SeqDepPitch_AcrossBirds.birds{birdnum}.birdname;
    exptname=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.ExptID;
    
    tmp=1; % start passing
    
    % Check 1) FIGURE OUT IF KEEP THIS EXPERIMENTS BASED ON DURATION REQUIREMENTS
    if OnlyAnalyzeIfHasAllDaysData==1;
        % throw out if not enough WN days
        if length(SeqDepPitch_AcrossBirds.AllExpt(i).Target.DayVals_DurWN_zscore_targdirsign{1}) < max(DaysToPlot);
            tmp=0;
        elseif any(isnan(SeqDepPitch_AcrossBirds.AllExpt(i).Target.DayVals_DurWN_zscore_targdirsign{1}(DaysToPlot)))
            tmp=0;
        end
    end
    
    % Check 2) goes into multidir?
    if isfield(SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.DATES, 'MultiDir_OneDayBeforeStart_Ind');
        OneDayBeforeMultidirStart=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.DATES.MultiDir_OneDayBeforeStart_Ind;
        NumBaselineDays=max(SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Data_PlotLearning.Params.SeqFilter.BaselineDays);
        NumDaysNeeded=max(DaysToPlot);
        
        if NumBaselineDays+NumDaysNeeded>OneDayBeforeMultidirStart;
            % then bad, not enough single dir days
            tmp=0;
            disp(['experiment ' num2str(i) ' miltidir days started too soon']);
        end
        
    end
    
         % Check 3) goes into samedir?
    if isfield(SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.DATES, 'SameDir_OneDayBeforeStart_Ind');
        OneDayBeforeMultidirStart=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.DATES.SameDir_OneDayBeforeStart_Ind;
        NumBaselineDays=max(SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Data_PlotLearning.Params.SeqFilter.BaselineDays);
        NumDaysNeeded=max(DaysToPlot);
        
        if NumBaselineDays+NumDaysNeeded>OneDayBeforeMultidirStart;
            % then bad, not enough single dir days
            tmp=0;
            disp(['experiment ' birdname '-' exptname ' miltidir days started too soon']);
        end
        
    end

    
    % 2) Throw out specific experiments that "drove more" pitch
    % change during the days window.
    if ThrowOutIfDriveMore==1;
        birdname=SeqDepPitch_AcrossBirds.birds{birdnum}.birdname;
        exptname=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.ExptID;
        
        % check to see whether this bird/expt is on the list
        inds1=find(strcmp(DriveMoreList, birdname));
        for kk=1:length(inds1);
            if strcmp(DriveMoreList{inds1(kk)+1},exptname);
                % THEN CHECK TO SEE IF DAYS IS PAST START OF DRIVE
                % MORE
                LastWNDayNeeded=max(DaysToPlot);
                LastDayBeforeDriveMore_WNdayInds=DriveMoreList{inds1(kk)+2};
                
                if LastWNDayNeeded>LastDayBeforeDriveMore_WNdayInds;
                    disp(['SKIPPING ' birdname '-' exptname ' because starts "driving more" before the last desired day to plot']);
                    tmp=0;
                end
            end
        end
    end
    
    if tmp==0;
        disp(['threw out expt: ' birdname '-' exptname]);
        continue;
        
    end
    
    
    
    % ====== TARGET
    % -- INITIATE THE array for all days = [-1 1 2 ...]
    Data_days=nan(1,max(DaysToPlot)+1);
    % --- COLLECT DATA
    Data_days(1)=SeqDepPitch_AcrossBirds.AllExpt(i).Target.DayVals_Baseline_zscore_targdirsign{1}(end); % baseline
    if length(SeqDepPitch_AcrossBirds.AllExpt(i).Target.DayVals_DurWN_zscore_targdirsign{1})<max(DaysToPlot); % not enough days, take max days
        tmpind=length(SeqDepPitch_AcrossBirds.AllExpt(i).Target.DayVals_DurWN_zscore_targdirsign{1});
        Data_days(DaysToPlot(1)+1:tmpind+1)=SeqDepPitch_AcrossBirds.AllExpt(i).Target.DayVals_DurWN_zscore_targdirsign{1}(DaysToPlot(1):tmpind);
    else
        Data_days(DaysToPlot+1)=SeqDepPitch_AcrossBirds.AllExpt(i).Target.DayVals_DurWN_zscore_targdirsign{1}(DaysToPlot); % WN
    end
    
    % -- plot
    X=[-1 DaysToPlot];
    plot(X, Data_days, '-k');
    
    % -- save, to plot mean across experiments
    STRUCT_OUT.Targets=[STRUCT_OUT.Targets; Data_days];
    
    
    % ====== NONTARGETS (1 mean value for each experiments);
    % -- INITIATE THE array for all days = [-1 1 2 ...]
    Data_days=nan(1,max(DaysToPlot)+1);
    % --- Baseline
    func=@(x) x{1}(end);
    mean_val_today=mean(arrayfun(func, SeqDepPitch_AcrossBirds.AllExpt(i).Nontargets.DayVals_Baseline_zscore_targdirsign));
    
    Data_days(1)=mean_val_today;
    
    % --- WN days
    for j=1:length(DaysToPlot);
        daynum=DaysToPlot(j);
        
        func=@(x) x{1}(daynum);
        mean_val_today=mean(arrayfun(func, SeqDepPitch_AcrossBirds.AllExpt(i).Nontargets.DayVals_DurWN_zscore_targdirsign));
        
        Data_days(j+1)=mean_val_today;
    end
    
    % --- plot
    X=[-1 DaysToPlot];
    plot(X, Data_days, '-','Color',[0.6 0.6 0.6]);
    
    % --- save
    STRUCT_OUT.Nontargets=[STRUCT_OUT.Nontargets; Data_days];
    
    
    % ======= SIMILAR
    % -- INITIATE THE array for all days = [-1 1 2 ...]
    Data_days=nan(1,max(DaysToPlot)+1);
    % --- Baseline
    func=@(x) x{1}(end);
    mean_val_today=mean(arrayfun(func, SeqDepPitch_AcrossBirds.AllExpt(i).Similar.DayVals_Baseline_zscore_targdirsign));
    
    Data_days(1)=mean_val_today;
    
    % --- WN days
    for j=1:length(DaysToPlot);
        daynum=DaysToPlot(j);
        
        func=@(x) x{1}(daynum);
        mean_val_today=mean(arrayfun(func, SeqDepPitch_AcrossBirds.AllExpt(i).Similar.DayVals_DurWN_zscore_targdirsign));
        
        Data_days(j+1)=mean_val_today;
    end
    
    % --- plot
    X=[-1 DaysToPlot];
    plot(X, Data_days, '-', 'Color',[0.2 0.2 0.8]);
    
    % --- save
    STRUCT_OUT.Similar=[STRUCT_OUT.Similar; Data_days];
    
    
    % ========== DIFFERENT
    % -- INITIATE THE array for all days = [-1 1 2 ...]
    Data_days=nan(1,max(DaysToPlot)+1);
    % --- Baseline
    func=@(x) x{1}(end);
    mean_val_today=mean(arrayfun(func, SeqDepPitch_AcrossBirds.AllExpt(i).Different.DayVals_Baseline_zscore_targdirsign));
    
    Data_days(1)=mean_val_today;
    
    % --- WN days
    for j=1:length(DaysToPlot);
        daynum=DaysToPlot(j);
        
        func=@(x) x{1}(daynum);
        mean_val_today=mean(arrayfun(func, SeqDepPitch_AcrossBirds.AllExpt(i).Different.DayVals_DurWN_zscore_targdirsign));
        
        Data_days(j+1)=mean_val_today;
    end
    
    % --- plot
    X=[-1 DaysToPlot];
    plot(X, Data_days, '-r', 'Color',[0.8 0.2 0.2]);
    
    % --- save
    STRUCT_OUT.Diff=[STRUCT_OUT.Diff; Data_days];
end
grid on


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ==== PLOT ONLY MEANS
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('means (across experiments)');

% -- Targets
Y=nanmean(STRUCT_OUT.Targets,1);
Ysem=lt_sem(STRUCT_OUT.Targets);

errorbar(X, Y, Ysem, 'o-k','MarkerFaceColor','k','MarkerSize',7)

% -- Nontarg
Y=nanmean(STRUCT_OUT.Nontargets,1);
Ysem=lt_sem(STRUCT_OUT.Nontargets);

errorbar(X, Y, Ysem, 'o-k','MarkerFaceColor','w','MarkerSize',7)

% -- Similar
Y=nanmean(STRUCT_OUT.Similar,1);
Ysem=lt_sem(STRUCT_OUT.Similar);

errorbar(X, Y, Ysem, 'o-k','MarkerFaceColor','b','MarkerSize',7)

% -- Diff
Y=nanmean(STRUCT_OUT.Diff,1);
Ysem=lt_sem(STRUCT_OUT.Diff);

errorbar(X, Y, Ysem, 'o-k','MarkerFaceColor','r','MarkerSize',7)

grid on




% ===== ZOOM INTO NONTARGETS (raw data)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('all nontargets');

% -- Nontarg
Y=STRUCT_OUT.Nontargets;

plot(X, Y, 'Color',[0.6 0.6 0.6])
%     errorbar(X, Y, Ysem, 'o-k','MarkerFaceColor','w','MarkerSize',7)

grid on

% ===== ZOOM INTO SIMILAR (raw data)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('similar ');

% -- Nontarg
Y=STRUCT_OUT.Similar;

plot(X, Y, 'Color',[0.2 0.2 0.8])
%     errorbar(X, Y, Ysem, 'o-k','MarkerFaceColor','w','MarkerSize',7)

grid on


% ===== ZOOM INTO DIFFERENT (raw data)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('different');

% -- Nontarg
Y=STRUCT_OUT.Diff;

plot(X, Y, 'Color',[0.8 0.2 0.2])
%     errorbar(X, Y, Ysem, 'o-k','MarkerFaceColor','w','MarkerSize',7)

grid on





% ====== ZOOM INTO THE NONTARGETS (means)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('zooming in');

% -- Nontarg
Y=nanmean(STRUCT_OUT.Nontargets,1);
Ysem=lt_sem(STRUCT_OUT.Nontargets);

errorbar(X, Y, Ysem, 'o-k','MarkerFaceColor','w','MarkerSize',7)

% -- Similar
Y=nanmean(STRUCT_OUT.Similar,1);
Ysem=lt_sem(STRUCT_OUT.Similar);

errorbar(X, Y, Ysem, 'o-k','MarkerFaceColor','b','MarkerSize',7)

% -- Diff
Y=nanmean(STRUCT_OUT.Diff,1);
Ysem=lt_sem(STRUCT_OUT.Diff);

errorbar(X, Y, Ysem, 'o-k','MarkerFaceColor','r','MarkerSize',7)

grid on




% ======== PLOT RATIO RELATIVE TO THE TARGETS
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('divided by target learning');
xlabel('day post WN on');
Targ_mean=nanmean(STRUCT_OUT.Targets,1);
Targ_mean=repmat(Targ_mean, size(STRUCT_OUT.Targets,1),1); % one row for each expt

% -- Nontarg
tmp=STRUCT_OUT.Nontargets./Targ_mean;
Y=nanmean(tmp,1);
Ysem=lt_sem(tmp);

errorbar(X(2:end), Y(2:end), Ysem(2:end), 'o-k','MarkerFaceColor','w','MarkerSize',7)

% -- Similar
tmp=STRUCT_OUT.Similar./Targ_mean;
Y=nanmean(tmp,1);
Ysem=lt_sem(tmp);

errorbar(X(2:end), Y(2:end), Ysem(2:end), 'o-b','MarkerFaceColor','w','MarkerSize',7)

% -- Diff
tmp=STRUCT_OUT.Diff./Targ_mean;
Y=nanmean(tmp,1);
Ysem=lt_sem(tmp);

errorbar(X(2:end), Y(2:end), Ysem(2:end), 'o-r','MarkerFaceColor','w','MarkerSize',7)

grid on

% ======== PLOT DIFFERENCE FROM YESTERDAY
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('difference rel to yesterday');
xlabel('today (day post WN on)');
ylabel('differnece from yesterday (diff of z-score)');
grid on;

% -- Targets
Ydiff=diff(STRUCT_OUT.Targets,1,2);

Ymean=nanmean(Ydiff,1);
Ysem=lt_sem(Ydiff);

errorbar(X(2:end), Ymean, Ysem, '-ok','Color','k','MarkerSize',7);

lt_plot_zeroline

% -- Nontarg
Ydiff=diff(STRUCT_OUT.Nontargets,1,2);

Ymean=nanmean(Ydiff,1);
Ysem=lt_sem(Ydiff);

errorbar(X(2:end), Ymean, Ysem, '-o','Color','k','MarkerSize',7);

lt_plot_zeroline

% -- Sim
Ydiff=diff(STRUCT_OUT.Similar,1,2);

Ymean=nanmean(Ydiff,1);
Ysem=lt_sem(Ydiff);

errorbar(X(2:end), Ymean, Ysem, '-o','Color','b','MarkerSize',7);

lt_plot_zeroline

% -- Diff
Ydiff=diff(STRUCT_OUT.Diff,1,2);

Ymean=nanmean(Ydiff,1);
Ysem=lt_sem(Ydiff);

errorbar(X(2:end), Ymean, Ysem, '-o','Color','r','MarkerSize',7);

lt_plot_zeroline



% ===== DIVIDE PITCH DIFFERENCE BY THAT OF TARGET
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('difference rel to yesterday (divide target diff)');
xlabel('today (day post WN on)');
ylabel('diff: norm to global targ mean')
grid on;

% tareget
Y=nanmean(STRUCT_OUT.Targets,1);
Ytarg=diff(Y);
Ytarg=repmat(Ytarg,size(STRUCT_OUT.Targets,1) ,1);

% -- Nontarg
Ydiff=diff(STRUCT_OUT.Nontargets,1,2);
Ydiff=Ydiff./Ytarg;

Ymean=nanmean(Ydiff,1);
Ysem=lt_sem(Ydiff);

errorbar(X(2:end), Ymean, Ysem, '-o','Color','k','MarkerSize',7);

lt_plot_zeroline

% -- Similar
Ydiff=diff(STRUCT_OUT.Similar,1,2);
Ydiff=Ydiff./Ytarg;

Ymean=nanmean(Ydiff,1);
Ysem=lt_sem(Ydiff);

errorbar(X(2:end), Ymean, Ysem, '-o','Color','b','MarkerSize',7);

lt_plot_zeroline

% -- Diff
Ydiff=diff(STRUCT_OUT.Diff,1,2);
Ydiff=Ydiff./Ytarg;

Ymean=nanmean(Ydiff,1);
Ysem=lt_sem(Ydiff);

errorbar(X(2:end), Ymean, Ysem, '-o','Color','r','MarkerSize',7);

lt_plot_zeroline


% ===== DIVIDE PITCH DIFFERENCE BY THAT OF TARGET [zooming into
% region where targ is learning]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('difference rel to yesterday (divide target diff)');
xlabel('today (day post WN on)');
ylabel('diff: norm to global targ mean')
grid on;

% tareget
Y=nanmean(STRUCT_OUT.Targets,1);
Ytarg=diff(Y);
Ytarg=repmat(Ytarg,size(STRUCT_OUT.Targets,1) ,1);

% -- Nontarg
Ydiff=diff(STRUCT_OUT.Nontargets,1,2);
Ydiff=Ydiff./Ytarg;

Ymean=nanmean(Ydiff,1);
Ysem=lt_sem(Ydiff);

errorbar(X(2:end), Ymean, Ysem, '-o','Color','k','MarkerSize',7);

lt_plot_zeroline

% -- Similar
Ydiff=diff(STRUCT_OUT.Similar,1,2);
Ydiff=Ydiff./Ytarg;

Ymean=nanmean(Ydiff,1);
Ysem=lt_sem(Ydiff);

errorbar(X(2:end), Ymean, Ysem, '-o','Color','b','MarkerSize',7);

lt_plot_zeroline

% -- Diff
Ydiff=diff(STRUCT_OUT.Diff,1,2);
Ydiff=Ydiff./Ytarg;

Ymean=nanmean(Ydiff,1);
Ysem=lt_sem(Ydiff);

errorbar(X(2:end), Ymean, Ysem, '-o','Color','r','MarkerSize',7);

lt_plot_zeroline

xlim([0 5]);
ylim([-1 1]);




% ===== DIVIDE PITCH DIFFERENCE BY THAT OF TARGET
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('difference rel to yesterday (divide target diff)');
xlabel('today (day post WN on)');
ylabel('diff: norm to within expt targ')
grid on;

% --- tareget
Ytarg=diff(STRUCT_OUT.Targets,1,2);

% -- Nontarg
Ydiff=diff(STRUCT_OUT.Nontargets,1,2);
Ydiff=Ydiff./Ytarg;

Ymean=nanmean(Ydiff,1);
Ysem=lt_sem(Ydiff);

errorbar(X(2:end), Ymean, Ysem, '-o','Color','k','MarkerSize',7);

lt_plot_zeroline

% -- Similar
Ydiff=diff(STRUCT_OUT.Similar,1,2);
Ydiff=Ydiff./Ytarg;

Ymean=nanmean(Ydiff,1);
Ysem=lt_sem(Ydiff);

errorbar(X(2:end), Ymean, Ysem, '-o','Color','b','MarkerSize',7);

lt_plot_zeroline

% -- Diff
Ydiff=diff(STRUCT_OUT.Diff,1,2);
Ydiff=Ydiff./Ytarg;

Ymean=nanmean(Ydiff,1);
Ysem=lt_sem(Ydiff);

errorbar(X(2:end), Ymean, Ysem, '-o','Color','r','MarkerSize',7);

lt_plot_zeroline



% ===== PLOT EACH TARGET ON ITS OWN, with experiment label
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('each expt"s target');
grid on;
tmp=[SeqDepPitch_AcrossBirds.AllExpt(:).Target];
for j=1:length(tmp);
    try
        plot(tmp(j).DayVals_DurWN_zscore_targdirsign{1}(DaysToPlot), '-','Color',[rand rand rand]);
        birdname=SeqDepPitch_AcrossBirds.birds{tmp(j).BirdNum_all}.birdname;
        exptname=SeqDepPitch_AcrossBirds.birds{tmp(j).BirdNum_all}.experiment{tmp(j).ExptNum_all}.ExptID;
        
        lt_plot_text(DaysToPlot(end)+0.2, tmp(j).DayVals_DurWN_zscore_targdirsign{1}(DaysToPlot(end)), [birdname '-' exptname])
    catch err
        
    end
    
end


lt_subtitle('Each datapoint is mean within one experiment');




%% ====== regressions (nontarg to targ)
clear STRUCT_OUT2;
STRUCT_OUT2.Targets=[];
STRUCT_OUT2.Nontargets=[];
STRUCT_OUT2.Similar=[];
STRUCT_OUT2.Diff=[];

% === COLLECT DATA
ExptsToKeep=[];

for i=1:NumExperiments;
    birdnum=SeqDepPitch_AcrossBirds.AllExpt(i).birdnum;
    exptnum=SeqDepPitch_AcrossBirds.AllExpt(i).exptnum;
    
    birdname=SeqDepPitch_AcrossBirds.birds{birdnum}.birdname;
    exptname=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.ExptID;
    
    tmp=1; % start passing
    
    % Check 1) FIGURE OUT IF KEEP THIS EXPERIMENTS BASED ON DURATION REQUIREMENTS
    if OnlyAnalyzeIfHasAllDaysData==1;
        % throw out if not enough WN days
        if length(SeqDepPitch_AcrossBirds.AllExpt(i).Target.DayVals_DurWN_zscore_targdirsign{1}) < max(DaysToPlot);
            disp(['not enough WN days ' birdname '-' exptname]);
            tmp=0;
        elseif any(isnan(SeqDepPitch_AcrossBirds.AllExpt(i).Target.DayVals_DurWN_zscore_targdirsign{1}(DaysToPlot)))
            disp(['some WN days nan ' birdname '-' exptname]);
            tmp=0;
        end
    end
    
    % Check 2) goes into multidir?
    if isfield(SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.DATES, 'MultiDir_OneDayBeforeStart_Ind');
        OneDayBeforeMultidirStart=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.DATES.MultiDir_OneDayBeforeStart_Ind;
        NumBaselineDays=max(SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Data_PlotLearning.Params.SeqFilter.BaselineDays);
        NumDaysNeeded=max(DaysToPlot);
        
        if NumBaselineDays+NumDaysNeeded>OneDayBeforeMultidirStart;
            % then bad, not enough single dir days
            tmp=0;
            disp(['experiment ' birdname '-' exptname ' miltidir days started too soon']);
        end
        
    end
    
     % Check 3) goes into samedir?
    if isfield(SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.DATES, 'SameDir_OneDayBeforeStart_Ind');
        OneDayBeforeMultidirStart=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.DATES.SameDir_OneDayBeforeStart_Ind;
        NumBaselineDays=max(SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Data_PlotLearning.Params.SeqFilter.BaselineDays);
        NumDaysNeeded=max(DaysToPlot);
        
        if NumBaselineDays+NumDaysNeeded>OneDayBeforeMultidirStart;
            % then bad, not enough single dir days
            tmp=0;
            disp(['experiment ' birdname '-' exptname ' miltidir days started too soon']);
        end
        
    end
   
    
    % 2) Throw out specific experiments that "drove more" pitch
    % change during the days window.
    if ThrowOutIfDriveMore==1;
        birdname=SeqDepPitch_AcrossBirds.birds{birdnum}.birdname;
        exptname=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.ExptID;
        
        % check to see whether this bird/expt is on the list
        inds1=find(strcmp(DriveMoreList, birdname));
        for kk=1:length(inds1);
            if strcmp(DriveMoreList{inds1(kk)+1},exptname);
                % THEN CHECK TO SEE IF DAYS IS PAST START OF DRIVE
                % MORE
                LastWNDayNeeded=max(DaysToPlot);
                LastDayBeforeDriveMore_WNdayInds=DriveMoreList{inds1(kk)+2};
                
                if LastWNDayNeeded>LastDayBeforeDriveMore_WNdayInds;
                    disp(['SKIPPING ' birdname '-' exptname ' because starts "driving more" before the last desired day to plot']);
                    tmp=0;
                end
            end
        end
    end
    
    if tmp==0;
        disp(['threw out expt: ' birdname '-' exptname]);
        continue;
        
    end
    
    ExptsToKeep=[ExptsToKeep i];
    
     
end



%% =========================== PLOT THINGS (each syl relative to its own target)
% NOTE: for this to work well need to only do this on expts with full data
% for these days.

syltypefield='Similar';
    NontargArray=[];
    TargArray=[];
% ===== for each day, plot regression (target shift vs. nontarg shift)
for ee=1:length(ExptsToKeep)
    
    exptnum=ExptsToKeep(ee);
    
    % --- NONtarget syl (either same or diff type)
    NumSyls=length(SeqDepPitch_AcrossBirds.AllExpt(exptnum).(syltypefield).LearningRelTarg_all);
    
    for kk=1:NumSyls
        
    % --- nontarget
    Data_days=nan(1,max(DaysToPlot)+1);
    
    Data_days(1)=SeqDepPitch_AcrossBirds.AllExpt(exptnum).(syltypefield).DayVals_Baseline_zscore_targdirsign{kk}(end); % baseline
    Data_days(DaysToPlot+1)=SeqDepPitch_AcrossBirds.AllExpt(exptnum).(syltypefield).DayVals_DurWN_zscore_targdirsign{kk}(DaysToPlot); % WN

    NontargArray=[NontargArray; Data_days];
    
    
    % --- target
    Data_days=nan(1,max(DaysToPlot)+1);
    
    Data_days(1)=SeqDepPitch_AcrossBirds.AllExpt(exptnum).Target.DayVals_Baseline_zscore_targdirsign{1}(end); % baseline
    Data_days(DaysToPlot+1)=SeqDepPitch_AcrossBirds.AllExpt(exptnum).Target.DayVals_DurWN_zscore_targdirsign{1}(DaysToPlot); % WN
    
    TargArray=[TargArray; Data_days];
    
    end
end

% ========= 1) regression for each day, targ vs nontarg (cumulative shift)
lt_figure; hold on;
numrows=ceil(size(NontargArray, 2)/3);
hplots=[];
for i=1:size(NontargArray, 2);
    % this day
    hplot=lt_subplot(numrows, 3, i); hold on;
    title(['WN day: ' num2str(i-1)]);
    xlabel('targ cum shift');
    ylabel('nontarg cum shift');
    
    X=TargArray(:, i);
    Y=NontargArray(:, i);
    
    lt_regress(Y, X, 1, 0, 1, 1, 'k');
    hplots=[hplots hplot];
    line(xlim, [0 0]);
    line([0 0], ylim);
end
linkaxes(hplots, 'xy');




% ========= 2) regression for diff from previous day, targ vs nontarg (cumulative shift)
lt_figure; hold on;
numrows=ceil(size(NontargArray, 2)/3);
hplots=[];
for i=2:size(NontargArray, 2);
    % this day
    hplot=lt_subplot(numrows, 3, i-1); hold on;
    title(['WN day ' num2str(i-2) ' minus ' num2str(i-1)]);
    xlabel('targ day diff');
    ylabel('nontarg day diff');
    
    X=TargArray(:, i)-TargArray(:, i-1);
    Y=NontargArray(:, i)-NontargArray(:, i-1);
    
    lt_regress(Y, X, 1, 0, 1, 1, 'k');
    hplots=[hplots hplot];
    
    line(xlim, [0 0]);
    line([0 0], ylim);
end
linkaxes(hplots, 'xy');



% ======== 3) DAY TO DAY absolute shift by nontarget (forget about target)
lt_figure; hold on;
title('mean day to day shift by nontarget (not normalized to target');
xlabel('day (1 means day 1 minus last bline day');
ylabel('pitch diff');

Y=diff(NontargArray, 1, 2);
Ymean=mean(Y, 1);
Ysem=lt_sem(Y);
X=1:size(Y, 2);

plot(X, Y, 'ob');
lt_plot_bar(X, Ymean, {'Errors', Ysem});

lt_plot_zeroline;


% ============= 3.b - repeated measures anova to test whether day has
% effect
[p, table]=anova_rm(Y)

[~, p]=ttest(Y(:,1), Y(:,4))


%% BELOW: OLD VERSION, STILL GOOD FOR ALL SYLS ANALYSIS
% %% === FIRST, Figure out what datapoints are eligible for this analysis
% % === LIMIT ANALYSIS ONLY TO SYLLABLES WITH DATA FOR ALL DAYS IN
% % DaysToPlot, (i.e. enoguh days, and no nan)
% if OnlyAnalyzeIfHasAllDaysData==1;
%     Inds_HasData=[];
%     for i=1:length(SeqDepPitch_AcrossBirds.AllSyllables.AcousticDist_all);
%         tmp=1; % start as pass. if fail anuything, then goes to 0
%         % failure: not enough days
%         if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{i})<length(DaysToPlot);
%             tmp=0;
%         elseif any(isnan(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{i}(DaysToPlot)));
%             % failure: enough days, but some are nantmp=0;
%             tmp=0;
%         end
%         Inds_HasData=[Inds_HasData tmp];
%     end
%
%     Inds_HasData=logical(Inds_HasData);
% else
%     Inds_HasData=logical(ones(1,length(SeqDepPitch_AcrossBirds.AllSyllables.AcousticDist_all)));
% end
%
% numthrownout=sum(~Inds_HasData);
% numkept=sum(Inds_HasData);
% disp(' --- ');
% disp(['Threw out ' num2str(numthrownout) ' out of ' num2str(numthrownout+numkept) ' total syllables']);
%
%
% %% EACH EXPERIMENT ONE DATAPOINT [OBSOLETE USE ABOVE]
% % each experiment gives target, nontargets, similar, diff
% lt_figure; hold on;
%
% % === FIRST, go through each experiment. for each, collect the inds for that expt and plot.
% % Each experiment gets its own set of inds
% for i=1:NumBirds;
%     numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%
%     for ii=1:numexperiments;
%
%         birdnum=i;
%         exptnum=ii;
%
%         % ===== PLOT AND COLLECT DATA JUST FOR THIS EXPERIMENT
%         inds_expt = SeqDepPitch_AcrossBirds.AllSyllables.BirdNum_all==birdnum & ...
%             SeqDepPitch_AcrossBirds.AllSyllables.ExptNum_all==exptnum;
%
%
%         % -- Targets
%         inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1 & Inds_HasData & inds_expt);
%         if ~isempty(inds)>0;
%             % ==== Collect across days
%             Data_across_days=nan(1,max(DaysToPlot)+1);
%             % 1) last baseline day
%             tmp=[];
%             for j=inds;
%                 tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{j}(end)];
%             end
%             Data_across_days(1)=mean(tmp); % baseline last day
%
%             % 2) WN days
%             for j=DaysToPlot;
%                 daynum=j;
%
%                 tmp=[];
%                 for jj=inds;
%                     if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{jj})>=daynum
%                         tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{jj}(daynum)];
%                     else
%                         tmp=[tmp nan];
%                     end
%                 end
%
%                 Data_across_days(daynum+1) = mean(tmp);
%             end
%
%         % plot and save
%         X=[-1 DaysToPlot];
%         plot(X, Data_across_days, '-k');
%
%         end
%
%
%
% %         -- Nontargets
%         inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & Inds_HasData & inds_expt);
%         if ~isempty(inds)>0;
%             % ==== Collect across days
%             Data_across_days=nan(1,max(DaysToPlot)+1);
%             % 1) last baseline day
%             tmp=[];
%             for j=inds;
%                 tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{j}(end)];
%             end
%             Data_across_days(1) = mean(tmp);
%
%             % 2) WN days
%             for j=DaysToPlot;
%                 daynum=j;
%
%                 tmp=[];
%                 for jj=inds;
%                     if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{jj})>=daynum
%                         tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{jj}(daynum)];
%                     else
%                         tmp=[tmp nan];
%                     end
%                 end
%
%                 Data_across_days(daynum+1) = mean(tmp);
%             end
%         X=[-1 DaysToPlot];
%         plot(X, Data_across_days, '--k');
%         end
%
%
%
% % %       STRUCTURE.
% %
% %         % -- Similar
% %         inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData & inds_expt);
% %         if ~isempty(inds)>0;
% %             % ==== Collect across days
% %             Data_across_days=[];
% %             % 1) last baseline day
% %             tmp=[];
% %             for j=inds;
% %                 tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{j}(end)];
% %             end
% %             Data_across_days=[Data_across_days mean(tmp)];
% %
% %             % 2) WN days
% %             for j=DaysToPlot;
% %                 daynum=j;
% %
% %                 tmp=[];
% %                 for jj=inds;
% %                     if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{jj})>=daynum
% %                         tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{jj}(daynum)];
% %                     else
% %                         tmp=[tmp nan];
% %                     end
% %                 end
% %
% %                 Data_across_days=[Data_across_days mean(tmp)];
% %             end
% %         end
% %
% %         X=[-1 DaysToPlot];
% %         plot(X, Data_across_days, '-b');
% %
% %         % -- Diff
% %         inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData & inds_expt);
% %         if ~isempty(inds)>0;
% %             % ==== Collect across days
% %             Data_across_days=[];
% %             % 1) last baseline day
% %             tmp=[];
% %             for j=inds;
% %                 tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{j}(end)];
% %             end
% %             Data_across_days=[Data_across_days mean(tmp)];
% %
% %             % 2) WN days
% %             for j=DaysToPlot;
% %                 daynum=j;
% %
% %                 tmp=[];
% %                 for jj=inds;
% %                     if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{jj})>=daynum
% %                         tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{jj}(daynum)];
% %                     else
% %                         tmp=[tmp nan];
% %                     end
% %                 end
% %
% %                 Data_across_days=[Data_across_days mean(tmp)];
% %             end
% %         end
% %
% %         X=[-1 DaysToPlot];
% %         plot(X, Data_across_days, '-r');
%
%
%     end
% end
%
%
%




%% ACROSS ALL EXPERIMENTS - COMBINE ALL SYLS

% ++++++++++++++++++++++++++++++++  FIRST, Figure out what datapoints are eligible for this analysis
% === LIMIT ANALYSIS ONLY TO SYLLABLES WITH DATA FOR ALL DAYS IN
% DaysToPlot, (i.e. enoguh days, and no nan)
Inds_HasData=[];
for i=1:length(SeqDepPitch_AcrossBirds.AllSyllables.AcousticDist_all);
    tmp=1; % start as pass. if fail anuything, then goes to 0
    % failure: not enough days
    
    birdnum=SeqDepPitch_AcrossBirds.AllSyllables.BirdNum_all(i);
    exptnum=SeqDepPitch_AcrossBirds.AllSyllables.ExptNum_all(i);
    
    birdname=SeqDepPitch_AcrossBirds.birds{birdnum}.birdname;
    exptname=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.ExptID;
    
    % Check 1) complete days
    if OnlyAnalyzeIfHasAllDaysData==1;
        
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{i})<length(DaysToPlot);
            tmp=0;
            disp(['lacking enough days ' birdname '-' exptname])
        elseif any(isnan(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{i}(DaysToPlot)));
            % failure: enough days, but some are nantmp=0;
            tmp=0;
             disp(['some days are nan ' birdname '-' exptname])
        end
    end
    
    % --- Check 2) goes into multidir days
    % does day window go into days when multidir learning has
    % started?
    if isfield(SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.DATES, 'MultiDir_OneDayBeforeStart_Ind');
        OneDayBeforeMultidirStart=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.DATES.MultiDir_OneDayBeforeStart_Ind;
        NumBaselineDays=max(SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Data_PlotLearning.Params.SeqFilter.BaselineDays);
        NumDaysNeeded=max(DaysToPlot);
        
        if NumBaselineDays+NumDaysNeeded>OneDayBeforeMultidirStart;
            % then bad, not enough single dir days
            tmp=0;
            disp(['syllable num ' num2str(i) ' miltidir days started too soon']);
                        disp(['multidir start too soon ' birdname '-' exptname])

        end
    end
    
    % ==== Check 2.5) goes into samedir phase
    if isfield(SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.DATES, 'SameDir_OneDayBeforeStart_Ind');
        OneDayBeforeMultidirStart=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.DATES.SameDir_OneDayBeforeStart_Ind;
        NumBaselineDays=max(SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Data_PlotLearning.Params.SeqFilter.BaselineDays);
        NumDaysNeeded=max(DaysToPlot);
        
        if NumBaselineDays+NumDaysNeeded>OneDayBeforeMultidirStart;
            % then bad, not enough single dir days
            tmp=0;
            disp(['syllable num ' num2str(i) ' miltidir days started too soon']);
                        disp(['samedir start too soon ' birdname '-' exptname])

        end
    end
    
    
    
    % Check 3) Goes into "drive more" days
    % Throw out specific experiments that "drove more" pitch
    % change during the days window.
    if ThrowOutIfDriveMore==1;
        
        % check to see whether this bird/expt is on the list
        inds1=find(strcmp(DriveMoreList, birdname));
        for kk=1:length(inds1);
            if strcmp(DriveMoreList{inds1(kk)+1},exptname);
                % THEN CHECK TO SEE IF DAYS IS PAST START OF DRIVE
                % MORE
                LastWNDayNeeded=max(DaysToPlot);
                LastDayBeforeDriveMore_WNdayInds=DriveMoreList{inds1(kk)+2};
                
                if LastWNDayNeeded>LastDayBeforeDriveMore_WNdayInds;
                    disp(['SKIPPING ' birdname '-' exptname ' because starts "driving more" before the last desired day to plot']);
                    tmp=0;
                end
            end
        end
    end
    
    % --- DONE
    Inds_HasData=[Inds_HasData tmp];
end

Inds_HasData=logical(Inds_HasData);
% else
%     Inds_HasData=logical(ones(1,length(SeqDepPitch_AcrossBirds.AllSyllables.AcousticDist_all)));

% -----
numthrownout=sum(~Inds_HasData);
numkept=sum(Inds_HasData);
disp(' --- ');
disp(['Threw out ' num2str(numthrownout) ' out of ' num2str(numthrownout+numkept) ' total syllables']);



%% ++++++++++++++++++++++++++++++++++ START PLOTS - FIRST GENERAL LEARNING
count=1;
SubplotsPerFig=8;
subplotrows=2;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];


% ==================== SUBPLOT 1 - for each day plot one dot for mean across all targets, flipping sign if
% learning in other dir
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('mean across all experiments');
xlabel('days, post WN on');
ylabel('z-score rel baseline');

% -- Targets
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1 & Inds_HasData);
% plot last bline day
tmp=[];
for ii=inds;
    tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end)];
end

% plot(-1, tmp, '.k')
errorbar(-1, nanmean(tmp), lt_sem(tmp),'-ok','MarkerSize',7);
% plot WN days
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
%     plot(daynum, tmp, '.k')
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-ok','MarkerSize',7);
end

% -- NonTarget
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & Inds_HasData);
% plot last bline day
tmp=[];
for ii=inds;
    tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end)];
end

errorbar(-1, nanmean(tmp), lt_sem(tmp),'-ok','MarkerSize',7);
% plot wn days
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-ok','MarkerSize',7);
end

% -- Similar
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData);
% plot last bline day
tmp=[];
for ii=inds;
    tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end)];
end

errorbar(-1, nanmean(tmp), lt_sem(tmp),'-ob','MarkerSize',7);
% plot wn days
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-ob','MarkerSize',7);
end

% -- Different
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData);
% plot last bline day
tmp=[];
for ii=inds;
    tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end)];
end

errorbar(-1, nanmean(tmp), lt_sem(tmp),'-or','MarkerSize',7);
% plot wn days
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-or','MarkerSize',7);
end

lt_plot_zeroline
line([0 0], ylim , 'Color','k');
legend('bk:targets','wh:nontargets','bu:similar','rd:diff')



% ==================== SUBPLOT 1b - ONLY PLOTTING NONTARGETS
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('][ONLY NONTARGETS] day means (bk, wh, bu, rd: targ, nontarg, sim, diff)');
xlabel('days, post WN on');
ylabel('z-score rel baseline');

% -- NonTarget
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & Inds_HasData);
% plot last bline day
tmp=[];
for ii=inds;
    tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end)];
end

errorbar(-1, nanmean(tmp), lt_sem(tmp),'-ok','MarkerSize',7);
% plot wn days
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-ok','MarkerSize',7);
end

% -- Similar
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData);
% plot last bline day
tmp=[];
for ii=inds;
    tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end)];
end

errorbar(-1, nanmean(tmp), lt_sem(tmp),'-ob','MarkerSize',7);
% plot wn days
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-ob','MarkerSize',7);
end

% -- Different
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData);
% plot last bline day
tmp=[];
for ii=inds;
    tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end)];
end

errorbar(-1, nanmean(tmp), lt_sem(tmp),'-or','MarkerSize',7);
% plot wn days
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-or','MarkerSize',7);
end

lt_plot_zeroline
line([0 0], ylim , 'Color','k');



%% =========== FRACTION OF WN LEARNING

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('bk(nontargs), bl(sim), rd(diff)');
xlabel('day');
ylabel('norm to global targ mean');

% -- Targets
Target_Means=[];
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    % == Collect each days vals
    Target_Means(daynum)=nanmean(tmp);
end


% -- All nontargets
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & Inds_HasData);
    
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    % get relative to target
    tmp=tmp./Target_Means(daynum);
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-ok','MarkerSize',7, 'MarkerFaceColor','k');
end




% -- Similar
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData);
    
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    % get relative to target
    tmp=tmp./Target_Means(daynum);
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-ob','MarkerSize',7);
end


% -- diff
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData);
    
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    % get relative to target
    tmp=tmp./Target_Means(daynum);
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-or','MarkerSize',7);
end


lt_plot_zeroline


%% ===================== fraction of cumulative learning, as a regression



%% +++++++++++++++++++ PRESYL PLOTS

% ==================== SUBPLOT 1b - ONLY PLOTTING NONTARGETS
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[PRESYL SIM(bk), DIFF(open)]');
xlabel('days, post WN on');
ylabel('z-score rel baseline');

% -- Similar [PRESIM]
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1);
% plot last bline day
tmp=[];
for ii=inds;
    tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end)];
end

errorbar(-1, nanmean(tmp), lt_sem(tmp),'-ob','MarkerSize',7, 'MarkerFaceColor','k');
% plot wn days
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-ob','MarkerSize',7, 'MarkerFaceColor','k');
end


% -- Similar [PREDIFF]
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0);
% plot last bline day
tmp=[];
for ii=inds;
    tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end)];
end

errorbar(-1, nanmean(tmp), lt_sem(tmp),'-ob','MarkerSize',7);
% plot wn days
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-ob','MarkerSize',7);
end

% -- DIFFERENT [PRESIM]
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1);
% plot last bline day
tmp=[];
for ii=inds;
    tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end)];
end

errorbar(-1, nanmean(tmp), lt_sem(tmp),'-or','MarkerSize',7, 'MarkerFaceColor','k');
% plot wn days
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-or','MarkerSize',7, 'MarkerFaceColor','k');
end


% -- DIFFERENT [PREDIFF]
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0);
% plot last bline day
tmp=[];
for ii=inds;
    tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end)];
end

errorbar(-1, nanmean(tmp), lt_sem(tmp),'-or','MarkerSize',7);
% plot wn days
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-or','MarkerSize',7);
end


lt_plot_zeroline
line([0 0], ylim , 'Color','k');





% ==================== SUBPLOT 2 - for each day plot ratio of pitch shift
% relative to target
% learning in other dir
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[PRESYL SIM (BLACK), DIFF (OPEN)');
xlabel('day');
ylabel('norm to global targ mean');

% -- Targets
Target_Means=[];
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    % == Collect each days vals
    Target_Means(daynum)=nanmean(tmp);
end

% -- Similar [presimilar]
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1);
    
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    % get relative to target
    tmp=tmp./Target_Means(daynum);
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-ob','MarkerSize',7, 'MarkerFaceColor','k');
end


% -- Similar [prediff]
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0);
    
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    % get relative to target
    tmp=tmp./Target_Means(daynum);
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-ob','MarkerSize',7);
end

% -- Diff [presimilar]
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1);
    
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    % get relative to target
    tmp=tmp./Target_Means(daynum);
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-or','MarkerSize',7, 'MarkerFaceColor','k');
end


% -- Diff [prediff]
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0);
    
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            tmp=[tmp SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)];
        else
            tmp=[tmp nan];
        end
    end
    
    % get relative to target
    tmp=tmp./Target_Means(daynum);
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-or','MarkerSize',7);
end



lt_plot_zeroline;





% ============================== SUBPLOT 3) FOR EACH DAY PLOT THE DIFFERENCE BETWEEN
% TODAY AND YESTERDAY
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[PRESYL] Pitch shift rel yesterday')
xlabel('day');
ylabel('shift from yesterday (z-score diff)');

% -- Targets
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'ok','MarkerFaceColor','k', 'MarkerSize',7);
    
end


% -- Similar [PRESIM]
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'ob', 'MarkerSize',7, 'MarkerFaceColor','k');
end

% -- Similar [PREDIFF]
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'ob', 'MarkerSize',7);
end


% -- diff [PRESIM]
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'or', 'MarkerSize',7, 'MarkerFaceColor','k');
end



% -- diff [PREDIFF]
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'or', 'MarkerSize',7);
end






% ============================== SUBPLOT 4) DIFFERENCES, BUT NORMALIZING TO
% TARGET
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[SIMILAR, PRESYL shaded] Pitch shift rel yesterday')
xlabel('day');
ylabel('shift from yesterday (first norm by targ, then mean)');

% -- Targets
TargetMeans=[];
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    % === collect values for this day
    TargetMeans(daynum)=nanmean(tmp);
    
end

% -- Similar [PRESIM]
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1);
Ymeans=[];
Ysems=[];
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    % normalize to target
    tmp=tmp./TargetMeans(daynum);
    
    % ---- COLLECT
    Ymeans=[Ymeans nanmean(tmp)];
    Ysems=[Ysems lt_sem(tmp)];
end    
    
    hbar=lt_plot_bar(DaysToPlot-0.15, Ymeans, {'Errors', Ysems, 'BarWidth',0.3});
    set(hbar, 'FaceColor', 'b');
    hold on
    
% -- Similar [PREDIFF]
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0);
Ymeans=[];
Ysems=[];
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    % normalize to target
    tmp=tmp./TargetMeans(daynum);
    
    % ---- COLLECT
    Ymeans=[Ymeans nanmean(tmp)];
    Ysems=[Ysems lt_sem(tmp)];
end    
    
    hbar=lt_plot_bar(DaysToPlot+0.15, Ymeans, {'Errors', Ysems, 'BarWidth',0.3});
%     set(hbar, 'FaceColor', '');
        ylim([-0.4 1])

    
% ============================== SUBPLOT 4) DIFFERENCES, BUT NORMALIZING TO
% TARGET
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[DIFF, PRESYL shaded] Pitch shift rel yesterday')
xlabel('day');
ylabel('shift from yesterday (first norm by targ, then mean)');

% -- Targets
TargetMeans=[];
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    % === collect values for this day
    TargetMeans(daynum)=nanmean(tmp);
    
end

% -- Different [PRESIM]
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1);
Ymeans=[];
Ysems=[];
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    % normalize to target
    tmp=tmp./TargetMeans(daynum);
    
    % ---- COLLECT
    Ymeans=[Ymeans nanmean(tmp)];
    Ysems=[Ysems lt_sem(tmp)];
end    
    
    hbar=lt_plot_bar(DaysToPlot-0.15, Ymeans, {'Errors', Ysems, 'BarWidth',0.3});
    set(hbar, 'FaceColor', 'r');
    hold on
    
    
% -- DIFF [PREDIFF]
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData ...
    & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0);
Ymeans=[];
Ysems=[];
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    % normalize to target
    tmp=tmp./TargetMeans(daynum);
    
    % ---- COLLECT
    Ymeans=[Ymeans nanmean(tmp)];
    Ysems=[Ysems lt_sem(tmp)];
end    
    
    hbar=lt_plot_bar(DaysToPlot+0.15, Ymeans, {'Errors', Ysems, 'BarWidth',0.3});

    
    ylim([-0.4 1])
    
    
    
    % -- ALL NONTARGETS
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Pitch shift rel yesterday [all nontargs]')
xlabel('day');
ylabel('shift from yesterday (first norm by targ, then mean)');

inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & Inds_HasData);

Ymeans=[];
Ysems=[];
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    % normalize to target
    tmp=tmp./TargetMeans(daynum);
    
    % ---- COLLECT
    Ymeans=[Ymeans nanmean(tmp)];
    Ysems=[Ysems lt_sem(tmp)];
end    
    
    hbar=lt_plot_bar(DaysToPlot+0.15, Ymeans, {'Errors', Ysems, 'BarWidth',0.3});

    
    ylim([-0.4 1])



%% TEST WHETHER EARLY DAYS SHOW MORE LOCKING TO TARGET SHIFT THAN LATER DAYS
% DAY SHIFTS

% ========= INPUTS
Earlydays=1:2;
Latedays=3:4;


% --- ALL NONTARGS
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 ...
    & Inds_HasData);
lt_seq_dep_pitch_ACROSSBIRDS_Separation_sub2;
title(['[NONTARGS] Early=' num2str(Earlydays) '; Late=' num2str(Latedays)])

% ------- SIMILAR
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 ...
    & Inds_HasData);
lt_seq_dep_pitch_ACROSSBIRDS_Separation_sub2;
title(['[SIM] Early=' num2str(Earlydays) '; Late=' num2str(Latedays)])

% ------- DIFFERENT
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
    & Inds_HasData);
lt_seq_dep_pitch_ACROSSBIRDS_Separation_sub2;
title(['[DIFF] Early=' num2str(Earlydays) '; Late=' num2str(Latedays)])



% % inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData ...
% %     & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0);
% 
% inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 ...
%     & Inds_HasData);
% 



    
    
%% TEST WHETHER EARLY DAYS HAVE HIGHER GENERALIZATION THAN LATER DAYS [cumulative learning at targ]

% ======================= INPUTS
Earlydays1=6:8;
Latedays1=9:11;
% inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData ...
%     & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0);

% ============== ALL NONTARGETS
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & Inds_HasData);
    
lt_seq_dep_pitch_ACROSSBIRDS_Separation_sub1;
title(['[All nontargs] Early=' num2str(Earlydays1) '; Late=' num2str(Latedays1)])



% ============== SIMILAR
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & ...
    SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData);
    
lt_seq_dep_pitch_ACROSSBIRDS_Separation_sub1;
title(['[SIMILAR] Early=' num2str(Earlydays1) '; Late=' num2str(Latedays1)])



% ============== DIFFERENT
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & ...
    SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData);
    
lt_seq_dep_pitch_ACROSSBIRDS_Separation_sub1;
title(['[DIFF] Early=' num2str(Earlydays1) '; Late=' num2str(Latedays1)])



% ============== 
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & ...
    SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==1 ...
    & Inds_HasData);
    
lt_seq_dep_pitch_ACROSSBIRDS_Separation_sub1;
title(['[SIMILAR-presim] Early=' num2str(Earlydays1) '; Late=' num2str(Latedays1)])



% ============== 
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & ...
    SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all==0 ...
    & Inds_HasData);
    
lt_seq_dep_pitch_ACROSSBIRDS_Separation_sub1;
title(['[SIMILAR-prediff] Early=' num2str(Earlydays1) '; Late=' num2str(Latedays1)])




%% CONTINUE PLOTTING +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% ============================== SUBPLOT 3) FOR EACH DAY PLOT THE DIFFERENCE BETWEEN
% TODAY AND YESTERDAY
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Pitch shift relative to yesterday')
xlabel('day');
ylabel('shift from yesterday (z-score diff)');

% -- Targets
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'ok','MarkerFaceColor','k', 'MarkerSize',7);
    
end

% -- nonTargets
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'ok','MarkerSize',7);
    
end

% -- Similar
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'ob', 'MarkerSize',7);
    
end

% -- Diff
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'or', 'MarkerSize',7);
    
end

lt_plot_zeroline;


% ============================== SUBPLOT 4) DIFFERENCES, BUT NORMALIZING TO
% TARGET
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Pitch shift relative to yesterday (divided by targ)')
xlabel('day');
ylabel('[diff] first norm by targ, then take mean');

% -- Targets
TargetMeans=[];
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    % === collect values for this day
    TargetMeans(daynum)=nanmean(tmp);
    
end

% -- nonTargets
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    % normalize to target
    tmp=tmp./TargetMeans(daynum);
    
%     plot(daynum, tmp, '.k')
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-ok', 'MarkerSize',7);
end

% -- Similar
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    % normalize to target
    tmp=tmp./TargetMeans(daynum);
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-ob','MarkerSize',7);
    
end

% -- Diff
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    % normalize to target
    tmp=tmp./TargetMeans(daynum);
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-or','MarkerSize',7);
    
end

lt_plot_zeroline;



% ============================== SUBPLOT 5) [SAME AS SUBPLOT 4, BUT ZOOMING IN TO WHEN TARG IS LAERNING
% DIFFERENCES, BUT NORMALIZING TO TARGET
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[ZOOMING IN Y] Pitch shift relative to yesterday (divided by targ)')
xlabel('day');
ylabel('z-score diff, divided by target)');

% -- Targets
TargetMeans=[];
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    
    % === collect values for this day
    TargetMeans(daynum)=nanmean(tmp);
    
end

% -- nonTargets
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    % normalize to target
    tmp=tmp./TargetMeans(daynum);
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-ok','MarkerFaceColor','w', 'MarkerSize',7);
    
end

% -- Similar
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    % normalize to target
    tmp=tmp./TargetMeans(daynum);
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-ob','MarkerFaceColor','b', 'MarkerSize',7);
    
end

% -- Diff
inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData);
for i=DaysToPlot;
    daynum=i;
    
    tmp=[];
    for ii=inds;
        if length(SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii})>=daynum
            % get difference in pitch of today minus yesterday
            if daynum==1;
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)-...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign{ii}(end);
            else
                val=SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum)- ...
                    SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign{ii}(daynum-1);
            end
            
            tmp=[tmp val];
            
        else
            tmp=[tmp nan];
        end
    end
    
    % normalize to target
    tmp=tmp./TargetMeans(daynum);
    
    errorbar(daynum, nanmean(tmp), lt_sem(tmp),'-or','MarkerFaceColor','r', 'MarkerSize',7);
    
end

lt_plot_zeroline;
xlim([0 5]);
ylim([-1 1]);

%% TROUBLESHOOTING - verifying that trajecvtory plots give same learning value at learning day as calculated previously
if (0)
    lt_figure; hold on;
    
    % targ
    inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1 & Inds_HasData);
    ffvals=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);
    
    lt_plot(1, mean(ffvals), {'Color','k'});
    
    
    % nontarg
    inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & Inds_HasData);
    ffvals=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);
    
    lt_plot(1, mean(ffvals), {'Color','k' , 'MarkerFaceColor','w'});
    
    % sim
    inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & Inds_HasData);
    ffvals=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);
    
    lt_plot(1, mean(ffvals), {'Color','b'});
    
    % diff
    inds=find(SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & Inds_HasData);
    ffvals=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All(inds);
    
    lt_plot(1, mean(ffvals), {'Color','r'});
end

