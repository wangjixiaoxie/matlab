function SeqDepPitch_AcrossBirds=lt_seq_dep_pitch_ACROSSBIRDS_RecalcBaseline(SeqDepPitch_AcrossBirds, Params,OverwriteLearningMetric, ConvertToHz)
%% LT 12/17/15 - recalcualte learning using last 2 baseline days. allows to compare with distrubtion
NumBirds=length(SeqDepPitch_AcrossBirds.birds);


if ConvertToHz==1  
    % since is not default, let use know
    figure; 
    lt_plot_text(0, 0.5, 'Convert to HZ, sure? (type anything)');
    pause;
end


%%  FIRST COLLECT LEARNING FOR ALL BIRDS (EVEN THOSE WITH SHORT BASELINES) - ALL BIRDS MUST HAVE BASELINE LAST 2 DAYS)

count=1;
subplotrows=4;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];

disp(' ');

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % baseline
        BaselineDayInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.BaselineDays;
       
        
        % learning
        WN_start_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        Learning_day_inds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.LearningMetric.dayIndsUsed;
        LearningWindow_days=length(Learning_day_inds);
        NumWNDays=Learning_day_inds(end)-WN_start_ind+1; % to end of epoch used to quantify learning
        
        [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title([birdname '-' exptname]);
        xlabel('old (all bline days)');
        ylabel('new (2 days)');
        
        % ===== COLLECT DATA
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        % ==== DISPLAY
        disp(' ');
        disp([birdname ' - ' exptname]);
        disp('Sample sizes: ');
        CCcount=0;
        
        plotcols=lt_make_plot_colors(length(syls_unique), 0,0);
        
        % 
        
        
        % ==== SKIP IF DON"T HAVE LAERNING VALS (e.g. expt not long enough)
        if isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syls_unique{1}).LEARNING.learning_metric.mean)
            disp(' ');
            disp(['SKIPPED!!!!  ' birdname '-' exptname ': learning metric is already nan'])
            
            continue
        end
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            % =============== 2) ONLY ONE BASELINE WINDOW (UP TO END OF BASELINE) Pitch diff for baseline window: up to last WN day,
            % with equal number of days to learning. If WN uses 2 day mean,
            % then this must match. [CONFIRMED WORKS WITH ANY BINSIZE OF
            % DAYS]
            Learning_start_inds=BaselineDayInds(end)-LearningWindow_days+1:BaselineDayInds(end); % during end of baseline
            Learning_end_inds=Learning_day_inds;
            
            if CCcount==0;
                
                disp(['Learning start - end: ' num2str([Learning_start_inds Learning_end_inds])]);
                CCcount=1;
            end
            
            % Collect those pitch differences
            learning_start_ffvals=[];
            for k=Learning_start_inds;
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                learning_start_ffvals=[learning_start_ffvals SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{k}];
                else
                learning_start_ffvals=[learning_start_ffvals cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})];
                end
            end
            
            learning_end_ffvals=[];
            if ~isnan(Learning_end_inds)
            for k=Learning_end_inds;
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                learning_end_ffvals=[learning_end_ffvals SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{k}];
                else                
                learning_end_ffvals=[learning_end_ffvals cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})];
                end
            end
            end
            
            % ===== DISPLAY VARIOUS THINGS ABOUT THIS SYL
            disp([syl ' - ' num2str([length(learning_start_ffvals), length(learning_end_ffvals)])]);
            
            %
            %             % =============== 3) SUBSAMPLE to lowest sample size (if multiple days in a
            %             % bin, first mix those days, then subsample)
            %             N_subsample = min([length(learning_start_ffvals), length(learning_end_ffvals), length(baseline_end_ffvals), length(baseline_start_ffvals)]);
            %
            %             inds=randperm(N_subsample);
            %             learning_start_ffvals=learning_start_ffvals(inds);
            %
            %             inds=randperm(N_subsample);
            %             learning_end_ffvals=learning_end_ffvals(inds);
            %
            %             inds=randperm(N_subsample);
            %             baseline_end_ffvals=baseline_end_ffvals(inds);
            %
            %             inds=randperm(N_subsample);
            %             baseline_start_ffvals=baseline_start_ffvals(inds);
            
            
            % =============== 4) CALCULATE STATISTICS
            % ----------- difference (of means)
            % learning
            mu1=mean(learning_start_ffvals);
            mu2=mean(learning_end_ffvals);
            
            mudiff_learning=mu2-mu1;
            mudiff_sem=lt_sem(learning_end_ffvals-mu1);
            
            
            %             % ----------- dprime
            %             % learning
            %             std1=std(learning_start_ffvals);
            %             std2=std(learning_end_ffvals);
            %             denominator=sqrt((1/2)*(std1^2 + std2^2));
            %
            %             dprime_learning=mudiff_learning/denominator;
            %
            %             % baseline
            %             std1=std(baseline_start_ffvals);
            %             std2=std(baseline_end_ffvals);
            %             denominator=sqrt((1/2)*(std1^2 + std2^2));
            %
            %             dprime_baseline=mudiff_baseline/denominator;
            
            % ----------- zscore (using std of early)
            % learning
            std1=std(learning_start_ffvals);
            zscore_learning=mudiff_learning/std1;
            zscore_learning_sem=lt_sem((learning_end_ffvals-mean(learning_start_ffvals))/std(learning_start_ffvals));
            
            
            %             % ----------- ranksum p values
            %             p_ranksum_learning=ranksum(learning_start_ffvals, learning_end_ffvals);
            %             p_ranksum_baseline=ranksum(baseline_start_ffvals, baseline_end_ffvals);
            %
            %
            
            % ===== plot
            try
            lt_plot_45degScatter(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean, zscore_learning, plotcols{j});
            catch err
            end
            
            % ======== OVERWRITE LEARNING STATS FOR THIS BIRD
            if OverwriteLearningMetric==1;
                if ConvertToHz==1
                % use hz diff
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean=mudiff_learning;
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.sem=mudiff_sem;
                
                else
                    % use zscore
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean=zscore_learning;
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.sem=zscore_learning_sem;
                end
            end
            
        end
        
        
        % == SHIFT REL TO TARGET
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            if OverwriteLearningMetric==1;
                
                learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
                
                plot(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ, ...
                    learning/SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean, 'o', 'Color',plotcols{j});
                
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ=learning/...
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
            end
            
        end
        
    end
end

lt_subtitle('closed: zscore, open: generalization');

%% BASELINE COLLECTION - window (day1 to day2) and 1 learning (day3 to day4), with day2=day3, and equal difference of days.
disp(' ');
disp(' +++++++++++++++ BASELINE +++++++++++++++ ');
NumBaseDaysALL=[];
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        % ==== SKIP IF LACKING ENOUGH LEARNING DAYS
        if isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syls_unique{1}).LEARNING.learning_metric.mean)
            disp(' ');
            disp(['SKIPPED!!!!  ' birdname '-' exptname ': learning metric is already nan'])
            
            continue
        end
        
        
        % baseline
        BaselineDayInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.BaselineDays;
        NumBaselineDays=length(BaselineDayInds);
        
        NumBaseDaysALL=[NumBaseDaysALL NumBaselineDays];
        % learning
        WN_start_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        Learning_day_inds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.LearningMetric.dayIndsUsed;
        LearningWindow_days=length(Learning_day_inds);
        NumWNDays=Learning_day_inds(end)-WN_start_ind+1; % to end of epoch used to quantify learning
        
        if NumWNDays+LearningWindow_days>NumBaselineDays
            % then cannot do analysis
            disp(' ');
            disp(['BASELINE ANALYSIS - Threw out: ' birdname ' - ' exptname ', since not enough baseline days']);
            continue;
        end
        
        % ===== COLLECT DATA
        
        % ==== DISPLAY
        disp(' ');
        disp([birdname ' - ' exptname]);
        disp('Sample sizes: ');
        CCcount=0;
        
        
                
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            % =============== 2) ONLY ONE BASELINE WINDOW (UP TO END OF BASELINE) Pitch diff for baseline window: up to last WN day,
            % with equal number of days to learning. If WN uses 2 day mean,
            % then this must match. [CONFIRMED WORKS WITH ANY BINSIZE OF
            % DAYS]
            Learning_start_inds=BaselineDayInds(end)-LearningWindow_days+1:BaselineDayInds(end); % during end of baseline
            Learning_end_inds=Learning_day_inds;
            
            Baseline_end_inds_Shared=Learning_start_inds; % shared with Learning start
            Baseline_start_inds_Shared=Learning_start_inds(end)-LearningWindow_days+1-NumWNDays:Learning_start_inds(end)-NumWNDays;
            
            if CCcount==0;
                
                disp(['Learning start - end: ' num2str([Learning_start_inds Learning_end_inds])]);
                disp(['Baseline start - end: ' num2str([Baseline_start_inds_Shared Baseline_end_inds_Shared])]);
                
                CCcount=1;
            end
            
            %             if strcmp(birdname, 'gr41gr90') & strcmp(exptname,'SeqDepPitchShift');
            %                 keyboard;
            %             end
            %               TROUBLESHOOTING ONLY
            
            % Collect those pitch differences
            baseline_end_ffvals=[];
            for k=Baseline_end_inds_Shared;
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                baseline_end_ffvals=[baseline_end_ffvals SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{k}];
                else
                baseline_end_ffvals=[baseline_end_ffvals cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})];
                end
            end
            
            baseline_start_ffvals=[];
            for k=Baseline_start_inds_Shared;
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                    baseline_start_ffvals=[baseline_start_ffvals SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{k}];
                else
                    baseline_start_ffvals=[baseline_start_ffvals cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})];
                end
            end
            
            % ===== DISPLAY VARIOUS THINGS ABOUT THIS SYL
            disp([syl ' - ' num2str([length(baseline_end_ffvals), length(baseline_start_ffvals)])]);
            
            %
            %             % =============== 3) SUBSAMPLE to lowest sample size (if multiple days in a
            %             % bin, first mix those days, then subsample)
            %             N_subsample = min([length(learning_start_ffvals), length(learning_end_ffvals), length(baseline_end_ffvals), length(baseline_start_ffvals)]);
            %
            %             inds=randperm(N_subsample);
            %             learning_start_ffvals=learning_start_ffvals(inds);
            %
            %             inds=randperm(N_subsample);
            %             learning_end_ffvals=learning_end_ffvals(inds);
            %
            %             inds=randperm(N_subsample);
            %             baseline_end_ffvals=baseline_end_ffvals(inds);
            %
            %             inds=randperm(N_subsample);
            %             baseline_start_ffvals=baseline_start_ffvals(inds);
            
            
            % =============== 4) CALCULATE STATISTICS
            % ----------- difference (of means)
            
            % baseline
            mu1=mean(baseline_start_ffvals);
            mu2=mean(baseline_end_ffvals);
            
            mudiff_baseline=mu2-mu1;
            
            mudiff_sem=lt_sem(baseline_end_ffvals-mu1);
            
            %             % ----------- dprime
            %             % learning
            %             std1=std(learning_start_ffvals);
            %             std2=std(learning_end_ffvals);
            %             denominator=sqrt((1/2)*(std1^2 + std2^2));
            %
            %             dprime_learning=mudiff_learning/denominator;
            %
            %             % baseline
            %             std1=std(baseline_start_ffvals);
            %             std2=std(baseline_end_ffvals);
            %             denominator=sqrt((1/2)*(std1^2 + std2^2));
            %
            %             dprime_baseline=mudiff_baseline/denominator;
            
            % ----------- zscore (using std of early)
            
            % baseline
            std1=std(baseline_start_ffvals);
            zscore_baseline=mudiff_baseline/std1;
            zscore_baseline_sem=lt_sem((baseline_end_ffvals-mean(baseline_start_ffvals))/std(baseline_start_ffvals));
            
            
            %             % ----------- ranksum p values
            %             p_ranksum_learning=ranksum(learning_start_ffvals, learning_end_ffvals);
            %             p_ranksum_baseline=ranksum(baseline_start_ffvals, baseline_end_ffvals);
            %
            %
            
            % ======= SAVE BASELINE DRIFT AS WELL
            if ConvertToHz==1
                % use hz
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.zscore_mean=mudiff_baseline;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.zscore_sem=mudiff_sem;
SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.NOTE='stats are hz, not score';
            
            else
                % use zscore
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.zscore_mean=zscore_baseline;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.zscore_sem=zscore_baseline_sem;
            disp(['baseline zscore (sem): ' num2str(zscore_baseline) '(' num2str(zscore_baseline_sem) ')']);
            end
        end
        
        
        % == SHIFT REL TO TARGET
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            baselinedrift=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.zscore_mean;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.zscore_rel_targ=baselinedrift/...
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
        end
        
    end
end


%% ==== PLOT DISTRUTION OF NUMBER OF BASELINE DAYS PER EXPERIMENT.

lt_figure; hold on;
title('dist of num baseline days across experiments');
xlabel('num base days');
ylabel('prob');

lt_plot_cdf(NumBaseDaysALL);

