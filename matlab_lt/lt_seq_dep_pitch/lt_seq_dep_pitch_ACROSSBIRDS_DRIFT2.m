function [SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_DRIFT2(SeqDepPitch_AcrossBirds, PARAMS, plot_text)
%% LT 9/7/15 - MATCH NUMBER OF DAYS OF WN WITH BASELINE DAYS, AND COMPARE DRIFT
%  Method 1: If have N WN days, then there are N+1 days from last baseline day to WN
%  day. Take all N+1 day windows during baseline and plot distribution of
%  those

% Method 2: Take day to day, and plot empirical distribution of day to day
% shifts

% NOTE: uses learning metric

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

%% DISPLAY NUMBER OF LEARNING DAYS AND BASELINE DAYS FOR ALL EXPERIMENTS

disp(' ');
disp('Bird - Experiment - NumBaselineDays - NumLearningDays');

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;

        BaselineDayInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.BaselineDays;
        
        WN_start_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        Learning_day_inds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.LearningMetric.dayIndsUsed;
        NumWNDays=Learning_day_inds(end)-WN_start_ind+1; % to end of epoch used to quantify learning
        
        disp([birdname ' - ' exptname ' - ' num2str(length(BaselineDayInds)) ' - ' num2str(NumWNDays)]);
    
    end
end

%% METHOD: 1 baseline window (day1 to day2) and 1 learning (day3 to day4), with day2=day3, and equal difference of days.
% === COLLECT BASELINE AND LEARNING STATS FOR ALL EXPERIMENTS THAT PASS CRITERIA TO PERFORM THIS ANALYSIS 
% i.e. need enough baseline days to match learning days + 1;
disp(' ');

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % baseline
        BaselineDayInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.BaselineDays;
        NumBaselineDays=length(BaselineDayInds);
        
        
        % learning
        WN_start_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        Learning_day_inds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.LearningMetric.dayIndsUsed;
        LearningWindow_days=length(Learning_day_inds);
        NumWNDays=Learning_day_inds(end)-WN_start_ind+1; % to end of epoch used to quantify learning
        
        if NumWNDays+LearningWindow_days>NumBaselineDays
            % then cannot do analysis
            disp(' ');
            disp(['Threw out: ' birdname ' - ' exptname ', since not enough baseline days']);
            continue;
        end
        
        % ===== COLLECT DATA
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        % ==== DISPLAY
        disp(' ');
        disp([birdname ' - ' exptname]);
        disp('Sample sizes: ');
        count=0;
        
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
            
            if count==0;
            
                disp(['Learning start - end: ' num2str([Learning_start_inds Learning_end_inds])]);
                disp(['Baseline start - end: ' num2str([Baseline_start_inds_Shared Baseline_end_inds_Shared])]);
                
                count=1;
            end
            
%             if strcmp(birdname, 'gr41gr90') & strcmp(exptname,'SeqDepPitchShift');
%                 keyboard;
%             end
%               TROUBLESHOOTING ONLY
            
            % Collect those pitch differences
            learning_start_ffvals=[];
            for k=Learning_start_inds;
                learning_start_ffvals=[learning_start_ffvals cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})];
            end
            
            learning_end_ffvals=[];
            for k=Learning_end_inds;
                learning_end_ffvals=[learning_end_ffvals cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})];
            end
            
            baseline_end_ffvals=[];
            for k=Baseline_end_inds_Shared;
                baseline_end_ffvals=[baseline_end_ffvals cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})];
            end
           
            baseline_start_ffvals=[];
            for k=Baseline_start_inds_Shared;
                baseline_start_ffvals=[baseline_start_ffvals cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})];
            end
        
            % ===== DISPLAY VARIOUS THINGS ABOUT THIS SYL
            disp([syl ' - ' num2str([length(learning_start_ffvals), length(learning_end_ffvals), length(baseline_end_ffvals), length(baseline_start_ffvals)])]);

            
            % =============== 3) SUBSAMPLE to lowest sample size (if multiple days in a
            % bin, first mix those days, then subsample)
            N_subsample = min([length(learning_start_ffvals), length(learning_end_ffvals), length(baseline_end_ffvals), length(baseline_start_ffvals)]);
            
            inds=randperm(N_subsample);
            learning_start_ffvals=learning_start_ffvals(inds);
            
            inds=randperm(N_subsample);
            learning_end_ffvals=learning_end_ffvals(inds);

            inds=randperm(N_subsample);
            baseline_end_ffvals=baseline_end_ffvals(inds);
            
            inds=randperm(N_subsample);
            baseline_start_ffvals=baseline_start_ffvals(inds);
            
            
            % =============== 4) CALCULATE STATISTICS
            % ----------- difference (of means)
            % learning
            mu1=mean(learning_start_ffvals);
            mu2=mean(learning_end_ffvals);
            
            mudiff_learning=mu2-mu1;
            
            % baseline
            mu1=mean(baseline_start_ffvals);
            mu2=mean(baseline_end_ffvals);
            
            mudiff_baseline=mu2-mu1;
            

            % ----------- dprime
            % learning
            std1=std(learning_start_ffvals);
            std2=std(learning_end_ffvals);
            denominator=sqrt((1/2)*(std1^2 + std2^2));
               
            dprime_learning=mudiff_learning/denominator;
            
            % baseline
            std1=std(baseline_start_ffvals);
            std2=std(baseline_end_ffvals);
            denominator=sqrt((1/2)*(std1^2 + std2^2));
               
            dprime_baseline=mudiff_baseline/denominator;
            
            % ----------- zscore (using std of early)
            % learning
            std1=std(learning_start_ffvals);
            zscore_learning=mudiff_learning/std1;
            
            % baseline
            std1=std(baseline_start_ffvals);
            zscore_baseline=mudiff_baseline/std1;
            
                        
            % ----------- ranksum p values
            p_ranksum_learning=ranksum(learning_start_ffvals, learning_end_ffvals);
            p_ranksum_baseline=ranksum(baseline_start_ffvals, baseline_end_ffvals);
            
            
            
            % ======= SAVE ALL THOSE STATS FOR THIS SYL
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.Learning_start_inds=Learning_start_inds;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.Learning_end_inds=Learning_end_inds;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.Baseline_start_inds=Baseline_start_inds_Shared;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.Baseline_end_inds=Baseline_end_inds_Shared;
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.DATA.(syl).learning_start_ffvals=learning_start_ffvals;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.DATA.(syl).learning_end_ffvals=learning_end_ffvals;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.DATA.(syl).baseline_start_ffvals=baseline_start_ffvals;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.DATA.(syl).baseline_end_ffvals=baseline_end_ffvals;
           
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.STATS.(syl).LEARNING.ffmean_diff=mudiff_learning;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.STATS.(syl).LEARNING.dprime=dprime_learning;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.STATS.(syl).LEARNING.zscore=zscore_learning;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.STATS.(syl).LEARNING.p_ranksum=p_ranksum_learning;
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.STATS.(syl).BASELINE.ffmean_diff=mudiff_baseline;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.STATS.(syl).BASELINE.dprime=dprime_baseline;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.STATS.(syl).BASELINE.zscore=zscore_baseline;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.STATS.(syl).BASELINE.p_ranksum=p_ranksum_baseline;
            
            
            
            
            
        end
    end
end

        

%% PLOT ACROSS ALL SYLS, ALL BIRDS
% plot_text=1; % birdname, expt, syl

% figure;
count=1;
SubplotsPerFig=8;
subplotrows=2;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];

Dprime_learning_All=[];
Dprime_base_All=[];
Target_All=[];
Similar_All=[];
Pval_learn_All=[];
Pval_base_All=[];
BirdName_All={};
ExptName_All={};
SylName_All={};
FFdiff_learn_All=[];
FFdiff_base_All=[];

PreSimilar_All=[];
TargLearnDir_All=[];

disp(' ----------------------- ');
disp(' ');

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        % skip if no DataDrift
        if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}, 'Data_DRIFT2');
            disp(['Skipped ' birdname ' - '  exptname ' (no DataDRIFT2, likely too few base days)']);
            continue;
        end
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            % =============== 1) Various things about syl: target? similar?
            istarget=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
            issimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            
                        % ==== targ learning dir
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            targ_learndir=sign(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(targsyl).ConsolEnd_meanFF);

                        % === preceding syllable similar?
            presyl_similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;

            % ====== dprime
            dprime_learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.STATS.(syl).LEARNING.dprime;
            dprime_base=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.STATS.(syl).BASELINE.dprime;
            
            % ===== pval (ranksum)
            pval_learn=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.STATS.(syl).LEARNING.p_ranksum;
            pval_base=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.STATS.(syl).BASELINE.p_ranksum;
            
           % ==== Pitch diff 
            ffdiff_learn=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.STATS.(syl).LEARNING.ffmean_diff;
            ffdiff_base=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_DRIFT2.OneBaselineComparison_shared.STATS.(syl).BASELINE.ffmean_diff;
           
% ===== COLLECT ACROSS ALL DATA
            Dprime_learning_All=[Dprime_learning_All dprime_learning ];
            Dprime_base_All=[Dprime_base_All dprime_base ];
            Pval_learn_All=[Pval_learn_All pval_learn];
            Pval_base_All=[Pval_base_All pval_base];
            FFdiff_learn_All=[FFdiff_learn_All ffdiff_learn];
            FFdiff_base_All=[FFdiff_base_All ffdiff_base];
            
                        TargLearnDir_All=[TargLearnDir_All targ_learndir];
            PreSimilar_All=[PreSimilar_All presyl_similar];

            
            Target_All=[Target_All istarget];
            Similar_All=[Similar_All issimilar];
            BirdName_All=[BirdName_All birdname];
            ExptName_All=[ExptName_All exptname];
            SylName_All=[SylName_All syl];
            
        end
        
    end
end

% +++++++++++++++++++++++++ SUBPLOTS
% === 1) dprime (dot-dot)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('dprimes (all nontargs)');
xlabel('baseline');
ylabel('WN');

% non-targets
inds=Target_All==0;
lt_plot(Dprime_base_All(inds), Dprime_learning_All(inds));

line([0 0] , ylim);
line(xlim, [0 0]);

% === 2) dprime (absolute) (dot-dot)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('abs(dprimes) (all nontargs)');
xlabel('baseline');
ylabel('WN');

% non-targets
inds=Target_All==0;
lt_plot(abs(Dprime_base_All(inds)), abs(Dprime_learning_All(inds)));

line(xlim, ylim)

% === 3) p-val [dot-dot]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('pval(ranksum) (nontargs:bk, targs:yel)');
xlabel('log10(p val) (baseline)');
ylabel('log10(p val) (WN)');

% non-targets
inds=Target_All==0;
lt_plot(log10(Pval_base_All(inds)), log10(Pval_learn_All(inds)));

% targets
inds=Target_All==1;
lt_plot(log10(Pval_base_All(inds)), log10(Pval_learn_All(inds)), {'Color','y'});
        
line([log10(0.05) log10(0.05)], ylim);
line(xlim, [log10(0.05) log10(0.05)]);

if plot_text==1;
    for i=1:length(BirdName_All);
        birdname=BirdName_All{i};
        expt=ExptName_All{i};
        string=[birdname(1:4) '-' expt(end-5:end) '-' SylName_All{i}];
        text(log10(Pval_base_All(i)), log10(Pval_learn_All(i)), string, 'Color','r');
    end
        
end


% ============ DISTRIBUTIONS OF PVALS
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[nontarget] baseline (bk) and learning (rd) pitch shift p-values')
xlabel('log10(ranksum pval)');

% get xcenters
Y=[log10(Pval_base_All(inds)), log10(Pval_learn_All(inds))];
[~, Xcenters, ~]=lt_plot_histogram(Y, '', 0, 0, 1, 1);

% -- nontargets
inds=Target_All==0;

[~,~,hbar]=lt_plot_histogram(log10(Pval_base_All(inds)), Xcenters,1,1,0,1);
set(hbar, 'Color', 'k');
[~,~,hbar]=lt_plot_histogram(log10(Pval_learn_All(inds)), Xcenters,1,1,0,1);
set(hbar, 'Color', 'r');


line([log10(0.05) log10(0.05)], ylim, 'LineWidth', 2);


lt_plot(log10(Pval_base_All(inds)), log10(Pval_learn_All(inds)));




% === 3) ff diff [dot-dot]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('pitch diff (nontargs:bk, targs:yel)');
xlabel('mean ff diff (baseline)');
ylabel('mean ff diff (WN)');

% non-targets
inds=Target_All==0;
lt_plot(FFdiff_base_All(inds), FFdiff_learn_All(inds));

% targets
inds=Target_All==1;
lt_plot(FFdiff_base_All(inds), FFdiff_learn_All(inds), {'Color','y'});
        
if plot_text==1;
    for i=1:length(BirdName_All);
        birdname=BirdName_All{i};
        expt=ExptName_All{i};
        string=[birdname(1:4) '-' expt(end-5:end) '-' SylName_All{i}];
        text(FFdiff_base_All(i), FFdiff_learn_All(i), string, 'Color','r');
    end
end

% === 4) pval [1-2]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('pval(ranksum) (nontargs:bk, targs:red)');
xlabel('baseline - WN');
ylabel('log10(pval) (ranksum)');

% nontargs
inds=Target_All==0;
Y=[log10(Pval_base_All(inds))' log10(Pval_learn_All(inds))'];
X=[1 2];

plot(X, Y, 'o-', 'Color','k');

% targs
inds=Target_All==1;
Y=[log10(Pval_base_All(inds))' log10(Pval_learn_All(inds))'];
X=[1 2];

plot(X, Y, 'o-', 'Color','r');

xlim([0 3]);

% ==== 5) d-prime (abs) [1 2];
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('dprime(absolute) (nontargs:bk, targs:red) [means]');
xlabel('baseline - WN');
ylabel('abs(dprime)');

% nontargs
inds=Target_All==0;
Y=[abs(Dprime_base_All(inds))' abs(Dprime_learning_All(inds))'];
X=[1 2];

plot(X, Y, 'o-', 'Color','k');
errorbar([2.5 2.8], mean(Y), lt_sem(Y), 's-k', 'MarkerSize',9,'MarkerFaceColor','k');

p=signrank(Y(:,1), Y(:,2), 'tail', 'left');
htext=lt_plot_pvalue(p);

% targs
inds=Target_All==1;
Y=[abs(Dprime_base_All(inds))' abs(Dprime_learning_All(inds))'];
X=[1 2];

plot(X, Y, 'o-', 'Color','r');
errorbar([2.5 2.8], mean(Y), lt_sem(Y), 's-r', 'MarkerSize',9,'MarkerFaceColor','r');

p=signrank(Y(:,1), Y(:,2), 'tail', 'left');
htext=lt_plot_pvalue(p);
set(htext, 'Color','r');

xlim([0 3.5]);


% ==== 6) d-prime [1 2];
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('dprime (nontargs:bk, targs:red)');
xlabel('baseline - WN');
ylabel('dprime');

% nontargs
inds=Target_All==0;
Y=[Dprime_base_All(inds)' Dprime_learning_All(inds)'];
X=[1 2];

plot(X, Y, 'o-', 'Color','k');
errorbar([2.5 2.8], mean(Y), lt_sem(Y), 's-k', 'MarkerSize',9,'MarkerFaceColor','k');


% targs
inds=Target_All==1;
Y=[Dprime_base_All(inds)' Dprime_learning_All(inds)'];
X=[1 2];

plot(X, Y, 'o-', 'Color','r');
errorbar([2.5 2.8], mean(Y), lt_sem(Y), 's-r', 'MarkerSize',9,'MarkerFaceColor','r');

xlim([0 3.5]);


% ===== 6) dprime (distributions, not paired);
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('dprime histograms');
xlabel('baseline - WN');
ylabel('dprime');

% nontargs
inds=Target_All==0;
Y=[Dprime_base_All(inds)' Dprime_learning_All(inds)'];
X=[1 2];

distributionPlot(Y,'histOpt',0, 'showMM', 6);

% ===== 6) dprime (distributions, abs);
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('dprime(abs) histograms');
xlabel('baseline - WN');
ylabel('dprime');

% nontargs
inds=Target_All==0;
Y=[abs(Dprime_base_All(inds))' abs(Dprime_learning_All(inds))'];
X=[1 2];

distributionPlot(Y,'histOpt',0, 'showMM', 6);

% ===== 6) mean diff (distributions);
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('mean pitch diff');
xlabel('baseline - WN');
ylabel('pitch diff (late - early)');

% nontargs
inds=Target_All==0;
Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];
X=[1 2];

distributionPlot(Y,'histOpt',0, 'showMM', 6);

lt_plot_zeroline;

% ===== 6) mean diff absolute (distributions);
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('mean pitch diff (abs)');
xlabel('baseline - WN');
ylabel('pitch diff (late - early)');

% nontargs
inds=Target_All==0;
Y=[abs(FFdiff_base_All(inds))' abs(FFdiff_learn_All(inds))'];
X=[1 2];

distributionPlot(Y,'histOpt',0, 'showMM', 6);



% ===== 7) TEST NULL HYPOTHESIS - diff of abs mean pitch shift should be 0
% as a population
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('abs mean Pitch diff');
xlabel('Learning diff minus baseline diff (paired)');
ylabel('count');
grid on;

% nontargs
inds=Target_All==0;
y1=abs(FFdiff_base_All(inds));
y2=abs(FFdiff_learn_All(inds));

Y=y2-y1;

Y=y2-y1;
lt_plot_histogram(Y,'',1,1)
lt_plot_cdf(Y)

% ===== 8) TEST NULL HYPOTHESIS - diff of abs mean pitch shift should be 0
% as a population
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('abs dprime');
xlabel('Learning dprime minus baseline dprime (paired)');
ylabel('count');
grid on;

% nontargs
inds=Target_All==0;
y1=abs(Dprime_base_All(inds));
y2=abs(Dprime_learning_All(inds));

Y=y2-y1;
lt_plot_histogram(Y,'',1,1)
lt_plot_cdf(Y)

       
%% ======= 9a) HISTOGRAM
COLLECTED_STATS=struct;

% nontargs
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
lt_plot_text(0,0.5, 'upcoming plots: ALL NONTARGS');

inds=Target_All==0;
Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];

OUTPUT=struct;
lt_seq_dep_pitch_ACROSSBIRDS_DRIFT3_sub1
COLLECTED_STATS.nontargs=OUTPUT;


% % --- SIMILAR
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% lt_plot_text(0,0.5, 'upcoming plots: SIMILAR');
% 
% inds=Target_All==0 & Similar_All==1;
% Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];
% 
% OUTPUT=struct;
% lt_seq_dep_pitch_ACROSSBIRDS_DRIFT3_sub1
% COLLECTED_STATS.similar=OUTPUT;


% % --- DIFF
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% lt_plot_text(0,0.5, 'upcoming plots: DIFFERENT');
% 
% inds=Target_All==0 & Similar_All==0;
% Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];
% 
% OUTPUT=struct;
% lt_seq_dep_pitch_ACROSSBIRDS_DRIFT3_sub1
% COLLECTED_STATS.different=OUTPUT;


% % --- PRE SIMILAR
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% lt_plot_text(0,0.5, 'upcoming plots: PRESYL SIMILAR');
% 
% inds=Target_All==0 & PreSimilar_All==1;
% Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];
% 
% OUTPUT=struct;
% lt_seq_dep_pitch_ACROSSBIRDS_DRIFT3_sub1
% COLLECTED_STATS.pre_sim=OUTPUT;
% 

% --- PRE SIMILAR [SIMILAR]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
lt_plot_text(0,0.5, 'upcoming plots: PRESIM [SIMILAR]');

inds=Target_All==0 & PreSimilar_All==1 & Similar_All==1;
Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];

OUTPUT=struct;
lt_seq_dep_pitch_ACROSSBIRDS_DRIFT3_sub1
COLLECTED_STATS.similar_presim=OUTPUT;


% --- PRE DIFF [SIMILAR]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
lt_plot_text(0,0.5, 'upcoming plots: PREDIFF [SIMILAR]');

inds=Target_All==0 & PreSimilar_All==0 & Similar_All==1;
Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];

OUTPUT=struct;
lt_seq_dep_pitch_ACROSSBIRDS_DRIFT3_sub1
COLLECTED_STATS.similar_prediff=OUTPUT;


% --- DIFF [presim]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
lt_plot_text(0,0.5, 'upcoming plots: PRESIM [DIFFERENT]');

inds=Target_All==0 & PreSimilar_All==1 & Similar_All==0;
Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];

OUTPUT=struct;
lt_seq_dep_pitch_ACROSSBIRDS_DRIFT3_sub1
COLLECTED_STATS.different_presim=OUTPUT;

% --- DIFF [prediff]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
lt_plot_text(0,0.5, 'upcoming plots: PREDIFF [DIFFERENT]');

inds=Target_All==0 & PreSimilar_All==0 & Similar_All==0;
Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];

OUTPUT=struct;
lt_seq_dep_pitch_ACROSSBIRDS_DRIFT3_sub1
COLLECTED_STATS.different_prediff=OUTPUT;


% 
% % --- PRE SIMILAR
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% lt_plot_text(0,0.5, 'upcoming plots: PRESYL DIFF');
% 
% inds=Target_All==0 & PreSimilar_All==0;
% Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];
% 
% OUTPUT=struct;
% lt_seq_dep_pitch_ACROSSBIRDS_DRIFT3_sub1
% COLLECTED_STATS.pre_diff=OUTPUT;
% 



%% =============== SUMMARY PLOT OF ALL EFFECT SIZES/pvalues
lt_figure; hold on;
syl_type_fieldnames=fieldnames(COLLECTED_STATS);
EffectsToPlot_List={'abs_val','positive_vs_basematched','negative_vs_basematched'};
ylabel('Learning minus baseline (hz)');
title('Learning pitch shift minus duration-matched baseline drift (mean/SEM)');

Y_all=[];
Ysem_all=[];
Yp_all=[];
Y_p_twotail=[];

Xlabels={};
cc=1;
X=[];

for i=1:length(syl_type_fieldnames);
    syl_type_field=syl_type_fieldnames{i};
    
    
    for ii=1:length(EffectsToPlot_List);
        effect_type_field=EffectsToPlot_List{ii};
        
        if isfield(COLLECTED_STATS.(syl_type_field).Paired_analyses, effect_type_field);
            
            FFdiff_mean=COLLECTED_STATS.(syl_type_field).Paired_analyses.(effect_type_field).ff_diff_mean;
            FFdiff_sem=COLLECTED_STATS.(syl_type_field).Paired_analyses.(effect_type_field).ff_diff_sem;
            
            Y_all=[Y_all FFdiff_mean];
            Ysem_all=[Ysem_all FFdiff_sem];
            
            p_twotail=COLLECTED_STATS.(syl_type_field).Paired_analyses.(effect_type_field).ff_diff_p_notail;
            Y_p_twotail=[Y_p_twotail p_twotail];
            
        else
            % then doesn't have data for this type, e.g. no negative
            % generalizers
            Y_all=[Y_all nan];
            Ysem_all=[Ysem_all nan];
            Y_p_twotail=[Y_p_twotail nan];
        end
        
        tmpfield=effect_type_field;
        if strcmp(effect_type_field, 'positive_vs_basematched');
            tmpfield='Positive';
        elseif strcmp(effect_type_field, 'negative_vs_basematched');
            tmpfield='Negative';
        elseif strcmp(effect_type_field, 'abs_val');
            tmpfield='AbsVal';
        end
            
            Xlabels=[Xlabels [upper(syl_type_field) '-' tmpfield]];
        
        X=[X cc];
        cc=cc+1;
        
    end
    
    cc=cc+1; % a gap between syl types
    
end

lt_plot_bar(X, Y_all, {'Errors', Ysem_all});
Ylim=ylim;
for i=1:length(Y_p_twotail);
    if Y_p_twotail(i)<0.05;
        color='r';
    else
        color='k';
    end
    
    lt_plot_text(X(i)-0.3, Y_all(i)+10, [num2str(Y_p_twotail(i), '%3.2g')], color); % plot pvalue
    
end
lt_plot_text(X(1), Ylim(2)-10, 'two-tailed sign rank p values', 'b');



set(gca, 'XTick', X);
set(gca, 'XTickLabel', Xlabels);
rotateXLabels(gca, 45)



%% QUANTIFY WHAT PERCENT OF SYLS SHOW POSITIVE AND NEG GENERALIZATION

% --- all nontargets
sign_desired=-1; % 1 or -1
inds=Target_All==0;
Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];
cycles=1000;

lt_seq_dep_pitch_ACROSSBIRDS_DRIFT3_sub2

% --- similar
sign_desired=1; % 1 or -1
inds=Target_All==0 & Similar_All==1;
Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];
cycles=1000;

lt_seq_dep_pitch_ACROSSBIRDS_DRIFT3_sub2

% --- different
sign_desired=-1; % 1 or -1
inds=Target_All==0 & Similar_All==0;
Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];
cycles=1000;

lt_seq_dep_pitch_ACROSSBIRDS_DRIFT3_sub2

        
