function [SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_DRIFT2(SeqDepPitch_AcrossBirds, PARAMS)
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
        
        if NumWNDays+1>NumBaselineDays
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
plot_text=1; % birdname, expt, syl

% figure;
count=1;
SubplotsPerFig=6;
subplotrows=2;
subplotcols=3;
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

       
%% 


        
        
        
        
        
%         
%         for j=1:length(syls_unique);
%             syl=syls_unique{j};
%             
%             % syl count
%             syl_count=syl_count+1;
%             
%             % =============== STATS (ACROSS ALL BASELINE DAYS)
%             Baseline_Mean_UsingVals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF;
%             Baseline_STD_UsingVals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).stdFF;
%             
%             
%             % ====================== STATS (EACH BASELINE DAY VS OTHER
%             % BASELINE DAYS)
%             Deviation_VsOtherDaysMean_ff=[];
%             Deviation_VsOtherDaysMean_z=[];
%             BaselineDayMeans=[];
%             BaselineDaySTD=[];
%             BaselineDaySEM=[];
%             Pvals_Ranksum=[];
%             for k=BaselineDayInds;
%                 
%                 % skip today if lacking data
%                 if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})
%                     continue
%                 end
%                 
%                 % ----- 1) GET TODAY'S DATA
%                 ffvals_today=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k});
%                 
%                 % ---- 2) GET EMPIRICAL PROBABILITY DISTRIBUTION OF ALL OTHER
%                 % DAYS
%                 % -- go through all the other baseline days + learning day and collect
%                 % ffvals
%                 ffvals_otherDays=[];
%                 for kk=BaselineDayInds;
%                     if kk==k;
%                         % skip today
%                         continue;
%                     end
%                     
%                     % skip today if lacking data
%                     if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{kk})
%                         continue
%                     end
%                     
%                     ffvals_otherDays=[ffvals_otherDays cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{kk})];
%                 end
%                 
%                 % +++++++++++ STATS FOR TODAY
%                 p_ranksum=ranksum(ffvals_today, ffvals_otherDays);
%                 Pvals_Ranksum=[Pvals_Ranksum p_ranksum];
%                 Deviation_VsOtherDaysMean_ff=[Deviation_VsOtherDaysMean_ff mean(ffvals_today)-mean(ffvals_otherDays)];
%                 Deviation_VsOtherDaysMean_z=[Deviation_VsOtherDaysMean_z (mean(ffvals_today)-mean(ffvals_otherDays))/std(ffvals_otherDays)];
%                 BaselineDayMeans=[BaselineDayMeans mean(ffvals_today)];
%                 BaselineDaySTD=[BaselineDaySTD std(ffvals_today)];
%                 BaselineDaySEM=[BaselineDaySEM lt_sem(ffvals_today)];
%             end
%             
%             
%             % ============== STATS (ONE FOR EACH SYL, SUMMARIZING BASELINE
%             % DAYS)
%             NumBaselineDays=length(BaselineDayMeans);
%             
%             % Day means as zscores
%            BaselineDayMeans_zscore=(BaselineDayMeans-Baseline_Mean_UsingVals)./Baseline_STD_UsingVals;
%             
%             
%             
%             
%             % --- FF
%             MaxDev_VsOtherDay_pos_FF=max(Deviation_VsOtherDaysMean_ff);
%             MaxDev_VsOtherDay_neg_FF=min(Deviation_VsOtherDaysMean_ff);
%             MaxDev_VsOtherDay_abs_FF=max(abs(Deviation_VsOtherDaysMean_ff));
%             MeanDev_VsOtherDay_abs_FF=mean(abs(Deviation_VsOtherDaysMean_ff));
%             
%             % zscore
%             MaxDev_VsOtherDay_pos_Z=max(Deviation_VsOtherDaysMean_z);
%             MaxDev_VsOtherDay_neg_Z=min(Deviation_VsOtherDaysMean_z);
%             MaxDev_VsOtherDay_abs_Z=max(abs(Deviation_VsOtherDaysMean_z));
%             MeanDev_VsOtherDay_abs_Z=mean(abs(Deviation_VsOtherDaysMean_z));
%             
%             % std over days (deviations)
%             STD_DevOverDays_ff=std(Deviation_VsOtherDaysMean_ff);
%             STD_DevOverDays_z=std(Deviation_VsOtherDaysMean_z);
%             
%             % std of day means
%             STD_DayMeans_ff=std(BaselineDayMeans);
%             STD_DayMeans_z=std(BaselineDayMeans_zscore);
%             
%             
%             
%             % =================== STORE VALS FOR THIS SYL
% %             try numel(BaselineStats.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_pos_FF);
% %             catch err % field does not exist yet. make it
% %                             BaselineStats.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_pos_FF(1)=MaxDev_VsOtherDay_pos_FF;
% %             BaselineStats.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_neg_FF(1)=MaxDev_VsOtherDay_neg_FF;
% %             BaselineStats.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_abs_FF(1)=MaxDev_VsOtherDay_abs_FF;
% %             BaselineStats.DEVIATIONS_EachDayVsOtherDays.MeanDev_VsOtherDay_abs_FF(1)=MeanDev_VsOtherDay_abs_FF;
% %             BaselineStats.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_pos_Z(1)=MaxDev_VsOtherDay_pos_Z;
% %             BaselineStats.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_neg_Z(1)=MaxDev_VsOtherDay_neg_Z;
% %              BaselineStats.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_abs_Z(1)=MaxDev_VsOtherDay_abs_Z;
% %               BaselineStats.DEVIATIONS_EachDayVsOtherDays.MeanDev_VsOtherDay_abs_Z(1)=MeanDev_VsOtherDay_abs_Z;
% %         
% %             BaselineStats.NumBaselineDays(syl_count)=NumBaselineDays;
% %             
% %               BaselineStats.BaselineDayMeans{syl_count}=BaselineDayMeans;
% %               BaselineStats.BaselineDaySTD{syl_count}=BaselineDaySTD;
% %               BaselineStats.BaselineDaySEM{syl_count}=BaselineDaySEM;
% %               BaselineStats.Deviation_VsOtherDaysMean_ff{syl_count}=Deviation_VsOtherDaysMean_ff;
% %               BaselineStats.Deviation_VsOtherDaysMean_z{syl_count}=Deviation_VsOtherDaysMean_z;
% % 
% %                 
% %             end
%             BaselineStats.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_pos_FF(syl_count)=MaxDev_VsOtherDay_pos_FF;
%             BaselineStats.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_neg_FF(syl_count)=MaxDev_VsOtherDay_neg_FF;
%             BaselineStats.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_abs_FF(syl_count)=MaxDev_VsOtherDay_abs_FF;
%             BaselineStats.DEVIATIONS_EachDayVsOtherDays.MeanDev_VsOtherDay_abs_FF(syl_count)=MeanDev_VsOtherDay_abs_FF;
%             BaselineStats.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_pos_Z(syl_count)=MaxDev_VsOtherDay_pos_Z;
%             BaselineStats.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_neg_Z(syl_count)=MaxDev_VsOtherDay_neg_Z;
%              BaselineStats.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_abs_Z(syl_count)=MaxDev_VsOtherDay_abs_Z;
%               BaselineStats.DEVIATIONS_EachDayVsOtherDays.MeanDev_VsOtherDay_abs_Z(syl_count)=MeanDev_VsOtherDay_abs_Z;
%         
%             BaselineStats.NumBaselineDays(syl_count)=NumBaselineDays;
%             
%               BaselineStats.BaselineDayMeans{syl_count}=BaselineDayMeans;
%               BaselineStats.BaselineDayMeans_zscore{syl_count}=BaselineDayMeans_zscore;
%               
%               BaselineStats.BaselineDaySTD{syl_count}=BaselineDaySTD;
%               BaselineStats.BaselineDaySEM{syl_count}=BaselineDaySEM;
%               BaselineStats.Deviation_VsOtherDaysMean_ff{syl_count}=Deviation_VsOtherDaysMean_ff;
%               BaselineStats.Deviation_VsOtherDaysMean_z{syl_count}=Deviation_VsOtherDaysMean_z;
%              BaselineStats.Pvals_ranksum_EachDayVsOthers{syl_count}=Pvals_Ranksum;
%               
%               BaselineStats.STD_DevOverDays_ff(syl_count)=STD_DevOverDays_ff;
%               BaselineStats.STD_DevOverDays_z(syl_count)=STD_DevOverDays_z;
%               BaselineStats.STD_DayMeans_ff(syl_count)=STD_DayMeans_ff;
%               BaselineStats.STD_DayMeans_z(syl_count)=STD_DayMeans_z;
%               
%         end
%     end
% end
% 
% SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS=BaselineStats;


% %% FOR ALL EXPERIMENTS, PLOT LINES OF LEARNING FOR ALL SYLLABLES, OVERLAYING BASELINE WITH SAME DURATION DURING LEARNING 
% % [ALSO OVERLAY DURATION UP TO THE DAYS USED TO GET LEARNING METRIC
% 
% 
% 
% %% FOR ALL SYLS IN ALL EXPERIMENTS PLOT 1) HISTOGRAM OF DAY TO DAY CHANGES AND 2) HISTOGRAM OF ALL BASELINE RENDITIONS
% 
% if (0) % skipping for now, a lot fo information
% norm_pdf=1;
% 
% for i=1:NumBirds;
%     birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
%     numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%     
%     for ii=1:numexperiments;
%         exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
%         
%         count=1;
%         SubplotsPerFig=20;
%         subplotrows=4;
%         subplotcols=5;
%         fignums_alreadyused=[];
%         hfigs=[];
% 
% 
%         syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
%         
%         for j=1:length(syls_unique);
%             syl=syls_unique{j};
%             
%             
%             % ===== GET VECTOR OF ALL FF DURING LEARNING (DAYS NUM SAME AS
%             % BASELINE]
%             numdays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.BaselineDays);
%             firstWNday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
%             Days=firstWNday:firstWNday+numdays-1;
%             
%             ffvals=[];
%             tvals=[];
%             for k=Days;
%                 try % because might not have enough WN days to match baseline days
%                 ffvals=[ffvals cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})];
%                 tvals=[tvals cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals{k})];
%                 catch err
%                 end
%             end
%             
%             % ---- Put into output structure
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).DRIFT.DuringLearning.SameNumDaysAsBaseline.ffvals=ffvals;
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).DRIFT.DuringLearning.SameNumDaysAsBaseline.tvals=tvals;
%             
%             
%             % ========== GET VECTOR OF LEARNING METRIC DAYS
%             firstday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolStartInd;
%             numdays=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.Consol_DayBins;
%             Days=firstday:firstday+numdays-1;
%             
%             ffvals=[];
%             tvals=[];
%             for k=Days;
%                 ffvals=[ffvals cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})];
%                 tvals=[tvals cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals{k})];
%             end
%             
%             % ---- Put into output structure
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).DRIFT.DuringLearning.ConsolEarlyDays.ffvals=ffvals;
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).DRIFT.DuringLearning.ConsolEarlyDays.tvals=tvals;
%             
%             
%             
%             % ==== GET VECTOR OF BASELINE FF AND TIMES
%             tvals_base=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).Tvals;
%             ffvals_base=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF;
%             
%             
%             % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%             % PLOTS, THIS SYL
%             
%             % =========================
%             % PLOT BASELINE (one subplot for each syl)
%             [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
%             title(syl);
%             
%             [nffvals, xcenters]=hist(ffvals_base);
%             
%             % normalize to get pdf
%             if norm_pdf==1;
%             nffvals=nffvals./sum(nffvals);
%             end
%             
%             % plot   
% %             lt_plot_bar(xcenters-2, nffvals, {'Color', 'b'});
%             lt_plot_area(xcenters, nffvals, 'b', 0.3)
%             
%             % ============ PLOT FFVALS OVER LEARNING (SAME NUM DAYS AS
%             % BASELINE)
%             ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).DRIFT.DuringLearning.SameNumDaysAsBaseline.ffvals;            
%             nffvals=hist(ffvals, xcenters);
%             
%             % normalize to get pdf
%             if norm_pdf==1;
%             nffvals=nffvals./sum(nffvals);
%             end
%             
%             % plot            
% %             lt_plot_bar(xcenters, nffvals, {'Color', 'r'});
%             lt_plot_area(xcenters, nffvals, 'r', 0.3)
%             hold on;
%            
%             % ============= PLOT FFVALS CONSOLIDATION EARLY
%             ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).DRIFT.DuringLearning.ConsolEarlyDays.ffvals;
%             
%             nffvals=hist(ffvals, xcenters);
%             
%             % normalize to get pdf
%             if norm_pdf==1;
%             nffvals=nffvals./sum(nffvals);
%             end
%             
%             % plot            
% %             lt_plot_bar(xcenters+2, nffvals, {'Color', 'k'});
%             lt_plot_area(xcenters, nffvals, 'k', 0.3)
%             hold on;
%             
%             
%             
%         end
%         
%         lt_subtitle(['baseline(blue); learning(red, days equal baseline); consol early(black); ' birdname '-' exptname]);
%     end
% end
% 
% end
% 
% 
% 
% 
% 
% %% PLOT ALL LEARNING DISTRIBUTION, COMPARISON TO BASELINE MAGNITUDE
% % 1) Plot one for each bird (across all experiments), and 2) Plot one
% % across all birds
% 
% % FIRST, MAGNITUDE OF LEARNING (ZSCORE)
% % 1) PDF, all syls
% count=1;
% SubplotsPerFig=9;
% subplotrows=3;
% subplotcols=3;
% fignums_alreadyused=[];
% hfigs=[];
% 
% %  PDF learning
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('all syls and targs');
% xlabel('learning, targ dir');
% Learning_All=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All;
% 
% % nontargets
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% Y=Learning_All(inds);
% [~, Xcenters]=lt_plot_histogram(Y);
% 
% % targets
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1;
% Y=Learning_All(inds);
% [~,~,hbar]=lt_plot_histogram(Y,Xcenters);
% set(hbar, 'FaceColor','y');
% 
% % 2) Separate similar and diff
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('similar(blue), diff(red), targ(yel)');
% xlabel('learning, targ dir');
% 
% % nontargets/similar
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
% Y=Learning_All(inds);
% [~,~,hbar]=lt_plot_histogram(Y,Xcenters);
% set(hbar, 'FaceColor','b');
% 
% % nontargets/different
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
% Y=Learning_All(inds);
% [~,~,hbar]=lt_plot_histogram(Y,Xcenters);
% set(hbar, 'FaceColor','r');
% 
% % targets
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1;
% Y=Learning_All(inds);
% [~,~,hbar]=lt_plot_histogram(Y,Xcenters)
% set(hbar, 'FaceColor','y');
% 
% 
% % 3) cdf learning all syls, similar, diff, targets
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('sim, diff, all, targs');
% xlabel('learning (zscore, direction of targ)');
% grid on;
% 
% % all nontargs
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% Y=Learning_All(inds);
% 
% [F, X]=ecdf(Y);
% plot(X, F, 'k-', 'LineWidth',2);
% 
% % nontargs/similar
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
% Y=Learning_All(inds);
% 
% [F, X]=ecdf(Y);
% plot(X, F, 'b-', 'LineWidth',2);
% 
% % nontargs/different
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
% Y=Learning_All(inds);
% 
% [F, X]=ecdf(Y);
% plot(X, F, 'r-', 'LineWidth',2);
% 
% % targets
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1;
% Y=Learning_All(inds);
% 
% [F, X]=ecdf(Y);
% plot(X, F, 'g-', 'LineWidth',2);
% 
% 
% % 3) histogram of various baseline measures
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('baseline, [abs deviation (blue=max), (red=mean) ] (gr: std)');
% 
% % max abs deviation for all syls (mirror around zero)
% MaxBaselineDev_z=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_abs_Z;
% [~,~,hbar]=lt_plot_histogram([MaxBaselineDev_z -MaxBaselineDev_z], Xcenters);
% set(hbar, 'FaceColor','b');
% 
% % mean abs dev
% MeanBaselineDev_z=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.DEVIATIONS_EachDayVsOtherDays.MeanDev_VsOtherDay_abs_Z;
% [~,~,hbar]=lt_plot_histogram([MeanBaselineDev_z -MeanBaselineDev_z], Xcenters);
% set(hbar, 'FaceColor','r');
% 
% % std over days
% Y=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.STD_DevOverDays_z;
% [~,~,hbar]=lt_plot_histogram([Y -Y], Xcenters);
% set(hbar, 'FaceColor','g');
% 
% % distribution of all baseline days (deviations)
% Y=[];
% for i=1:length(SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.Deviation_VsOtherDaysMean_z);
%     Y=[Y SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.Deviation_VsOtherDaysMean_z{i}];
% end
% [~,~,hbar]=lt_plot_histogram(Y, Xcenters);
% set(hbar, 'FaceColor','y');
% 
% % 4) Overlay learning historgram (each syl) and all baseline days (one set
% % for each syl).
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('shift by nontargs and all baseline days (vs. other baseline days)');
% xlabel('zscore');
% 
% % -- LEARNING
% Learning_All=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All;
% % nontargets
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% Y=Learning_All(inds);
% [~, Xcenters, hbar]=lt_plot_histogram(Y, '', 1, 1);
% set(hbar, 'FaceColor', 'r')
% % --- BASELINE DAYS
% Y=[];
% for i=1:length(SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.Deviation_VsOtherDaysMean_z);
%     Y=[Y SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.Deviation_VsOtherDaysMean_z{i}];
% end
% [~,~,hbar]=lt_plot_histogram(Y, Xcenters, 1,1);
% set(hbar, 'FaceColor','k');
% 
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('shift by nontargs and all baseline days (vs. same baseline mean)');
% xlabel('zscore');
% 
% % -- LEARNING
% Learning_All=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All;
% % nontargets
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% Y=Learning_All(inds);
% [~, Xcenters, hbar]=lt_plot_histogram(Y, '', 1, 1);
% set(hbar, 'FaceColor', 'r')
% % --- BASELINE DAYS
% Y=[];
% for i=1:length(SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.BaselineDayMeans_zscore);
%     Y=[Y SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.BaselineDayMeans_zscore{i}];
% end
% [~,~,hbar]=lt_plot_histogram(Y, Xcenters, 1,1);
% set(hbar, 'FaceColor','k');
% 
% 
% 
% 
% 
% % 4) cdf of all + various baseline measures
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('black: nontargets; bl: max base abs dev, cyan: all baseline days');
% xlabel('learning (zscore, direction of targ');
% grid on;
% 
% % all nontargs
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% Y=Learning_All(inds);
% [F, X]=ecdf(Y);
% plot(X, F, 'k-', 'LineWidth',2);
% 
% % baseline stats (using max z dev)
% MaxBaselineDev_z=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_abs_Z;
% Y=[MaxBaselineDev_z -MaxBaselineDev_z];
% [F, X]=ecdf(Y);
% plot(X, F, 'b-', 'LineWidth',2);
% % draw a line for 95%iles 
% X=prctile(Y, [2.5 97.5]);
% line([X(1) X(1)], ylim, 'Color','b','LineStyle','--');
% line([X(2) X(2)], ylim, 'Color','b','LineStyle','--');
% 
% % baseline stats (using mean z dev)
% % MeanBaselineDev_z=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.DEVIATIONS_EachDayVsOtherDays.MeanDev_VsOtherDay_abs_Z;
% % Y=[MeanBaselineDev_z -MeanBaselineDev_z];
% % [F, X]=ecdf(Y);
% % plot(X, F, 'r-', 'LineWidth',2);
% % % draw a line for 95%iles 
% % % draw a line for 95%iles 
% % X=prctile(Y, [2.5 97.5]);
% % line([X(1) X(1)], ylim, 'Color','r','LineStyle','--');
% % line([X(2) X(2)], ylim, 'Color','r','LineStyle','--');
% 
% % % against std of baseline days
% % Y=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.STD_DevOverDays_z;
% % Y=[Y -Y];
% % [F, X]=ecdf(Y);
% % plot(X, F, 'g-', 'LineWidth',2);
% % % draw a line for 95%iles 
% % % draw a line for 95%iles 
% % X=prctile(Y, [2.5 97.5]);
% % line([X(1) X(1)], ylim, 'Color','g','LineStyle','--');
% % line([X(2) X(2)], ylim, 'Color','g','LineStyle','--');
% 
% % all baseline deviations concatenated
% Y=[];
% for i=1:length(SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.Deviation_VsOtherDaysMean_z);
%     Y=[Y SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.Deviation_VsOtherDaysMean_z{i}];
% end
% [F, X]=ecdf(Y);
% plot(X, F, 'c-', 'LineWidth',2);
% % draw a line for 95%iles 
% % draw a line for 95%iles 
% X=prctile(Y, [2.5 97.5]);
% line([X(1) X(1)], ylim, 'Color','c','LineStyle','--');
% line([X(2) X(2)], ylim, 'Color','c','LineStyle','--');
% 
% 
% 
% % --- 5) Normalize all learning for each syl vs. its own baseline
% % all nontargets, using baseline stats
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('All nontargets, (blue: vs. max abs dev) (gr: vs. 2*sigma)');
% xlabel('learning, norm to baselined');
% grid on;
% Learning_TargDir_All=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All;
% 
% % mean absolute dev (red)
% % inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% % Mean_abs_dev=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.DEVIATIONS_EachDayVsOtherDays.MeanDev_VsOtherDay_abs_Z;
% % 
% % Learning_norm_baseline=Learning_TargDir_All./Mean_abs_dev;
% % [~,Xcenters,hbar]=lt_plot_histogram(Learning_norm_baseline(inds));
% % set(hbar, 'FaceColor','r');
% 
% % max absolute dev (blue)
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% Max_Abs_Dev=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_pos_Z;
% 
% Learning_norm_baseline=Learning_TargDir_All./Max_Abs_Dev;
% [~,Xcenters,hbar]=lt_plot_histogram(Learning_norm_baseline(inds));
% set(hbar, 'FaceColor','b');
% 
% % 2*std (assuming gaussian, 95% within 2sigma)
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% STD_base=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.STD_DevOverDays_z;
% STD_base_double=2*STD_base;
% 
% Learning_norm_baseline=Learning_TargDir_All./STD_base_double;
% [~,~,hbar]=lt_plot_histogram(Learning_norm_baseline(inds),Xcenters);
% set(hbar, 'FaceColor','g');
% 
% 
% % dray line for 1
% line([-1 -1], ylim, 'Color','k')
% line([1 1], ylim, 'Color','k');
% 
% 
% 
% 
% % --- 6) cdf of learning norm to baseline (all)
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('all nontargets, (bl: max dev) (gr: 2*sigma_dev) (mag:2*sigma_means)');
% grid on;
% xlabel('learning, norm baseline');
% Learning_TargDir_All=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All;
% 
% % mean absolute dev (red)
% % inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% % Mean_abs_dev=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.DEVIATIONS_EachDayVsOtherDays.MeanDev_VsOtherDay_abs_Z;
% % 
% % Y=Learning_TargDir_All./Mean_abs_dev;
% % Y=Y(inds);
% % [F, X]=ecdf(Y);
% % plot(X, F, 'r-', 'LineWidth',3);
% 
% % max absolute dev (bl)
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% Learning_TargDir_All=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All;
% Max_Abs_Dev=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_pos_Z;
% 
% Y=Learning_TargDir_All./Max_Abs_Dev;
% Y=Y(inds);
% [F, X]=ecdf(Y);
% plot(X, F, 'b-', 'LineWidth',2);
% 
% % 2*std 
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% Learning_TargDir_All=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All;
% STD_base=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.STD_DevOverDays_z;
% STD_base_double=2*STD_base;
% 
% Y=Learning_TargDir_All./STD_base_double;
% Y=Y(inds);
% [F, X]=ecdf(Y);
% plot(X, F, 'g-', 'LineWidth',2);
% 
% 
% % 2* STD (of day means)
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% STD_vals=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.STD_DayMeans_z;
% STD_vals_double=2*STD_vals;
% 
% Y=Learning_TargDir_All./STD_vals_double;
% Y=Y(inds);
% [F, X]=ecdf(Y);
% plot(X, F, 'm-', 'LineWidth',2);
% 
% 
% % dray line for 1
% line([-1 -1], ylim, 'Color','k')
% line([1 1], ylim, 'Color','k');
% 
% 
% % --- 6) cdf of learning norm to baseline (similar)
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('similar, (bl: max dev) (gr: 2*sigma_dev) (mag:2*sigma_means)');
% grid on;
% xlabel('learning, norm baseline');
% Learning_TargDir_All=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All;
% 
% % mean absolute dev (red)
% % inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
% % Mean_abs_dev=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.DEVIATIONS_EachDayVsOtherDays.MeanDev_VsOtherDay_abs_Z;
% % 
% % Y=Learning_TargDir_All./Mean_abs_dev;
% % Y=Y(inds);
% % [F, X]=ecdf(Y);
% % plot(X, F, 'r-', 'LineWidth',3);
% 
% % max absolute dev (bl)
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
% Learning_TargDir_All=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All;
% Max_Abs_Dev=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_pos_Z;
% 
% Y=Learning_TargDir_All./Max_Abs_Dev;
% Y=Y(inds);
% [F, X]=ecdf(Y);
% plot(X, F, 'b-', 'LineWidth',2);
% 
% % 2*std 
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
% Learning_TargDir_All=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All;
% STD_base=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.STD_DevOverDays_z;
% STD_base_double=2*STD_base;
% 
% Y=Learning_TargDir_All./STD_base_double;
% Y=Y(inds);
% [F, X]=ecdf(Y);
% plot(X, F, 'g-', 'LineWidth',2);
% 
% % 2* STD (of day means)
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
% STD_vals=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.STD_DayMeans_z;
% STD_vals_double=2*STD_vals;
% 
% Y=Learning_TargDir_All./STD_vals_double;
% Y=Y(inds);
% [F, X]=ecdf(Y);
% plot(X, F, 'm-', 'LineWidth',2);
% 
% 
% % dray line for 1
% line([-1 -1], ylim, 'Color','k')
% line([1 1], ylim, 'Color','k');
% 
% 
% % --- 8) cdf of learning norm to baseline (different)
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('(different syls, (bl: max dev) (gr: 2*sigma) (mag:2*sigma_means)');
% grid on;
% xlabel('learning, norm baseline');
% Learning_TargDir_All=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All;
% 
% % mean absolute dev (red)
% % inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
% % Mean_abs_dev=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.DEVIATIONS_EachDayVsOtherDays.MeanDev_VsOtherDay_abs_Z;
% % 
% % Y=Learning_TargDir_All./Mean_abs_dev;
% % Y=Y(inds);
% % [F, X]=ecdf(Y);
% % plot(X, F, 'r-', 'LineWidth',3);
% 
% % max absolute dev (bl)
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
% Learning_TargDir_All=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All;
% Max_Abs_Dev=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.DEVIATIONS_EachDayVsOtherDays.MaxDev_VsOtherDay_pos_Z;
% 
% Y=Learning_TargDir_All./Max_Abs_Dev;
% Y=Y(inds);
% [F, X]=ecdf(Y);
% plot(X, F, 'b-', 'LineWidth',2);
% 
% % 2*std(of deviations)
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
% Learning_TargDir_All=SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All;
% STD_base=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.STD_DevOverDays_z;
% STD_base_double=2*STD_base;
% 
% Y=Learning_TargDir_All./STD_base_double;
% Y=Y(inds);
% [F, X]=ecdf(Y);
% plot(X, F, 'g-', 'LineWidth',2);
% 
% % 2* STD (of day means)
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
% STD_vals=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.STD_DayMeans_z;
% STD_vals_double=2*STD_vals;
% 
% Y=Learning_TargDir_All./STD_vals_double;
% Y=Y(inds);
% [F, X]=ecdf(Y);
% plot(X, F, 'm-', 'LineWidth',2);
% 
% 
% % dray line for 1
% line([-1 -1], ylim, 'Color','k')
% line([1 1], ylim, 'Color','k');
% 
% 
% %% ====== PLOT LEARNING NORM TARGET
% 
% count=1;
% SubplotsPerFig=6;
% subplotrows=2;
% subplotcols=3;
% fignums_alreadyused=[];
% hfigs=[];
% 
% 
% % ---- 9) LEARNING, NORM TARGET
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('all syls');
% xlabel('learning (norm target)');
% 
% 
% % all nontargets
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% Learning_reltarg=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
% [~,Xcenters,hbar]=lt_plot_histogram(Learning_reltarg(inds));
% % set(hbar, 'FaceColor','b');
% Ymean=mean(Learning_reltarg(inds));
% plot(Ymean,0, '^k', 'MarkerFaceColor','k','MarkerSize',9);
% [~, p]=ttest(Learning_reltarg(inds));
% lt_plot_pvalue(p);
% 
% 
% % ------ 10) learning, norm targ, similar and diff
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('similar');
% xlabel('learning (norm target)');
% 
% 
% % similar 
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
% Learning_reltarg=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
% [~,~,hbar]=lt_plot_histogram(Learning_reltarg(inds), Xcenters);
% set(hbar, 'FaceColor','b');
% 
% Ymean=mean(Learning_reltarg(inds));
% plot(Ymean,0, '^b', 'MarkerFaceColor','b','MarkerSize',9);
% [~, p]=ttest(Learning_reltarg(inds));
% lt_plot_pvalue(p);
% 
% 
% % ------ 10) learning, norm targ, similar and diff
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('different');
% xlabel('learning (norm target)');
% 
% % different 
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
% Learning_reltarg=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
% [~,~,hbar]=lt_plot_histogram(Learning_reltarg(inds), Xcenters);
% set(hbar, 'FaceColor','r');
% Ymean=mean(Learning_reltarg(inds));
% plot(Ymean,0, '^r', 'MarkerFaceColor','r','MarkerSize',9);
% [~, p]=ttest(Learning_reltarg(inds));
% lt_plot_pvalue(p);
% 
% 
% % ------ 11) Similar vs. different
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('similar vs. different');
% 
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
% Learning_reltarg=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
% Y1=Learning_reltarg(inds);
% 
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
% Learning_reltarg=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
% Y2=Learning_reltarg(inds);
% 
% distributionPlot({Y1, Y2}, 'xNames',{'similar','diff'}, 'showMM',4, 'addSpread',1);
% [h, p]=ttest2(Y1, Y2);
% lt_plot_pvalue(p);
% 
% line(xlim, [0 0]);
% 
% 
% 
% 
% % ---- cdf of learning, norm target
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('all, similar, and diff');
% xlabel('learning (norm target)');
% grid on
% 
% % all nontargs
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% Learning_reltarg=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
% 
% Y=Learning_reltarg(inds);
% [F, X]=ecdf(Y);
% plot(X, F, 'k-', 'LineWidth',2);
% 
% 
% % similar
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1;
% Learning_reltarg=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
% 
% Y=Learning_reltarg(inds);
% [F, X]=ecdf(Y);
% plot(X, F, 'b-', 'LineWidth',2);
% 
% % diff
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0;
% Learning_reltarg=SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all;
% 
% Y=Learning_reltarg(inds);
% [F, X]=ecdf(Y);
% plot(X, F, 'r-', 'LineWidth',2);
% 
% 
% 
% %% HOW MANY ARE SIGNIFICANTLY DIFFERENT FROM BASELINE? (ranksum: WN day vs. all baseline lumped)
% 
% Pvals_ranksum_all=[];
% Pvals_ttest_all=[];
% Pvals_ttest_greater_all=[];
% Pvals_ttest_less_all=[];
% 
% Pvals_ttest_baseline_all=[];
% 
% for i=1:length(SeqDepPitch_AcrossBirds.AllSyllables.AcousticDist_all);
%     
%     birdnum=SeqDepPitch_AcrossBirds.AllSyllables.BirdNum_all(i);
%     exptnum=SeqDepPitch_AcrossBirds.AllSyllables.ExptNum_all(i);
%     syl=SeqDepPitch_AcrossBirds.AllSyllables.Syl_all{i};
%     
%     % --- extract baseline Ffvals
%     baseline_FFvals=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF;
% 
%     % --- WN days FF vals
%     if PARAMS.global.learning_metric=='zscore';
%         WN_FFvals=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Data_ZSCORE.Epoch.TargPassLearnThresh.DATA.(syl).FFvals;
%     elseif PARAMS.global.learning_metric=='ff_consolstart';
%         WN_FFvals=SeqDepPitch_AcrossBirds.birds{birdnum}.experiment{exptnum}.Data_PostProcessed.(syl).ConsolStart_FFvals_minbase+mean(baseline_FFvals);
%     end
%    
%     % ====== ranksum
%     p=ranksum(baseline_FFvals, WN_FFvals);
%     
%     Pvals_ranksum_all=[Pvals_ranksum_all p];
%     
%     % ==== ttest
%     [~,p]=ttest2(baseline_FFvals, WN_FFvals);
%     
%     Pvals_ttest_all=[Pvals_ttest_all p];
%     
%     % -- One sided
%     % WN>baseline
%     [~,p]=ttest2(WN_FFvals, baseline_FFvals, 'Tail','right');
%     Pvals_ttest_greater_all=[Pvals_ttest_greater_all p];
%     
%     % WN<baseline
%     [~,p]=ttest2(WN_FFvals, baseline_FFvals, 'Tail','left');
%     Pvals_ttest_less_all=[Pvals_ttest_less_all p];
%     
% end
% 
% % === STORE IN STRUCTURE
% SeqDepPitch_AcrossBirds.AllSyllables.WN_vs_BASE.Pvals_ranksum_all=Pvals_ranksum_all;
% SeqDepPitch_AcrossBirds.AllSyllables.WN_vs_BASE.Pvals_ttest_all=Pvals_ttest_all;
% SeqDepPitch_AcrossBirds.AllSyllables.WN_vs_BASE.Pvals_ttest_greater_all=Pvals_ttest_greater_all;
% SeqDepPitch_AcrossBirds.AllSyllables.WN_vs_BASE.Pvals_ttest_less_all=Pvals_ttest_less_all;
% 
% 
% 
% % ==== PLOT 
% count=1;
% SubplotsPerFig=9;
% subplotrows=3;
% subplotcols=3;
% fignums_alreadyused=[];
% hfigs=[];
% 
% 
% % --1) distribution of rank sum/ttest p vals
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('WN days FF vs. all baseline');
% xlabel('log10(pvalue (rank sum(bk), ttest(rd)), BASELINEranksum(bl)');
% inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% 
% % rank sum
% lt_plot_cdf(log10(Pvals_ranksum_all(inds)));
% % ttest
% lt_plot_cdf(log10(Pvals_ttest_all(inds)),'r');
% 
% line([log10(0.05) log10(0.05)], ylim);
% xlim([-10 0]);
% 
% % same, but for all baseline days
% Y=[];
% for i=1:length(SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.Pvals_ranksum_EachDayVsOthers);
%     
%     if inds(i)==1;
%     tmp=SeqDepPitch_AcrossBirds.AllSyllables.BASELINE_STATS.Pvals_ranksum_EachDayVsOthers{i};
%     Y=[Y tmp];
%     end
% end
% lt_plot_cdf(log10(Y),'b');
% 
% 
% 
% % --- 2) One sided - all syls
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('All nontargs');
% xlabel('log10(pvalue) (cyan: WN>base),(gr:WN<base)');
% inds = SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% 
% % WN>base
% lt_plot_cdf(log10(Pvals_ttest_greater_all(inds)), 'c');
% 
% % WN<base
% lt_plot_cdf(log10(Pvals_ttest_less_all(inds)), 'g');
% 
% line([log10(0.05) log10(0.05)], ylim);
% xlim([-10 0]);
% 
% 
% 
% 
% 
% % --- 2) one sided tests (similar)
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('Similar');
% xlabel('log10(pvalue) (cyan: WN>base),(gr:WN<base)');
% inds = SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% 
% % WN>base
% lt_plot_cdf(log10(Pvals_ttest_greater_all(inds)), 'c');
% 
% % WN>base
% lt_plot_cdf(log10(Pvals_ttest_less_all(inds)), 'g');
% 
% line([log10(0.05) log10(0.05)], ylim);
% xlim([-10 0]);
% 
% % --- 2) one sided tests (diff)
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('Different');
% xlabel('log10(pvalue) (cyan: WN>base),(gr:WN<base)');
% inds = SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
% 
% % WN>base
% lt_plot_cdf(log10(Pvals_ttest_greater_all(inds)), 'c');
% 
% % WN>base
% lt_plot_cdf(log10(Pvals_ttest_less_all(inds)), 'g');
% 
% line([log10(0.05) log10(0.05)], ylim);
% xlim([-10 0]);
% 
% 
% pause;
% 
% %% WHAT IS THE MAGNITUDE OF BASELINE DRIFT? (Linear regression, and extrapolate to the number of days used, and plot taht distribution vs. the actual learning distribution)
% 
% 
% %% FLUCTUATION (part 1) - HOW LIKEY IS MEASURED LEARNING GIVEN BASELINE DAY TO DAY FLUCTUATION? (NOT DRIFT, AS NOT TAKING INTO ACCOUNT ACROSS DAYS STRUCTURE)
% % for each day, have multiple datapoints. ask, given all the other baseline
% % datapoints + the WN day, likelighood of picking out those values in that day. do that
% % for all baseline days and the WN day
% 
% % NOTE: 8/13 - DONE A LOT, but still in progress -. Summary: have done:
% % lump all days (baseline + 1 WN day) and for each day for all values get
% % z-score rel to shuffled baseline (i.e. other days). Then have a
% % distribution for baseline (combine all baseline days) and for WN day.  I
% % have computed "effect size" which is WN day z-scored against the new
% % baseline distribution (of z-scores) and KS test (of Wn distribution vs
% % baseline combined. Problem with KS test is that many baseline days will
% % also be significantly different from combined baseline days.  Z-score
% % seems like a reasonable metric of learning.
% 
% temporary_test_on_baseline_day=0; % sanity check, if choose blien day as WN day, how many significnat?
% 
% for i=1:NumBirds;
%     
%     numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%     
%     for ii=1:numexperiments;
%         
%         syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
%         BaselineDayInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.BaselineDays;
%         % -- WN day
%         % NOTE: Use the last learning day, since should only use one
%         % day to be consistent with taking each baseline day separately
%         LearningDayInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DayInds;
%         LastLearningDayInd=LearningDayInds(end);
%         
%         if temporary_test_on_baseline_day==1;
%         LastLearningDayInd=1;
%         disp('WARNING - uising basekine day for WN, sanity cheeck');
%         end
%         
%         DRIFT=struct;
%         for j=1:length(syls_unique);
%             syl=syls_unique{j};
%             
%             % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%             % For each baseline day, ask how likely that day's
%             % distribution of data is given the other baseline days + the
%             % learning day
%             Zscores_base_all=[]; % will collect across all baseline days
%             
%             for k=BaselineDayInds;
%                 
%                 % skip today if lacking data
%                 if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})
%                     continue
%                 end
%                 
%                 % ----- 1) GET TODAY'S DATA
%                 ffvals_today=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k});
%                 
%                 % ---- 2) GET EMPIRICAL PROBABILITY DISTRIBUTION OF ALL OTHER
%                 % DAYS
%                 % -- go through all the other baseline days + learning day and collect
%                 % ffvals
%                 ffvals_otherDays=[];
%                 for kk=BaselineDayInds;
%                     if kk==k;
%                         % skip today
%                         continue;
%                     end
%                     
%                     % skip today if lacking data
%                     if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{kk})
%                         continue
%                     end
%                     
%                     ffvals_otherDays=[ffvals_otherDays cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{kk})];
%                 end
%                 
%                 % --- go to WN day and add on those vals
%                 ffvals_otherDays=[ffvals_otherDays cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{LastLearningDayInd})];
%                 
%                 
%                 
%                 % ======= CALCULATE DISTANCE DIFFERENT WAYS
%                 % 1) For each datapoint, get a z-score
%                 meanFF_otherdays=mean(ffvals_otherDays);
%                 stdFF_otherdays=std(ffvals_otherDays);
%                 
%                 zscores_today=(ffvals_today-meanFF_otherdays)./stdFF_otherdays;
%                                 
%                 % save zscores for each baseline day, to do null plotting
%                 DRIFT.Fluctuation.EachRendZscored.Baseline.Zscores_EachDay{k}=zscores_today;
% 
%                 % -------- COLLECT ACROSS BASELINE DAYS
%                 Zscores_base_all=[Zscores_base_all zscores_today];
%                 
%                 
% %                 2) For each datapoint, get a "likelyhood" by fitting a
% %                 gaussian to the "other data" empirical distribution and
% %                 estimating probability (then take the WN day mean and ask
% %                 whether the likelihood that of that mean is within the
% %                 distribution of likelihoods?
%                 
%                 
%                 
%                 % 3) Distance of probability distributions
%                 
%                 
%             end
%             
%             % ====== OUTPUT METRICS OF BASELINE VARIABILITY
%             DRIFT.Fluctuation.EachRendZscored.Baseline.Zscores_all=Zscores_base_all;
%             
%             
%             % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%             % For learning days, ask how likely the datapoints are given
%             % baseline distribution
%             
%             % ----- 1) GET TODAY'S DATA
%             ffvals_today=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{LastLearningDayInd});
%             
%             % ---- 2) GET BASELINE DISTRIBUTION
%             ffvals_Baseline=[];
%             for k=BaselineDayInds;
%                 % skip today if lacking data
%                 if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})
%                     continue
%                 end
%                 
%                 % ----- 1) GET TODAY'S DATA
%                 ffvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k});
%                 
%                 % --- 2) COMPILE ACROSS DAYS
%                 ffvals_Baseline=[ffvals_Baseline ffvals];
%             end
%             
%             % ======= CALCULATE DISTANCE DIFFERENT WAYS
%             % 1) For each datapoint, get a z-score
%             meanFF_baseline=mean(ffvals_Baseline);
%             stdFF_baseline=std(ffvals_Baseline);
%             
%             zscores_today=(ffvals_today-meanFF_baseline)./stdFF_baseline;
%             
%             
%             % 2) Distance of probability distributions
%             
%             % ====== OUTPUT 
%             DRIFT.Fluctuation.EachRendZscored.LearningDay.Zscores_all=zscores_today;
%             
%             % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%             % ======== OUTPUT FOR THIS SYL
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).Fluctuation=DRIFT.Fluctuation;
%             
%             
%         end
%     end
% end
% 
% 
% %% FLUCTUATION (part 2) - CALCUALATE SIGNIFICANCE, AND GET DISTANCE MEASURED USING THESE DISTRIBUTIONS
% 
% for i=1:NumBirds;
%     birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
%     numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%     
%     for ii=1:numexperiments;
%         exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
%         
%         count=1;
%         SubplotsPerFig=16;
%         subplotrows=4;
%         subplotcols=4;
%         fignums_alreadyused=[];
%         hfigs=[];
%         
%         syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
%         
%         pvals_KS_all=[];
%         effect_sizes_zscore_all=[];
%         for j=1:length(syls_unique);
%             syl=syls_unique{j};
%             
%             % --- note down if is target
%             if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
%                 targind=j;
%             end
%             
%             Baseline_zscores=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).Fluctuation.EachRendZscored.Baseline.Zscores_all;
%             WNday_zscores=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).Fluctuation.EachRendZscored.LearningDay.Zscores_all;
%             
%             % ===================== Use KS test to ask whether the two
%             % distributions are sign. diff
%             % diff
%             [h, p]=kstest2(Baseline_zscores, WNday_zscores);
%             
%             % aside, check whether sample sizes are reasonable
%             n1=numel(WNday_zscores);
%             n2=numel(Baseline_zscores);
%             
%             tmp=(n1*n2)/(n1+n2);
%             if tmp<4;
%                 disp(['K-S test sample size kinda low; ' birdname '-' exptname '-' syl]);
%             end
%             
%             % ---- SAVE p-val
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).Fluctuation.EachRendZscored.Comparing_LearningVsBaseline.KS_div_NotAbsVal.p=p;
%             pvals_KS_all=[pvals_KS_all p];
%             
%             % =========================== CALCULATE effect size (zscore of
%             % learning day mean rel baseline distrib);
%             baseline_mean=mean(Baseline_zscores);
%             baseline_std=std(Baseline_zscores);
%             
%             EffectSize_zscore_mean=(mean(WNday_zscores)-baseline_mean)/baseline_std;
%             EffectSize_zscore_sem=lt_sem((WNday_zscores-baseline_mean)./baseline_std);
%             
%             % ---- SAVE
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).Fluctuation.EachRendZscored.Comparing_LearningVsBaseline.EffectSize_zscore.mean=EffectSize_zscore_mean;
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).Fluctuation.EachRendZscored.Comparing_LearningVsBaseline.EffectSize_zscore.sem=EffectSize_zscore_sem;
%             effect_sizes_zscore_all=[effect_sizes_zscore_all EffectSize_zscore_mean];
%             
%             % ============== PLOT DISTRIBUTIONS OF zscores (not absolute val)
%             [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
%             title([syl '; p=' num2str(p)]);
%             
%             
%             [~, Xcenters, hbar]=lt_plot_histogram(Baseline_zscores);
%             set(hbar, 'FaceColor','w')
%             [~,~,hbar]=lt_plot_histogram(WNday_zscores);
%             
%             set(hbar, 'FaceColor', 'r');
%             
%             % plot effect size
%             lt_plot(EffectSize_zscore_mean, 0, {'Marker', '^'});
%             
%         end
%         
%         % ==== PLOT ACROSS ALL SYLS
%         % 1) ------------- p-values vs. effect size (including target)
%         [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
%         title('p-values(KS test) vs. (zscores, rel to baseline zscore distrib)');
%         ylabel('p val (log10)'); xlabel('zscore');
%         
%         lt_plot(effect_sizes_zscore_all, log10(pvals_KS_all), {'Color','r'});
%         
%         line(xlim, [log10(0.05) log10(0.05)], 'Color','b');
%         
%         % plot syl names
%         for j=1:length(syls_unique);
%             lt_plot_text(effect_sizes_zscore_all(j)+0.1, log10(pvals_KS_all(j)), syls_unique{j});
%         end
%         
%         % 2) ------------ p-values vs. effect size (excluding target)
%         [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
%         title('p-values(KS test) vs. (zscores, rel to baseline zscore distrib)');
%         ylabel('p val (log10)'); xlabel('zscore');
%         
%         % get inds excluding targ
%         inds=1:length(syls_unique);
%         inds(targind)=[];
%         
%         lt_plot(effect_sizes_zscore_all(inds), log10(pvals_KS_all(inds)), {'Color','r'});
%         
%         line(xlim, [log10(0.05) log10(0.05)], 'Color','b');
%         
%         % plot syl names
%         for j=1:length(syls_unique);
%             if j~=targind
%                 lt_plot_text(effect_sizes_zscore_all(j)+0.1, log10(pvals_KS_all(j)), syls_unique{j});
%             end
%         end
%         
%         
%         lt_subtitle([birdname '-' exptname '; baseline(bk) and WNday(red) z-scores rel to other days']);
%     end
% end
% 
% 
% % ======== PLOT ACROSS BIRDS
% pvals_KS_all=[];
% effect_sizes_zscore_all=[];
% TargStatus_all=[];
% SimilarStatus_all=[];
% 
% baseline_pvals=[];
% baseline_effectsizes=[];
% baseline_SimilarStatus_all=[];
% baseline_TargStatus_all=[];
% 
% 
% for i=1:NumBirds;
%     birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
%     numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%     
%     for ii=1:numexperiments;
%         exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
%         
%         syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
%         
%         targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
%         sign_of_targ_learning=sign(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean);
%         %
%         %         sign_of_targ_learning=sign(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Learning_by_target.WN_start_daybins); % to flip if targ is neg learning
%         
%         for j=1:length(syls_unique);
%             syl=syls_unique{j};
%             
%             % -- is target?
%             istarget=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
%             TargStatus_all=[TargStatus_all istarget];
%             
%             % -- is similar?
%             issimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
%             SimilarStatus_all=[SimilarStatus_all issimilar];
%             
%             % ---- Get pval
%             pval=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).Fluctuation.EachRendZscored.Comparing_LearningVsBaseline.KS_div_NotAbsVal.p;
%             pvals_KS_all=[pvals_KS_all pval];
%             
%             
%             % ---- Get effect size
%             effectsize=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).Fluctuation.EachRendZscored.Comparing_LearningVsBaseline.EffectSize_zscore.mean;
%             effectsize=sign_of_targ_learning*effectsize;
%             effect_sizes_zscore_all=[effect_sizes_zscore_all effectsize];
%             
%             % ---- "Null data" - For all baseline days, get effect size and
%             % a pval
%             Baseline_zscores=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).Fluctuation.EachRendZscored.Baseline.Zscores_all; % across all baseline days
%             baseline_mean=mean(Baseline_zscores);
%             baseline_std=std(Baseline_zscores);
%             
%             for k=1:length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).Fluctuation.EachRendZscored.Baseline.Zscores_EachDay); % for each baseline day
%                 
%                 if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).Fluctuation.EachRendZscored.Baseline.Zscores_EachDay{k});
%                     continue; 
%                 end
%                 
%                 zscores_today=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).Fluctuation.EachRendZscored.Baseline.Zscores_EachDay{k};
%                 
%                 % -- KS test
%                 [h, p]=kstest2(Baseline_zscores, zscores_today);
%                 
%                 SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).Fluctuation.EachRendZscored.Comparing_BaselineDayVsBaseline{k}.KS_div_NotAbsVal.p=p;
%                 
%                 
%                 % -- Effect size
%                 EffectSize_zscore_mean=(mean(zscores_today)-baseline_mean)/baseline_std;
%                 EffectSize_zscore_sem=lt_sem((zscores_today-baseline_mean)./baseline_std);
%                 
%                 SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).Fluctuation.EachRendZscored.Comparing_BaselineDayVsBaseline{k}.EffectSize_zscore.mean=EffectSize_zscore_mean;
%                 SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).Fluctuation.EachRendZscored.Comparing_BaselineDayVsBaseline{k}.EffectSize_zscore.sem=EffectSize_zscore_sem;
%                 
%                 
%                 % ==== Collect across all baseline days, all syls, all
%                 % expts, for plotting
%                 baseline_pvals=[baseline_pvals p];
%                 baseline_effectsizes=[baseline_effectsizes EffectSize_zscore_mean];
%                 baseline_SimilarStatus_all=[baseline_SimilarStatus_all issimilar];
%                 baseline_TargStatus_all=[baseline_TargStatus_all istarget];
%                 
%                 
%             end
%             
%         end
%     end
% end
% 
% 
% lt_figure; hold on;
% hsplot=[];
% % === ALL SYLS WN DAY
% hsplot(1)=lt_subplot(3,2,1); hold on
% title('all expts: pvals (KS) vs. effect sizes (zscore, accounting for baseline)');
% ylabel('pvals, log10');
% xlabel('effectsize (z)');
% 
% % similar
% inds=SimilarStatus_all==1;
% lt_plot(effect_sizes_zscore_all(inds), log10(pvals_KS_all(inds)), {'Color','b'});
% 
% % different
% inds=SimilarStatus_all==0;
% lt_plot(effect_sizes_zscore_all(inds), log10(pvals_KS_all(inds)), {'Color','r'});
% 
% 
% % targets
% inds=TargStatus_all==1;
% lt_plot(effect_sizes_zscore_all(inds), log10(pvals_KS_all(inds)), {'Color','k'});
% 
% 
% 
% % === REPLOT, but using baseline days
% hsplot(2)=lt_subplot(3,2,2); hold on;
% title('All baseline days, all syls, all expts');
% xlabel('effectsize (z)');
% 
% % similar
% inds=baseline_SimilarStatus_all==1;
% lt_plot(baseline_effectsizes(inds), log10(baseline_pvals(inds)), {'Color','b'});
% 
% % different
% inds=baseline_SimilarStatus_all==0;
% lt_plot(baseline_effectsizes(inds), log10(baseline_pvals(inds)), {'Color','r'});
% 
% 
% % targets
% inds=baseline_TargStatus_all==1;
% lt_plot(baseline_effectsizes(inds), log10(baseline_pvals(inds)), {'Color','k'});
% 
% 
% % ===== PLOT 3 - histogram of effect sizes
% lt_subplot(3,2,3); hold on;
% title('PDF of effect sizes (WN days, simialr(bl), diff(red); baseline days, all(bk) [all syls/expts]');
% 
% % ---- WN DAY - nontargets (Similar)
% inds=SimilarStatus_all==1 & TargStatus_all==0;
% [Ybins, Xcenters]=lt_plot_histogram(effect_sizes_zscore_all(inds), '', 0, 1);
% lt_plot_area(Xcenters, Ybins, 'b')
% 
% % ---- WN DAY - non targets (different)
% inds=SimilarStatus_all==0 & TargStatus_all==0;
% [Ybins, Xcenters]=lt_plot_histogram(effect_sizes_zscore_all(inds), Xcenters, 0, 1);
% lt_plot_area(Xcenters, Ybins, 'r')
% 
% 
% % ---- ALL BASELINE DAYS (all syls)
% [Ybins, Xcenters]=lt_plot_histogram(baseline_effectsizes, Xcenters, 0, 1);
% lt_plot_area(Xcenters, Ybins, 'k')
% 
% % --- WN DAY (target)
% inds= TargStatus_all==1;
% [Ybins, Xcenters]=lt_plot_histogram(effect_sizes_zscore_all(inds), '', 0, 1);
% lt_plot_area(Xcenters, Ybins, 'y')
% 
% 
% ylabel('prob density');
% xlabel('effect size');
% 
% 
% % ===== PLOT 4 - histogram of effect sizes (absolute values)
% lt_subplot(3,2,4); hold on;
% title('PDF of abs. effect sizes (WN days, simialr(bl), diff(red); baseline days, all(bk) [all syls/expts]');
% 
% % ---- WN DAY - nontargets (Similar)
% inds=SimilarStatus_all==1 & TargStatus_all==0;
% [Ybins, Xcenters]=lt_plot_histogram(abs(effect_sizes_zscore_all(inds)), '', 0, 1);
% lt_plot_area(Xcenters, Ybins, 'b')
% 
% % ---- WN DAY - non targets (different)
% inds=SimilarStatus_all==0 & TargStatus_all==0;
% [Ybins, Xcenters]=lt_plot_histogram(abs(effect_sizes_zscore_all(inds)), Xcenters, 0, 1);
% lt_plot_area(Xcenters, Ybins, 'r')
% 
% % ---- ALL BASELINE DAYS (all syls)
% [Ybins, Xcenters]=lt_plot_histogram(abs(baseline_effectsizes), Xcenters, 0, 1);
% lt_plot_area(Xcenters, Ybins, 'k')
% 
% % --- WN DAY (target)
% inds= TargStatus_all==1;
% [Ybins, Xcenters]=lt_plot_histogram(abs(effect_sizes_zscore_all(inds)), '', 0, 1);
% lt_plot_area(Xcenters, Ybins, 'y')
% 
% line(xlim, [0.95 0.96], 'LineStyle', '--')
% line(xlim, [0.05 0.06], 'LineStyle', '--')
% 
% % ===== PLOT 5 - histogram of effect sizes (CDF)
% lt_subplot(3,2,5); hold on; grid on
% title('CDF of effect sizes (WN days, simialr(bl), diff(red); baseline days, all(bk) [all syls/expts]');
% 
% % ---- WN DAY - nontargets (Similar)
% inds=SimilarStatus_all==1 & TargStatus_all==0;
% [Ybins, Xcenters]=lt_plot_histogram(effect_sizes_zscore_all(inds), '', 0, 1);
% Ycdf=cumsum(Ybins); % get cdf
% lt_plot(Xcenters, Ycdf, {'LineStyle' ,'-', 'Marker' ,'none', 'Color','b'})
% 
% % ---- WN DAY - nontargets (Similar)
% inds=SimilarStatus_all==0 & TargStatus_all==0;
% [Ybins, Xcenters]=lt_plot_histogram(effect_sizes_zscore_all(inds), '', 0, 1);
% Ycdf=cumsum(Ybins); % get cdf
% lt_plot(Xcenters, Ycdf, {'LineStyle' ,'-', 'Marker' ,'none', 'Color','r'})
% 
% % ---- ALL BASELINE DAYS (all syls)
% [Ybins, Xcenters]=lt_plot_histogram(baseline_effectsizes, Xcenters, 0, 1);
% Ycdf=cumsum(Ybins); % get cdf
% lt_plot(Xcenters, Ycdf, {'LineStyle' ,'-', 'Marker' ,'none', 'Color','k'})
% 
% % ===== PLOT 5 - histogram of effect sizes (absolute val CDF)
% lt_subplot(3,2,6); hold on; grid on
% title('CDF of effect sizes (WN days, simialr(bl), diff(red); baseline days, all(bk) [all syls/expts]');
% 
% % ---- WN DAY - nontargets (Similar)
% inds=SimilarStatus_all==1 & TargStatus_all==0;
% [Ybins, Xcenters]=lt_plot_histogram(abs(effect_sizes_zscore_all(inds)), '', 0, 1);
% Ycdf=cumsum(Ybins); % get cdf
% lt_plot(Xcenters, Ycdf, {'LineStyle' ,'-', 'Marker' ,'none', 'Color','b'})
% 
% % ---- WN DAY - nontargets (Similar)
% inds=SimilarStatus_all==0 & TargStatus_all==0;
% [Ybins, Xcenters]=lt_plot_histogram(abs(effect_sizes_zscore_all(inds)), '', 0, 1);
% Ycdf=cumsum(Ybins); % get cdf
% lt_plot(Xcenters, Ycdf, {'LineStyle' ,'-', 'Marker' ,'none', 'Color','r'})
% 
% % ---- ALL BASELINE DAYS (all syls)
% [Ybins, Xcenters]=lt_plot_histogram(abs(baseline_effectsizes), Xcenters, 0, 1);
% Ycdf=cumsum(Ybins); % get cdf
% lt_plot(Xcenters, Ycdf, {'LineStyle' ,'-', 'Marker' ,'none', 'Color','k'})
% 
% line(xlim, [0.95 0.96], 'LineStyle', '--')
% line(xlim, [0.05 0.06], 'LineStyle', '--')
% 
% % ===
% linkaxes(hsplot, 'xy')
% 
% % === NOTE TO SELF
% lt_figure; 
% lt_plot_text(0.1,0.5,'as function of cutoff based on baseline days, how many sim/diff syls are sign different?')
% 
% 
% 
% %% SHUFFLING METHOD (method 1) - Permutation test [WN data + baseline data combined]
% temporary_test_on_baseline_day=0; % keep 0, use 1 to do sanity check*( baseline day as WN day)
% % TO DO: 
% % 1( ONE PLOT FOR EACH BIRD. ONE FOR EACH EXPERIMENT - P VALS
% % DISTRIBUTION
% % 2) Use not just that one learning day, but other days as well.
% 
% % ===== 1)  IS THERE DAY EFFECT IN BASELINE? (ANOVA). If so, then shuffling all
% % baseline days is flawed.
% 
% 
% % ====== 2) PERFORM SHUFFLE ANALYSIS (PERMUTATION TEST)
% Zscore_Actual_all=[];
% Pvals_all=[];
% Similar_all=[];
% Target_all=[];
% Birdnum_all=[];
% Exptnum_all=[];
% Pvals_Baseline_all=[];
% Pvals_Baseline_all_ZscoreToAllVals=[];
%                 Target_baseline_all=[];
%                 Similar_baseline_all=[];
% 
% 
% for i=1:NumBirds;
%     birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
%     numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%     
%     for ii=1:numexperiments;
%         exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
%         
%         count=1;
% SubplotsPerFig=16;
% subplotrows=4;
% subplotcols=4;
% fignums_alreadyused=[];
% hfigs=[];
% 
%         syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
%         BaselineDayInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.BaselineDays;
%         
%         % -- WN day
%         % NOTE: Use the last learning day, since should only use one
%         % day to be consistent with taking each baseline day separately
%         LearningDayInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DayInds;
%         LastLearningDayInd=LearningDayInds(end);
%         
%         if temporary_test_on_baseline_day==1;
%         LastLearningDayInd=1;
%         disp('WARNING - uising basekine day for WN, sanity cheeck');
%         end
%         
%         
%         for j=1:length(syls_unique);
%             syl=syls_unique{j};
%             
%             % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%             % ======== Option 1) SHUFFLE ALL BASELINE DATA into one vector
%             Baseline_AllFFvals=[]; % will collect across all baseline days
%             
%             for k=BaselineDayInds;
%                 
%                 % skip today if lacking data
%                 if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})
%                     continue
%                 end
%                 
%                 % ----- 1) GET TODAY'S DATA (ffvals)
%                 ffvals_today=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k});
%                 
%                 % --- compile across days
%                 Baseline_AllFFvals=[Baseline_AllFFvals ffvals_today];
%             end
%             
%             % +++++++++++++++++++++ SIGNIFICANCE TESTS
%             
%             
%             % ========= Option 1) Premutation test using all data (combine
%             % baseline and WNday)
%             ActualFFvals_WN=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{LastLearningDayInd});
%             ActualFFvals_baseline=Baseline_AllFFvals;
%             Z=[ActualFFvals_baseline ActualFFvals_WN]; % combine data into one vector
%             Ntest=length(ActualFFvals_WN);
%             
%             % -- Calculate actual distance (zscore)
%             Xbase=ActualFFvals_baseline;
%             Xtest=ActualFFvals_WN;
%             
%             mean_base=mean(Xbase);
%             std_base=std(Xbase);
%             mean_test=mean(Xtest);
%             
%             D_Actual=(mean_test-mean_base)/std_base;
%             
%             % ---- Get distribution of D_shuffles
%             NumPerms=10000;
%             D_shuffle_all=[];
%             for k=1:NumPerms;
%             
%                 % shuffle data vector
%                 Inds=randperm(length(Z));
%                 
%                 Zshuff=Z(Inds);
%                 
%                 Xtest=Zshuff(Inds(1:Ntest));
%                 Xbase=Zshuff(Inds(Ntest+1:end));
%                 
%                 % -- Calculate distance (zscore)
%                 mean_base=mean(Xbase);
%                 std_base=std(Xbase);
%                 mean_test=mean(Xtest);
% 
%                 Dshuff=(mean_test-mean_base)/std_base;
%                
%                 % --- collect D vals
%                 D_shuffle_all=[D_shuffle_all Dshuff];
%             end
%                 
%             % Get p-value of actual difference
%             % (two sided)
%             p=length(find(abs(D_shuffle_all)>abs(D_Actual)))/length(D_shuffle_all);
%             
%             % if p is zero, make it 1/*N
%             if p==0;
%                 p=1/NumPerms;
%             end
%                 
%             
%             
%             
%             % --- OUTPUT
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).PermutationTest.Zscore_Actual=D_Actual;
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).PermutationTest.Zscore_ShuffleAll=D_shuffle_all;
%             SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).PermutationTest.p_twosided=p;
%             
%             
%             % ==== NULL VALUES - ALSO GET p-values for same analysis using baseline days
%             % (subsample baseline days to match sample size of learning day
%             for k=BaselineDayInds;
%                 ffvals_today=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k});
%                 
%                 % get zscore rel to actual baseline (over other baseline
%                 % days)
%                 Baselineallvals=[]; % will collect across all baseline days
%                 for kk=BaselineDayInds;
%                     
%                     % skip today if lacking data
%                     if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{kk})
%                         continue
%                     end
%                     
%                     if kk==k;
%                         continue;
%                     end
%                     
%                     % ----- 1) GET TODAY'S DATA (ffvals)
%                     ffvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{kk});
%                     
%                     % --- compile across days
%                     Baselineallvals=[Baselineallvals ffvals];
%                 end
%                 
%                 Zscore=(mean(ffvals_today)-mean(Baselineallvals))/std(Baselineallvals);
%                 
%                 % pvalue for that zscore
%                 p_base=length(find(abs(D_shuffle_all)>abs(Zscore)))/length(D_shuffle_all);
%                 
%                 % if p is zero, make it 1/*N
%                 if p_base==0;
%                     p_base=1/NumPerms;
%                 end
%                 
%                 
%                 % ==== SAVE OUTPUT
%                 SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).PermutationTest.BaselineDay_vs_Baseline{k}.p_twosided=p_base;
%             
%                 % --- COLLECT TO PLOT
%                 Target_baseline_all=[Target_baseline_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
%                 Similar_baseline_all=[Similar_baseline_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
%             Pvals_Baseline_all=[Pvals_Baseline_all p_base];
%             
%             
%             % ====== DO SAME, but 1) zscoring to all vals (including WN
%             % day), and 2) resampling to match sampel size of WN day
%             BaselinePlusWNvals=[Baselineallvals cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{LastLearningDayInd})];
%             
%             % subsample
%               if numel(ffvals_today)>Ntest;
%                   inds=randperm(Ntest);
%                   
%                   ffvals_today_subsample=ffvals_today(inds);
%               end
%             
%             % -- zscore
%             Zscore2=(mean(ffvals_today_subsample)-mean(BaselinePlusWNvals))/std(BaselinePlusWNvals);
%             
%             % -- get pval
%                 p_base=length(find(abs(D_shuffle_all)>abs(Zscore2)))/length(D_shuffle_all);
%                 
%                 % if p is zero, make it 1/*N
%                 if p_base==0;
%                     p_base=1/NumPerms;
%                 end
%             
%                 % ==== SAVE OUTPUT
%                 SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Drift.(syl).PermutationTest.BaselineDay_vs_Baseline{k}.p_twosided_WNinBaseline=p_base;
% 
%             Pvals_Baseline_all_ZscoreToAllVals=[Pvals_Baseline_all_ZscoreToAllVals p_base];
%             end
%             
%             
%             
%             
%         % --- COLLECT for all birds
%             Zscore_Actual_all=[Zscore_Actual_all D_Actual];
%             Pvals_all=[Pvals_all p];
%             Similar_all=[Similar_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
%             Target_all=[Target_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
%             Birdnum_all=[Birdnum_all i];
%             Exptnum_all=[Exptnum_all ii];
% 
%             
%             % ==== PLOT 
%             % Dvals distribution 
%             [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
%             title([syl '; p=' num2str(p)]);
%             
%             lt_plot_histogram(D_shuffle_all);
%             
%             % actual D val
%             line([D_Actual D_Actual], ylim);
%             
%         end
%         
% 
% 
%         % ==== PLOT all syls this bird
%         [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('All syls for this bird');
% 
% 
% 
% 
% 
% lt_subtitle(['Permutation test; ' birdname '-' exptname]);
%     end
% end
% 
% % ==== PLOT ACROSS ALL SYLS
% count=1;
% SubplotsPerFig=6;
% subplotrows=2;
% subplotcols=3;
% fignums_alreadyused=[];
% hfigs=[];
% 
% % ==== 1) PVALS vs. effect size
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('pval vs. learning');
% ylabel('log10(p-val)');
% xlabel('actual zscore');
% 
% % - similar (nontargs)
% inds=Similar_all==1 & Target_all==0;
% lt_plot(Zscore_Actual_all(inds), log10(Pvals_all(inds)), {'Color','b'});
% 
% % - diff (nontargs)
% inds=Similar_all==0 & Target_all==0;
% lt_plot(Zscore_Actual_all(inds), log10(Pvals_all(inds)), {'Color','r'});
% 
% % - targets
% inds=Target_all==1;
% lt_plot(Zscore_Actual_all(inds), log10(Pvals_all(inds)), {'Color','y'});
% 
% 
% % ==== 2) PVALS vs. effect size (abs)
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('pval vs. learning (abs)');
% ylabel('log10(p-val)');
% xlabel('actual zscore');
% 
% % - similar (nontargs)
% inds=Similar_all==1 & Target_all==0;
% lt_plot(abs(Zscore_Actual_all(inds)), log10(Pvals_all(inds)), {'Color','b'});
% 
% % - diff (nontargs)
% inds=Similar_all==0 & Target_all==0;
% lt_plot(abs(Zscore_Actual_all(inds)), log10(Pvals_all(inds)), {'Color','r'});
% 
% % - targets
% inds=Target_all==1;
% lt_plot(abs(Zscore_Actual_all(inds)), log10(Pvals_all(inds)), {'Color','y'});
% 
% 
% % ==== 3) PDF OF P-values
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('PDF of p valules');
% 
% % similar
% inds=Similar_all==1 & Target_all==0;
% [~, Xcenters, hbar]=lt_plot_histogram(log10(Pvals_all(inds)), '', 1, 1);
% set(hbar, 'FaceColor', 'b');
% 
% % diff
% inds=Similar_all==0 & Target_all==0;
% [~, Xcenters, hbar]=lt_plot_histogram(log10(Pvals_all(inds)), Xcenters, 1, 1);
% set(hbar, 'FaceColor', 'r');
% 
% line([log10(0.05) log10(0.05)], ylim)
% 
% % ==== 4) PDF of p values (using baseline days instead of WN days)
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('null: baseline days as WN');
% xlabel('p-values (using same shuffled dat');
% 
% inds=Target_baseline_all==0;
% % not including WN day 
% [~, ~, hbar]=lt_plot_histogram(log10(Pvals_Baseline_all(inds)), Xcenters, 1, 1);
% set(hbar, 'FaceColor', 'k');
% 
% % including WN day
% [~, ~, hbar]=lt_plot_histogram(log10(Pvals_Baseline_all_ZscoreToAllVals(inds)), Xcenters, 1, 1);
% set(hbar, 'FaceColor', 'w');
% 
% lt_plot_text(-5, 0.2, 'bk: vs. only baseline days; wh: vs including WN day (and also subsample to match N)');
% line([log10(0.05) log10(0.05)], ylim)
% 
% 
% % ==== 5) CDF OF P-values
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('CDF of p valules (nontargs,bk, sim/diff,bl/rd, targets,yel');
% 
% % nontargets
% inds=Target_all==0;
% [Ybins, Xcenters]=lt_plot_histogram(log10(Pvals_all(inds)), '', 0, 1);
% Ycdf=cumsum(Ybins);
% plot(Xcenters, Ycdf, '-k');
% 
% % similar
% inds=Similar_all==1 & Target_all==0;
% [Ybins, Xcenters]=lt_plot_histogram(log10(Pvals_all(inds)), Xcenters, 0, 1);
% Ycdf=cumsum(Ybins);
% plot(Xcenters, Ycdf, '-b');
% 
% % diff
% inds=Similar_all==0 & Target_all==0;
% [Ybins, Xcenters]=lt_plot_histogram(log10(Pvals_all(inds)), Xcenters, 0, 1);
% Ycdf=cumsum(Ybins);
% plot(Xcenters, Ycdf, '-r');
% 
% % targ
% inds=Target_all==1;
% [Ybins, Xcenters]=lt_plot_histogram(log10(Pvals_all(inds)), Xcenters, 0, 1);
% Ycdf=cumsum(Ybins);
% plot(Xcenters, Ycdf, '-y');
% 
% % baseline days (vs just baseline)
% inds=Target_baseline_all==0;
% 
% [Ybins, Xcenters]=lt_plot_histogram(log10(Pvals_Baseline_all(inds)), Xcenters, 0, 1);
% Ycdf=cumsum(Ybins);
% plot(Xcenters, Ycdf, '--k');
% 
% % baseline days (vs. [baseline WN]);
% [Ybins, Xcenters]=lt_plot_histogram(log10(Pvals_Baseline_all_ZscoreToAllVals(inds)), Xcenters, 0, 1);
% Ycdf=cumsum(Ybins);
% plot(Xcenters, Ycdf, '--r');
% 
% lt_plot_text(-4, 1.1, 'dashed: baseline days; bk: vs. other baseline days; rd: vs [baseline WN]')
% 
% line([log10(0.05) log10(0.05)], ylim, 'Color','k')
% 
% % 7) ==== NOTE TO SELF
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% lt_plot_text(0,0.3,'is expected for baseline p vals to be comparable to actual WN (since H0 is all days are baseline)');
% 
% 
% %% SHUFFLING METHOD (method 2) - Resample baseline data to simulate WN - ask about significance of WN day given that baseline
% % == LOOK JUST AT BASELINE. MODEL VARIABILITY AS SUM OF day to day
% % fluctuation + 2) residuals. Each time resample pick randomly and
% % indepednently: 1) day mean and 2) residuals (from entire pool)
% % --- NOTE: Thought about it - cannot keep Wn day in the shuffle - that
% % will certainly lead to nothing being significnat (i.e. effectively will
% % have ~1/5 prob to having that magnitude diff, if there are 4 baseline
% % days)
% % --- SHOULD NOT WEIGHT PROB OF DAY BASED ON SAMPLE SIZE - sem in some
% % sense accounts for uncertainty in mean of that day
% % ---- WILL do both subtraction of actual baseline mean and resampled
% % baseline mean
% 
% for i=1:NumBirds;
%     birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
%     numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%     
%     for ii=1:numexperiments;
%         exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
%         
% %         count=1;
% % SubplotsPerFig=16;
% % subplotrows=4;
% % subplotcols=4;
% % fignums_alreadyused=[];
% % hfigs=[];
% % 
%         syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
%         BaselineDayInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.BaselineDays;
%         
%         % -- WN day
%         % NOTE: Use the last learning day, since should only use one
%         % day to be consistent with taking each baseline day separately
%         LearningDayInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DayInds;
%         LastLearningDayInd=LearningDayInds(end);
%         
%         
%         if temporary_test_on_baseline_day==1;
%         LastLearningDayInd=1;
%         disp('WARNING - uising basekine day for WN, sanity cheeck');
%         end
%         
%         
%         for j=1:length(syls_unique);
%             syl=syls_unique{j};
%             
%             % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%             % ======== COLLECT ALL BASELINE DATA
%             Baseline_AllFFvals=[]; % will collect across all baseline days
%             Baseline_AllFFResiduals=[];
%             Baseline_AllDayFFMeans=[];
%             Baseline_AllDayFFsem=[];
%             
%             for k=BaselineDayInds;
%                 
%                 % skip today if lacking data
%                 if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})
%                     continue
%                 end
%                 
%                 % ----- 1) GET TODAY'S DATA (ffvals)
%                 ffvals_today=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k});
%                 
%                 ffmean=mean(ffvals_today);
%                 ffresiduals=ffvals_today-ffmean;
%                 ffsem=lt_sem(ffvals_today);
%                 
%                 % --- compile across days
%                 Baseline_AllFFvals=[Baseline_AllFFvals ffvals_today];
%                 Baseline_AllFFResiduals=[Baseline_AllFFResiduals ffresiduals];
%                 Baseline_AllDayFFMeans=[Baseline_AllDayFFMeans ffmean];
%                 Baseline_AllDayFFsem=[Baseline_AllDayFFsem ffsem];
%             end
%             
%             % ++++++++++++++++++++++++++++++++++++++++++
%             % ======= GENERATE NULL DISTRIBUTION OF ZSCORES
%             % sample size is that of WN day
%             ffvals_WNday=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{LastLearningDayInd});
%             Nsamps=numel(ffvals_WNday);
%             
%             Nperms=10000;
%             for k=1:Nperms;
%                % 1)  Pick a day mean (inlcuding noise, which is gaussian
%                % using sem of actual data for taht day)
%                ind=randi(length(Baseline_AllDayFFMeans), 1);
%                mu=Baseline_AllDayFFMeans(ind);
%                sigma=Baseline_AllDayFFsem(ind);
%                
%                FF_daymean_withnoise=normrnd(mu, sigma);
%                
%                % 2) Add residuals, (from pool shuffled over all baseline days and pick
%                % (with replacement)
%                inds=randperm(length(Baseline_AllFFResiduals));
%                Residuals_to_add=Baseline_AllFFResiduals(inds(1:Nsamps));
%                
%                FFvals_resampled=Residuals_to_add+FF_daymean_with_Noise_Resid;
%                 
%                % ================ NOTE: PAUSED HERE!!!!!!!!!!!!!!!!!!!!!
%                 
%                 
%                 
%                 
%             end
%             
%             
%             
%             
%             
%             
%         end
%         
%         
%     end
% end
% 
% 
% 
% 
% 
% 
% 
% %% HOW LIKLEY IS MEASURED LEARNING GIVEN BASELINE ACROSS DAYS DRIFT.
% 
% 
