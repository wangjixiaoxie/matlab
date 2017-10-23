function [SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_HITRATE(SeqDepPitch_AcrossBirds, PARAMS);
%% Hit rate does not explain pattern of genearlization?

%% PARAMS

NumBirds=length(SeqDepPitch_AcrossBirds.birds);



%% FOR EACH EXPERIMENT, PLOT LEARNING AS FUNCTION OF MOTIF DISTANCE
count=1;
SubplotsPerFig=6;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

% collect stats over all expts
        SylList_AllExpts={};
        Similar_AllExpts=[];
        HitRate_AllExpts=[];
        LearningRelTarg_AllExpts=[];
        IsTarg_AllExpts = [];
        LearningHz_AllExpts = [];
        SampleSize_AllExpts = [];
        
for i=1:NumBirds;
    
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        SylList_all={};
        Similar_all=[];
        HitRate_all=[];
        LearningRelTarg_all=[];
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            % ==== EXTRACT DATA
            
            % -- Collect syls
            SylList_all=[SylList_all syl];
            SylList_AllExpts=[SylList_AllExpts syl];
            
            % -- is targ?
            istarg = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
            IsTarg_AllExpts = [IsTarg_AllExpts istarg];
            
            % --- note if similar/diff
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            Similar_all=[Similar_all similar];
            Similar_AllExpts=[Similar_AllExpts similar];
            
            % --- Collect hit rate at that syllable (mean up to start
            % of consolid day (i.e. time used for quantifying learning)
            
            day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd; % 1st WN day
            if (0) % actually just taking first 4 days, since that is what I used for learning
                if PARAMS.global.learning_metric=='zscore';
                    day2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DayInds(end);
                elseif PARAMS.global.learning_metric=='ff_consolstart';
                    day2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolStartInd+SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.Consol_DayBins-1; % last day in quahntified period (e.g. consol start).
                end
            else
                numemptydays = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode;
            day2 = day1+3+numemptydays;
            end
            
            numhits_all=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).HIT_RATE_overdays.NumHits(day1:day2));
            numtotal_all=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).HIT_RATE_overdays.NumTotal(day1:day2));
            
            hitrate=numhits_all/numtotal_all;
            HitRate_all=[HitRate_all hitrate];
            HitRate_AllExpts=[HitRate_AllExpts hitrate];
            SampleSize_AllExpts = [SampleSize_AllExpts numtotal_all];
            
            % +++ ASIDE - PUT INTO STRUCTURE
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).HITRATE.HitRate_MeanOverDaysToConsolStart=hitrate;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.Inds_WNon_to_ConsolStartBin=day1:day2;
            
            % --- Collect learning
            if (0)
            learning_rel_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.consolid_start_rel_targ;
            else
                learning_rel_targ = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ;
                learning_hz = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
                learning_hz = learning_hz * SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
            end
            LearningRelTarg_all=[LearningRelTarg_all learning_rel_targ];
            LearningRelTarg_AllExpts=[LearningRelTarg_AllExpts learning_rel_targ];
            LearningHz_AllExpts = [LearningHz_AllExpts learning_hz];
            
        end
        % === PLOT (scatter for each expt)
        
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title([birdname '-' exptname]);
        xlabel('hit rate (%)');
        ylabel('learning (rel targ)');
        
        % learning vs. hit rate
        Y=LearningRelTarg_all;
        X=HitRate_all;
        
        % similar
        inds=Similar_all==1;
        
        Xtmp=X(inds);
        Ytmp=Y(inds);
        
        lt_plot(100*Xtmp, Ytmp, {'Color' ,'b'});
        
        % different
        inds=Similar_all==0;
        
        Xtmp=X(inds);
        Ytmp=Y(inds);
        
        lt_plot(100*Xtmp, Ytmp, {'Color' ,'r'});
        
        % lines
        lt_plot_zeroline
        xlim([-5 100]);
    end
end

%% --- PLOT ACROSS ALL EXPERIMENTS
lt_figure; hold on;

lt_subplot(1,3,1); hold on;
        title('All expts (learning vs. hit rate)');
        xlabel('hit rate (%)');
        ylabel('learning (rel targ)');
        
        % learning vs. hit rate
        Y=LearningRelTarg_AllExpts;
        X=HitRate_AllExpts;
        
        % similar
        inds=Similar_AllExpts==1;
        
        Xtmp=X(inds);
        Ytmp=Y(inds);
        
        plot(100*Xtmp, Ytmp, 'ob');
        
        % different
        inds=Similar_AllExpts==0;
        
        Xtmp=X(inds);
        Ytmp=Y(inds);
        
        plot(100*Xtmp, Ytmp, 'or');
        
        % lines
        lt_plot_zeroline
        xlim([-5 100]);

        
% --- PLOT ABSOLUTE LEARNING VS. HIT RATE
lt_subplot(1,3,2); hold on;
        title('All expts (absolute learning vs. hit rate)');
        xlabel('hit rate (%)');
        ylabel('absolute learning (rel targ)');
        
        % learning vs. hit rate
        Y=LearningRelTarg_AllExpts;
        Y=abs(Y);
        X=HitRate_AllExpts;
        
        % similar
        inds=Similar_AllExpts==1;
        
        Xtmp=X(inds);
        Ytmp=Y(inds);
        
        plot(100*Xtmp, Ytmp, 'ob');
        
        % different
        inds=Similar_AllExpts==0;
        
        Xtmp=X(inds);
        Ytmp=Y(inds);
        
        plot(100*Xtmp, Ytmp, 'or');
        
        % lines
        lt_plot_zeroline
        xlim([-5 100]);

% === PLOT HISTOGRAM OF HIT RATES
lt_subplot(2,3,3); hold on;
        title('Histogram of hit rates: nontargets');
        xlabel('hit rate (%)');
        ylabel('count');
        
                X=HitRate_AllExpts;
        Y=LearningRelTarg_AllExpts;
        
        % nontargets
        inds=Y~=1;
        
        [nbins, centers]=hist(100*X(inds));
        bar(centers, nbins, 'FaceColor', 'k');
        
        % targets
lt_subplot(2,3,6); hold on;
title('targets')
          xlabel('hit rate (%)');
        ylabel('count');
      
        inds=Y==1;
        [nbins, centers]=hist(100*X(inds));
        bar(centers, nbins, 'FaceColor','y');
        
        

% ---- Divide absolute learning into bins, are hit rates different?
% half half

% Remove targets
inds=Y==1;
X(inds)=[];
Y(inds)=[];

% sort out 1st and second half
[~, inds]=sort(Y);
    
Y_sorted=Y(inds);
X_sorted=X(inds);

halfind=ceil(length(Y_sorted)/2);

lt_figure; hold on;
title('divided into two bins based on learning');
xlabel('Mean hit rate');
ylabel('Learning');

errorbar(mean(X_sorted(1:halfind)), mean(Y_sorted(1:halfind)), lt_sem(Y_sorted(1:halfind)), 'ok');
errorbar(mean(X_sorted(halfind+1:end)), mean(Y_sorted(halfind+1:end)), lt_sem(Y_sorted(halfind+1:end)), 'ok');



%% ========= NEW FIGS

lt_figure; hold on;

% ======== 1) cdf of hit rate on targ
lt_subplot(3,2,1); hold on;

% - targ
inds = IsTarg_AllExpts==1;
plotcol = 'k';

Y = HitRate_AllExpts(inds);
if (0)
    lt_plot_cdf(Y, plotcol, 0);
end

% - pdf
lt_plot_histogram(100*Y, '', 1, 1, '', '', plotcol)
lt_plot_text(50, 0.2, ['N=' num2str(length(Y))]);


% ========== 2) pdf only
lt_subplot(3,2,2); hold on;
Xcenters = 0.0025:0.005:(0.02+max(HitRate_AllExpts(IsTarg_AllExpts==0)));

% === SAME
inds = IsTarg_AllExpts==0 & Similar_AllExpts==1;
plotcol = 'b';

Y = HitRate_AllExpts(inds);
% - pdf
lt_plot_histogram(100*Y, 100*Xcenters, 1, 1, '', '', plotcol)
lt_plot_text(1, 0.5, ['N=' num2str(length(Y))], plotcol);

% ---------- difftype
inds = Similar_AllExpts==0;
plotcol = 'r';

Y = HitRate_AllExpts(inds);
% - pdf
lt_plot_histogram(100*Y, 100*Xcenters, 1, 1, '', '', plotcol)
lt_plot_text(1, 0.5, ['N=' num2str(length(Y))], plotcol);



% ===================== CORRELATION BETWEEN GENERALIZATION AND HIT RATE?
% -- targ
lt_subplot(3,2,3); hold on;
title('targ');
xlabel('Hit rate');
ylabel('Learning (hz)');
inds = IsTarg_AllExpts==1;
plotcol = 'k';

X = HitRate_AllExpts(inds);
Y = LearningHz_AllExpts(inds);

lt_regress(Y, X, 1, 0, 1, 1, plotcol, 0);


% === nontarg
lt_subplot(3,2,4); hold on;
xlabel('Hit rate');
ylabel('Learning (hz)');

% - same
inds = IsTarg_AllExpts==0 & Similar_AllExpts==1;
plotcol = 'b';

X = HitRate_AllExpts(inds);
Y = LearningHz_AllExpts(inds);

lt_regress(Y, X, 1, 0, 1, 1, plotcol, 0);

% diff
inds = Similar_AllExpts==0;
plotcol = 'r';

X = HitRate_AllExpts(inds);
Y = LearningHz_AllExpts(inds);

lt_regress(Y, X, 1, 0, 1, 1, plotcol, 0);

% === nontarg [SEPARATE INTO NO HITS AND >0 HITS]
lt_subplot(3,2,5); hold on;
title('same'); 
xlabel('no hit - hit');
ylabel('ff shift, targ dir');

% - same
inds = IsTarg_AllExpts==0 & Similar_AllExpts==1;
plotcol = 'b';

X = HitRate_AllExpts(inds)>0;
Y = LearningHz_AllExpts(inds);
plot(X,Y, 'o');

Ymean = [];
Ysem = [];
Ymean(1) = mean(Y(X==0)); % mean learning for those ina  hit rate class
Ymean(2) = mean(Y(X==1)); % mean learning for those ina  hit rate class
Ysem(1) = lt_sem(Y(X==0));
Ysem(2) = lt_sem(Y(X==1));

lt_plot_bar([0 1], Ymean, {'Errors', Ysem});
p = ranksum(Y(X==0), Y(X==1));
lt_plot_pvalue(p);

% each one sig from 0?
p1 = signrank(Y(X==0));
p2 = signrank(Y(X==1));
lt_plot_text(0, Ymean(1)*1.1, ['p=' num2str(p1)], 'r');
lt_plot_text(1, Ymean(2)*1.1, ['p=' num2str(p2)], 'r');



% ----------- diff
lt_subplot(3,2,6); hold on;
title('diff');

inds = Similar_AllExpts==0;
plotcol = 'r';

X = HitRate_AllExpts(inds)>0;
Y = LearningHz_AllExpts(inds);
plot(X,Y, 'o');

Ymean = [];
Ysem = [];
Ymean(1) = mean(Y(X==0)); % mean learning for those ina  hit rate class
Ymean(2) = mean(Y(X==1)); % mean learning for those ina  hit rate class
Ysem(1) = lt_sem(Y(X==0));
Ysem(2) = lt_sem(Y(X==1));

lt_plot_bar([0 1], Ymean, {'Errors', Ysem});
p = ranksum(Y(X==0), Y(X==1));
lt_plot_pvalue(p);
% each one sig from 0?
p1 = signrank(Y(X==0));
p2 = signrank(Y(X==1));
lt_plot_text(0, Ymean(1)*1.1, ['p=' num2str(p1)], 'r');
lt_plot_text(1, Ymean(2)*1.1, ['p=' num2str(p2)], 'r');


%% ========= FIGURES SAME AS ABOVE, BUT PLOTTED FOR ILLUSTRATOR

lt_figure; hold on;

% ======== 1) cdf of hit rate on targ
lt_subplot(3,2,1); hold on;

% - targ
inds = IsTarg_AllExpts==1;
plotcol = 'k';

Y = HitRate_AllExpts(inds);
if (0)
    lt_plot_cdf(Y, plotcol, 0);
end

% - pdf
[~, xcenters] = lt_plot_histogram(100*Y, '', 1, 1, '', '', plotcol)
lt_plot_text(50, 0.2, ['N=' num2str(length(Y))]);


% ========== SAME (BUT SAME SCALE AS TARG)
lt_subplot(3,2,5); hold on;

% === SAME
inds = IsTarg_AllExpts==0 & Similar_AllExpts==1;
plotcol = 'b';

Y = HitRate_AllExpts(inds);
% - pdf
lt_plot_histogram(100*Y, xcenters, 1, 1, '', '', plotcol)
lt_plot_text(1, 0.5, ['N=' num2str(length(Y))], plotcol);

% ---------- difftype
inds = Similar_AllExpts==0;
plotcol = 'r';

Y = HitRate_AllExpts(inds);
% - pdf
lt_plot_histogram(100*Y, xcenters, 1, 1, '', '', plotcol)
lt_plot_text(1, 0.5, ['N=' num2str(length(Y))], plotcol);



% ========== 2) pdf only
lt_subplot(3,2,2); hold on;
Xcenters = 0.0025:0.005:(0.02+max(HitRate_AllExpts(IsTarg_AllExpts==0)));

% === SAME
inds = IsTarg_AllExpts==0 & Similar_AllExpts==1;
plotcol = 'b';

Y = HitRate_AllExpts(inds);
% - pdf
lt_plot_histogram(100*Y, 100*Xcenters, 1, 1, '', '', plotcol)
lt_plot_text(1, 0.5, ['N=' num2str(length(Y))], plotcol);

% ---------- difftype
inds = Similar_AllExpts==0;
plotcol = 'r';

Y = HitRate_AllExpts(inds);
% - pdf
lt_plot_histogram(100*Y, 100*Xcenters, 1, 1, '', '', plotcol)
lt_plot_text(1, 0.5, ['N=' num2str(length(Y))], plotcol);



% ===================== CORRELATION BETWEEN GENERALIZATION AND HIT RATE?
% - same
lt_subplot(3,2,3); hold on;
xlabel('Hit rate');
ylabel('Learning (hz)');

inds = IsTarg_AllExpts==0 & Similar_AllExpts==1;
plotcol = 'b';

X = HitRate_AllExpts(inds);
Y = LearningHz_AllExpts(inds);

plot(100*X, Y, 'o', 'Color', plotcol);

lt_regress(Y, 100*X, 0, 0, 1, 1, plotcol, 1);
xlim([-1 7]);
lt_plot_zeroline;
% -- figure out range and median for those with >0
lt_plot_annotation(1, ['median ' num2str(median(X(X>0))) ', range ' ...
    num2str(min(X(X>0))) '-' num2str(max(X))], 'r');


% diff
lt_subplot(3,2,4); hold on;
xlabel('Hit rate');
ylabel('Learning (hz)');

inds = Similar_AllExpts==0;
plotcol = 'r';

X = HitRate_AllExpts(inds);
Y = LearningHz_AllExpts(inds);

plot(100*X, Y, 'o', 'Color', plotcol);
lt_regress(Y, 100*X, 0, 0, 1, 1, plotcol, 1);

% --
xlim([-1 7]);
lt_plot_zeroline;
% -- figure out range and median for those with >0
lt_plot_annotation(1, ['median ' num2str(median(X(X>0))) ', range ' ...
    num2str(min(X(X>0))) '-' num2str(max(X))], 'r');


