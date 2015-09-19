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
            
            % --- note if similar/diff
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            Similar_all=[Similar_all similar];
            Similar_AllExpts=[Similar_AllExpts similar];
            
            % --- Collect hit rate at that syllable (mean up to start
            % of consolid day (i.e. time used for quantifying learning)
            
            day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd; % 1st WN day
            if PARAMS.global.learning_metric=='zscore';
                day2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DayInds(end);
            elseif PARAMS.global.learning_metric=='ff_consolstart';
            day2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolStartInd+SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.Consol_DayBins-1; % last day in quahntified period (e.g. consol start).
            end
                
            numhits_all=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).HIT_RATE_overdays.NumHits(day1:day2));
            numtotal_all=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).HIT_RATE_overdays.NumTotal(day1:day2));
            
            hitrate=numhits_all/numtotal_all;
            HitRate_all=[HitRate_all hitrate];
            HitRate_AllExpts=[HitRate_AllExpts hitrate];
            
            % +++ ASIDE - PUT INTO STRUCTURE
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).HITRATE.HitRate_MeanOverDaysToConsolStart=hitrate;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.Inds_WNon_to_ConsolStartBin=day1:day2;
            
            % --- Collect learning
            learning_rel_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.consolid_start_rel_targ;
            LearningRelTarg_all=[LearningRelTarg_all learning_rel_targ];
            LearningRelTarg_AllExpts=[LearningRelTarg_AllExpts learning_rel_targ];
            
            
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

% --- PLOT ACROSS ALL EXPERIMENTS
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






