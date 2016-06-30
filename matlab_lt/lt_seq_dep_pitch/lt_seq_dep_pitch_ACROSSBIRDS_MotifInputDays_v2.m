function [SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_MotifInputDays_v2(SeqDepPitch_AcrossBirds, PARAMS, DaysToPlot, plot_hit_rate, equal_y_scale, AnnotateSylPositions, plotSizeForIllustrator)
% === INPUT ARRAY OF DAYS TO PLOT PITCH SHIFT FOR
% e.g.
% DaysToPlot=[5 6]; % WN day 1 = 1; % leave EMPTY to use learning metric
% plot_hit_rate=1; plots as text over bars
% equal_y_scale=1; all with same y axis (max)

%% PARAMS

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

if ~exist('plot_hit_rate' ,'var');
    plot_hit_rate=1;
end

if ~exist('equal_y_scale' ,'var');
    equal_y_scale=0;
end

if ~exist('DaysToPlot','var');
    use_learning_metric=1;
elseif isempty(DaysToPlot);
    use_learning_metric=1;
else
    use_learning_metric=0;
end


%% FOR EACH EXPERIMENT, PLOT A BUNCH OF THINGS
% 1) learning by position
% 2) acoustic dist from targ by position (with error bars)
% 3) Pos in acoustic PCA space (connect in order by lines)
% 4) Learning vs. acoustic distance


count=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

if plotSizeForIllustrator==1;
    subplotrows=2;
end

if equal_y_scale==1;
    axis_handles_all=[];
    ylim_all=[];
end

for i=1:NumBirds;
    
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:numexperiments;
        
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        motif_list=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder;
        
        SylList_all={};
        Learning_all=[];
        Learning_sem_all=[];
        Similar_all=[];
        motif_boundaries=[];
        HitRate_all=[];
        MotifDistance_all=[];
        IsTarget_all=[];
        DistFromTarg_all=[];
        AcoustDist_all=[];
        FVpcaMean_all={};
        FVpcaCov_all={};
        SylNameFull_all={};
        PosOnItsOwnMotif_all=[];
        
        
        
        for j=1:length(motif_list);
            motif=motif_list{j};
            
            motif_boundaries=[motif_boundaries length(motif)]; % i.e. when new motifs are started;
            
            for jj=1:length(motif);
                syl=motif{jj};
                
                if ~any(strcmp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique, syl));
                    continue;
                end
                
                
                if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl), 'LEARNING')
                    continue;
                end
                
                
                if use_learning_metric==1;
                    learning_mean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
                    learning_sem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.sem;
                else
                    % ---- use the days
                        WNday1_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
                    if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                        % 1) collect all ffvals
                        ffvals_all=[];
                        baseline_ffmean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
                        baseline_ffstd=std(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_WithinTimeWindow);
                        
                        for k=DaysToPlot;
                            daynum=WNday1_ind+k-1;
                            ffvals_tmp=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{daynum};
                            
                            ffvals_all=[ffvals_all ffvals_tmp];
                        end
                        
                    else
                        % 1) collect all ffvals
                        ffvals_all=[];
                        baseline_ffmean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF;
                        baseline_ffstd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).stdFF;
                        
                        for k=DaysToPlot;
                            daynum=WNday1_ind+k-1;
                            ffvals_tmp=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{daynum});
                            
                            ffvals_all=[ffvals_all ffvals_tmp];
                        end
                    end
                    
                    % -- zscore
                    ffvals_z=(ffvals_all-baseline_ffmean)./baseline_ffstd;
                    
                    % --- mean + std of zscore
                    learning_mean=mean(ffvals_z);
                    learning_sem=lt_sem(ffvals_z);
                end
                
                
                % -- each syl, position on its own motif
                if (0)
                    pos_own=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).regexp_PosInMotif_thissyl;
                    PosOnItsOwnMotif_all=[PosOnItsOwnMotif_all pos_own];
                end
                
                SylNameFull_all=[SylNameFull_all syl];
                
                
                % -- Collect FF across syls across all motifs
                Learning_all=[Learning_all learning_mean];
                Learning_sem_all=[Learning_sem_all learning_sem];
                
                % --- is it target?
                istarget=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
                IsTarget_all=[IsTarget_all istarget];
                
                % -- Note Down distance from target
                motif_distance=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).distance_from_targ;
                MotifDistance_all=[MotifDistance_all motif_distance];
                
                
                % -- Collect syls
                syl_lower=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl;
                SylList_all=[SylList_all syl_lower];
                
                % --- note if similar/diff
                similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
                Similar_all=[Similar_all similar];
                
                % --- Collect hit rate at that syllable (mean up to start
                % of consolid day (i.e. time used for quantifying learning)
                try
                    day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd; % 1st WN day
                    day2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.dayIndsUsed(end);
                    
                    numhits_all=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).HIT_RATE_overdays.NumHits(day1:day2));
                    numtotal_all=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).HIT_RATE_overdays.NumTotal(day1:day2));
                    
                    hitrate=numhits_all/numtotal_all;
                    HitRate_all=[HitRate_all hitrate];
                catch err
                end
                
                % ==== collect annnotated position from target
                DistFromTarg=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).distance_from_targ;
                DistFromTarg_all=[DistFromTarg_all DistFromTarg];
                
                % ==== acoustic dist from target
                acoustdist=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore;
                AcoustDist_all=[AcoustDist_all acoustdist];
                
                % ==== coordinate in PC2 space
                PCAmean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_mean(1:2);
                PCAcov=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_cov(1:2, 1:2);
                
                FVpcaMean_all=[FVpcaMean_all PCAmean];
                FVpcaCov_all=[FVpcaCov_all PCAcov];
                
                
            end
        end
        
        
        
        
        % +++++++++++++++++++++++ PLOTS
        % 1) ==== LEARNING REL POSITION
        if isempty(Learning_all)
            disp(['PROBLEM, learning empty for one expt.bird']);
            continue;
        end
        
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title([birdname '-' exptname]);
        ylabel('Pitch shift (z-score)');
        
        X=1:length(Learning_all);
        Y=Learning_all;
        Ysem=Learning_sem_all;
        
        % similar
        inds=Similar_all==1 & IsTarget_all==0;
        
        Xtmp=X;
        Ytmp=Learning_all;
        Ysemtemp=Ysem;
        
        Ytmp(~inds)=nan;
        Ysemtemp(~inds)=nan;
        
        if ~isempty(Xtmp);
            hbar=lt_plot_bar(Xtmp, Ytmp, {'Errors', Ysemtemp, 'Color','b', 'BarWidth',0.9});
        end
        hold on;
        
        % diff
        inds=Similar_all==0;
        
        Xtmp=X;
        Ytmp=Learning_all;
        Ysemtemp=Ysem;
        
        Ytmp(~inds)=nan;
        Ysemtemp(~inds)=nan;
        if ~isempty(Xtmp);
            hbar=lt_plot_bar(Xtmp, Ytmp, {'Errors', Ysemtemp, 'Color','r', 'BarWidth',0.9});
            
            hold on;
        end
        
        % target
        inds=IsTarget_all==1;
        
        Xtmp=X;
        Ytmp=Learning_all;
        Ysemtemp=Ysem;
        
        Ytmp(~inds)=nan;
        Ysemtemp(~inds)=nan;
        if ~isempty(Xtmp);
            hbar=lt_plot_bar(Xtmp, Ytmp, {'Errors', Ysemtemp, 'Color','y', 'BarWidth',0.9});
            set(gca, 'XTick', 1:length(Learning_all), 'XTickLabel', SylList_all);
            set(gca, 'XTickLabel', SylList_all);
        end
        
        
        
        % vertical line separating motifs
        for j=1:length(motif_boundaries)
            X=sum(motif_boundaries(1:j));
            
            line([X+0.5 X+0.5] , ylim, 'Color', 'r', 'LineWidth', 2.5);
        end
        
        % hit rate
        if plot_hit_rate==1;
            try
                Ylim=ylim;
                
                if use_learning_metric==1;
                    y_add=0.5;
                else
                    y_add=20;
                end
                for j=1:length(HitRate_all);
                    hrate_str=[num2str(100*HitRate_all(j), '%2.1f')];
                    if HitRate_all(j)>0;
                        text(j-0.2, heaviside(Y(j))*Y(j)+y_add, hrate_str, 'FontSize',13, 'FontWeight', 'bold','Color','r');
                    else
                        text(j-0.2, heaviside(Y(j))*Y(j)+y_add, hrate_str, 'FontSize',13, 'FontWeight', 'bold');
                    end
                    
                end
            catch err
            end
        end
        
        
        % --- annotate syl positions
        if AnnotateSylPositions==1;
            Ylim=ylim;
            
            if use_learning_metric==1;
                y_add=0.8;
            else
                y_add=40;
            end
            for j=1:length(DistFromTarg_all);
                plotstring=[num2str(DistFromTarg_all(j), '%2.1f')];
                text(j-0.2, heaviside(Y(j))*Y(j)+y_add, plotstring, 'FontSize',13, 'FontWeight', 'bold');
                
            end
        end
        
        
        % If equalize y acis, collect all y scales
        if equal_y_scale==1;
            axis_handles_all=[axis_handles_all gca];
            
            ylim_all=[ylim_all; ylim];
        end
        
        
        %         % ============ 2) pos in acoustic space (vector)
        %         [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        %         title([birdname '-' exptname]);
        %         xlabel('pc1');
        %         ylabel('pc2');
        %
        %         % -- targ
        %         inds=IsTarget_all==1;
        % %         basecol=[0 1 0];
        %
        %         fvmeans=FVpcaMean_all(inds);
        %         fvcov=FVpcaCov_all(inds);
        %         syllist=SylList_all(inds);
        %         posownmotif=PosOnItsOwnMotif_all(inds);
        %
        %         plotcols={[0.1 0.8 0.1]};
        %         for j=1:length(fvmeans);
        %             color=plotcols{j};
        %             mu=fvmeans{j};
        %             cov=fvcov{j};
        %             syl=syllist{j};
        %             pos=posownmotif(j);
        %
        %             lt_plot_text(mu(1), mu(2), [syl num2str(pos)], color, 10)
        %
        %             h=error_ellipse(cov, mu, 'style','-');
        %             set(h, 'Color',color);
        %         end
        %
        %
        %         % -- similar
        %         inds=Similar_all==1 & IsTarget_all==0;
        %         basecol=[0 0 1];
        %
        %                 if any(inds);
        %
        %         fvmeans=FVpcaMean_all(inds);
        %         fvcov=FVpcaCov_all(inds);
        %         syllist=SylList_all(inds);
        %         posownmotif=PosOnItsOwnMotif_all(inds);
        %
        %         plotcols=lt_make_plot_colors(length(fvmeans), 1, basecol);
        %         for j=1:length(fvmeans);
        %             color=plotcols{j};
        %             mu=fvmeans{j};
        %             cov=fvcov{j};
        %             syl=syllist{j};
        %             pos=posownmotif(j);
        %
        %             lt_plot_text(mu(1), mu(2), [syl num2str(pos)], color, 10)
        %
        %             h=error_ellipse(cov, mu, 'style','-');
        %             set(h, 'Color',color);
        %         end
        %                 end
        %
        %         % -- diff
        %         inds=Similar_all==0;
        %         basecol=[1 0 0];
        %
        %                 if any(inds);
        %
        %         fvmeans=FVpcaMean_all(inds);
        %         fvcov=FVpcaCov_all(inds);
        %         syllist=SylList_all(inds);
        %         posownmotif=PosOnItsOwnMotif_all(inds);
        %
        %         plotcols=lt_make_plot_colors(length(fvmeans), 1, basecol);
        %         for j=1:length(fvmeans);
        %             color=plotcols{j};
        %             mu=fvmeans{j};
        %             cov=fvcov{j};
        %             syl=syllist{j};
        %             pos=posownmotif(j);
        %
        %             lt_plot_text(mu(1), mu(2), [syl num2str(pos)], color, 10)
        %
        %             h=error_ellipse(cov, mu, 'style','-');
        %             set(h, 'Color',color);
        %         end
        %                 end
        %
        
        % ============ 2) pos in acoustic space (vector)
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title([birdname '-' exptname]);
        xlabel('pc1');
        ylabel('pc2');
        
        % -- targ
        inds=IsTarget_all==1;
        %         basecol=[0 1 0];
        
        fvmeans=FVpcaMean_all(inds);
        fvcov=FVpcaCov_all(inds);
        syllist=SylNameFull_all(inds);
        
        plotcols={[0.1 0.8 0.1]};
        for j=1:length(fvmeans);
            color=plotcols{j};
            mu=fvmeans{j};
            cov=fvcov{j};
            syl=syllist{j};
            
            lt_plot_text(mu(1), mu(2), [syl], color, 8)
            
            h=error_ellipse(cov, mu, 'style','-');
            set(h, 'Color',color);
        end
        
        
        % -- similar
        inds=Similar_all==1 & IsTarget_all==0;
        basecol=[0 0 1];
        
        if any(inds);
            
            fvmeans=FVpcaMean_all(inds);
            fvcov=FVpcaCov_all(inds);
            syllist=SylNameFull_all(inds);
            
            plotcols=lt_make_plot_colors(length(fvmeans), 1, basecol);
            for j=1:length(fvmeans);
                color=plotcols{j};
                mu=fvmeans{j};
                cov=fvcov{j};
                syl=syllist{j};
                
                lt_plot_text(mu(1), mu(2), [syl], color, 8)
                
                h=error_ellipse(cov, mu, 'style','-');
                set(h, 'Color',color);
            end
        end
        
        % -- diff
        inds=Similar_all==0;
        basecol=[1 0 0];
        
        if any(inds);
            
            fvmeans=FVpcaMean_all(inds);
            fvcov=FVpcaCov_all(inds);
            syllist=SylNameFull_all(inds);
            
            plotcols=lt_make_plot_colors(length(fvmeans), 1, basecol);
            for j=1:length(fvmeans);
                color=plotcols{j};
                mu=fvmeans{j};
                cov=fvcov{j};
                syl=syllist{j};
                
                lt_plot_text(mu(1), mu(2), [syl], color, 8)
                
                h=error_ellipse(cov, mu, 'style','-');
                set(h, 'Color',color);
            end
        end
        
        
        % ============ 3) Acoustic dist(to targ) rel position
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title([birdname '-' exptname]);
        ylabel('acoust dist from targ (eucl of z)');
        X=1:length(AcoustDist_all);
        
        lt_plot_bar(X, AcoustDist_all);
        
        % vertical line separating motifs
        for j=1:length(motif_boundaries)
            X=sum(motif_boundaries(1:j));
            
            line([X+0.5 X+0.5] , ylim, 'Color', 'r', 'LineWidth', 2.5);
        end
        % horizontal line for acoustic thr
        line(xlim, [PARAMS.SylClassify.SylDistCutoff PARAMS.SylClassify.SylDistCutoff]);
        
        set(gca, 'XTick', 1:length(AcoustDist_all), 'XTickLabel', SylList_all);
        
        
        % ============ 4) Learning rel acoustic dist.
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title([birdname '-' exptname]);
        ylabel('pitch shift');
        xlabel('acoust dist from targ');
        
        % -- targ
        inds=IsTarget_all==1;
        color='k';
        
        X=AcoustDist_all(inds);
        Y=Learning_all(inds);
        lt_plot(X, Y, {'Color', color});
        
        % -- similar
        inds=~IsTarget_all & Similar_all==1;
        color='b';
        
        X=AcoustDist_all(inds);
        Y=Learning_all(inds);
        lt_plot(X, Y, {'Color', color});
        
        % -- diff
        inds=Similar_all==0;
        color='r';
        
        X=AcoustDist_all(inds);
        Y=Learning_all(inds);
        lt_plot(X, Y, {'Color', color});
        
        xlim([0 10]);
        lt_plot_zeroline;
        
    end
end

% ---- equalize y axis
if equal_y_scale==1;
    
    ylim_min=min(ylim_all(:,1));
    ylim_max=max(ylim_all(:,2));
    
    for i=1:length(axis_handles_all);
        h=axis_handles_all(i);
        
        set(h, 'Ylim', [ylim_min ylim_max]);
    end
    
    linkaxes(axis_handles_all, 'y');
end


%% ====== SEPARATE INTO DIFF REPEAT CLASSES, PLOT MEAN SHIFT OVERLAYED

disp('NOTE: plots all regexp motifs for each expt - doesnt try to figure out which one has target');

use_raw_ff=0; % then won't subtract base
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

numdaystoextend=4;
DaysToPlot_extended=[DaysToPlot DaysToPlot(end)+1:DaysToPlot(end)+numdaystoextend]; % add days to get moe data

for i=1:NumBirds
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        NumMotifs=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.expressions);
        for j=1:NumMotifs
            motifname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.expressions{j};
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title({[birdname '-' exptname '-' motifname], ['targ: ' targsyl]});
            ylabel('pitch');
            
            numsubclasses=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.subexpressions{j});
            plotcols=lt_make_plot_colors(numsubclasses, 0); % plot cols
            SubClassNames={};
            Hfigs=[];
            
            for jj=1:numsubclasses
                
                subclassname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.subexpressions{j}{jj};
                numnotes=length(subclassname);
                
                % ============ COLLECT AND PLOT FOR THIS SUBCLASS
                % -- baseline
                ffmean_base=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.data_ParsedIntoSubclasses{j}.sub_class{jj}.FFmean;
                % only continue with this subclass if baseline has >8
                % renditions
                numrends=size(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.data_ParsedIntoSubclasses{j}.sub_class{jj}.FFvals, 1);
                if numrends<9;
                    disp(['SKIPPED - baseline too few rends: (' num2str(numrends) '): ' subclassname '-' birdname '-' exptname]);
                    continue
                end
                
                % -- learning
                FFvalsAll=[];
                for d=DaysToPlot_extended;
                    day=WNday1+d-1;
                    % make sure does not encroach into bidir learning
                    if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart_Ind');
                        if day>SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind;
                            continue
                        end
                    end
                    if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'SameDir_OneDayBeforeStart_Ind');
                        if day>SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.SameDir_OneDayBeforeStart_Ind;
                            continue
                        end
                    end
                    
                    if length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.day_data)>=day
                        if ~isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.day_data{day});
                        if ~isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.day_data{day}.data_ParsedIntoSubclasses{j}.sub_class{jj}.FFvals);
                            ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.day_data{day}.data_ParsedIntoSubclasses{j}.sub_class{jj}.REL_TO_BASELINE.FFvals_minus_base;
                            FFvalsAll=[FFvalsAll; ffvals];
                        end
                        end
                    end
                end
                
                if isempty(FFvalsAll)
                    % then no data for this class,
                    continue
                end
                
                ffmean=nanmean(FFvalsAll, 1);
                ffsem=lt_sem(FFvalsAll);
                if use_raw_ff==1;
                    ffmean=ffmean+ffmean_base;
                end
                
                % ============ PLOT FOR THIS SUBCLASS
                hfig=lt_plot(1:size(FFvalsAll, 2), ffmean, {'Errors', ffsem, 'LineStyle', '-', 'Color', plotcols{jj}});
                
                % =========== COLLECT SUBCLASS (WITHIN ONE CLASS)
                
                % ==== collect names for legend
                SubClassNames=[SubClassNames subclassname];
                Hfigs=[Hfigs hfig];
            end
            
            lt_plot_zeroline
            
            % plot legend
%             legend(Hfigs, SubClassNames);
            XLABEL={};
            for kk=1:length(SubClassNames{end});
                XLABEL=[XLABEL SubClassNames{end}(kk)];
            end
            set(gca, 'XTick', 1:length(SubClassNames{end}), 'XTickLabel', XLABEL);
            
        end
    end
end

%% =========== pairwise correlations (first diff from others, and from last?)

% OLD WAY - overlayed lines in each plot.
if (0)
    disp('NOTE: plots all regexp motifs for each expt - doesnt try to figure out which one has target');

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


for i=1:NumBirds
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        NumMotifs=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.expressions);
        for j=1:NumMotifs
            motifname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.expressions{j};
            numsubclasses=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.subexpressions{j});
            SubClassNames={};
            Hfigs=[];
            
            for jj=2:numsubclasses-1 % skip first and last, as usually not much data.
                
                subclassname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.subexpressions{j}{jj};
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title({[birdname '-' exptname '-' motifname], ['targ: ' targsyl]});
                ylabel('corr (motif)');
                xlabel(subclassname);
                
                % ============ COLLECT AND PLOT FOR THIS SUBCLASS
                % -- baseline
                ffvals_base=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.data_ParsedIntoSubclasses{j}.sub_class{jj}.FFvals;
                if any(isnan(ffvals_base))
                    disp('WHANTTT - nan in base? shouldnt be problem')
                end
                
                % only continue with this subclass if baseline has >8
                % renditions
                numrends=size(ffvals_base, 1);
                if numrends<20;
                    disp(['SKIPPED - baseline too few rends: (' num2str(numrends) '): ' subclassname '-' birdname '-' exptname]);
                    continue
                end
                
                % === calculate all pairwise corrs
                corrvals=corr(ffvals_base);
                
                % === plot line for corr, rel to each syl
                plotcols=lt_make_plot_colors(size(corrvals,1), 0); % plot cols
                for jjj=1:size(corrvals,1);
                    corrtoplot=corrvals(jjj, :);
                    
                    % ============ PLOT FOR THIS SUBCLASS
                    hfig=plot(1:length(corrtoplot), corrtoplot, 'o-', 'Color', plotcols{jjj});
                    
                    % =========== COLLECT SUBCLASS (WITHIN ONE CLASS)
                    
                    % ==== collect names for legend
                    SubClassNames=[SubClassNames subclassname];
                    Hfigs=[Hfigs hfig];
                end
                
                % plot legend
                lt_plot_zeroline
                %             legend(Hfigs, SubClassNames);
                XLABEL={};
                for kk=1:length(SubClassNames{end});
                    XLABEL=[XLABEL SubClassNames{end}(kk)];
                end
                set(gca, 'XTick', 1:length(SubClassNames{end}), 'XTickLabel', XLABEL);
            end
            
        end
    end
end
end


% NEW WAY - one subplot for each line, but choose the middle subclass
disp('NOTE: plots all regexp motifs for each expt - doesnt try to figure out which one has target');


for i=1:NumBirds
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        NumMotifs=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.expressions);
        for j=1:NumMotifs
            figcount=1;
            subplotrows=4;
            subplotcols=2;
            fignums_alreadyused=[];
            hfigs=[];
            
            motifname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.expressions{j};
            numsubclasses=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.subexpressions{j});
            SubClassNames={};
            Hfigs=[];
            
            jj=ceil((3/5)*numsubclasses); % choose about 2/3 in.
            subclassname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.subexpressions{j}{jj};
            
            
            % ============ COLLECT AND PLOT FOR THIS SUBCLASS
            % -- baseline
            ffvals_base=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.data_ParsedIntoSubclasses{j}.sub_class{jj}.FFvals;
            if any(isnan(ffvals_base))
                disp('WHANTTT - nan in base? shouldnt be problem')
            end
            
            % only continue with this subclass if baseline has >8
            % renditions
            numrends=size(ffvals_base, 1);
            if numrends<20;
                disp(['SKIPPED - baseline too few rends: (' num2str(numrends) '): ' subclassname '-' birdname '-' exptname]);
                continue
            end
            
            % === calculate all pairwise corrs
            corrvals=corr(ffvals_base);
            
            % === plot line for corr, rel to each syl
            plotcols=lt_make_plot_colors(size(corrvals,1), 0); % plot cols
            Hsplots=[];
            
            for jjj=1:size(corrvals,1);
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title({[birdname '-' exptname '-' motifname], ['targ: ' targsyl]});
                ylabel('corr (motif)');
                xlabel(subclassname);
                
                corrtoplot=corrvals(jjj, :);
                
                % ============ PLOT FOR THIS SUBCLASS
                hfig=plot(1:length(corrtoplot), corrtoplot, 'o-', 'Color', plotcols{jjj});
                
                % =========== COLLECT SUBCLASS (WITHIN ONE CLASS)
                
                % ==== collect names for legend
                SubClassNames=[SubClassNames subclassname];
                Hfigs=[Hfigs hfig];
                Hsplots=[Hsplots hsplot];
                % plot legend
                lt_plot_zeroline
                %             legend(Hfigs, SubClassNames);
                XLABEL={};
                for kk=1:length(SubClassNames{end});
                    XLABEL=[XLABEL SubClassNames{end}(kk)];
                end
                set(gca, 'XTick', 1:length(SubClassNames{end}), 'XTickLabel', XLABEL);
            end
            
        end
        linkaxes( Hsplots, 'xy');
        ylim([-0.2 1.2]);
        
    end
end


