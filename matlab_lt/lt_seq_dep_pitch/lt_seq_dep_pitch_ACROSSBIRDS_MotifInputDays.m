function [SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_MotifInputDays(SeqDepPitch_AcrossBirds, PARAMS, DaysToPlot, plot_hit_rate, equal_y_scale, AnnotateSylPositions)
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


%% FOR EACH EXPERIMENT, PLOT LEARNING AS FUNCTION OF MOTIF DISTANCE [USING LEARNING METRIC]

count=1;
SubplotsPerFig=6;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

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
                
                                        WNday1_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;

                if use_learning_metric==1;
                    learning_mean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
                    learning_sem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.sem;
                else
                    % ---- use the days
                    if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                        % 1) collect all ffvals
                        ffvals_all=[];
                        baseline_ffmean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
                        baseline_ffstd=std(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_WithinTimeWindow);
                        
                        for k=DaysToPlot;
                            daynum=WNday1_ind+k-1;
                            ffvals_tmp=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{daynum});
                            
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
                
            end
        end
        
        
        
        
        
        % ==== PLOT (bar)
        if isempty(Learning_all)
            disp(['PROBLEM, learning empty for one expt.bird']);
            continue;
        end
        
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
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

