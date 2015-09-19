function [SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_MOTIFS(SeqDepPitch_AcrossBirds, PARAMS, use_learning_metric)
%% Ask whether genearlziation is stronger to same motif and if closer in motif?
%% Do syls in other motifs tend to generalization together?

%% PARAMS

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

if ~exist('use_learning_metric' ,'var');
    use_learning_metric=0;
end



%% FOR EACH EXPERIMENT, PLOT LEARNING AS FUNCTION OF MOTIF DISTANCE

if use_learning_metric==0;
count=1;
SubplotsPerFig=6;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];


for i=1:NumBirds;
   
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        motif_list=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder;
        
        SylList_all={};
        FFmean_all=[];
        FFsem_all=[];
        Similar_all=[];
        motif_boundaries=[];
        HitRate_all=[];
        MotifDistance_all=[];
        IsTarget_all=[];
        
        for j=1:length(motif_list);
            motif=motif_list{j};
            
            motif_boundaries=[motif_boundaries length(motif)]; % i.e. when new motifs are started;
            
            for jj=1:length(motif);
                syl=motif{jj};
                
                ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).ConsolStart_FFvals_minbase;
                ffmean=mean(ffvals);
                ffsem=lt_sem(ffvals);
                
                % -- Collect FF across syls across all motifs
                FFmean_all=[FFmean_all ffmean];
                FFsem_all=[FFsem_all ffsem];
                
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
                day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd; % 1st WN day
                day2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolStartInd+SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.Consol_DayBins-1; % last day in quahntified period (e.g. consol start).
                
                numhits_all=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).HIT_RATE_overdays.NumHits(day1:day2));
                numtotal_all=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).HIT_RATE_overdays.NumTotal(day1:day2));
                
                hitrate=numhits_all/numtotal_all;
                HitRate_all=[HitRate_all hitrate];
                
            end
        end
                
        % ==== PLOT (bar)
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title([birdname '-' exptname]);
        
        X=1:length(FFmean_all);
        Y=FFmean_all;
        Ysem=FFsem_all;
        
        % similar
        inds=Similar_all==1;
        
        Xtmp=X(inds);
        Ytmp=Y(inds);
        Ysemtemp=Ysem(inds);
        if ~isempty(Xtmp);
            [hbar, ~]=barwitherr(Ysemtemp, Xtmp,  Ytmp);
            set(hbar, 'FaceColor', 'b');
        end
        hold on;

        % diff
        inds=Similar_all==0;
        
        Xtmp=X(inds);
        Ytmp=Y(inds);
        Ysemtemp=Ysem(inds);
        
        if ~isempty(Xtmp);
            [hbar, ~]=barwitherr(Ysemtemp, Xtmp,  Ytmp);
            set(hbar, 'FaceColor', 'r');
            set(gca, 'XTick', 1:length(FFmean_all), 'XTickLabel', SylList_all);
            set(gca, 'XTickLabel', SylList_all);
            hold on;
        end
        
        % target
        inds=IsTarget_all==1;
        
        Xtmp=X(inds);
        Ytmp=Y(inds);
        Ysemtemp=Ysem(inds);
            [hbar, ~]=barwitherr(Ysemtemp, Xtmp,  Ytmp);
            set(hbar, 'FaceColor', 'y');
        
        % Other things
        ylabel('FF (hz, minus baseline)');
        
        % vertical line separating motifs
        for j=1:length(motif_boundaries)
            X=sum(motif_boundaries(1:j));
            
            line([X+0.5 X+0.5] , ylim, 'Color', 'r', 'LineWidth', 4);
        end
        
        % hit rate
        Ylim=ylim;
        for j=1:length(HitRate_all);
            hrate_str=[num2str(100*HitRate_all(j), '%2.1f')];
            if HitRate_all(j)>0;
                text(j-0.2, heaviside(Y(j))*Y(j)+20, hrate_str, 'FontSize',13, 'FontWeight', 'bold','Color','r');
            else
                text(j-0.2, heaviside(Y(j))*Y(j)+20, hrate_str, 'FontSize',13, 'FontWeight', 'bold');
            end
                
        end
        
    end
end

end
%% FOR EACH EXPERIMENT, PLOT LEARNING AS FUNCTION OF MOTIF DISTANCE [USING LEARNING METRIC]

if use_learning_metric==1;

count=1;
SubplotsPerFig=6;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];


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
        
        for j=1:length(motif_list);
            motif=motif_list{j};
            
            motif_boundaries=[motif_boundaries length(motif)]; % i.e. when new motifs are started;
            
            for jj=1:length(motif);
                syl=motif{jj};
                
                learning_mean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
                learning_sem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.sem;
                
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
                day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd; % 1st WN day
                day2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolStartInd+SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.Consol_DayBins-1; % last day in quahntified period (e.g. consol start).
                
                numhits_all=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).HIT_RATE_overdays.NumHits(day1:day2));
                numtotal_all=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).HIT_RATE_overdays.NumTotal(day1:day2));
                
                hitrate=numhits_all/numtotal_all;
                HitRate_all=[HitRate_all hitrate];
                
            end
        end
                
        % ==== PLOT (bar)
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title([birdname '-' exptname]);
        
        X=1:length(Learning_all);
        Y=Learning_all;
        Ysem=Learning_sem_all;
        
        % similar
        inds=Similar_all==1;
        
        Xtmp=X(inds);
        Ytmp=Y(inds);
        Ysemtemp=Ysem(inds);
        if ~isempty(Xtmp);
            [hbar, ~]=barwitherr(Ysemtemp, Xtmp,  Ytmp);
            set(hbar, 'FaceColor', 'b');
        end
        hold on;

        % diff
        inds=Similar_all==0;
        
        Xtmp=X(inds);
        Ytmp=Y(inds);
        Ysemtemp=Ysem(inds);
        
        if ~isempty(Xtmp);
            [hbar, ~]=barwitherr(Ysemtemp, Xtmp,  Ytmp);
            set(hbar, 'FaceColor', 'r');
            set(gca, 'XTick', 1:length(Learning_all), 'XTickLabel', SylList_all);
            set(gca, 'XTickLabel', SylList_all);
            hold on;
        end
        
        % target
        inds=IsTarget_all==1;
        
        Xtmp=X(inds);
        Ytmp=Y(inds);
        Ysemtemp=Ysem(inds);
            [hbar, ~]=barwitherr(Ysemtemp, Xtmp,  Ytmp);
            set(hbar, 'FaceColor', 'y');
        
        % Other things
        if use_learning_metric==1;
            ylabel('z-score learning');
        else
        ylabel('FF (hz, minus baseline)');
        end
        
        % vertical line separating motifs
        for j=1:length(motif_boundaries)
            X=sum(motif_boundaries(1:j));
            
            line([X+0.5 X+0.5] , ylim, 'Color', 'r', 'LineWidth', 4);
        end
        
        % hit rate
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
        
    end
end
end

