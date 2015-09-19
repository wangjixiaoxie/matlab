function [PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANlearning(SeqDepPitch_AcrossBirds, PARAMS, norm_by_targsyl)
%% LT 7/21/15 - Plots LMAN learning
% norm_by_targsyl = 1; for across birds plots will normalize all metrics
% (learning, AFP, MP) by target syl learning. Default: 1

close all;

%% SORT OUT ONLY THE THE EXPEIRMENTS THAT HAVE LMAN INACTIVATION DATA

% copy strcuture, save backup.
SeqDepPitch_AcrossBirds_ORIGINAL=SeqDepPitch_AcrossBirds;

filter = 'LMAN';
[SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);



%% PARAMS

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

if ~exist('norm_by_targsyl','var');
    norm_by_targsyl=1;
end


%% REMOVE SYLS THAT SHOULD NOT BE ANALYZED (I.E O/L WITH WN, for non-catch analyses)

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        % Check whether this birdname expt pair has syl that should be removed.
        % If so, remove them from unique syls
        
        inds=find(strcmp(PARAMS.global.LMAN.SylsToRemove_SingleDir, birdname));
        
        for j=1:length(inds);
            
            expt_toremove=PARAMS.global.LMAN.SylsToRemove_SingleDir{inds(j)+1};
            syls_toremove=PARAMS.global.LMAN.SylsToRemove_SingleDir{inds(j)+2};
            
            % IF CURRENT EXPERIMENT IS THE ONE TO REMOVE, THEN DO SO
            if strcmp(exptname, expt_toremove);
                
                for k=1:length(syls_toremove);
                    
                    tmp_sylremove=syls_toremove{k};
                    
                    syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
                    
                    ind_to_remove=strcmp(syls_unique, tmp_sylremove);
                    
                    syls_unique(ind_to_remove)=[];
                    
                    % Put back into structure
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique=syls_unique;
                    
                    % tell user this is done
                    disp([birdname '-' exptname ': removed ' tmp_sylremove]);
                end
            end
        end
    end
end



%% PLOT BAR PLOTS OF LEARNING FOR ALL EXPERIMENTS (BAR PLOTS FOR CONSOLIDATION START, END, BIDIR START, END)
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        % ==== ONE FIGURE PER EXPT
        lt_figure; hold on;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        % ========== COLLECT AND PLOT ONE BY ONE EACH EPOCH (E.G. CONSOLID START...)
        EpochNameList={'days_consolid_early', 'days_consolid_late', 'days_bidir_early', 'days_bidir_late'};
        for k=1:length(EpochNameList);
            epochfield=EpochNameList{k};
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, epochfield);
                
                Y_FFmean_pbs=[];
                Y_FFmean_musc=[];
                Y_FFsem_pbs=[];
                Y_FFsem_musc=[];
                Y_syls={};
                Y_similar_diff=[];
                Y_istarg=[];
                Y_AFP_bias=[];
                
                for j=1:length(SylsUnique);
                    syl=SylsUnique{j};
                    
                    % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                    % MUSC)
                    FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                    FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                    
                    FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                    FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                    
                    
                    % ===== OUTPUT DATA
                    Y_FFmean_pbs=[Y_FFmean_pbs FF_PBS];
                    Y_FFmean_musc=[Y_FFmean_musc FF_MUSC];
                    Y_FFsem_pbs=[Y_FFsem_pbs FFsem_PBS];
                    Y_FFsem_musc=[Y_FFsem_musc FFsem_MUSC];
                    Y_syls=[Y_syls, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
                    Y_similar_diff=[Y_similar_diff SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
                    Y_istarg=[Y_istarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
                    Y_AFP_bias=[Y_AFP_bias FF_PBS-FF_MUSC];
                end
                
                % ====== PLOT
                lt_subplot(2,2,k); hold on;
                title(epochfield);
                
                % ++++++++++++++++++++++++++++++++++ SIMILAR SYLS
                Inds=Y_similar_diff==1;
                X=find(Inds==1);
                
                % --- PBS
                Y=Y_FFmean_pbs(Inds);
                Ysem=Y_FFsem_pbs(Inds);
                
                [hBar hErrorbar]=barwitherr(Ysem, X-0.25, Y);
                set(hBar, 'FaceColor', 'b');
                set(hBar, 'BarWidth', 0.25);
                set(hBar, 'EdgeColor', 'k');
                set(hBar, 'LineWidth', 3);
                hold on;
                
                % --- MUSC
                Y=Y_FFmean_musc(Inds);
                Ysem=Y_FFsem_musc(Inds);
                
                [hBar hErrorbar]=barwitherr(Ysem, X, Y);
                set(hBar, 'FaceColor', 'r');
                set(hBar, 'BarWidth', 0.25);
                set(hBar, 'EdgeColor', 'k');
                set(hBar, 'LineWidth', 3);
                hold on;
                
                % ---- AFP BIAS
                Y=Y_AFP_bias(Inds);
                
                hBar=lt_plot_bar(X+0.25, Y, {'Color','w', 'BarWidth',0.25});
                set(hBar, 'EdgeColor', 'k');
                set(hBar, 'LineWidth', 3);
                
                
                % ++++++++++++++++++++++++++++++++++ DIFF SYLS
                % --- PBS
                Inds=Y_similar_diff==0;
                X=find(Inds==1);
                
                % --- PBS
                Y=Y_FFmean_pbs(Inds);
                Ysem=Y_FFsem_pbs(Inds);
                
                [hBar hErrorbar]=barwitherr(Ysem, X-0.25, Y);
                set(hBar, 'FaceColor', 'b');
                set(hBar, 'BarWidth', 0.25);
                set(hBar, 'EdgeColor', 'k');
                hold on;
                
                % --- MUSC
                Y=Y_FFmean_musc(Inds);
                Ysem=Y_FFsem_musc(Inds);
                
                [hBar hErrorbar]=barwitherr(Ysem, X, Y);
                set(hBar, 'FaceColor', 'r');
                set(hBar, 'BarWidth', 0.25);
                set(hBar, 'EdgeColor', 'k');
                hold on;
                
                % ---- AFP BIAS
                Y=Y_AFP_bias(Inds);
                
                hBar=lt_plot_bar(X+0.25, Y, {'Color','w', 'BarWidth',0.25});
                set(hBar, 'EdgeColor', 'k');
                
                
                
                % ++++++++++++++++++++++++++++++++++ TARGET SYL
                Inds=Y_istarg==1;
                X=find(Inds==1);
                
                % --- PBS
                Y=Y_FFmean_pbs(Inds);
                Ysem=Y_FFsem_pbs(Inds);
                
                [hBar hErrorbar]=barwitherr(Ysem, X-0.25, Y);
                set(hBar, 'FaceColor', 'b');
                set(hBar, 'BarWidth', 0.25);
                set(hBar, 'EdgeColor', 'y');
                set(hBar, 'LineWidth', 3);
                hold on;
                
                % --- MUSC
                Y=Y_FFmean_musc(Inds);
                Ysem=Y_FFsem_musc(Inds);
                
                [hBar hErrorbar]=barwitherr(Ysem, X, Y);
                set(hBar, 'FaceColor', 'r');
                set(hBar, 'BarWidth', 0.25);
                set(hBar, 'EdgeColor', 'y');
                set(hBar, 'LineWidth', 3);
                hold on;
                
                % ---- AFP BIAS
                Y=Y_AFP_bias(Inds);
                
                hBar=lt_plot_bar(X+0.25, Y, {'Color','w', 'BarWidth',0.25});
                set(hBar, 'EdgeColor', 'y');
                set(hBar, 'LineWidth', 3);
                
                
                % ++++++++++++++++++++++++++++++++++++++++++++++ GLOBAL
                set(gca, 'XTick', 1:length(Y_syls));
                set(gca, 'XTickLabel', Y_syls)
            end
        end
        
        % ===== GLOBAL
        lt_subtitle([birdname '-' exptname]);
    end
end


%% PLOT BAR PLOTS OF LEARNING FOR ALL EXPERIMENTS [ARRANGE IN ORDER OF ACOUSTIC SIMILARITY]
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        % ==== ONE FIGURE PER EXPT
        lt_figure; hold on;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        % ========== COLLECT AND PLOT ONE BY ONE EACH EPOCH (E.G. CONSOLID START...)
        EpochNameList={'days_consolid_early', 'days_consolid_late', 'days_bidir_early', 'days_bidir_late'};
        for k=1:length(EpochNameList);
            epochfield=EpochNameList{k};
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, epochfield);
                
                Y_FFmean_pbs=[];
                Y_FFmean_musc=[];
                Y_FFsem_pbs=[];
                Y_FFsem_musc=[];
                Y_syls={};
                Y_similar_diff=[];
                Y_istarg=[];
                Y_AFP_bias=[];
                Y_AcousticDist=[];
                
                for j=1:length(SylsUnique);
                    syl=SylsUnique{j};
                    
                    % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                    % MUSC)
                    FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                    FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                    
                    FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                    FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                    
                    
                    % ===== OUTPUT DATA
                    Y_FFmean_pbs=[Y_FFmean_pbs FF_PBS];
                    Y_FFmean_musc=[Y_FFmean_musc FF_MUSC];
                    Y_FFsem_pbs=[Y_FFsem_pbs FFsem_PBS];
                    Y_FFsem_musc=[Y_FFsem_musc FFsem_MUSC];
                    Y_syls=[Y_syls, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
                    Y_similar_diff=[Y_similar_diff SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
                    Y_istarg=[Y_istarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
                    Y_AFP_bias=[Y_AFP_bias FF_PBS-FF_MUSC];
                    Y_AcousticDist=[Y_AcousticDist SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];
                end
                
                
                % ====== FIRST REARRANGE ALL DATA IN ORDER OF ACOUSTIC SIM;
                [~, IndsSorted]=sort(Y_AcousticDist); % from most similar to least
                Y_FFmean_pbs=Y_FFmean_pbs(IndsSorted);
                Y_FFmean_musc=Y_FFmean_musc(IndsSorted);
                Y_FFsem_pbs=Y_FFsem_pbs(IndsSorted);
                Y_FFsem_musc=Y_FFsem_musc(IndsSorted);
                Y_syls=Y_syls(IndsSorted);
                Y_similar_diff=Y_similar_diff(IndsSorted);
                Y_istarg=Y_istarg(IndsSorted);
                Y_AFP_bias=Y_AFP_bias(IndsSorted);
                Y_AcousticDist=Y_AcousticDist(IndsSorted);
                               
                
                % ====== PLOT
                lt_subplot(2,2,k); hold on;
                title(epochfield);
                
                % ++++++++++++++++++++++++++++++++++ SIMILAR SYLS
                Inds=Y_similar_diff==1;
                X=find(Inds==1);
                
                % --- PBS
                Y=Y_FFmean_pbs(Inds);
                Ysem=Y_FFsem_pbs(Inds);
                
                [hBar hErrorbar]=barwitherr(Ysem, X-0.25, Y);
                set(hBar, 'FaceColor', 'b');
                set(hBar, 'BarWidth', 0.25);
                set(hBar, 'EdgeColor', 'k');
                set(hBar, 'LineWidth', 3);
                hold on;
                
                % --- MUSC
                Y=Y_FFmean_musc(Inds);
                Ysem=Y_FFsem_musc(Inds);
                
                [hBar hErrorbar]=barwitherr(Ysem, X, Y);
                set(hBar, 'FaceColor', 'r');
                set(hBar, 'BarWidth', 0.25);
                set(hBar, 'EdgeColor', 'k');
                set(hBar, 'LineWidth', 3);
                hold on;
                
                % ---- AFP BIAS
                Y=Y_AFP_bias(Inds);
                
                hBar=lt_plot_bar(X+0.25, Y, {'Color','w', 'BarWidth',0.25});
                set(hBar, 'EdgeColor', 'k');
                set(hBar, 'LineWidth', 3);
                
                
                % ++++++++++++++++++++++++++++++++++ DIFF SYLS
                % --- PBS
                Inds=Y_similar_diff==0;
                X=find(Inds==1);
                
                % --- PBS
                Y=Y_FFmean_pbs(Inds);
                Ysem=Y_FFsem_pbs(Inds);
                
                [hBar hErrorbar]=barwitherr(Ysem, X-0.25, Y);
                set(hBar, 'FaceColor', 'b');
                set(hBar, 'BarWidth', 0.25);
                set(hBar, 'EdgeColor', 'k');
                hold on;
                
                % --- MUSC
                Y=Y_FFmean_musc(Inds);
                Ysem=Y_FFsem_musc(Inds);
                
                [hBar hErrorbar]=barwitherr(Ysem, X, Y);
                set(hBar, 'FaceColor', 'r');
                set(hBar, 'BarWidth', 0.25);
                set(hBar, 'EdgeColor', 'k');
                hold on;
                
                % ---- AFP BIAS
                Y=Y_AFP_bias(Inds);
                
                hBar=lt_plot_bar(X+0.25, Y, {'Color','w', 'BarWidth',0.25});
                set(hBar, 'EdgeColor', 'k');
                
                
                
                % ++++++++++++++++++++++++++++++++++ TARGET SYL
                Inds=Y_istarg==1;
                X=find(Inds==1);
                
                % --- PBS
                Y=Y_FFmean_pbs(Inds);
                Ysem=Y_FFsem_pbs(Inds);
                
                [hBar hErrorbar]=barwitherr(Ysem, X-0.25, Y);
                set(hBar, 'FaceColor', 'b');
                set(hBar, 'BarWidth', 0.25);
                set(hBar, 'EdgeColor', 'y');
                set(hBar, 'LineWidth', 3);
                hold on;
                
                % --- MUSC
                Y=Y_FFmean_musc(Inds);
                Ysem=Y_FFsem_musc(Inds);
                
                [hBar hErrorbar]=barwitherr(Ysem, X, Y);
                set(hBar, 'FaceColor', 'r');
                set(hBar, 'BarWidth', 0.25);
                set(hBar, 'EdgeColor', 'y');
                set(hBar, 'LineWidth', 3);
                hold on;
                
                % ---- AFP BIAS
                Y=Y_AFP_bias(Inds);
                
                hBar=lt_plot_bar(X+0.25, Y, {'Color','w', 'BarWidth',0.25});
                set(hBar, 'EdgeColor', 'y');
                set(hBar, 'LineWidth', 3);
                
                
                % ++++++++++++++++++++++++++++++++++++++++++++++ GLOBAL
                set(gca, 'XTick', 1:length(Y_syls));
                set(gca, 'XTickLabel', Y_syls)
            end
        end
        
        % ===== GLOBAL
        lt_subtitle(['Sorted in decreasing acoustic similarity;' birdname '-' exptname]);
    end
end


%% PLOT CORRELATION BETWEEN MP GENEARALIZATION VS ACOUSTIC DISTANCE
Y_all=struct;
X_all=struct;
SimDiff_all=struct;
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
%     if strcmp(birdname, 'gr41gr90');
%         continue
%     end
%     
    for ii=1:numexperiments;
        
        % ==== ONE FIGURE PER EXPT
        lt_figure; hold on;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        % ========== COLLECT AND PLOT ONE BY ONE EACH EPOCH (E.G. CONSOLID START...)
        EpochNameList={'days_consolid_early', 'days_consolid_late', 'days_bidir_early', 'days_bidir_late'};
        for k=1:length(EpochNameList);
            epochfield=EpochNameList{k};
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, epochfield);
                
                Y_FFmean_pbs=[];
                Y_FFmean_musc=[];
                Y_FFsem_pbs=[];
                Y_FFsem_musc=[];
                Y_syls={};
                Y_similar_diff=[];
                Y_istarg=[];
                Y_AFP_bias=[];
                Y_AcousticDist=[];
                
                for j=1:length(SylsUnique);
                    syl=SylsUnique{j};
                    
                    % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                    % MUSC)
                    FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                    FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                    
                    FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                    FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                    
                    
                    % ===== OUTPUT DATA
                    Y_FFmean_pbs=[Y_FFmean_pbs FF_PBS];
                    Y_FFmean_musc=[Y_FFmean_musc FF_MUSC];
                    Y_FFsem_pbs=[Y_FFsem_pbs FFsem_PBS];
                    Y_FFsem_musc=[Y_FFsem_musc FFsem_MUSC];
                    Y_syls=[Y_syls, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
                    Y_similar_diff=[Y_similar_diff SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
                    Y_istarg=[Y_istarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
                    Y_AFP_bias=[Y_AFP_bias FF_PBS-FF_MUSC];
                    Y_AcousticDist=[Y_AcousticDist SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];
                end
                
                
                % ====== PLOT
                lt_subplot(2,2,k); hold on;
                
                X=Y_AcousticDist;
                Y=Y_FFmean_musc;
                
                % Convert FF to rel target
                Y_FF_targ=Y_FFmean_musc(Y_istarg==1);
                Y=Y./Y_FF_targ;
                
                % Keep only non-targets
                inds=Y_istarg==0;
                X=X(inds);
                Y=Y(inds);
                Y_similar_diff_plot=Y_similar_diff(inds);
                
                % ++++++++++++++++++++++++++++++++++++++ ALL SYLS
                lt_regress(Y, X, 1);
                title(epochfield);
                xlabel('acoustic distance');
                ylabel('MP generalization');
                
                % +++++++++++++++++++++++++++++++++++ OVERLAY SIMILAR SYLS
                Inds=Y_similar_diff_plot==1;
                Xplot=X(Inds);
                Yplot=Y(Inds);
                
                lt_plot(Xplot, Yplot, {'Color','b'});
                
                % +++++++++++++++++++++++++++++++++++ OVERLAY DIFF SYLS
                Inds=Y_similar_diff_plot==0;
                Xplot=X(Inds);
                Yplot=Y(Inds);
                
                lt_plot(Xplot, Yplot, {'Color','r'});
                
                
                % ============================ COLLECT DATA TO PLOT FOR ALL
                % EXPERIMENTS
                if ~isfield(Y_all, epochfield);
                    Y_all.(epochfield)=[]; % ready to plot
                    X_all.(epochfield)=[];
                    SimDiff_all.(epochfield)=[];
                end
                
                Y_all.(epochfield)=[Y_all.(epochfield) Y]; % ready to plot
                X_all.(epochfield)=[X_all.(epochfield) X];
                SimDiff_all.(epochfield)=[SimDiff_all.(epochfield) Y_similar_diff_plot];
                
            end
        end
        
        % ===== GLOBAL
        lt_subtitle(['MP generalization(FF(musc)/FF(musc) vs. acoustic dist;' birdname '-' exptname]);
        

    end
end

% ++++++++++++++++++++++++ PLOT ACROSS ALL BIRDS/EXPTS
lt_figure; hold on;
title('MP generalization(FF(musc)/FF(musc) vs. acoustic dist; all experiments (early consolid)');
Y=Y_all.days_consolid_early;
X=X_all.days_consolid_early;

lt_regress(Y, X, 1);
xlabel('acoustic distance');
ylabel('MP generalization');



%% PLOT CORRELATION BETWEEN AFP BIAS VS ACOUSTIC DISTANCE
Y_all=struct;
X_all=struct;
SimDiff_all=struct;

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
%         if strcmp(birdname, 'gr41gr90');
%         continue
%     end
%     
    
    for ii=1:numexperiments;
        
        % ==== ONE FIGURE PER EXPT
        lt_figure; hold on;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        % ========== COLLECT AND PLOT ONE BY ONE EACH EPOCH (E.G. CONSOLID START...)
        EpochNameList={'days_consolid_early', 'days_consolid_late', 'days_bidir_early', 'days_bidir_late'};
        for k=1:length(EpochNameList);
            epochfield=EpochNameList{k};
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, epochfield);
                
                Y_FFmean_pbs=[];
                Y_FFmean_musc=[];
                Y_FFsem_pbs=[];
                Y_FFsem_musc=[];
                Y_syls={};
                Y_similar_diff=[];
                Y_istarg=[];
                Y_AFP_bias=[];
                Y_AcousticDist=[];
                
                for j=1:length(SylsUnique);
                    syl=SylsUnique{j};
                    
                    % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                    % MUSC)
                    FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                    FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                    
                    FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                    FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                    
                    
                    % ===== OUTPUT DATA
                    Y_FFmean_pbs=[Y_FFmean_pbs FF_PBS];
                    Y_FFmean_musc=[Y_FFmean_musc FF_MUSC];
                    Y_FFsem_pbs=[Y_FFsem_pbs FFsem_PBS];
                    Y_FFsem_musc=[Y_FFsem_musc FFsem_MUSC];
                    Y_syls=[Y_syls, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
                    Y_similar_diff=[Y_similar_diff SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
                    Y_istarg=[Y_istarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
                    Y_AFP_bias=[Y_AFP_bias FF_PBS-FF_MUSC];
                    Y_AcousticDist=[Y_AcousticDist SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];
                end
                
                
                % ====== PLOT
                lt_subplot(2,2,k); hold on;
                
                X=Y_AcousticDist;
                Y=Y_AFP_bias;
                
                % Convert FF to rel target
                Y_FF_targ=Y_AFP_bias(Y_istarg==1);
                Y=Y./Y_FF_targ;
                
                % Keep only non-targets
                inds=Y_istarg==0;
                X=X(inds);
                Y=Y(inds);
                Y_similar_diff_plot=Y_similar_diff(inds);
                
                % ++++++++++++++++++++++++++++++++++++++ ALL SYLS
                lt_regress(Y, X, 1);
                title(epochfield);
                xlabel('acoustic distance');
                ylabel('AFP bias (hz)');
                
                % +++++++++++++++++++++++++++++++++++ OVERLAY SIMILAR SYLS
                Inds=Y_similar_diff_plot==1;
                Xplot=X(Inds);
                Yplot=Y(Inds);
                
                lt_plot(Xplot, Yplot, {'Color','b'});
                
                % +++++++++++++++++++++++++++++++++++ OVERLAY DIFF SYLS
                Inds=Y_similar_diff_plot==0;
                Xplot=X(Inds);
                Yplot=Y(Inds);
                
                lt_plot(Xplot, Yplot, {'Color','r'});
                
                 % ============================ COLLECT DATA TO PLOT FOR ALL
                % EXPERIMENTS
                if ~isfield(Y_all, epochfield);
                    Y_all.(epochfield)=[]; % ready to plot
                    X_all.(epochfield)=[];
                    SimDiff_all.(epochfield)=[];
                end
                
                Y_all.(epochfield)=[Y_all.(epochfield) Y]; % ready to plot
                X_all.(epochfield)=[X_all.(epochfield) X];
                SimDiff_all.(epochfield)=[SimDiff_all.(epochfield) Y_similar_diff_plot];

                
            end
        end
        
        % ===== GLOBAL
        lt_subtitle(['AFP Bias (FF(pbs)-FF(musc)) vs. acoustic dist;' birdname '-' exptname]);
    end
end


% ++++++++++++++++++++++++ PLOT ACROSS ALL BIRDS/EXPTS
lt_figure; hold on;
title('AFP Bias (FF(pbs)-FF(musc)) vs. acoustic dist; all experiments (early consolid)');
Y=Y_all.days_consolid_early;
X=X_all.days_consolid_early;

lt_regress(Y, X, 1);
xlabel('acoustic distance');
ylabel('AFP Bias');


%% PLOT CORRELATION BETWEEN MP GENEARALIZATION VS Correlation (baseline, PBS)
Y_all=struct;
X_all=struct;
SimDiff_all=struct;
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
%     if strcmp(birdname, 'gr41gr90');
%         continue
%     end
%     
    for ii=1:numexperiments;
        
        % ==== ONE FIGURE PER EXPT
        lt_figure; hold on;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        % ========== COLLECT AND PLOT ONE BY ONE EACH EPOCH (E.G. CONSOLID START...)
        EpochNameList={'days_consolid_early', 'days_consolid_late', 'days_bidir_early', 'days_bidir_late'};
        for k=1:length(EpochNameList);
            epochfield=EpochNameList{k};
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, epochfield);
                
                Y_FFmean_pbs=[];
                Y_FFmean_musc=[];
                Y_FFsem_pbs=[];
                Y_FFsem_musc=[];
                Y_syls={};
                Y_similar_diff=[];
                Y_istarg=[];
                Y_AFP_bias=[];
                Y_AcousticDist=[];
                Y_Corr=[];
                
                for j=1:length(SylsUnique);
                    syl=SylsUnique{j};
                    
                    % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                    % MUSC)
                    FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                    FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                    
                    FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                    FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                    
                    
                    % ===== OUTPUT DATA
                    Y_FFmean_pbs=[Y_FFmean_pbs FF_PBS];
                    Y_FFmean_musc=[Y_FFmean_musc FF_MUSC];
                    Y_FFsem_pbs=[Y_FFsem_pbs FFsem_PBS];
                    Y_FFsem_musc=[Y_FFsem_musc FFsem_MUSC];
                    Y_syls=[Y_syls, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
                    Y_similar_diff=[Y_similar_diff SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
                    Y_istarg=[Y_istarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
                    Y_AFP_bias=[Y_AFP_bias FF_PBS-FF_MUSC];
                    Y_AcousticDist=[Y_AcousticDist SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];
                    targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
                    Y_Corr=[Y_Corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl)];
                end
                
                
                % ====== PLOT
                lt_subplot(2,2,k); hold on;
                
                X=Y_Corr;
                Y=Y_FFmean_musc;
                
                % Convert FF to rel target
                Y_FF_targ=Y_FFmean_musc(Y_istarg==1);
                Y=Y./Y_FF_targ;
                
                % Keep only non-targets
                inds=Y_istarg==0;
                X=X(inds);
                Y=Y(inds);
                Y_similar_diff_plot=Y_similar_diff(inds);
                
                % ++++++++++++++++++++++++++++++++++++++ ALL SYLS
                lt_regress(Y, X, 1);
                title(epochfield);
                xlabel('corr (song)');
                ylabel('MP generalization');
                
                % +++++++++++++++++++++++++++++++++++ OVERLAY SIMILAR SYLS
                Inds=Y_similar_diff_plot==1;
                Xplot=X(Inds);
                Yplot=Y(Inds);
                
                lt_plot(Xplot, Yplot, {'Color','b'});
                
                % +++++++++++++++++++++++++++++++++++ OVERLAY DIFF SYLS
                Inds=Y_similar_diff_plot==0;
                Xplot=X(Inds);
                Yplot=Y(Inds);
                
                lt_plot(Xplot, Yplot, {'Color','r'});
                
                
                % ============================ COLLECT DATA TO PLOT FOR ALL
                % EXPERIMENTS
                if ~isfield(Y_all, epochfield);
                    Y_all.(epochfield)=[]; % ready to plot
                    X_all.(epochfield)=[];
                    SimDiff_all.(epochfield)=[];
                end
                
                Y_all.(epochfield)=[Y_all.(epochfield) Y]; % ready to plot
                X_all.(epochfield)=[X_all.(epochfield) X];
                SimDiff_all.(epochfield)=[SimDiff_all.(epochfield) Y_similar_diff_plot];
                
            end
        end
        
        % ===== GLOBAL
        lt_subtitle(['MP generalization(FF(musc)/FF(musc) vs. corr (MUSC, song-by-song);' birdname '-' exptname]);
        

    end
end

% ++++++++++++++++++++++++ PLOT ACROSS ALL BIRDS/EXPTS
lt_figure; hold on;
title('MP generalization(FF(musc)/FF(musc) vs. corr (MUSC, song by song); all experiments (early consolid)');
Y=Y_all.days_consolid_early;
X=X_all.days_consolid_early;

lt_regress(Y, X, 1);
xlabel('corr (song)');
ylabel('MP generalization');


%% PLOT CORRELATION BETWEEN MP GENEARALIZATION VS Correlation (baseline, MUSC)
Y_all=struct;
X_all=struct;
SimDiff_all=struct;
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
%     if strcmp(birdname, 'gr41gr90');
%         continue
%     end
%     
    for ii=1:numexperiments;
        
        % ==== ONE FIGURE PER EXPT
        lt_figure; hold on;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        % ========== COLLECT AND PLOT ONE BY ONE EACH EPOCH (E.G. CONSOLID START...)
        EpochNameList={'days_consolid_early', 'days_consolid_late', 'days_bidir_early', 'days_bidir_late'};
        for k=1:length(EpochNameList);
            epochfield=EpochNameList{k};
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, epochfield);
                
                Y_FFmean_pbs=[];
                Y_FFmean_musc=[];
                Y_FFsem_pbs=[];
                Y_FFsem_musc=[];
                Y_syls={};
                Y_similar_diff=[];
                Y_istarg=[];
                Y_AFP_bias=[];
                Y_AcousticDist=[];
                Y_Corr=[];
                
                for j=1:length(SylsUnique);
                    syl=SylsUnique{j};
                    
                    % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                    % MUSC)
                    FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                    FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                    
                    FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                    FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                    
                    
                    % ===== OUTPUT DATA
                    Y_FFmean_pbs=[Y_FFmean_pbs FF_PBS];
                    Y_FFmean_musc=[Y_FFmean_musc FF_MUSC];
                    Y_FFsem_pbs=[Y_FFsem_pbs FFsem_PBS];
                    Y_FFsem_musc=[Y_FFsem_musc FFsem_MUSC];
                    Y_syls=[Y_syls, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
                    Y_similar_diff=[Y_similar_diff SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
                    Y_istarg=[Y_istarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
                    Y_AFP_bias=[Y_AFP_bias FF_PBS-FF_MUSC];
                    Y_AcousticDist=[Y_AcousticDist SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];
                    targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
                    Y_Corr=[Y_Corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl)];
                end
                
                
                % ====== PLOT
                lt_subplot(2,2,k); hold on;
                
                X=Y_Corr;
                Y=Y_FFmean_musc;
                
                % Convert FF to rel target
                Y_FF_targ=Y_FFmean_musc(Y_istarg==1);
                Y=Y./Y_FF_targ;
                
                % Keep only non-targets
                inds=Y_istarg==0;
                X=X(inds);
                Y=Y(inds);
                Y_similar_diff_plot=Y_similar_diff(inds);
                
                % ++++++++++++++++++++++++++++++++++++++ ALL SYLS
                lt_regress(Y, X, 1);
                title(epochfield);
                xlabel('corr (song)');
                ylabel('MP generalization');
                
                % +++++++++++++++++++++++++++++++++++ OVERLAY SIMILAR SYLS
                Inds=Y_similar_diff_plot==1;
                Xplot=X(Inds);
                Yplot=Y(Inds);
                
                lt_plot(Xplot, Yplot, {'Color','b'});
                
                % +++++++++++++++++++++++++++++++++++ OVERLAY DIFF SYLS
                Inds=Y_similar_diff_plot==0;
                Xplot=X(Inds);
                Yplot=Y(Inds);
                
                lt_plot(Xplot, Yplot, {'Color','r'});
                
                
                % ============================ COLLECT DATA TO PLOT FOR ALL
                % EXPERIMENTS
                if ~isfield(Y_all, epochfield);
                    Y_all.(epochfield)=[]; % ready to plot
                    X_all.(epochfield)=[];
                    SimDiff_all.(epochfield)=[];
                end
                
                Y_all.(epochfield)=[Y_all.(epochfield) Y]; % ready to plot
                X_all.(epochfield)=[X_all.(epochfield) X];
                SimDiff_all.(epochfield)=[SimDiff_all.(epochfield) Y_similar_diff_plot];
                
            end
        end
        
        % ===== GLOBAL
        lt_subtitle(['MP generalization(FF(musc)/FF(musc) vs. corr (MUSC, song-by-song);' birdname '-' exptname]);
        

    end
end

% ++++++++++++++++++++++++ PLOT ACROSS ALL BIRDS/EXPTS
lt_figure; hold on;
title('MP generalization(FF(musc)/FF(musc) vs. corr (MUSC, song by song); all experiments (early consolid)');
Y=Y_all.days_consolid_early;
X=X_all.days_consolid_early;

lt_regress(Y, X, 1);
xlabel('corr (song)');
ylabel('MP generalization');



%% PLOT CORRELATION BETWEEN AFP BIAS VS CORR (baseline, PBS)
Y_all=struct;
X_all=struct;
SimDiff_all=struct;

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
%         if strcmp(birdname, 'gr41gr90');
%         continue
%     end
%     
    
    for ii=1:numexperiments;
        
        % ==== ONE FIGURE PER EXPT
        lt_figure; hold on;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        % ========== COLLECT AND PLOT ONE BY ONE EACH EPOCH (E.G. CONSOLID START...)
        EpochNameList={'days_consolid_early', 'days_consolid_late', 'days_bidir_early', 'days_bidir_late'};
        for k=1:length(EpochNameList);
            epochfield=EpochNameList{k};
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, epochfield);
                
                Y_FFmean_pbs=[];
                Y_FFmean_musc=[];
                Y_FFsem_pbs=[];
                Y_FFsem_musc=[];
                Y_syls={};
                Y_similar_diff=[];
                Y_istarg=[];
                Y_AFP_bias=[];
                Y_AcousticDist=[];
                Y_Corr=[];
                
                for j=1:length(SylsUnique);
                    syl=SylsUnique{j};
                    
                    % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                    % MUSC)
                    FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                    FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                    
                    FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                    FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                    
                    
                    % ===== OUTPUT DATA
                    Y_FFmean_pbs=[Y_FFmean_pbs FF_PBS];
                    Y_FFmean_musc=[Y_FFmean_musc FF_MUSC];
                    Y_FFsem_pbs=[Y_FFsem_pbs FFsem_PBS];
                    Y_FFsem_musc=[Y_FFsem_musc FFsem_MUSC];
                    Y_syls=[Y_syls, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
                    Y_similar_diff=[Y_similar_diff SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
                    Y_istarg=[Y_istarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
                    Y_AFP_bias=[Y_AFP_bias FF_PBS-FF_MUSC];
                    Y_AcousticDist=[Y_AcousticDist SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];
                    targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
                    Y_Corr=[Y_Corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl)];

                end
                
                
                % ====== PLOT
                lt_subplot(2,2,k); hold on;
                
                X=Y_Corr;
                Y=Y_AFP_bias;
                
                % Convert FF to rel target
                Y_FF_targ=Y_AFP_bias(Y_istarg==1);
                Y=Y./Y_FF_targ;
                
                % Keep only non-targets
                inds=Y_istarg==0;
                X=X(inds);
                Y=Y(inds);
                Y_similar_diff_plot=Y_similar_diff(inds);
                
                % ++++++++++++++++++++++++++++++++++++++ ALL SYLS
                lt_regress(Y, X, 1);
                title(epochfield);
                xlabel('corr(song)');
                ylabel('AFP bias (hz)');
                
                % +++++++++++++++++++++++++++++++++++ OVERLAY SIMILAR SYLS
                Inds=Y_similar_diff_plot==1;
                Xplot=X(Inds);
                Yplot=Y(Inds);
                
                lt_plot(Xplot, Yplot, {'Color','b'});
                
                % +++++++++++++++++++++++++++++++++++ OVERLAY DIFF SYLS
                Inds=Y_similar_diff_plot==0;
                Xplot=X(Inds);
                Yplot=Y(Inds);
                
                lt_plot(Xplot, Yplot, {'Color','r'});
                
                 % ============================ COLLECT DATA TO PLOT FOR ALL
                % EXPERIMENTS
                if ~isfield(Y_all, epochfield);
                    Y_all.(epochfield)=[]; % ready to plot
                    X_all.(epochfield)=[];
                    SimDiff_all.(epochfield)=[];
                end
                
                Y_all.(epochfield)=[Y_all.(epochfield) Y]; % ready to plot
                X_all.(epochfield)=[X_all.(epochfield) X];
                SimDiff_all.(epochfield)=[SimDiff_all.(epochfield) Y_similar_diff_plot];

                
            end
        end
        
        % ===== GLOBAL
        lt_subtitle(['AFP Bias (FF(pbs)-FF(musc)) vs. corr(PBS, song);' birdname '-' exptname]);
    end
end


% ++++++++++++++++++++++++ PLOT ACROSS ALL BIRDS/EXPTS
lt_figure; hold on;
title('AFP Bias (FF(pbs)-FF(musc)) vs. corr(PBS, song); all experiments (early consolid)');
Y=Y_all.days_consolid_early;
X=X_all.days_consolid_early;

lt_regress(Y, X, 1);
xlabel('corr(song)');
ylabel('AFP Bias');



%% PLOT CORRELATION BETWEEN AFP BIAS VS CORR (baseline, MUSC)
Y_all=struct;
X_all=struct;
SimDiff_all=struct;

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
%         if strcmp(birdname, 'gr41gr90');
%         continue
%     end
%     
    
    for ii=1:numexperiments;
        
        % ==== ONE FIGURE PER EXPT
        lt_figure; hold on;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        % ========== COLLECT AND PLOT ONE BY ONE EACH EPOCH (E.G. CONSOLID START...)
        EpochNameList={'days_consolid_early', 'days_consolid_late', 'days_bidir_early', 'days_bidir_late'};
        for k=1:length(EpochNameList);
            epochfield=EpochNameList{k};
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, epochfield);
                
                Y_FFmean_pbs=[];
                Y_FFmean_musc=[];
                Y_FFsem_pbs=[];
                Y_FFsem_musc=[];
                Y_syls={};
                Y_similar_diff=[];
                Y_istarg=[];
                Y_AFP_bias=[];
                Y_AcousticDist=[];
                Y_Corr=[];
                
                for j=1:length(SylsUnique);
                    syl=SylsUnique{j};
                    
                    % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                    % MUSC)
                    FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                    FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                    
                    FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                    FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                    
                    
                    % ===== OUTPUT DATA
                    Y_FFmean_pbs=[Y_FFmean_pbs FF_PBS];
                    Y_FFmean_musc=[Y_FFmean_musc FF_MUSC];
                    Y_FFsem_pbs=[Y_FFsem_pbs FFsem_PBS];
                    Y_FFsem_musc=[Y_FFsem_musc FFsem_MUSC];
                    Y_syls=[Y_syls, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
                    Y_similar_diff=[Y_similar_diff SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
                    Y_istarg=[Y_istarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
                    Y_AFP_bias=[Y_AFP_bias FF_PBS-FF_MUSC];
                    Y_AcousticDist=[Y_AcousticDist SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];
                    targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
                    Y_Corr=[Y_Corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl)];

                end
                
                
                % ====== PLOT
                lt_subplot(2,2,k); hold on;
                
                X=Y_Corr;
                Y=Y_AFP_bias;
                
                % Convert FF to rel target
                Y_FF_targ=Y_AFP_bias(Y_istarg==1);
                Y=Y./Y_FF_targ;
                
                % Keep only non-targets
                inds=Y_istarg==0;
                X=X(inds);
                Y=Y(inds);
                Y_similar_diff_plot=Y_similar_diff(inds);
                
                % ++++++++++++++++++++++++++++++++++++++ ALL SYLS
                lt_regress(Y, X, 1);
                title(epochfield);
                xlabel('corr(song)');
                ylabel('AFP bias (hz)');
                
                % +++++++++++++++++++++++++++++++++++ OVERLAY SIMILAR SYLS
                Inds=Y_similar_diff_plot==1;
                Xplot=X(Inds);
                Yplot=Y(Inds);
                
                lt_plot(Xplot, Yplot, {'Color','b'});
                
                % +++++++++++++++++++++++++++++++++++ OVERLAY DIFF SYLS
                Inds=Y_similar_diff_plot==0;
                Xplot=X(Inds);
                Yplot=Y(Inds);
                
                lt_plot(Xplot, Yplot, {'Color','r'});
                
                 % ============================ COLLECT DATA TO PLOT FOR ALL
                % EXPERIMENTS
                if ~isfield(Y_all, epochfield);
                    Y_all.(epochfield)=[]; % ready to plot
                    X_all.(epochfield)=[];
                    SimDiff_all.(epochfield)=[];
                end
                
                Y_all.(epochfield)=[Y_all.(epochfield) Y]; % ready to plot
                X_all.(epochfield)=[X_all.(epochfield) X];
                SimDiff_all.(epochfield)=[SimDiff_all.(epochfield) Y_similar_diff_plot];

                
            end
        end
        
        % ===== GLOBAL
        lt_subtitle(['AFP Bias (FF(pbs)-FF(musc)) vs. corr(MUSC, song);' birdname '-' exptname]);
    end
end


% ++++++++++++++++++++++++ PLOT ACROSS ALL BIRDS/EXPTS
lt_figure; hold on;
title('AFP Bias (FF(pbs)-FF(musc)) vs. corr(MUSC, song); all experiments (early consolid)');
Y=Y_all.days_consolid_early;
X=X_all.days_consolid_early;

lt_regress(Y, X, 1);
xlabel('corr(song)');
ylabel('AFP Bias');


pause

%% PLOT AFP BIAS VS. MP LEARNING, DISTRIBUTIONS, AND OTHER THINGS
count=1;
SubplotsPerFig=9;
subplotrows=3;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];


                    Learning_all=[];
                    MPlearning_all=[];
                    AFPbias_all=[];
                    SimDiff_all=[];
                    TargStatus_all=[];
                    PreSylSimilar_all=[];

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
%     if birdname=='gr41gr90';
%         continue;
%     end
%     
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
    for ii=1:numexperiments;
        
        % ==== ONE FIGURE PER EXPT
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        % ========== ONLY PLOT CONSOLID START 
            epochfield='days_consolid_early';
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, epochfield);
                
                
                Y_FFmean_pbs=[];
                Y_FFmean_musc=[];
                Y_FFsem_pbs=[];
                Y_FFsem_musc=[];
                Y_syls={};
                Y_similar_diff=[];
                Y_istarg=[];
                Y_AFP_bias=[];
                Y_AcousticDist=[];
                Y_Corr=[];
                Y_presimilar=[];
                for j=1:length(SylsUnique);
                    syl=SylsUnique{j};
                    
                    % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                    % MUSC)
                    FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                    FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                    
                    FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                    FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                    
                    
                    % ===== OUTPUT DATA
                    Y_FFmean_pbs=[Y_FFmean_pbs FF_PBS];
                    Y_FFmean_musc=[Y_FFmean_musc FF_MUSC];
                    Y_FFsem_pbs=[Y_FFsem_pbs FFsem_PBS];
                    Y_FFsem_musc=[Y_FFsem_musc FFsem_MUSC];
                    Y_syls=[Y_syls, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
                    Y_similar_diff=[Y_similar_diff SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
                    Y_istarg=[Y_istarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
                    Y_AFP_bias=[Y_AFP_bias FF_PBS-FF_MUSC];
                    Y_AcousticDist=[Y_AcousticDist SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];
                    targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
                    Y_Corr=[Y_Corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl)];

                    Y_presimilar=[Y_presimilar SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl];
                end
                
                
                % ================= Flip sign if learning at targsyl is negative
                if Y_FFmean_pbs(Y_istarg==1)<0;
                    Y_FFmean_pbs=-1.*Y_FFmean_pbs;
                    Y_FFmean_musc=-1.*Y_FFmean_musc;
                    Y_AFP_bias=-1.*Y_AFP_bias;
                end
                
                % ========= Normalize by targsyl if desired (PBS learning
                % by taergsyl)
                if norm_by_targsyl==1;
                    inds=Y_istarg==1;
                    learning_by_targ=Y_FFmean_pbs(inds);
                    
                    Y_FFmean_pbs=Y_FFmean_pbs./learning_by_targ;
                    Y_FFmean_musc=Y_FFmean_musc./learning_by_targ;
                    Y_AFP_bias=Y_AFP_bias./learning_by_targ;
                end
                
                % ====== PLOT
                [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
                title([birdname '-' exptname]);
                xlabel('MP learning');
                ylabel('AFP bias');
                
                % ++++++++++++++++++++++++++++++++++++++ SIMILAR SYLS
                inds=Y_similar_diff==1;
                X=Y_FFmean_musc(inds);
                Y=Y_AFP_bias(inds);
                
                lt_plot(X,Y, {'Color','b'});
                
                % +++++++++++++++++++++++++++++++++++ OVERLAY DIFF SYLS
                inds=Y_similar_diff==0;
                X=Y_FFmean_musc(inds);
                Y=Y_AFP_bias(inds);
                
                lt_plot(X,Y, {'Color','r'});
                
                % +++++++++++++++++++++++++++++++++++ OVERLAY TARGET SYLS
                inds=Y_istarg==1;
                X=Y_FFmean_musc(inds);
                Y=Y_AFP_bias(inds);
                
                lt_plot(X,Y, {'Color','k'});
                
                % +++++ global stuff
                lt_plot_line_xy;
                if norm_by_targsyl==1;
                xlim([-1 1]);
                ylim([-1 1]);
                end
                
                 % ============================ COLLECT DATA TO PLOT FOR ALL
                % EXPERIMENTS
                if any(~isnan(Y_FFmean_pbs)); % if any are not nan.
                    Learning_all=[Learning_all Y_FFmean_pbs];
                    MPlearning_all=[MPlearning_all Y_FFmean_musc];
                    AFPbias_all=[AFPbias_all Y_AFP_bias];
                    SimDiff_all=[SimDiff_all Y_similar_diff];
                    TargStatus_all=[TargStatus_all Y_istarg];
                    PreSylSimilar_all=[PreSylSimilar_all Y_presimilar];
                end
            end
        

    end
end


% ++++++++++++++++++++++++++++ PLOTS FOR ALL EXPERIMENTS
% ======= 1) AFP bias vs. MP learning [SCATTER]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments');
xlabel('MP learning');
ylabel('AFP bias');

                
% --- SIMILAR
inds=SimDiff_all==1;
X=MPlearning_all(inds);
Y=AFPbias_all(inds);

lt_plot(X,Y, {'Color','b'});

% --- DIFF
inds=SimDiff_all==0;
X=MPlearning_all(inds);
Y=AFPbias_all(inds);

lt_plot(X,Y, {'Color','r'});

% --- TARG
inds=TargStatus_all==1;
X=MPlearning_all(inds);
Y=AFPbias_all(inds);

lt_plot(X,Y, {'Color','k'});

% --- all
lt_plot_line_xy;


% ======= 2) AFP bias vs. learning [SCATTER]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments');
xlabel('Learning(PBS)');
ylabel('AFP bias');

% -- regression of all (minus targets)
inds=TargStatus_all==0;
X=Learning_all(inds);
Y=AFPbias_all(inds);

lt_regress(Y, X, 1);
                
% --- SIMILAR
inds=SimDiff_all==1;
X=Learning_all(inds);
Y=AFPbias_all(inds);

lt_plot(X,Y, {'Color','b'});

% --- DIFF
inds=SimDiff_all==0;
X=Learning_all(inds);
Y=AFPbias_all(inds);

lt_plot(X,Y, {'Color','r'});

% --- TARG
inds=TargStatus_all==1;
X=Learning_all(inds);
Y=AFPbias_all(inds);

lt_plot(X,Y, {'Color','k'});

% --- all
lt_plot_line_xy;
line([-1 1], [-1 1]);


% ======= 3) MP learning vs. learning [SCATTER]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments');
xlabel('Learning (PBS)');
ylabel('MP learning (MUSC)');

% --- regression of all (minus targs)
inds=TargStatus_all==0;
X=Learning_all(inds);
Y=MPlearning_all(inds);

lt_regress(Y, X, 1);

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0;
X=Learning_all(inds);
Y=MPlearning_all(inds);

lt_plot(X,Y, {'Color','b'});

% --- DIFF
inds=SimDiff_all==0;
X=Learning_all(inds);
Y=MPlearning_all(inds);

lt_plot(X,Y, {'Color','r'});

% --- TARG
inds=TargStatus_all==1;
X=Learning_all(inds);
Y=MPlearning_all(inds);

lt_plot(X,Y, {'Color','k'});

% --- all
lt_plot_line_xy;
line([-1 1], [-1 1]);

% ++++++++++++++++++++++++++++++++++++++++++++++++
% ================================ 4) MP-Learning-AFP [DOTS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments (MPlearning--Learning--AFPbias)');

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X,Y, 'o-b');

% --- DIFF
inds=SimDiff_all==0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X+4,Y, 'o-r');
legend({'diff'})

% --- TARG
inds=TargStatus_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X+8,Y, 'o-k');

xlim([0 13])
legend('blue=similar', 'red=diff', 'black=targ')
lt_plot_zeroline;


% ======= 4.1) MP-Learning-AFP [DOTS, MEANS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments (MEANS, MPlearning--Learning--AFPbias)');

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X, Ymean, Ysem, 's-b');
% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X, Ytmp_mean, 'o-b')


% --- DIFF
inds=SimDiff_all==0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+4, Ymean, Ysem, 's-r');

% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X+4, Ytmp_mean, 'o-r')


% --- TARG
inds=TargStatus_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+8, Ymean, Ysem, 's-k');

xlim([0 13])
lt_plot_zeroline;



% ================================ 4) MP-Learning-AFP [DOTS] [presyl
% similar]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('PRESYL SIMILAR (MPlearning--Learning--AFPbias)');

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0 & PreSylSimilar_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X,Y, 'o-b');

% --- DIFF
inds=SimDiff_all==0 & PreSylSimilar_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X+4,Y, 'o-r');
legend({'diff'})

% --- TARG
inds=TargStatus_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X+8,Y, 'o-k');

xlim([0 13])
legend('blue=similar', 'red=diff', 'black=targ')
lt_plot_zeroline;


% ======= 4.1) MP-Learning-AFP [DOTS, MEANS] [presyl sim]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('PRESYL SIM (MEANS, MPlearning--Learning--AFPbias)');

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0 & PreSylSimilar_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X, Ymean, Ysem, 's-b');
% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X, Ytmp_mean, 'o-b')


% --- DIFF
inds=SimDiff_all==0 & PreSylSimilar_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+4, Ymean, Ysem, 's-r');

% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X+4, Ytmp_mean, 'o-r')


% --- TARG
inds=TargStatus_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+8, Ymean, Ysem, 's-k');

xlim([0 13])
lt_plot_zeroline;



% ================================ 4) MP-Learning-AFP [DOTS] [presyl
% diff]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('PRESYL DIFF (MPlearning--Learning--AFPbias)');

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0 & PreSylSimilar_all==0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X,Y, 'o-b');

% --- DIFF
inds=SimDiff_all==0 & PreSylSimilar_all==0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X+4,Y, 'o-r');
legend({'diff'})

% --- TARG
inds=TargStatus_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X+8,Y, 'o-k');

xlim([0 13])
legend('blue=similar', 'red=diff', 'black=targ')
lt_plot_zeroline;


% ======= 4.1) MP-Learning-AFP [DOTS, MEANS] [presyl DIFF]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('PRESYL DIFF (MEANS, MPlearning--Learning--AFPbias)');

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0 & PreSylSimilar_all==0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X, Ymean, Ysem, 's-b');
% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X, Ytmp_mean, 'o-b')


% --- DIFF
inds=SimDiff_all==0 & PreSylSimilar_all==0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+4, Ymean, Ysem, 's-r');

% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X+4, Ytmp_mean, 'o-r')


% --- TARG
inds=TargStatus_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+8, Ymean, Ysem, 's-k');

xlim([0 13])
lt_plot_zeroline;







% ++++++++++++++++++++++++++++++++++++++++++++++++
% ============================== MP-Learning-AFP [DOTS, flip sign depending on learning sign];
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('SIGN FLIP IF LEARNING NEG (MPlearning--Learning--AFPbias)');

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0;
X=[1 2 3];

sign_to_mult=sign(Learning_all(inds)); % figure out sign of learning
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];
Y=repmat(sign_to_mult, 3,1).*Y;

plot(X,Y, 'o-b');

% --- DIFF
inds=SimDiff_all==0;
X=[1 2 3];
sign_to_mult=sign(Learning_all(inds)); % figure out sign of learning
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];
Y=repmat(sign_to_mult, 3,1).*Y;

plot(X+4,Y, 'o-r');
legend({'diff'})

% --- TARG
inds=TargStatus_all==1;
X=[1 2 3];
sign_to_mult=sign(Learning_all(inds)); % figure out sign of learning
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];
Y=repmat(sign_to_mult, 3,1).*Y;

plot(X+8,Y, 'o-k');

xlim([0 13])
legend('blue=similar', 'red=diff', 'black=targ')
lt_plot_zeroline;


% ====== MEANS
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('SIGN FLIP IF LEARNING NEG (MEANS, MPlearning--Learning--AFPbias)');

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0;
X=[1 2 3];
sign_to_mult=sign(Learning_all(inds)); % figure out sign of learning
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];
Y=repmat(sign_to_mult, 3,1).*Y;

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X, Ymean, Ysem, 's-b');
% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X, Ytmp_mean, 'o-b')


% --- DIFF
inds=SimDiff_all==0;
X=[1 2 3];
sign_to_mult=sign(Learning_all(inds)); % figure out sign of learning
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];
Y=repmat(sign_to_mult, 3,1).*Y;

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+4, Ymean, Ysem, 's-r');

% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X+4, Ytmp_mean, 'o-r')


% --- TARG
inds=TargStatus_all==1;
X=[1 2 3];
sign_to_mult=sign(Learning_all(inds)); % figure out sign of learning
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];
Y=repmat(sign_to_mult, 3,1).*Y;

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+8, Ymean, Ysem, 's-k');

xlim([0 13])
lt_plot_zeroline;




% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ [Only
% positive generalizers];
% ================================ 4) MP-Learning-AFP [DOTS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Only positive generalizers (MPlearning--Learning--AFPbias)');

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0 & Learning_all>0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X,Y, 'o-b');

% --- DIFF
inds=SimDiff_all==0 & Learning_all>0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X+4,Y, 'o-r');
legend({'diff'})

% --- TARG
inds=TargStatus_all==1 & Learning_all>0
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X+8,Y, 'o-k');

xlim([0 13])
legend('blue=similar', 'red=diff', 'black=targ')
lt_plot_zeroline;


% ======= 4.1) MP-Learning-AFP [DOTS, MEANS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('ONLY POSIT LEARNERS (MEANS, MPlearning--Learning--AFPbias)');

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0  & Learning_all>0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X, Ymean, Ysem, 's-b');
% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X, Ytmp_mean, 'o-b')


% --- DIFF
inds=SimDiff_all==0  & Learning_all>0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+4, Ymean, Ysem, 's-r');

% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X+4, Ytmp_mean, 'o-r')


% --- TARG
inds=TargStatus_all==1  & Learning_all>0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+8, Ymean, Ysem, 's-k');

xlim([0 13])
lt_plot_zeroline;


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ [Only
% LEARNING GREATER THAN X];
pos_gener_thresh=0.1;
% ================================ 4) MP-Learning-AFP [DOTS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['Only generalizers greater than ' num2str(pos_gener_thresh) ,'(MPlearning--Learning--AFPbias)']);

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0 & Learning_all>pos_gener_thresh;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X,Y, 'o-b');

% --- DIFF
inds=SimDiff_all==0 & Learning_all>pos_gener_thresh;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X+4,Y, 'o-r');
legend({'diff'})

% --- TARG
inds=TargStatus_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X+8,Y, 'o-k');

xlim([0 13])
legend('blue=similar', 'red=diff', 'black=targ')
lt_plot_zeroline;


% ======= 4.1) MP-Learning-AFP [DOTS, MEANS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['Only generalizers greater than ' num2str(pos_gener_thresh) ,'(MEANS, MPlearning--Learning--AFPbias)']);

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0  & Learning_all>pos_gener_thresh;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X, Ymean, Ysem, 's-b');
% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X, Ytmp_mean, 'o-b')


% --- DIFF
inds=SimDiff_all==0  & Learning_all>pos_gener_thresh;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+4, Ymean, Ysem, 's-r');

% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X+4, Ytmp_mean, 'o-r')


% --- TARG
inds=TargStatus_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+8, Ymean, Ysem, 's-k');

xlim([0 13])
lt_plot_zeroline;


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ [Only
% LEARNING LESS THAN X];
neg_gener_thresh=-0.1;
% ================================ 4) MP-Learning-AFP [DOTS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['Only generalizers less than ' num2str(neg_gener_thresh) ,'(MPlearning--Learning--AFPbias)']);

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0 & Learning_all<neg_gener_thresh;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X,Y, 'o-b');

% --- DIFF
inds=SimDiff_all==0 & Learning_all<neg_gener_thresh;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X+4,Y, 'o-r');
legend({'diff'})

% --- TARG
inds=TargStatus_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X+8,Y, 'o-k');

xlim([0 13])
legend('blue=similar', 'red=diff', 'black=targ')
lt_plot_zeroline;


% ======= 4.1) MP-Learning-AFP [DOTS, MEANS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['Only generalizers less than ' num2str(neg_gener_thresh) ,'(MEANS, MPlearning--Learning--AFPbias)']);

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0  & Learning_all<neg_gener_thresh;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X, Ymean, Ysem, 's-b');
% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X, Ytmp_mean, 'o-b')


% --- DIFF
inds=SimDiff_all==0  & Learning_all<neg_gener_thresh;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+4, Ymean, Ysem, 's-r');

% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X+4, Ytmp_mean, 'o-r')


% --- TARG
inds=TargStatus_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+8, Ymean, Ysem, 's-k');

xlim([0 13])
lt_plot_zeroline;

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ [Only
% NEGATIVE generalizers];
% ================================ 4) MP-Learning-AFP [DOTS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Only NEGATIVE generalizers (MPlearning--Learning--AFPbias)');

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0 & Learning_all<0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X,Y, 'o-b');

% --- DIFF
inds=SimDiff_all==0 & Learning_all<0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X+4,Y, 'o-r');
legend({'diff'})

% --- TARG
inds=TargStatus_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X+8,Y, 'o-k');

xlim([0 13])
legend('blue=similar', 'red=diff', 'black=targ')
lt_plot_zeroline;


% ======= 4.1) MP-Learning-AFP [DOTS, MEANS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('ONLY NEGATIVE LEARNERS (MEANS, MPlearning--Learning--AFPbias)');

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0  & Learning_all<0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X, Ymean, Ysem, 's-b');
% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X, Ytmp_mean, 'o-b')


% --- DIFF
inds=SimDiff_all==0  & Learning_all<0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+4, Ymean, Ysem, 's-r');

% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X+4, Ytmp_mean, 'o-r')


% --- TARG
inds=TargStatus_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+8, Ymean, Ysem, 's-k');

xlim([0 13])
lt_plot_zeroline;

% ++++++++++++++++++++++++++++++++++++++++++++++++
% ========================================= 4) MP-Learning-AFP [DOTS, ABSOLUTE VAL]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments (MPlearning--Learning--AFPbias, ABSOLUTE VALS)');

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X,abs(Y), 'o-b');

% --- DIFF
inds=SimDiff_all==0;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X+4,abs(Y), 'o-r');
legend({'diff'})

% --- TARG
inds=TargStatus_all==1;
X=[1 2 3];
Y=[MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)];

plot(X+8,abs(Y), 'o-k');

xlim([0 13])
legend('blue=similar', 'red=diff', 'black=targ')
lt_plot_zeroline;


% ======= 4.1) MP-Learning-AFP [DOTS, MEANS, ABSOLUTE VALS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments (MEANS, MPlearning--Learning--AFPbias)');

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0;
X=[1 2 3];
Y=abs([MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)]);

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X, Ymean, Ysem, 's-b');
% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X, Ytmp_mean, 'o-b')


% --- DIFF
inds=SimDiff_all==0;
X=[1 2 3];
Y=abs([MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)]);

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+4, Ymean, Ysem, 's-r');

% plot scaled up so learning is at 1
scalefactor=1/Ymean(2);
Ytmp_mean=scalefactor.*Ymean;
plot(X+4, Ytmp_mean, 'o-r')


% --- TARG
inds=TargStatus_all==1;
X=[1 2 3];
Y=abs([MPlearning_all(inds); Learning_all(inds); AFPbias_all(inds)]);

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X+8, Ymean, Ysem, 's-k');

xlim([0 13])
lt_plot_zeroline;


% ++++++++++++++++++++++++++++++++++++++++++++++ PLOT MP/LEARNING AS A
% FUNCTION OF GENERALIZATION
binsize = 10; % syls to run over;
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['binsize: ' num2str(binsize) '; MP/learning vs. generalization(PBS)']);

% --- all nontargets
inds=TargStatus_all==0;

X=Learning_all(inds); % PBS learning
Y=MPlearning_all(inds)./Learning_all(inds); % MP bias/PBS learning

lt_plot(X, Y, {'Color' ,'k'});

% running average
[X, inds2]=sort(X);
Y=Y(inds2);

Y_run=lt_running_stats(Y,binsize); 
X_run=lt_running_stats(X, binsize);

errorbar(X_run.Mean, Y_run.Mean, Y_run.SEM, 's-k');

% --- similar
inds=SimDiff_all==1 & TargStatus_all==0;

X=Learning_all(inds); % PBS learning
Y=MPlearning_all(inds)./Learning_all(inds); % MP bias/PBS learning

lt_plot(X, Y, {'Color' ,'b'});

% --- diff
inds=SimDiff_all==0 & TargStatus_all==0;

X=Learning_all(inds); % PBS learning
Y=MPlearning_all(inds)./Learning_all(inds); % MP bias/PBS learning

lt_plot(X, Y, {'Color' ,'r'});


% ---- targets
inds=TargStatus_all==1;

X=Learning_all(inds); % PBS learning
Y=MPlearning_all(inds)./Learning_all(inds); % MP bias/PBS learning

lt_plot(X, Y, {'Color' ,'g'});
errorbar(0.95, mean(Y), lt_sem(Y), 'sg');

% line for mean
line(xlim, [mean(Y) mean(Y)], 'Color','g');

% --- other lines
lt_plot_zeroline;
line(xlim, [1 1], 'Color','k','LineStyle','--');


% +++++++++++++++++++++++++++++++++++++++++++++++++++++
% ======= 5) MP-AFP [DOTS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments (MPlearning--AFPbias)');

% --- SIMILAR
inds=SimDiff_all==1 & TargStatus_all==0;
X=[1 2];
Y=[MPlearning_all(inds); AFPbias_all(inds)];

plot(X,Y, 'o-b');

% --- DIFF
inds=SimDiff_all==0;
X=[1 2];
Y=[MPlearning_all(inds); AFPbias_all(inds)];

plot(X,Y, 'o-r');

% --- TARG
inds=TargStatus_all==1;
X=[1 2];
Y=[MPlearning_all(inds); AFPbias_all(inds)];

plot(X,Y, 'o-k');

xlim([0 3])
lt_plot_zeroline;


% ======= 5.1) MP-AFP [DOTS, MEANS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments (MPlearning--AFPbias)');

% --- SIMILAR, not targ
inds=SimDiff_all==1 & TargStatus_all==0;
X=[1 2];
Y=[MPlearning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X, Ymean, Ysem, 's-b');


% --- DIFF
inds=SimDiff_all==0;
X=[1 2];
Y=[MPlearning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X, Ymean, Ysem, 's-r');

% --- TARG
inds=TargStatus_all==1;
X=[1 2];
Y=[MPlearning_all(inds); AFPbias_all(inds)];

Ymean=mean(Y, 2);
Ysem=lt_sem(Y');
errorbar(X, Ymean, Ysem, 's-k');

xlim([0 3])
lt_plot_zeroline;


% ======= 6) HISTOGRAM [NON-TARGETS, ALL]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments (NON-TARGETS)');
xlabel('hz rel. learning by target');
ylabel('count');

% --- ALL SYLS
% inds=TargStatus_all==0 & SimDiff_all==0;
inds=TargStatus_all==0;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);

if norm_by_targsyl==1;
    binsize=0.1;
Xcenters=-1:binsize:1;
xlim([-1 1]);

else
    binsize=25;
    Xcenters=-250:binsize:250;
end

Nbins_MP=hist(Dat_MP, Xcenters);
Nbins_AFP=hist(Dat_AFP, Xcenters);

% --- PLOT
lt_plot_bar(Xcenters-binsize/10, Nbins_MP, {'Color','r'});
lt_plot_bar(Xcenters, Nbins_AFP, {'Color','b'});

% plot means
MPmean=mean(Dat_MP);
AFPmean=mean(Dat_AFP);

line([MPmean MPmean], ylim, 'Color', 'r');
line([AFPmean AFPmean], ylim, 'Color', 'b');

legend({'MPlearning', 'AFPbias'})


% ======= HISTOGRAM - SIMILAR SYLS
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments (SIMILAR SYS)');
xlabel('learning by target');
ylabel('count');

% --- ALL SYLS
% --- similar'
inds=TargStatus_all==0 & SimDiff_all==1;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);

if norm_by_targsyl==1;
    binsize=0.1;
Xcenters=-1:binsize:1;
xlim([-1 1]);

else
    binsize=25;
    Xcenters=-250:binsize:250;
end

Nbins_MP=hist(Dat_MP, Xcenters);
Nbins_AFP=hist(Dat_AFP, Xcenters);

% --- PLOT
lt_plot_bar(Xcenters-binsize/10, Nbins_MP, {'Color','r'});
lt_plot_bar(Xcenters, Nbins_AFP, {'Color','b'});

% plot means
MPmean=mean(Dat_MP);
AFPmean=mean(Dat_AFP);

line([MPmean MPmean], ylim, 'Color', 'r');
line([AFPmean AFPmean], ylim, 'Color', 'b');

legend({'MPlearning', 'AFPbias'})

% ======= HISTOGRAM - DIFF SYLS
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments (DIFF SYS)');
xlabel('learning by target');
ylabel('count');

% --- ALL SYLS
% --- similar'
inds=TargStatus_all==0 & SimDiff_all==0;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);

if norm_by_targsyl==1;
    binsize=0.1;
Xcenters=-1:binsize:1;
xlim([-1 1]);

else
    binsize=25;
    Xcenters=-250:binsize:250;
end

Nbins_MP=hist(Dat_MP, Xcenters);
Nbins_AFP=hist(Dat_AFP, Xcenters);

% --- PLOT
lt_plot_bar(Xcenters-binsize/10, Nbins_MP, {'Color','r'});
lt_plot_bar(Xcenters, Nbins_AFP, {'Color','b'});

% plot means
MPmean=mean(Dat_MP);
AFPmean=mean(Dat_AFP);

line([MPmean MPmean], ylim, 'Color', 'r');
line([AFPmean AFPmean], ylim, 'Color', 'b');

legend({'MPlearning', 'AFPbias'})


% ++++++++++++++++++++++++++++++++++++++++++
% ======= 6) histogram [NON-TARGETS, POS LEARNERS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('POS LEARNERS (NON-TARGETS)');
xlabel('hz rel. learning by target');
ylabel('count');

% --- ALL SYLS
% inds=TargStatus_all==0 & SimDiff_all==0;
inds=TargStatus_all==0 & Learning_all>0;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);

if norm_by_targsyl==1;
    binsize=0.1;
Xcenters=-0.95:binsize:0.95;
xlim([-1 1]);

else
    binsize=25;
    Xcenters=-225:binsize:225;
end

Nbins_MP=hist(Dat_MP, Xcenters);
Nbins_AFP=hist(Dat_AFP, Xcenters);

% --- PLOT
lt_plot_bar(Xcenters, Nbins_MP, {'Color','r'});
lt_plot_bar(Xcenters, Nbins_AFP, {'Color','b'});

% plot means
MPmean=mean(Dat_MP);
AFPmean=mean(Dat_AFP);

line([MPmean MPmean], ylim, 'Color', 'r');
line([AFPmean AFPmean], ylim, 'Color', 'b');


legend({'MPlearning', 'AFPbias', 'Learning'})


% ++++++++++++++++++++++++++++++++++++++++++
% ======= 6) histogram [NON-TARGETS, NEG LEARNERS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('NEG LEARNERS (NON-TARGETS)');
xlabel('hz rel. learning by target');
ylabel('count');

% --- ALL SYLS
% inds=TargStatus_all==0 & SimDiff_all==0;
inds=TargStatus_all==0 & Learning_all<0;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);

if norm_by_targsyl==1;
    binsize=0.1;
Xcenters=-0.95:binsize:0.95;
xlim([-1 1]);

else
    binsize=25;
    Xcenters=-225:binsize:225;
end

Nbins_MP=hist(Dat_MP, Xcenters);
Nbins_AFP=hist(Dat_AFP, Xcenters);

% --- PLOT
lt_plot_bar(Xcenters, Nbins_MP, {'Color','r'});
lt_plot_bar(Xcenters, Nbins_AFP, {'Color','b'});

% plot means
MPmean=mean(Dat_MP);
AFPmean=mean(Dat_AFP);

line([MPmean MPmean], ylim, 'Color', 'r');
line([AFPmean AFPmean], ylim, 'Color', 'b');


legend({'MPlearning', 'AFPbias', 'Learning'})



% ++++++++++++++++++++++++++++++++++++++++++
% ======= 6) HISTOGRAM [NON-TARGETS, POS LEARNERS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('POS LEARNERS (NON-TARGETS)');
xlabel('hz rel. learning by target');
ylabel('count');

% --- ALL SYLS
% inds=TargStatus_all==0 & SimDiff_all==0;
inds=TargStatus_all==0 & Learning_all>0;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);

if norm_by_targsyl==1;
    binsize=0.1;
Xcenters=-0.95:binsize:0.95;
xlim([-1 1]);

else
    binsize=25;
    Xcenters=-225:binsize:225;
end

Nbins_MP=hist(Dat_MP, Xcenters);
Nbins_AFP=hist(Dat_AFP, Xcenters);

% --- PLOT
lt_plot_bar(Xcenters, Nbins_MP, {'Color','r'});
lt_plot_bar(Xcenters, Nbins_AFP, {'Color','b'});

% plot means
MPmean=mean(Dat_MP);
AFPmean=mean(Dat_AFP);

line([MPmean MPmean], ylim, 'Color', 'r');
line([AFPmean AFPmean], ylim, 'Color', 'b');


legend({'MPlearning', 'AFPbias', 'Learning'})



% ======= 7) HISTOGRAM [NON-TARGETS, ALL, ABSOLUTE VALUE]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments (NON-TARGETS, ABSOLUTE VALUE)');
xlabel('hz rel. learning by target');
ylabel('count');

% --- ALL SYLS
% inds=TargStatus_all==0 & SimDiff_all==0;
inds=TargStatus_all==0;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);

% -- take absolute value;
Dat_MP=abs(Dat_MP);
Dat_AFP=abs(Dat_AFP);

if norm_by_targsyl==1;
    binsize=0.1;
Xcenters=0:binsize/2:1;
xlim([-.1 1]);

else
    binsize=25;
    Xcenters=0:binsize/2:250;
end

Nbins_MP=hist(Dat_MP, Xcenters);
Nbins_AFP=hist(Dat_AFP, Xcenters);

% --- PLOT
lt_plot_bar(Xcenters-binsize/10, Nbins_MP, {'Color','r'});
lt_plot_bar(Xcenters, Nbins_AFP, {'Color','b'});

% plot means
MPmean=mean(Dat_MP);
AFPmean=mean(Dat_AFP);

line([MPmean MPmean], ylim, 'Color', 'r');
line([AFPmean AFPmean], ylim, 'Color', 'b');

legend({'MPlearning', 'AFPbias'});


% ======= 7.1) HISTOGRAM [NON-TARGETS, ALL, ABSOLUTE VALUE, INCLUDING
% LEARNING]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments (NON-TARGETS, ABSOLUTE VALUE)');
xlabel('hz rel. learning by target');
ylabel('count');

% --- ALL SYLS
% inds=TargStatus_all==0 & SimDiff_all==0;
inds=TargStatus_all==0;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);
Dat_learning=Learning_all(inds);

% -- take absolute value;
Dat_MP=abs(Dat_MP);
Dat_AFP=abs(Dat_AFP);
Dat_learning=abs(Dat_learning);

if norm_by_targsyl==1;
    binsize=0.1;
Xcenters=0:binsize/2:1;
xlim([-.1 1]);

else
    binsize=25;
    Xcenters=0:binsize/2:250;
end

Nbins_MP=hist(Dat_MP, Xcenters);
Nbins_AFP=hist(Dat_AFP, Xcenters);
Nbins_learning=hist(Dat_learning, Xcenters);

% --- PLOT
lt_plot_bar(Xcenters-binsize/10, Nbins_MP, {'Color','r'});
lt_plot_bar(Xcenters, Nbins_AFP, {'Color','b'});
lt_plot_bar(Xcenters+binsize/10, Nbins_learning, {'Color','w'});

% plot means
MPmean=mean(Dat_MP);
AFPmean=mean(Dat_AFP);
Learningmean=mean(Dat_learning);

line([MPmean MPmean], ylim, 'Color', 'r');
line([AFPmean AFPmean], ylim, 'Color', 'b');
line([Learningmean Learningmean], ylim, 'Color', 'k');

legend({'MPlearning', 'AFPbias', 'Learning'})


% ======= 8) HISTOGRAM [TARGETS, ALL]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments (Only targets)');
xlabel('hz rel. learning by target');
ylabel('count');

% --- ALL SYLS
% inds=TargStatus_all==0 & SimDiff_all==0;
inds=TargStatus_all==1;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);
% Dat_learning=Learning_all(inds);

if norm_by_targsyl==1;
    binsize=0.1;
Xcenters=0:binsize/2:1;
xlim([0 1.1]);

else
    binsize=25;
    Xcenters=0:binsize/2:250;
end

Nbins_MP=hist(Dat_MP, Xcenters);
Nbins_AFP=hist(Dat_AFP, Xcenters);

% --- PLOT
lt_plot_bar(Xcenters-binsize/10, Nbins_MP, {'Color','r'});
lt_plot_bar(Xcenters, Nbins_AFP, {'Color','b'});

% plot means
MPmean=mean(Dat_MP);
AFPmean=mean(Dat_AFP);

line([MPmean MPmean], ylim, 'Color', 'r');
line([AFPmean AFPmean], ylim, 'Color', 'b');

legend({'MPlearning', 'AFPbias'})


% =========== 9) Plot AFP BIAS, MP [sorted by MP]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments (NONTARGETS)');
xlabel('index, after sorting');
ylabel('hz, rel learning by targ');

% --- ALL SYLS
inds=TargStatus_all==0;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);
Dat_learning=Learning_all(inds);
X=1:length(Dat_MP);

% ---- SORT BY MP BIAS
[~, inds]=sort(Dat_MP);
Dat_MP_sorted=Dat_MP(inds);
Dat_AFP_sorted=Dat_AFP(inds);

lt_plot_bar(X-0.1, Dat_MP_sorted, {'Color', 'r'})
lt_plot_bar(X+0.1, Dat_AFP_sorted, {'Color', 'b'})

legend('MP learning', 'AFP learning')

% =========== 9.1) Plot AFP BIAS, MP [sorted by AFP]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments (NONTARGETS)');
xlabel('index, after sorting');
ylabel('hz, rel learning by targ');

% --- ALL SYLS
inds=TargStatus_all==0;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);
Dat_learning=Learning_all(inds);
X=1:length(Dat_MP);

% ---- SORT BY AFP BIAS
[~, inds]=sort(Dat_AFP);
Dat_MP_sorted=Dat_MP(inds);
Dat_AFP_sorted=Dat_AFP(inds);

lt_plot_bar(X-0.1, Dat_MP_sorted, {'Color', 'r'})
lt_plot_bar(X+0.1, Dat_AFP_sorted, {'Color', 'b'})

legend('MP learning', 'AFP learning')

% =========== 9.1) Plot AFP BIAS, MP [sorted by LEARNING]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments (NONTARGETS, SORTED BY LEARNING)');
xlabel('from most negative to most positive genearlization');
ylabel('hz, rel learning by targ');

% --- ALL SYLS
inds=TargStatus_all==0;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);
Dat_learning=Learning_all(inds);
X=1:length(Dat_MP);

% ---- SORT BY LEARNING
[Dat_learning_sorted, inds]=sort(Dat_learning);
Dat_MP_sorted=Dat_MP(inds);
Dat_AFP_sorted=Dat_AFP(inds);

lt_plot_bar(X-0.1, Dat_MP_sorted, {'Color', 'r'})
lt_plot_bar(X+0.1, Dat_AFP_sorted, {'Color', 'b'})
plot(Dat_learning_sorted, '-k')

legend('MP learning', 'AFP learning', 'Learning')

% =========== 9.1) Plot AFP BIAS, MP [sorted by LEARNING, SEPARATING MP, AFP, AND LEARNING]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('MP learning (NONTARGETS, SORTED BY LEARNING)');
xlabel('from most negative to most positive genearlization');
ylabel('hz, rel learning by targ');

% --- ALL SYLS
inds=TargStatus_all==0;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);
Dat_learning=Learning_all(inds);
X=1:length(Dat_MP);

% ---- SORT BY LEARNING
[Dat_learning_sorted, inds]=sort(Dat_learning);
Dat_MP_sorted=Dat_MP(inds);
Dat_AFP_sorted=Dat_AFP(inds);

lt_plot_bar(X, Dat_MP_sorted, {'Color', 'r'})
plot(Dat_learning_sorted, '-k')

legend('MP learning', 'Learning')

% =========== 9.1) Plot AFP BIAS, MP [sorted by LEARNING, SEPARATING MP, AFP, AND LEARNING]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('AFP learning (NONTARGETS, SORTED BY LEARNING)');
xlabel('from most negative to most positive genearlization');
ylabel('hz, rel learning by targ');

% --- ALL SYLS
inds=TargStatus_all==0;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);
Dat_learning=Learning_all(inds);
X=1:length(Dat_MP);

% ---- SORT BY LEARNING
[Dat_learning_sorted, inds]=sort(Dat_learning);
Dat_MP_sorted=Dat_MP(inds);
Dat_AFP_sorted=Dat_AFP(inds);

lt_plot_bar(X, Dat_AFP_sorted, {'Color', 'b'})
plot(Dat_learning_sorted, '-k')

legend('MP learning', 'Learning')

% =========== 9.2) Plot AFP BIAS, MP [sorted by LEARNING, SEPARATING MP,
% AFP, AND LEARNING] - PRESYL SIMILAR
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('AFP learning (PRESYL SIMILAR, SORTED BY LEARNING)');
xlabel('from most negative to most positive genearlization');
ylabel('hz, rel learning by targ');

% --- ALL SYLS
inds=TargStatus_all==0 & PreSylSimilar_all==1;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);
Dat_learning=Learning_all(inds);
X=1:length(Dat_MP);

% ---- SORT BY LEARNING
[Dat_learning_sorted, inds]=sort(Dat_learning);
Dat_MP_sorted=Dat_MP(inds);
Dat_AFP_sorted=Dat_AFP(inds);

lt_plot_bar(X-0.1, Dat_MP_sorted, {'Color', 'r'})
lt_plot_bar(X+0.1, Dat_AFP_sorted, {'Color', 'b'})
plot(Dat_learning_sorted, '-k')

legend('MP learning', 'Learning')

% =========== 9.2) [diff]Plot AFP BIAS, MP [sorted by LEARNING, SEPARATING MP,
% AFP, AND LEARNING] - PRESYL SIMILAR
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('AFP learning (DIFF SYLS, PRESYL SIMILAR, SORTED BY LEARNING)');
xlabel('from most negative to most positive genearlization');
ylabel('hz, rel learning by targ');

% --- ALL SYLS
inds=TargStatus_all==0 & PreSylSimilar_all==1 & SimDiff_all==0;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);
Dat_learning=Learning_all(inds);
X=1:length(Dat_MP);

% ---- SORT BY LEARNING
[Dat_learning_sorted, inds]=sort(Dat_learning);
Dat_MP_sorted=Dat_MP(inds);
Dat_AFP_sorted=Dat_AFP(inds);

lt_plot_bar(X-0.1, Dat_MP_sorted, {'Color', 'r'})
lt_plot_bar(X+0.1, Dat_AFP_sorted, {'Color', 'b'})
plot(Dat_learning_sorted, '-k')

legend('MP learning', 'Learning')

% =========== 9.2) [simialr]Plot AFP BIAS, MP [sorted by LEARNING, SEPARATING MP,
% AFP, AND LEARNING] - PRESYL SIMILAR
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('AFP learning (SIMILAR SYLS, PRESYL SIMILAR, SORTED BY LEARNING)');
xlabel('from most negative to most positive genearlization');
ylabel('hz, rel learning by targ');

% --- ALL SYLS
inds=TargStatus_all==0 & PreSylSimilar_all==1 & SimDiff_all==1;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);
Dat_learning=Learning_all(inds);
X=1:length(Dat_MP);

% ---- SORT BY LEARNING
[Dat_learning_sorted, inds]=sort(Dat_learning);
Dat_MP_sorted=Dat_MP(inds);
Dat_AFP_sorted=Dat_AFP(inds);

lt_plot_bar(X-0.1, Dat_MP_sorted, {'Color', 'r'})
lt_plot_bar(X+0.1, Dat_AFP_sorted, {'Color', 'b'})
plot(Dat_learning_sorted, '-k')

legend('MP learning', 'Learning')


% =========== 9.2) Plot AFP BIAS, MP [sorted by LEARNING, SEPARATING MP,
% AFP, AND LEARNING] - PRESYL DIFF
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('AFP learning (PRESYL DIFF, SORTED BY LEARNING)');
xlabel('from most negative to most positive genearlization');
ylabel('hz, rel learning by targ');

% --- ALL SYLS
inds=TargStatus_all==0 & PreSylSimilar_all==0;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);
Dat_learning=Learning_all(inds);
X=1:length(Dat_MP);

% ---- SORT BY LEARNING
[Dat_learning_sorted, inds]=sort(Dat_learning);
Dat_MP_sorted=Dat_MP(inds);
Dat_AFP_sorted=Dat_AFP(inds);

lt_plot_bar(X-0.1, Dat_MP_sorted, {'Color', 'r'})
lt_plot_bar(X+0.1, Dat_AFP_sorted, {'Color', 'b'})

plot(Dat_learning_sorted, '-k')

legend('MP learning', 'Learning')

% =========== 9.2) [SIMILAR] Plot AFP BIAS, MP [sorted by LEARNING, SEPARATING MP,
% AFP, AND LEARNING] - PRESYL DIFF
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[SIMILAR] (PRESYL DIFF, SORTED BY LEARNING)');
xlabel('from most negative to most positive genearlization');
ylabel('hz, rel learning by targ');

% --- ALL SYLS
inds=TargStatus_all==0 & PreSylSimilar_all==0 & SimDiff_all==1;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);
Dat_learning=Learning_all(inds);
X=1:length(Dat_MP);

% ---- SORT BY LEARNING
[Dat_learning_sorted, inds]=sort(Dat_learning);
Dat_MP_sorted=Dat_MP(inds);
Dat_AFP_sorted=Dat_AFP(inds);

lt_plot_bar(X-0.1, Dat_MP_sorted, {'Color', 'r'})
lt_plot_bar(X+0.1, Dat_AFP_sorted, {'Color', 'b'})

plot(Dat_learning_sorted, '-k')

legend('MP learning', 'Learning')

% =========== 9.2) [DIFF] Plot AFP BIAS, MP [sorted by LEARNING, SEPARATING MP,
% AFP, AND LEARNING] - PRESYL DIFF
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[DIFF] (PRESYL DIFF, SORTED BY LEARNING)');
xlabel('from most negative to most positive genearlization');
ylabel('hz, rel learning by targ');

% --- ALL SYLS
inds=TargStatus_all==0 & PreSylSimilar_all==0 & SimDiff_all==0;

Dat_MP=MPlearning_all(inds);
Dat_AFP=AFPbias_all(inds);
Dat_learning=Learning_all(inds);
X=1:length(Dat_MP);

% ---- SORT BY LEARNING
[Dat_learning_sorted, inds]=sort(Dat_learning);
Dat_MP_sorted=Dat_MP(inds);
Dat_AFP_sorted=Dat_AFP(inds);

lt_plot_bar(X-0.1, Dat_MP_sorted, {'Color', 'r'})
lt_plot_bar(X+0.1, Dat_AFP_sorted, {'Color', 'b'})

plot(Dat_learning_sorted, '-k')

legend('MP learning', 'Learning')


        % ===== GLOBAL
        lt_subtitle('All values norm to learning(PBS) by targ');


