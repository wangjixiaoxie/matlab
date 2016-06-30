function [PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANbidir(SeqDepPitch_AcrossBirds, PARAMS, use_final_extracted_windows, NormToTarg)
%% LT 8/4/15 - PLOTS bidir learning experiments for MUSCIMOL inactivation

%% SORT OUT ONLY THE THE EXPEIRMENTS THAT HAVE LMAN INACTIVATION DATA

% copy strcuture, save backup.
SeqDepPitch_AcrossBirds_ORIGINAL=SeqDepPitch_AcrossBirds;

filter = 'LMAN';
[SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);


%% PARAMS

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

TotalNumExperiments=0;
for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
TotalNumExperiments=TotalNumExperiments+numexperiments;
end


%% REMOVE SYLS THAT SHOULD NOT BE ANALYZED (I.E O/L WITH WN, for non-catch analyses)

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        % Check whether this birdname expt pair has syl that should be removed.
        % If so, remove them from unique syls
        
        inds=find(strcmp(PARAMS.global.LMAN.SylsToRemove_BiDir, birdname));
        
        for j=1:length(inds);
            
            expt_toremove=PARAMS.global.LMAN.SylsToRemove_BiDir{inds(j)+1};
            syls_toremove=PARAMS.global.LMAN.SylsToRemove_BiDir{inds(j)+2};
            
            % IF CURRENT EXPERIMENT IS THE ONE TO REMOVE, THEN DO SO
            if strcmp(exptname, expt_toremove);
                
                for k=1:length(syls_toremove);
                    
                    tmp_sylremove=syls_toremove{k};
                    
                    syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
                    
                    ind_to_remove=strcmp(syls_unique, tmp_sylremove);
                    
                    if ~any(ind_to_remove)
                        % then syl not found...
                        disp('Problem - syl not found for inds to remove');
                        continue; 
                    end
                    
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
        
        % ==== MAKE SURE THIS HAS BIDIR LEARNING
        if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'days_bidir_early');
            continue
        end
                
        % ==== ONE FIGURE PER EXPT
        lt_figure; hold on;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        tmp=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls;
        targsyl_other=tmp{~strcmp(tmp,targsyl)}; % find the other syl assuming that multidirsyls contains the target and the other target
        
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


%% ========== PLOT segragating into targ2, targ1, and other syls
flip_sign_if_neg=1; % if neg learning by targ

lt_figure; hold on;
count=1;
SubplotsPerFig=6;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        
        % ==== MAKE SURE THIS HAS BIDIR LEARNING
        if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'days_bidir_early');
            continue
        end
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        tmp=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls;
        targsyl_other=tmp{~strcmp(tmp,targsyl)}; % find the other syl assuming that multidirsyls contains the target and the other target
        
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title([birdname '-' exptname]);

        % -- to collect vals across syls
            Y_Learning_AllOthers=[];
                Y_MP_AllOthers=[];
                Y_AFP_AllOthers=[];
            
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            
            if use_final_extracted_windows==1;
                EpochNameList={'days_consolid_late', 'final_extracted_window_bidir'};
                EpochNameList={'final_extracted_window', 'final_extracted_window_bidir'};
            else
                % USE BIDIR LATE. If don't have, then use bidir early
                if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, 'days_bidir_late');
                    EpochNameList={'days_consolid_late', 'days_bidir_late'};
                else
                    EpochNameList={'days_consolid_late', 'days_bidir_early'};
                end
            end
            
            
            X=1:length(EpochNameList);
            Y_Learning=[]; % i.e. for this syl, one entry for each epoch.
            Y_MP=[];
            Y_AFP=[];
            
            
            for k=1:length(EpochNameList);
                epochfield=EpochNameList{k};
                
                % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                % MUSC)
                FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                
                FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                
                AFP_bias=FF_PBS-FF_MUSC;
                
                issimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
                istarget=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
                
                if strcmp(targsyl_other, syl);
                    istarg_other=1;
                else
                    istarg_other=0;
                end
                
                % ===== COLLECT FOR THIS SYL
                Y_Learning=[Y_Learning FF_PBS]; % i.e. for this syl, one entry for each epoch.
                Y_MP=[Y_MP FF_MUSC];
                Y_AFP=[Y_AFP AFP_bias];
                issimilar;
                istarget;
                istarg_other;
                
            end
            
                % ====== FLIP SIGN 
                if flip_sign_if_neg==1;
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.days_consolid_late.(targsyl).meanFF_pbs<0; % targ syl <0 learning
                    Y_Learning=-1.*Y_Learning;
                    Y_MP=-1.*Y_MP;
                    Y_AFP=-1.*Y_AFP;
                end
                end
            
            % ===== PLOT THIS SYL
            if istarget==1
%                 plot(X, Y_Learning, 'o-k')
%                 plot(X+length(X), Y_MP, 'o-r')
%                 plot(X+2*length(X), Y_AFP, 'o-b')
                
                lt_plot_bar(X, Y_Learning, {'Color','k'})
                lt_plot_bar(X+1+length(X), Y_MP, {'Color','r'})
                lt_plot_bar(X+2+2*length(X), Y_AFP, {'Color','b'})
                
            elseif istarg_other==1
%                 plot(X+9, Y_Learning, 'o-k')
%                 plot(X+10+length(X), Y_MP, 'o-r')
%                 plot(X+11+2*length(X), Y_AFP, 'o-b')
%                 
                lt_plot_bar(X+9, Y_Learning, {'Color','k'})
                lt_plot_bar(X+10+length(X), Y_MP, {'Color','r'})
                lt_plot_bar(X+11+2*length(X), Y_AFP, {'Color','b'})
                
            else % is all other syls
                plot(X+19, Y_Learning, 'o-k')
                plot(X+20+length(X), Y_MP, 'o-r')
                plot(X+21+2*length(X), Y_AFP, 'o-b')
                
                 % ----- COLLECT TO PLOT MEAN
                Y_Learning_AllOthers=[Y_Learning_AllOthers; Y_Learning];
                Y_MP_AllOthers=[Y_MP_AllOthers; Y_MP];
                Y_AFP_AllOthers=[Y_AFP_AllOthers; Y_AFP];
                
            end
        end
            % === PLOT MEAN (for other syls)
                lt_plot_bar(X+19, mean(Y_Learning_AllOthers), {'Color','k'})
                lt_plot_bar(X+20+length(X), mean(Y_MP_AllOthers), {'Color','r'})
                lt_plot_bar(X+21+2*length(X), mean(Y_AFP_AllOthers), {'Color','b'})

                % lines
                line([X(2)+2+2*length(X)+1, X(2)+2+2*length(X)+1], ylim);
                line([X(2)+11+2*length(X)+1, X(2)+11+2*length(X)+1], ylim);
                xlim([0 28]);
                
                % text as legend
                Ylim=ylim;
                Yrange=Ylim(2)-Ylim(1);
                
                text(X(2), Ylim(2)-Yrange/10, 'Target','FontSize',13,'FontWeight','bold')
                text(X(2)+9, Ylim(2)-Yrange/10, 'New Target','FontSize',13,'FontWeight','bold')
                text(X(2)+19, Ylim(2)-Yrange/10, 'All Other Syls','FontSize',13,'FontWeight','bold')
                
        % ===== GLOBAL
        lt_plot_zeroline
    end
end

lt_subtitle('Learning(black), MP(red), AFP(blue), "consolidation end" vs "bidir end (or early if no end)"')


%% ======================== PLOT segragating into targ2, targ1, and other syls
%==== ALSO PLOT ALL EXPERIMENTS IN ONE PLOT
lt_figure; hold on;
count=1;
SubplotsPerFig=6;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        
        % ==== MAKE SURE THIS HAS BIDIR LEARNING
        if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'days_bidir_early');
            continue
        end
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        tmp=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls;
        targsyl_other=tmp{~strcmp(tmp,targsyl)}; % find the other syl assuming that multidirsyls contains the target and the other target
        
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title([birdname '-' exptname]);
        ylabel('norm to learning at targ');
        
        % -- to collect vals across syls
        Y_Learning_AllOthers=[];
        Y_MP_AllOthers=[];
        Y_AFP_AllOthers=[];
        
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            
            if use_final_extracted_windows==1;
                EpochNameList={'final_extracted_window', 'final_extracted_window_bidir'};
            else
                % USE BIDIR LATE. If don't have, then use bidir early
                if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, 'days_bidir_late');
                    EpochNameList={'days_consolid_late', 'days_bidir_late'};
                else
                    EpochNameList={'days_consolid_late', 'days_bidir_early'};
                end
            end
            
                        
            X=1:length(EpochNameList);
            Y_Learning=[]; % i.e. for this syl, one entry for each epoch.
            Y_MP=[];
            Y_AFP=[];
            
            
            for k=1:length(EpochNameList);
                epochfield=EpochNameList{k};
                
                % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                % MUSC)
                FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                
                FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                
                AFP_bias=FF_PBS-FF_MUSC;
                
                issimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
                istarget=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
                if strcmp(targsyl_other, syl);
                    istarg_other=1;
                else
                    istarg_other=0;
                end
                
                % ===== COLLECT FOR THIS SYL
                Y_Learning=[Y_Learning FF_PBS]; % i.e. for this syl, one entry for each epoch.
                Y_MP=[Y_MP FF_MUSC];
                Y_AFP=[Y_AFP AFP_bias];
                issimilar;
                istarget;
                istarg_other;
                
            end
            
            % ========== Normalize to target learning (PBS, consolid end);
            learning_at_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.days_consolid_late.(targsyl).meanFF_pbs;
            
            Y_Learning=Y_Learning./learning_at_targ;
            Y_MP=Y_MP./learning_at_targ;
            Y_AFP=Y_AFP./learning_at_targ;
            
%             % ====== FLIP SIGN
%             if flip_sign_if_neg==1;
%                 if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.days_consolid_late.(targsyl).meanFF_pbs<0; % targ syl <0 learning
%                     Y_Learning=-1.*Y_Learning;
%                     Y_MP=-1.*Y_MP;
%                     Y_AFP=-1.*Y_AFP;
%                 end
%             end
            

            % ===== PLOT THIS SYL
            if istarget==1
                %                 plot(X, Y_Learning, 'o-k')
                %                 plot(X+length(X), Y_MP, 'o-r')
                %                 plot(X+2*length(X), Y_AFP, 'o-b')
                
                lt_plot_bar(X, Y_Learning, {'Color','k'})
                lt_plot_bar(X+1+length(X), Y_MP, {'Color','r'})
                lt_plot_bar(X+2+2*length(X), Y_AFP, {'Color','b'})
                
            elseif istarg_other==1
                %                 plot(X+9, Y_Learning, 'o-k')
                %                 plot(X+10+length(X), Y_MP, 'o-r')
                %                 plot(X+11+2*length(X), Y_AFP, 'o-b')
                %
                lt_plot_bar(X+9, Y_Learning, {'Color','k'})
                lt_plot_bar(X+10+length(X), Y_MP, {'Color','r'})
                lt_plot_bar(X+11+2*length(X), Y_AFP, {'Color','b'})
                
            else % is all other syls
                plot(X+19, Y_Learning, 'o-k')
                plot(X+20+length(X), Y_MP, 'o-r')
                plot(X+21+2*length(X), Y_AFP, 'o-b')
                
                % ----- COLLECT TO PLOT MEAN
                Y_Learning_AllOthers=[Y_Learning_AllOthers; Y_Learning];
                Y_MP_AllOthers=[Y_MP_AllOthers; Y_MP];
                Y_AFP_AllOthers=[Y_AFP_AllOthers; Y_AFP];
                
            end
        end
        % === PLOT MEAN (for other syls)
        lt_plot_bar(X+19, mean(Y_Learning_AllOthers), {'Color','k'})
        lt_plot_bar(X+20+length(X), mean(Y_MP_AllOthers), {'Color','r'})
        lt_plot_bar(X+21+2*length(X), mean(Y_AFP_AllOthers), {'Color','b'})
        
        % lines
        line([X(2)+2+2*length(X)+1, X(2)+2+2*length(X)+1], ylim);
        line([X(2)+11+2*length(X)+1, X(2)+11+2*length(X)+1], ylim);
        xlim([0 28]);
        
        % text as legend
        Ylim=ylim;
        Yrange=Ylim(2)-Ylim(1);
        
        text(X(2), Ylim(2)-Yrange/10, 'Target','FontSize',13,'FontWeight','bold')
        text(X(2)+9, Ylim(2)-Yrange/10, 'New Target','FontSize',13,'FontWeight','bold')
        text(X(2)+19, Ylim(2)-Yrange/10, 'All Other Syls','FontSize',13,'FontWeight','bold')
        
        % ===== GLOBAL
        lt_plot_zeroline
    end
end



%% ===== [COLLECT] PLOT ALL IN ONE PLOT
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('All experiments (norm. to target learning (consol end)');
% -- to collect vals across syls/expts
Y_Learning_AllOthers=[];
Y_MP_AllOthers=[];
Y_AFP_AllOthers=[];
Y_Learning_Target=[];
Y_MP_Target=[];
Y_AFP_Target=[];
Y_Learning_NewTarg=[];
Y_MP_NewTarg=[];
Y_AFP_NewTarg=[];

Y_ExptNum_Target=[];
Y_ExptNum_NewTarg=[];
Y_ExptNum_AllOthers=[];

Y_similar_NewTarg=[];
Y_presimilar_NewTarg=[];
Y_DistFromTarg_NewTarg=[];

plotcols=lt_make_plot_colors(TotalNumExperiments, 0, 0);
cc=1;

ExptCounter=1;

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        
        % ==== MAKE SURE THIS HAS BIDIR LEARNING
        if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'days_bidir_early');
            continue
        end
        
        
                
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        tmp=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls;
        targsyl_other=tmp{~strcmp(tmp,targsyl)}; % find the other syl assuming that multidirsyls contains the target and the other target
        rng('shuffle');
        plotcol=plotcols{cc};
        cc=cc+1;
        
            % == Only continue if both targsyl and targsyl_other are part of
            % SylsUnique
            if ~any(strcmp(SylsUnique, targsyl)) | ~any(strcmp(SylsUnique, targsyl_other));
                disp(['Skipping ' birdname '-' exptname ' - one or other targsyl not part of SylsUnique']);
            continue
            end
            
            
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            
            if use_final_extracted_windows==1;
                EpochNameList={'final_extracted_window', 'final_extracted_window_bidir'};
            else
                % USE BIDIR LATE. If don't have, then use bidir early
                if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, 'days_bidir_late');
                    EpochNameList={'days_consolid_late', 'days_bidir_late'};
                else
                    EpochNameList={'days_consolid_late', 'days_bidir_early'};
                end
            end
            
            
                        X=1:length(EpochNameList);
            Y_Learning=[]; % i.e. for this syl, one entry for each epoch.
            Y_MP=[];
            Y_AFP=[];
            
            
            for k=1:length(EpochNameList);
                epochfield=EpochNameList{k};
                
                % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                % MUSC)
                FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                
                FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                
                AFP_bias=FF_PBS-FF_MUSC;
                
                
                % ===== COLLECT FOR THIS SYL
                Y_Learning=[Y_Learning FF_PBS]; % i.e. for this syl, one entry for each epoch.
                Y_MP=[Y_MP FF_MUSC];
                Y_AFP=[Y_AFP AFP_bias];
%                 Y_similar=[Y_similar issimilar]
%                 Y_targetFirst=[Y_targetFirst istarget];
%                 Y_targetSecond=
%                 istarg_other;
                
            end
            
            % ====================
%             if strcmp(birdname, 'wh25pk77') & strcmp(exptname, 'SeqDepPitchLMAN')
%                 keyboard
%             end
            
                issimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
                presimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
                istarget=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
                if strcmp(targsyl_other, syl);
                    istarg_other=1;
                else
                    istarg_other=0;
                end
                
                distanceFromTarg=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).distance_from_targ;
                
                
                % ========== Normalize to target learning (PBS, consolid end);
                learning_at_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.days_consolid_late.(targsyl).meanFF_pbs;
                if NormToTarg==1;
                    Y_Learning=Y_Learning./learning_at_targ;
                    Y_MP=Y_MP./learning_at_targ;
                    Y_AFP=Y_AFP./learning_at_targ;
                else
                    % then don't norm, but flip sign if learn is neg at targ
                    
                    targLearnDir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
                    
                    Y_Learning=targLearnDir.*Y_Learning;
                    Y_MP=targLearnDir.*Y_MP;
                    Y_AFP=targLearnDir.*Y_AFP;
                    
                    
                end
            
            
            % ===== PLOT THIS SYL
            if istarget==1
                plot(X, Y_Learning, 'o-', 'Color',plotcol)
                plot(X+1+length(X), Y_MP, 'o-', 'Color',plotcol)
                plot(X+2+2*length(X), Y_AFP, 'o-', 'Color',plotcol)
                
                % ----- COLLECT TO PLOT MEAN
                Y_Learning_Target=[Y_Learning_Target; Y_Learning];
                Y_MP_Target=[Y_MP_Target; Y_MP];
                Y_AFP_Target=[Y_AFP_Target; Y_AFP];
                Y_ExptNum_Target=[Y_ExptNum_Target ExptCounter];
                
            elseif istarg_other==1
                plot(X+9, Y_Learning, 'o-', 'Color',plotcol)
                plot(X+10+length(X), Y_MP, 'o-', 'Color',plotcol)
                plot(X+11+2*length(X), Y_AFP, 'o-', 'Color',plotcol)
                
                % ----- COLLECT TO PLOT MEAN
                Y_Learning_NewTarg=[Y_Learning_NewTarg; Y_Learning];
                Y_MP_NewTarg=[Y_MP_NewTarg; Y_MP];
                Y_AFP_NewTarg=[Y_AFP_NewTarg; Y_AFP];
                Y_ExptNum_NewTarg=[Y_ExptNum_NewTarg ExptCounter];
                Y_similar_NewTarg=[Y_similar_NewTarg issimilar];
                Y_presimilar_NewTarg=[Y_presimilar_NewTarg presimilar];
                Y_DistFromTarg_NewTarg=[Y_DistFromTarg_NewTarg distanceFromTarg];

                
            else % is all other syls
                plot(X+19, Y_Learning, 'o-', 'Color',plotcol)
                plot(X+20+length(X), Y_MP, 'o-', 'Color',plotcol)
                plot(X+21+2*length(X), Y_AFP, 'o-', 'Color',plotcol)

                
                % ----- COLLECT TO PLOT MEAN
                Y_Learning_AllOthers=[Y_Learning_AllOthers; Y_Learning];
                Y_MP_AllOthers=[Y_MP_AllOthers; Y_MP];
                Y_AFP_AllOthers=[Y_AFP_AllOthers; Y_AFP];
                Y_ExptNum_AllOthers=[Y_ExptNum_AllOthers ExptCounter];
                
            end
            
            
            
        end
        ExptCounter=ExptCounter+1;
    end
end
% === PLOT MEAN
% -- target syl
lt_plot_bar(X, mean(Y_Learning_Target), {'Color','k'})
lt_plot_bar(X+1+length(X), mean(Y_MP_Target), {'Color','r'})
lt_plot_bar(X+2+2*length(X), mean(Y_AFP_Target), {'Color','b'})

% --- new target
lt_plot_bar(X+9, mean(Y_Learning_NewTarg), {'Color','k'})
lt_plot_bar(X+10+length(X), mean(Y_MP_NewTarg), {'Color','r'})
lt_plot_bar(X+11+2*length(X), mean(Y_AFP_NewTarg), {'Color','b'})

% ---- for other syls
lt_plot_bar(X+19, mean(Y_Learning_AllOthers), {'Color','k'})
lt_plot_bar(X+20+length(X), mean(Y_MP_AllOthers), {'Color','r'})
lt_plot_bar(X+21+2*length(X), mean(Y_AFP_AllOthers), {'Color','b'})

% lines
line([X(2)+2+2*length(X)+1, X(2)+2+2*length(X)+1], ylim);
line([X(2)+11+2*length(X)+1, X(2)+11+2*length(X)+1], ylim);
xlim([0 28]);

% text as legend
Ylim=ylim;
Yrange=Ylim(2)-Ylim(1);

text(X(2), Ylim(2)-Yrange/10, 'Target','FontSize',13,'FontWeight','bold')
text(X(2)+9, Ylim(2)-Yrange/10, 'New Target','FontSize',13,'FontWeight','bold')
text(X(2)+19, Ylim(2)-Yrange/10, 'All Other Syls','FontSize',13,'FontWeight','bold')

% ===== GLOBAL
lt_plot_zeroline



% ================================ PLOT SEPARATION (comparing just new
% target to the old target)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Old targ minus new targ (separation, all values normalized to target early learning)');

% all defined as target minus new target
Separation_Learning=Y_Learning_Target-Y_Learning_NewTarg;
Separation_MP=Y_MP_Target-Y_MP_NewTarg;
Separation_AFP=Y_AFP_Target-Y_AFP_NewTarg;

% --- PLOT
numexpts=size(Separation_Learning,2);
for i=1:numexpts;
    plot(Separation_Learning(i, :), 'o-k');
    plot(Separation_MP(i, :), 'o-r');
    plot(Separation_AFP(i, :), 'o-b');
end 
xlim([0 3]);
lt_plot_zeroline;

  
% ============= PLOT MEANS OF SEPARATIONS
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Old targ minus new targ (separation, Mean, all values normalized to target early learning)');
ylabel('separation, using values norm to learning by targ at start');

X=[1 2];

% learning
Ymean=mean(Separation_Learning,1);
Ysem=lt_sem(Separation_Learning);
errorbar(X, Ymean, Ysem, '-sk');

% MP
Ymean=mean(Separation_MP,1);
Ysem=lt_sem(Separation_MP);
errorbar(X, Ymean, Ysem, '-sr');

% AFP
Ymean=mean(Separation_AFP,1);
Ysem=lt_sem(Separation_AFP);
errorbar(X, Ymean, Ysem, '-sb');

xlim([0 3]);
lt_plot_zeroline

lt_subtitle('Learning(black), MP(red), AFP(blue), "consolidation end" vs "bidir end"')


% ================================ PLOT SEPARATION DIFFERENT WAY (comparing just new
% target to the old target)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Separation between new targ and old (consolid end vs. bidir end)');
ylabel('separation, using values norm to learning by targ at start');

% all defined as target minus new target
Separation_Learning=Y_Learning_Target-Y_Learning_NewTarg;
Separation_MP=Y_MP_Target-Y_MP_NewTarg;
Separation_AFP=Y_AFP_Target-Y_AFP_NewTarg;

% --- PLOT
numexpts=size(Separation_Learning,2);
for i=1:numexpts;
    
    % early
    plot(1, Separation_Learning(i, 1), 'o-k');
    plot(2, Separation_MP(i, 1), 'o-r');
    plot(3, Separation_AFP(i, 1), 'o-b');
    
    % late
    plot(5, Separation_Learning(i, 2), 'o-k');
    plot(6, Separation_MP(i, 2), 'o-r');
    plot(7, Separation_AFP(i, 2), 'o-b');
    
end 

% plot means
    lt_plot_bar(1, mean(Separation_Learning(:, 1)), {'Color','k'});
    lt_plot_bar(2, mean(Separation_MP(:, 1)), {'Color','r'});
    lt_plot_bar(3, mean(Separation_AFP(:, 1)), {'Color','b'});

    lt_plot_bar(5, mean(Separation_Learning(:, 2)), {'Color','k'});
    lt_plot_bar(6, mean(Separation_MP(:, 2)), {'Color','r'});
    lt_plot_bar(7, mean(Separation_AFP(:, 2)), {'Color','b'});
    
set(gca, 'XTick', [2 6]);
set(gca, 'XTickLabel', {'before bidir', 'end bidir'})
xlim([0 8]);
lt_plot_zeroline;


lt_subtitle('Learning(black), MP(red), AFP(blue), "consolidation end" vs "bidir end"')

%% +++++++++++++++ NEW PLOTS - focusing on AFP bias
%% % ====  [AFP BIAS]
count=1;
SubplotsPerFig=8;
subplotrows=2;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('AFP BIAS')
ylabel('pitch shift, norm to target learning (PBS)');

plotcols=lt_make_plot_colors(ExptCounter-1, 0, 0);
Yexptnums={Y_ExptNum_Target, Y_ExptNum_NewTarg};

% =================== SINGLE DIR PHASE
X=1; % location to plot
ind=1; % pre-start BIDIR

% --
Yraw={}; % each cell is one category of syllable
Ymean=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_AFP_Target;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_AFP_NewTarg;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));
% 
% % --- All nontargets
% x=3;
% DatArray=Y_AFP_AllOthers;
% 
% Yraw{x}=DatArray(:,ind);
% Ymean(x)=mean(DatArray(:,ind));
% Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
% 1) raw vals

for i=1:length(Yraw);
    for ii=1:ExptCounter-1;
    indsinds=Yexptnums{i}==ii;
    plot(i+0.2, Yraw{i}(indsinds), 'o' ,'Color',plotcols{ii});
    end
end

% 2) means
xx=X:X+1;
hbar=lt_plot_bar(xx, Ymean, {'Errors', Ysem, 'Color' 'k'});

Xlabel={'First Targ', 'Second Targ'};
set(gca, 'XTick', xx);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);




% =============================== DURING TWOP DIR PHASE
X=4;
ind=2; % pre-start BIDIR

Yraw={}; % each cell is one category of syllable
Ymean=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_AFP_Target;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_AFP_NewTarg;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));
% 
% % --- All nontargets
% x=3;
% DatArray=Y_AFP_AllOthers;
% 
% Yraw{x}=DatArray(:,ind);
% Ymean(x)=mean(DatArray(:,ind));
% Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
% 1) raw vals

for i=1:length(Yraw);
    for ii=1:ExptCounter-1;
    indsinds=Yexptnums{i}==ii;
    plot(i+X-1+0.2, Yraw{i}(indsinds), 'o' ,'Color',plotcols{ii});
    end
end

% 2) means
xx=X:X+1;
hbar=lt_plot_bar(xx, Ymean, {'Errors', Ysem, 'Color' 'k'});

Xlabel={'First Targ', 'Second Targ'};
set(gca, 'XTick', xx);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);





%% % ====  [LEARNING]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('LEARNING')
ylabel('pitch shift, norm to target learning (PBS)');

plotcols=lt_make_plot_colors(ExptCounter-1, 0, 0);
Yexptnums={Y_ExptNum_Target, Y_ExptNum_NewTarg};

% =================== SINGLE DIR PHASE
X=1; % location to plot
ind=1; % pre-start BIDIR

% --
Yraw={}; % each cell is one category of syllable
Ymean=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_Learning_Target;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_Learning_NewTarg;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));
% 
% % --- All nontargets
% x=3;
% DatArray=Y_AFP_AllOthers;
% 
% Yraw{x}=DatArray(:,ind);
% Ymean(x)=mean(DatArray(:,ind));
% Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
% 1) raw vals

for i=1:length(Yraw);
    for ii=1:ExptCounter-1;
    indsinds=Yexptnums{i}==ii;
    plot(i+0.2, Yraw{i}(indsinds), 'o' ,'Color',plotcols{ii});
    end
end

% 2) means
xx=X:X+1;
hbar=lt_plot_bar(xx, Ymean, {'Errors', Ysem, 'Color' 'k'});

Xlabel={'First Targ', 'Second Targ'};
set(gca, 'XTick', xx);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);




% =============================== DURING TWOP DIR PHASE
X=4;
ind=2; % pre-start BIDIR

Yraw={}; % each cell is one category of syllable
Ymean=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_Learning_Target;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_Learning_NewTarg;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));
% 
% % --- All nontargets
% x=3;
% DatArray=Y_AFP_AllOthers;
% 
% Yraw{x}=DatArray(:,ind);
% Ymean(x)=mean(DatArray(:,ind));
% Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
% 1) raw vals

for i=1:length(Yraw);
    for ii=1:ExptCounter-1;
    indsinds=Yexptnums{i}==ii;
    plot(i+X-1+0.2, Yraw{i}(indsinds), 'o' ,'Color',plotcols{ii});
    end
end

% 2) means
xx=X:X+1;
hbar=lt_plot_bar(xx, Ymean, {'Errors', Ysem, 'Color' 'k'});

Xlabel={'First Targ', 'Second Targ'};
set(gca, 'XTick', xx);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);


%% % ====  [MP BIAS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('MP BIAS')
ylabel('pitch shift, norm to target learning (PBS)');

plotcols=lt_make_plot_colors(ExptCounter-1, 0, 0);
Yexptnums={Y_ExptNum_Target, Y_ExptNum_NewTarg};

% =================== SINGLE DIR PHASE
X=1; % location to plot
ind=1; % pre-start BIDIR

% --
Yraw={}; % each cell is one category of syllable
Ymean=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_MP_Target;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_MP_NewTarg;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));
% 
% % --- All nontargets
% x=3;
% DatArray=Y_AFP_AllOthers;
% 
% Yraw{x}=DatArray(:,ind);
% Ymean(x)=mean(DatArray(:,ind));
% Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
% 1) raw vals

for i=1:length(Yraw);
    for ii=1:ExptCounter-1;
    indsinds=Yexptnums{i}==ii;
    plot(i+0.2, Yraw{i}(indsinds), 'o' ,'Color',plotcols{ii});
    end
end

% 2) means
xx=X:X+1;
hbar=lt_plot_bar(xx, Ymean, {'Errors', Ysem, 'Color' 'k'});

Xlabel={'First Targ', 'Second Targ'};
set(gca, 'XTick', xx);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);




% =============================== DURING TWOP DIR PHASE
X=4;
ind=2; % pre-start BIDIR

Yraw={}; % each cell is one category of syllable
Ymean=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_MP_Target;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_MP_NewTarg;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));
% 
% % --- All nontargets
% x=3;
% DatArray=Y_AFP_AllOthers;
% 
% Yraw{x}=DatArray(:,ind);
% Ymean(x)=mean(DatArray(:,ind));
% Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
% 1) raw vals

for i=1:length(Yraw);
    for ii=1:ExptCounter-1;
    indsinds=Yexptnums{i}==ii;
    plot(i+X-1+0.2, Yraw{i}(indsinds), 'o' ,'Color',plotcols{ii});
    end
end

% 2) means
xx=X:X+1;
hbar=lt_plot_bar(xx, Ymean, {'Errors', Ysem, 'Color' 'k'});

Xlabel={'First Targ', 'Second Targ'};
set(gca, 'XTick', xx);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);



%% like above, but more clunky, also plotted nontargets on same plot

if (0)
lt_seq_dep_pitch_ACROSSBIRDS_LMANbidir_junk
end


%% === LIKE ABOVE, BUT SAME TYPE ONLY
%% % ====  [AFP BIAS]
count=1;
SubplotsPerFig=8;
subplotrows=2;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[SAME] AFP BIAS')
ylabel('pitch shift, norm to target learning (PBS)');

plotcols=lt_make_plot_colors(ExptCounter-1, 0, 0);
Yexptnums={Y_ExptNum_Target, Y_ExptNum_NewTarg};

% =================== SINGLE DIR PHASE
X=1; % location to plot
ind=1; % pre-start BIDIR

% --
Yraw={}; % each cell is one category of syllable
Ymean=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_AFP_Target;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);


Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_AFP_NewTarg;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
% 1) raw vals
for i=1:length(Yraw{1});
    plot([X+0.2 X+1.2], [Yraw{1} Yraw{2}], 'o-');
end


% 2) means
xx=X:X+1;
hbar=lt_plot_bar(xx, Ymean, {'Errors', Ysem, 'Color' 'k'});

Xlabel={'First Targ', 'Second Targ'};
set(gca, 'XTick', xx);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);

% --- save Yraw for stats
YrawStats1=Yraw;



% =============================== DURING TWOP DIR PHASE
X=4;
ind=2; % pre-start BIDIR

% --
Yraw={}; % each cell is one category of syllable
Ymean=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_AFP_Target;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);


Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_AFP_NewTarg;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));



% ------------------ PLOT BARS
% 1) raw vals
for i=1:length(Yraw{1});
    plot([X+0.2 X+1.2], [Yraw{1} Yraw{2}], 'o-');
end


% 2) means
xx=X:X+1;
hbar=lt_plot_bar(xx, Ymean, {'Errors', Ysem, 'Color' 'k'});

% --- save Yraw for stats
YrawStats2=Yraw;

% ====== PERFORM STATS (change across phase?)
[~, p]=ttest(YrawStats1{1}, YrawStats2{1}); % first targ (one targ phase vs 2 targ phase)
disp(['First targ (one phase vs. two phase), p='  num2str(p)]);
[~, p]=ttest(YrawStats1{2}, YrawStats2{2}); % second targ (one targ phase vs 2 targ phase)
disp(['Second targ (one phase vs. two phase), p='  num2str(p)]);



Xlabel={'First Targ', 'Second Targ'};
set(gca, 'XTick', xx);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);



%% % ====  [LEARNING]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[SAME] LAERNING')
ylabel('pitch shift, norm to target learning (PBS)');

plotcols=lt_make_plot_colors(ExptCounter-1, 0, 0);
Yexptnums={Y_ExptNum_Target, Y_ExptNum_NewTarg};

% =================== SINGLE DIR PHASE
X=1; % location to plot
ind=1; % pre-start BIDIR

% --
Yraw={}; % each cell is one category of syllable
Ymean=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_Learning_Target;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);


Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_Learning_NewTarg;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
% 1) raw vals
for i=1:length(Yraw{1});
    plot([X+0.2 X+1.2], [Yraw{1} Yraw{2}], 'o-');
end


% 2) means
xx=X:X+1;
hbar=lt_plot_bar(xx, Ymean, {'Errors', Ysem, 'Color' 'k'});

Xlabel={'First Targ', 'Second Targ'};
set(gca, 'XTick', xx);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);

% --- save Yraw for stats
YrawStats1=Yraw;


% =============================== DURING TWOP DIR PHASE
X=4;
ind=2; % pre-start BIDIR

% --
Yraw={}; % each cell is one category of syllable
Ymean=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_Learning_Target;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);


Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_Learning_NewTarg;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
% 1) raw vals
for i=1:length(Yraw{1});
    plot([X+0.2 X+1.2], [Yraw{1} Yraw{2}], 'o-');
end


% 2) means
xx=X:X+1;
hbar=lt_plot_bar(xx, Ymean, {'Errors', Ysem, 'Color' 'k'});


% --- save Yraw for stats
YrawStats2=Yraw;

% ====== PERFORM STATS (change across phase?)
[~, p]=ttest(YrawStats1{1}, YrawStats2{1}); % first targ (one targ phase vs 2 targ phase)
disp(['First targ (one phase vs. two phase), p='  num2str(p)]);
[~, p]=ttest(YrawStats1{2}, YrawStats2{2}); % second targ (one targ phase vs 2 targ phase)
disp(['Second targ (one phase vs. two phase), p='  num2str(p)]);



Xlabel={'First Targ', 'Second Targ'};
set(gca, 'XTick', xx);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);




%% % ====  [MP BIAS]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[SAME] MP BIAS')
ylabel('pitch shift, norm to target learning (PBS)');

plotcols=lt_make_plot_colors(ExptCounter-1, 0, 0);
Yexptnums={Y_ExptNum_Target, Y_ExptNum_NewTarg};

% =================== SINGLE DIR PHASE
X=1; % location to plot
ind=1; % pre-start BIDIR

% --
Yraw={}; % each cell is one category of syllable
Ymean=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_MP_Target;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);


Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_MP_NewTarg;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
% 1) raw vals
for i=1:length(Yraw{1});
    plot([X+0.2 X+1.2], [Yraw{1} Yraw{2}], 'o-');
end


% 2) means
xx=X:X+1;
hbar=lt_plot_bar(xx, Ymean, {'Errors', Ysem, 'Color' 'k'});

Xlabel={'First Targ', 'Second Targ'};
set(gca, 'XTick', xx);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);

% --- save Yraw for stats
YrawStats1=Yraw;


% =============================== DURING TWOP DIR PHASE
X=4;
ind=2; % pre-start BIDIR

% --
Yraw={}; % each cell is one category of syllable
Ymean=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_MP_Target;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);


Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_MP_NewTarg;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
% 1) raw vals
for i=1:length(Yraw{1});
    plot([X+0.2 X+1.2], [Yraw{1} Yraw{2}], 'o-');
end


% 2) means
xx=X:X+1;
hbar=lt_plot_bar(xx, Ymean, {'Errors', Ysem, 'Color' 'k'});

% --- save Yraw for stats
YrawStats2=Yraw;

% ====== PERFORM STATS (change across phase?)
[~, p]=ttest(YrawStats1{1}, YrawStats2{1}); % first targ (one targ phase vs 2 targ phase)
disp(['First targ (one phase vs. two phase), p='  num2str(p)]);
[~, p]=ttest(YrawStats1{2}, YrawStats2{2}); % second targ (one targ phase vs 2 targ phase)
disp(['Second targ (one phase vs. two phase), p='  num2str(p)]);



Xlabel={'First Targ', 'Second Targ'};
set(gca, 'XTick', xx);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);





%% +++++++++++ PLOTS - PRE AND POST BIDIR ON SAME PLOT [AFP]
%% ====== SUBPLOT 6
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('AFP BIAS [(SingleDir WN) - (Bidir WN)]')
% ylabel('pitch shift, norm to target learning (PBS)');

% -- target 1
DatArray=Y_AFP_Target;
X=1:2;

for i=1:size(DatArray,1);
    plot(X, DatArray(i,:),'o-');
end

Ymean=mean(DatArray,1);
Ysem=lt_sem(DatArray);

lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' ,'k'});
hold on;
% - sign rank and ttest
% p1=signrank(DatArray(:,2), DatArray(:,1));
[h, p2]=ttest(DatArray(:,2), DatArray(:,1));
% if p1<0.1 | p2<0.1;
%      lt_plot_text(X(1), max(Ymean)+0.1, ['p=' num2str(p1, '%3.2g')], 'c');
%      lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
% end
if p2<0.1;
     lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
end


% --  target 2
DatArray=Y_AFP_NewTarg;
X=4:5;

for i=1:size(DatArray,1);
    plot(X, DatArray(i,:),'o-');
end

Ymean=mean(DatArray,1);
Ysem=lt_sem(DatArray);

lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' ,'k'});
hold on;
% - sign rank and ttest
[h, p2]=ttest(DatArray(:,2), DatArray(:,1));
% if p1<0.1 | p2<0.1;
%      lt_plot_text(X(1), max(Ymean)+0.1, ['p=' num2str(p1, '%3.2g')], 'c');
%      lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
% end
if p2<0.1;
     lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
end



% --  other nontargets
DatArray=Y_AFP_AllOthers;
X=7:8;

for i=1:size(DatArray,1);
    plot(X, DatArray(i,:),'o-');
end

Ymean=mean(DatArray,1);
Ysem=lt_sem(DatArray);

lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' ,'k'});
hold on;
% - sign rank and ttest
[h, p2]=ttest(DatArray(:,2), DatArray(:,1));
% if p1<0.1 | p2<0.1;
%      lt_plot_text(X(1), max(Ymean)+0.1, ['p=' num2str(p1, '%3.2g')], 'c');
%      lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
% end
if p2<0.1;
     lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
end



% --- labels

Xlabel={'First Targ', 'Second Targ', 'Nontargets (others)'};
set(gca, 'XTick', [2 5 8]);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);


% % ============ STAT - compare effect of bidir on diff pairs of syl types
% % target1 vs target 2
% Y1=diff(Y_AFP_Target, 1, 2);
% Y2=diff(Y_AFP_NewTarg, 1, 2);
% 
% [h,p]=ttest(Y1, Y2);
% 
% 
% 
% 
Xlim=xlim;
Ylim=ylim;
% lt_plot_text(Xlim(1)+0.1, Ylim(1)+0.2, {['cyan: signrank'],['mag: ttest']}, 'r');





%% ====== LEARNING 
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('LEARNING [(SingleDir WN) - (Bidir WN)]')
% ylabel('pitch shift, norm to target learning (PBS)');

% -- target 1
DatArray=Y_Learning_Target;
X=1:2;

for i=1:size(DatArray,1);
    plot(X, DatArray(i,:),'o-');
end

Ymean=mean(DatArray,1);
Ysem=lt_sem(DatArray);

lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' ,'k'});
hold on;
% - sign rank and ttest
[h, p2]=ttest(DatArray(:,2), DatArray(:,1));
% if p1<0.1 | p2<0.1;
%      lt_plot_text(X(1), max(Ymean)+0.1, ['p=' num2str(p1, '%3.2g')], 'c');
%      lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
% end
if p2<0.1;
     lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
end


% --  target 2
DatArray=Y_Learning_NewTarg;
X=4:5;

for i=1:size(DatArray,1);
    plot(X, DatArray(i,:),'o-');
end

Ymean=mean(DatArray,1);
Ysem=lt_sem(DatArray);

lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' ,'k'});
hold on;
% - sign rank and ttest
[h, p2]=ttest(DatArray(:,2), DatArray(:,1));
% if p1<0.1 | p2<0.1;
%      lt_plot_text(X(1), max(Ymean)+0.1, ['p=' num2str(p1, '%3.2g')], 'c');
%      lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
% end
if p2<0.1;
     lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
end


% --  other nontargets
DatArray=Y_Learning_AllOthers;
X=7:8;

for i=1:size(DatArray,1);
    plot(X, DatArray(i,:),'o-');
end

Ymean=mean(DatArray,1);
Ysem=lt_sem(DatArray);

lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' ,'k'});
hold on;
% - sign rank and ttest
[h, p2]=ttest(DatArray(:,2), DatArray(:,1));
% if p1<0.1 | p2<0.1;
%      lt_plot_text(X(1), max(Ymean)+0.1, ['p=' num2str(p1, '%3.2g')], 'c');
%      lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
% end
if p2<0.1;
     lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
end


% --- labels

Xlabel={'First Targ', 'Second Targ', 'Nontargets (others)'};
set(gca, 'XTick', [2 5 8]);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);


% lt_plot_text(Xlim(1)+0.1, Ylim(1)+0.2, {['cyan: signrank'],['mag: ttest']}, 'r');

%% MP 
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('MP BIAS [(SingleDir WN) - (Bidir WN)]')
% ylabel('pitch shift, norm to target learning (PBS)');

% -- target 1
DatArray=Y_MP_Target;
X=1:2;

for i=1:size(DatArray,1);
    plot(X, DatArray(i,:),'o-');
end

Ymean=mean(DatArray,1);
Ysem=lt_sem(DatArray);

lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' ,'k'});
hold on;
% - sign rank and ttest
[h, p2]=ttest(DatArray(:,2), DatArray(:,1));
% if p1<0.1 | p2<0.1;
%      lt_plot_text(X(1), max(Ymean)+0.1, ['p=' num2str(p1, '%3.2g')], 'c');
%      lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
% end
if p2<0.1;
     lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
end


% --  target 2
DatArray=Y_MP_NewTarg;
X=4:5;

for i=1:size(DatArray,1);
    plot(X, DatArray(i,:),'o-');
end

Ymean=mean(DatArray,1);
Ysem=lt_sem(DatArray);

lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' ,'k'});
hold on;
% - sign rank and ttest
[h, p2]=ttest(DatArray(:,2), DatArray(:,1));
% if p1<0.1 | p2<0.1;
%      lt_plot_text(X(1), max(Ymean)+0.1, ['p=' num2str(p1, '%3.2g')], 'c');
%      lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
% end
if p2<0.1;
     lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
end


% --  other nontargets
DatArray=Y_MP_AllOthers;
X=7:8;

for i=1:size(DatArray,1);
    plot(X, DatArray(i,:),'o-');
end

Ymean=mean(DatArray,1);
Ysem=lt_sem(DatArray);

lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' ,'k'});
hold on;
% - sign rank and ttest
[h, p2]=ttest(DatArray(:,2), DatArray(:,1));
% if p1<0.1 | p2<0.1;
%      lt_plot_text(X(1), max(Ymean)+0.1, ['p=' num2str(p1, '%3.2g')], 'c');
%      lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
% end
if p2<0.1;
     lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
end


% --- labels

Xlabel={'First Targ', 'Second Targ', 'Nontargets (others)'};
set(gca, 'XTick', [2 5 8]);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);


% lt_plot_text(Xlim(1)+0.1, Ylim(1)+0.2, {['cyan: signrank'],['mag: ttest']}, 'r');



%% ++++++++++++++++++++++++++++ [SAME TYPE ONLY]
%% +++++++++++ PLOTS - PRE AND POST BIDIR ON SAME PLOT [AFP]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[SAME TYPE] AFP BIAS [(SingleDir WN) - (Bidir WN)]')
% ylabel('pitch shift, norm to target learning (PBS)');

% -- target 1
DatArray=Y_AFP_Target;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);

X=1:2;

for i=1:size(DatArray,1);
    plot(X, DatArray(i,:),'o-');
end

Ymean=mean(DatArray,1);
Ysem=lt_sem(DatArray);

lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' ,'k'});
hold on;
% - sign rank and ttest
% p1=signrank(DatArray(:,2), DatArray(:,1));
[h, p2]=ttest(DatArray(:,2), DatArray(:,1));
% if p1<0.1 | p2<0.1;
%      lt_plot_text(X(1), max(Ymean)+0.1, ['p=' num2str(p1, '%3.2g')], 'c');
%      lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
% end
if p2<0.1;
     lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
end


% --  target 2
DatArray=Y_AFP_NewTarg;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);
X=4:5;

for i=1:size(DatArray,1);
    plot(X, DatArray(i,:),'o-');
end

Ymean=mean(DatArray,1);
Ysem=lt_sem(DatArray);

lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' ,'k'});
hold on;
% - sign rank and ttest
[h, p2]=ttest(DatArray(:,2), DatArray(:,1));
% if p1<0.1 | p2<0.1;
%      lt_plot_text(X(1), max(Ymean)+0.1, ['p=' num2str(p1, '%3.2g')], 'c');
%      lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
% end
if p2<0.1;
     lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
end




% --- labels

Xlabel={'First Targ', 'Second Targ'};
set(gca, 'XTick', [2 5]);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);


% % ============ STAT - compare effect of bidir on diff pairs of syl types
% % target1 vs target 2
% Y1=diff(Y_AFP_Target, 1, 2);
% Y2=diff(Y_AFP_NewTarg, 1, 2);
% 
% [h,p]=ttest(Y1, Y2);
% 
% 
% 
% 
Xlim=xlim;
Ylim=ylim;
% lt_plot_text(Xlim(1)+0.1, Ylim(1)+0.2, {['cyan: signrank'],['mag: ttest']}, 'r');





%% ====== LEARNING 
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[SAME] LEARNING [(SingleDir WN) - (Bidir WN)]')
% ylabel('pitch shift, norm to target learning (PBS)');

% -- target 1
DatArray=Y_Learning_Target;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);

X=1:2;

for i=1:size(DatArray,1);
    plot(X, DatArray(i,:),'o-');
end

Ymean=mean(DatArray,1);
Ysem=lt_sem(DatArray);

lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' ,'k'});
hold on;
% - sign rank and ttest
[h, p2]=ttest(DatArray(:,2), DatArray(:,1));
% if p1<0.1 | p2<0.1;
%      lt_plot_text(X(1), max(Ymean)+0.1, ['p=' num2str(p1, '%3.2g')], 'c');
%      lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
% end
if p2<0.1;
     lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
end


% --  target 2
DatArray=Y_Learning_NewTarg;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);
X=4:5;

for i=1:size(DatArray,1);
    plot(X, DatArray(i,:),'o-');
end

Ymean=mean(DatArray,1);
Ysem=lt_sem(DatArray);

lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' ,'k'});
hold on;
% - sign rank and ttest
[h, p2]=ttest(DatArray(:,2), DatArray(:,1));
% if p1<0.1 | p2<0.1;
%      lt_plot_text(X(1), max(Ymean)+0.1, ['p=' num2str(p1, '%3.2g')], 'c');
%      lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
% end
if p2<0.1;
     lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
end



% --- labels
Xlabel={'First Targ', 'Second Targ'};
set(gca, 'XTick', [2 5]);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);


% lt_plot_text(Xlim(1)+0.1, Ylim(1)+0.2, {['cyan: signrank'],['mag: ttest']}, 'r');

%% MP 
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[SAME] MP BIAS [(SingleDir WN) - (Bidir WN)]')
% ylabel('pitch shift, norm to target learning (PBS)');

% -- target 1
DatArray=Y_MP_Target;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);
X=1:2;

for i=1:size(DatArray,1);
    plot(X, DatArray(i,:),'o-');
end

Ymean=mean(DatArray,1);
Ysem=lt_sem(DatArray);

lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' ,'k'});
hold on;
% - sign rank and ttest
[h, p2]=ttest(DatArray(:,2), DatArray(:,1));
% if p1<0.1 | p2<0.1;
%      lt_plot_text(X(1), max(Ymean)+0.1, ['p=' num2str(p1, '%3.2g')], 'c');
%      lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
% end
if p2<0.1;
     lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
end


% --  target 2
DatArray=Y_MP_NewTarg;
inds=Y_similar_NewTarg==1;
DatArray=DatArray(inds, :);
X=4:5;

for i=1:size(DatArray,1);
    plot(X, DatArray(i,:),'o-');
end

Ymean=mean(DatArray,1);
Ysem=lt_sem(DatArray);

lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' ,'k'});
hold on;
% - sign rank and ttest
[h, p2]=ttest(DatArray(:,2), DatArray(:,1));
% if p1<0.1 | p2<0.1;
%      lt_plot_text(X(1), max(Ymean)+0.1, ['p=' num2str(p1, '%3.2g')], 'c');
%      lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
% end
if p2<0.1;
     lt_plot_text(X(1), max(Ymean)+0.2, ['p=' num2str(p2, '%3.2g')], 'm');
end



% --- labels

Xlabel={'First Targ', 'Second Targ'};
set(gca, 'XTick', [2 5]);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);


% lt_plot_text(Xlim(1)+0.1, Ylim(1)+0.2, {['cyan: signrank'],['mag: ttest']}, 'r');




%% ======================= NEW PLOTS VERSION 1






%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% ===== NEW PLOTS VERSION 2 - one plot for single dir (show reversion), one plot for bidir (show reversion)
count=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

hsplots=[];

%% ==== single dir

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Reversion (single dir phase)')
% ylabel('pitch shift, norm to target learning (PBS)');
hsplots=[hsplots hsplot];

column=1;

% === target
x=1;
Ylearn=Y_Learning_Target(:,column);
Ymp=Y_MP_Target(:,column);

plot([x-0.25 x+0.25], [Ylearn Ymp]', 'ko-');
lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
% -- ttest for reversion
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end

% === new targ
x=2;
Ylearn=Y_Learning_NewTarg(:,column);
Ymp=Y_MP_NewTarg(:,column);

plot([x-0.25 x+0.25], [Ylearn Ymp]', 'ko-');
lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end


% === other nontargs
x=3;
Ylearn=Y_Learning_AllOthers(:,column);
Ymp=Y_MP_AllOthers(:,column);

plot([x-0.25 x+0.25], [Ylearn Ymp]', 'ko-');
lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end


set(gca, 'Xtick', 1:x);
Xlabels={'First targ', 'Second targ', 'All other nontargs'};

set(gca, 'XtickLabel', Xlabels);

rotateXLabels(gca, 45);


%% ==== two dir

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Reversion (two dir phase)')
% ylabel('pitch shift, norm to target learning (PBS)');

column=2;
hsplots=[hsplots hsplot];

% === target
x=1;
Ylearn=Y_Learning_Target(:,column);
Ymp=Y_MP_Target(:,column);

plot([x-0.25 x+0.25], [Ylearn Ymp]', 'ko-');
lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
% -- ttest for reversion
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end


% === new targ
x=2;
Ylearn=Y_Learning_NewTarg(:,column);
Ymp=Y_MP_NewTarg(:,column);

plot([x-0.25 x+0.25], [Ylearn Ymp]', 'ko-');
lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end


% === other nontargs
x=3;
Ylearn=Y_Learning_AllOthers(:,column);
Ymp=Y_MP_AllOthers(:,column);

plot([x-0.25 x+0.25], [Ylearn Ymp]', 'ko-');
lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end


% ===
set(gca, 'Xtick', 1:x);
Xlabels={'First targ', 'Second targ', 'All other nontargs'};

set(gca, 'XtickLabel', Xlabels);

rotateXLabels(gca, 45);

linkaxes(hsplots, 'y');



%% ==== single dir [just means]

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Reversion (single dir phase)')
% ylabel('pitch shift, norm to target learning (PBS)');
hsplots=[hsplots hsplot];

column=1;

% === target
x=1;
Ylearn=Y_Learning_Target(:,column);
Ymp=Y_MP_Target(:,column);

lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
% -- ttest for reversion
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end

% === new targ
x=2;
Ylearn=Y_Learning_NewTarg(:,column);
Ymp=Y_MP_NewTarg(:,column);

lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end


% === other nontargs
x=3;
Ylearn=Y_Learning_AllOthers(:,column);
Ymp=Y_MP_AllOthers(:,column);

lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end


set(gca, 'Xtick', 1:x);
Xlabels={'First targ', 'Second targ', 'All other nontargs'};

set(gca, 'XtickLabel', Xlabels);

rotateXLabels(gca, 45);


%% ==== two dir [just means]

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Reversion (two dir phase)')
% ylabel('pitch shift, norm to target learning (PBS)');

column=2;
hsplots=[hsplots hsplot];

% === target
x=1;
Ylearn=Y_Learning_Target(:,column);
Ymp=Y_MP_Target(:,column);

lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
% -- ttest for reversion
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end


% === new targ
x=2;
Ylearn=Y_Learning_NewTarg(:,column);
Ymp=Y_MP_NewTarg(:,column);

lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end


% === other nontargs
x=3;
Ylearn=Y_Learning_AllOthers(:,column);
Ymp=Y_MP_AllOthers(:,column);

lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end


% ===
set(gca, 'Xtick', 1:x);
Xlabels={'First targ', 'Second targ', 'All other nontargs'};

set(gca, 'XtickLabel', Xlabels);

rotateXLabels(gca, 45);

linkaxes(hsplots, 'y');



%% +++++ plot color to represent syl type

%% ==== single dir

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Reversion (single dir phase)')
% ylabel('pitch shift, norm to target learning (PBS)');
hsplots=[hsplots hsplot];

column=1;

% === target
x=1;
Ylearn=Y_Learning_Target(:,column);
Ymp=Y_MP_Target(:,column);

% SS
inds=Y_similar_NewTarg & Y_presimilar_NewTarg;
color='b';
plot([x-0.25 x+0.25], [Ylearn(inds) Ymp(inds)]', 'o-', 'Color',color);
% SD
inds=Y_similar_NewTarg & ~Y_presimilar_NewTarg;
color='c';
plot([x-0.25 x+0.25], [Ylearn(inds) Ymp(inds)]', 'o-', 'Color',color);
% DS
inds=~Y_similar_NewTarg & Y_presimilar_NewTarg;
color='r';
if any(inds);
plot([x-0.25 x+0.25], [Ylearn(inds) Ymp(inds)]', 'o-', 'Color',color);
end
% DD
inds=~Y_similar_NewTarg & ~Y_presimilar_NewTarg;
color='m';
if any(inds);
plot([x-0.25 x+0.25], [Ylearn(inds) Ymp(inds)]', 'o-', 'Color',color);
end

lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
% -- ttest for reversion
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end

% === new targ
x=2;
Ylearn=Y_Learning_NewTarg(:,column);
Ymp=Y_MP_NewTarg(:,column);

% SS
inds=Y_similar_NewTarg & Y_presimilar_NewTarg;
color='b';
plot([x-0.25 x+0.25], [Ylearn(inds) Ymp(inds)]', 'o-', 'Color',color);
% SD
inds=Y_similar_NewTarg & ~Y_presimilar_NewTarg;
color='c';
plot([x-0.25 x+0.25], [Ylearn(inds) Ymp(inds)]', 'o-', 'Color',color);
% DS
inds=~Y_similar_NewTarg & Y_presimilar_NewTarg;
color='r';
if any(inds);
plot([x-0.25 x+0.25], [Ylearn(inds) Ymp(inds)]', 'o-', 'Color',color);
end
% DD
inds=~Y_similar_NewTarg & ~Y_presimilar_NewTarg;
color='m';
if any(inds);
plot([x-0.25 x+0.25], [Ylearn(inds) Ymp(inds)]', 'o-', 'Color',color);
end



lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end


% === other nontargs
x=3;
Ylearn=Y_Learning_AllOthers(:,column);
Ymp=Y_MP_AllOthers(:,column);

plot([x-0.25 x+0.25], [Ylearn Ymp]', 'ko-');
lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end


set(gca, 'Xtick', 1:x);
Xlabels={'First targ', 'Second targ', 'All other nontargs'};

set(gca, 'XtickLabel', Xlabels);

rotateXLabels(gca, 45);


%% ==== two dir

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Reversion (two dir phase)')
% ylabel('pitch shift, norm to target learning (PBS)');

column=2;
hsplots=[hsplots hsplot];

% === target
x=1;
Ylearn=Y_Learning_Target(:,column);
Ymp=Y_MP_Target(:,column);

% SS
inds=Y_similar_NewTarg & Y_presimilar_NewTarg;
color='b';
plot([x-0.25 x+0.25], [Ylearn(inds) Ymp(inds)]', 'o-', 'Color',color);
% SD
inds=Y_similar_NewTarg & ~Y_presimilar_NewTarg;
color='c';
plot([x-0.25 x+0.25], [Ylearn(inds) Ymp(inds)]', 'o-', 'Color',color);
% DS
inds=~Y_similar_NewTarg & Y_presimilar_NewTarg;
color='r';
if any(inds);
plot([x-0.25 x+0.25], [Ylearn(inds) Ymp(inds)]', 'o-', 'Color',color);
end
% DD
inds=~Y_similar_NewTarg & ~Y_presimilar_NewTarg;
color='m';
if any(inds);
plot([x-0.25 x+0.25], [Ylearn(inds) Ymp(inds)]', 'o-', 'Color',color);
end
lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
% -- ttest for reversion
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end


% === new targ
x=2;
Ylearn=Y_Learning_NewTarg(:,column);
Ymp=Y_MP_NewTarg(:,column);

% SS
inds=Y_similar_NewTarg & Y_presimilar_NewTarg;
color='b';
plot([x-0.25 x+0.25], [Ylearn(inds) Ymp(inds)]', 'o-', 'Color',color);
% SD
inds=Y_similar_NewTarg & ~Y_presimilar_NewTarg;
color='c';
plot([x-0.25 x+0.25], [Ylearn(inds) Ymp(inds)]', 'o-', 'Color',color);
% DS
inds=~Y_similar_NewTarg & Y_presimilar_NewTarg;
color='r';
if any(inds);
plot([x-0.25 x+0.25], [Ylearn(inds) Ymp(inds)]', 'o-', 'Color',color);
end
% DD
inds=~Y_similar_NewTarg & ~Y_presimilar_NewTarg;
color='m';
if any(inds);
plot([x-0.25 x+0.25], [Ylearn(inds) Ymp(inds)]', 'o-', 'Color',color);
end

% 
% plot([x-0.25 x+0.25], [Ylearn Ymp]', 'ko-');
lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end


% === other nontargs
x=3;
Ylearn=Y_Learning_AllOthers(:,column);
Ymp=Y_MP_AllOthers(:,column);

plot([x-0.25 x+0.25], [Ylearn Ymp]', 'ko-');
lt_plot_bar([x-0.2], mean([Ylearn]), {'Errors', [lt_sem(Ylearn)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean([Ymp]), {'Errors', [lt_sem(Ymp)], 'Color','r', 'BarWidth',0.3});
[h, p]=ttest(Ylearn, Ymp);
if p<0.1;
    lt_plot_text(x, max([Ylearn' Ymp']), num2str(p, '%3.2g'), 'b');
end


% ===
set(gca, 'Xtick', 1:x);
Xlabels={'First targ', 'Second targ', 'All other nontargs'};

set(gca, 'XtickLabel', Xlabels);

rotateXLabels(gca, 45);

linkaxes(hsplots, 'y');


%% +++++ ONLY SAME TYPE
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


%% ==== single dir

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[SAME] Reversion (single dir phase)')
% ylabel('pitch shift, norm to target learning (PBS)');
hsplots=[hsplots hsplot];

column=1;

% === target
x=1;
Ylearn=Y_Learning_Target(:,column);
Ymp=Y_MP_Target(:,column);

% ====
inds=Y_similar_NewTarg==1;
color='k';

Y1=Ylearn(inds);
Y2=Ymp(inds);
plot([x-0.25 x+0.25], [Y1 Y2]', '-', 'Color',color);

lt_plot_bar([x-0.2], mean(Y1), {'Errors', [lt_sem(Y1)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean(Y2), {'Errors', [lt_sem(Y2)], 'Color','r', 'BarWidth',0.3});
% -- ttest for reversion
[h, p]=ttest(Y1, Y2);
if p<0.1;
    lt_plot_text(x, max([Y1' Y2']), num2str(p, '%3.2g'), 'b');
end


% === new targ
x=2;
Ylearn=Y_Learning_NewTarg(:,column);
Ymp=Y_MP_NewTarg(:,column);

% ====
inds=Y_similar_NewTarg==1;
color='k';

Y1=Ylearn(inds);
Y2=Ymp(inds);
plot([x-0.25 x+0.25], [Y1 Y2]', '-', 'Color',color);

lt_plot_bar([x-0.2], mean(Y1), {'Errors', [lt_sem(Y1)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean(Y2), {'Errors', [lt_sem(Y2)], 'Color','r', 'BarWidth',0.3});
% -- ttest for reversion
[h, p]=ttest(Y1, Y2);
if p<0.1;
    lt_plot_text(x, max([Y1' Y2']), num2str(p, '%3.2g'), 'b');
end


set(gca, 'Xtick', 1:x);
Xlabels={'First targ', 'Second targ'};

set(gca, 'XtickLabel', Xlabels);

rotateXLabels(gca, 45);



%% ==== two dir

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[SAME] Reversion (two dir phase)')
% ylabel('pitch shift, norm to target learning (PBS)');
hsplots=[hsplots hsplot];

column=2;

% === target
x=1;
Ylearn=Y_Learning_Target(:,column);
Ymp=Y_MP_Target(:,column);

% ====
inds=Y_similar_NewTarg==1;
color='k';

Y1=Ylearn(inds);
Y2=Ymp(inds);
plot([x-0.25 x+0.25], [Y1 Y2]', '-', 'Color',color);

lt_plot_bar([x-0.2], mean(Y1), {'Errors', [lt_sem(Y1)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean(Y2), {'Errors', [lt_sem(Y2)], 'Color','r', 'BarWidth',0.3});
% -- ttest for reversion
[h, p]=ttest(Y1, Y2);
if p<0.1;
    lt_plot_text(x, max([Y1' Y2']), num2str(p, '%3.2g'), 'b');
end


% === new targ
x=2;
Ylearn=Y_Learning_NewTarg(:,column);
Ymp=Y_MP_NewTarg(:,column);

% ====
inds=Y_similar_NewTarg==1;
color='k';

Y1=Ylearn(inds);
Y2=Ymp(inds);
plot([x-0.25 x+0.25], [Y1 Y2]', '-', 'Color',color);

lt_plot_bar([x-0.2], mean(Y1), {'Errors', [lt_sem(Y1)], 'Color', 'k', 'BarWidth',0.3});
lt_plot_bar([x+0.2], mean(Y2), {'Errors', [lt_sem(Y2)], 'Color','r', 'BarWidth',0.3});
% -- ttest for reversion
[h, p]=ttest(Y1, Y2);
if p<0.1;
    lt_plot_text(x, max([Y1' Y2']), num2str(p, '%3.2g'), 'b');
end


set(gca, 'Xtick', 1:x);
Xlabels={'First targ', 'Second targ'};

set(gca, 'XtickLabel', Xlabels);

rotateXLabels(gca, 45);



%% ==== plot bars, but separate by adjacency




