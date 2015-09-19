function [PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANbidir(SeqDepPitch_AcrossBirds, PARAMS)
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
            
            EpochNameList={'days_consolid_late', 'days_bidir_late'};
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

lt_subtitle('Learning(black), MP(red), AFP(blue), "consolidation end" vs "bidir end"')


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
        
        % -- to collect vals across syls
        Y_Learning_AllOthers=[];
        Y_MP_AllOthers=[];
        Y_AFP_AllOthers=[];
        
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            EpochNameList={'days_consolid_late', 'days_bidir_late'};
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



% ===== PLOT ALL IN ONE PLOT
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
        
        plotcols=lt_make_plot_colors(TotalNumExperiments, 0, 0);
cc=1;

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
        
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            EpochNameList={'days_consolid_late', 'days_bidir_late'};
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
                plot(X, Y_Learning, 'o-', 'Color',plotcol)
                plot(X+1+length(X), Y_MP, 'o-', 'Color',plotcol)
                plot(X+2+2*length(X), Y_AFP, 'o-', 'Color',plotcol)
                
                % ----- COLLECT TO PLOT MEAN
                Y_Learning_Target=[Y_Learning_Target; Y_Learning];
                Y_MP_Target=[Y_MP_Target; Y_MP];
                Y_AFP_Target=[Y_AFP_Target; Y_AFP];
                
            elseif istarg_other==1
                plot(X+9, Y_Learning, 'o-', 'Color',plotcol)
                plot(X+10+length(X), Y_MP, 'o-', 'Color',plotcol)
                plot(X+11+2*length(X), Y_AFP, 'o-', 'Color',plotcol)
                
                % ----- COLLECT TO PLOT MEAN
                Y_Learning_NewTarg=[Y_Learning_NewTarg; Y_Learning];
                Y_MP_NewTarg=[Y_MP_NewTarg; Y_MP];
                Y_AFP_NewTarg=[Y_AFP_NewTarg; Y_AFP];
                
            else % is all other syls
                plot(X+19, Y_Learning, 'o-', 'Color',plotcol)
                plot(X+20+length(X), Y_MP, 'o-', 'Color',plotcol)
                plot(X+21+2*length(X), Y_AFP, 'o-', 'Color',plotcol)
                
                % ----- COLLECT TO PLOT MEAN
                Y_Learning_AllOthers=[Y_Learning_AllOthers; Y_Learning];
                Y_MP_AllOthers=[Y_MP_AllOthers; Y_MP];
                Y_AFP_AllOthers=[Y_AFP_AllOthers; Y_AFP];
                
            end
        end
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

