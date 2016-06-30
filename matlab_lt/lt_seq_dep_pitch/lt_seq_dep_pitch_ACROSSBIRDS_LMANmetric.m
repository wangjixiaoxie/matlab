function [PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANlearning(SeqDepPitch_AcrossBirds, PARAMS, norm_by_targsyl, epochfield_input)
%% LMAN bias localized to targ?  clean plots (final?) LT 11/21/15

%% PARAMS

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

if ~exist('norm_by_targsyl','var');
    norm_by_targsyl=1;
end


%% REMOVE SYLS THAT SHOULD NOT BE ANALYZED (I.E O/L WITH WN, for non-catch analyses)

disp('--');
disp('Removing syllables that should not be analyzed - i.e. WN overlap, since not using catch songs. REMOVED:');

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



%% [ COLLECT DATA]

epochfield=epochfield_input;

% LearningPBS_all=[];
%
% MPbias_all=[];
% AFPbias_all=[];
% SimDiff_all=[];
% TargStatus_all=[];
% PreSylSimilar_all=[];
% Expt_count_all=[];
% Yexpt_all={};
% Ybird_all={};
% Y_PosRelTarg_All=[];
%
% Generalization_MP_all=[];
% Generalization_AFP_all=[];
% Generalization_Learn_all=[];
%
% cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS=[];
% cvRatio_pvalue_UsingAllVals_ALLEXPTS=[];
%
% CVofAllDays_UsingValsDividedByDayMean_PBS_ALLEXPTS=[];
% CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLEXPTS=[];
%
% cvPBS_alldays_ALLEXPTS={};
% cvMUSC_alldays_ALLEXPTS={};
%
%
% expt_count=1;


Index1_AllSyls=[];
Index1_MP_AllSyls=[];
Index1_quotient_AllSyls=[];
Index3_AllSyls=[];
Index4_AllSyls=[];

SameType_AllSyls=[];


Index1_MeanOverSameType_AllExpt=[];
Index1_Quotient_MeanOverSameType_AllExpt=[];
Index4_MeanOverSameType_AllExpt=[];


for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        % ========== COLLECT DATA
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, epochfield);
            
            disp(['COLLLECTED DATA FOR : ' birdname '-' exptname]);
            
            
            % ===== STATS AT TARGET
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            FF_PBS_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl).meanFF_pbs;
            FF_MUSC_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl).meanFF_musc;
            FF_AFP_targ=FF_PBS_targ-FF_MUSC_targ;
            
            Index1_SameTypes=[];
            Index1_SameTypes_Quotient=[];
            Index4_SameTypes=[];
            
            for j=1:length(SylsUnique);
                syl=SylsUnique{j};
                
                istarget=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
                if istarget==1;
                    continue
                end
                
                
                % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                % MUSC)
                FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                FF_AFP=FF_PBS-FF_MUSC;
                
                
                % ==== FOR THIS SYL, CALCULATE VARIOUS INDICES OF
                % SPECIFICTIY RELATIVE TO TARGET
                % --- 1) calculate AFP bias at target relative to learning.
                % If add learning and AFP bias at the nontarget, how much
                % do you change the ratio? (+ is more target specific,
                % - is more nontarget specific)
                
                AFPdivideLearn_Targ=FF_AFP_targ/FF_PBS_targ;
                
                AFPdivideLearn_TargAndNontarg=(FF_AFP_targ + FF_AFP)/(FF_PBS_targ + FF_PBS);
                
                Index1=AFPdivideLearn_Targ-AFPdivideLearn_TargAndNontarg;
                Index1_quotient=AFPdivideLearn_Targ/AFPdivideLearn_TargAndNontarg;
                
                % ---- 2) SAME AS ABOVE, BUT USING MP BIAS
                MPdivideLearn_Targ=FF_MUSC_targ/FF_PBS_targ;
                MPdivideLearn_TargAndNontarg=(FF_MUSC_targ + FF_MUSC)/(FF_PBS_targ + FF_PBS);
                
                Index1_MP=MPdivideLearn_Targ-MPdivideLearn_TargAndNontarg;
                
                % ---- 3) Learn/learn vs. (learn + AFP)/(learn + AFP)
                % raeson, this might better (MORE HIGHLY) weigh nontargets that shifted
                % strongly
                Learn_Div_Learn=FF_PBS/FF_PBS_targ;
                tmp=(FF_PBS+FF_AFP)/(FF_PBS_targ + FF_AFP_targ);
                Index3=Learn_Div_Learn-tmp;
                
                % ---- 4) like above, but have sum of learn as denominator
                % (that way reduces error due to low learning)
                Learn_Div_Learn=FF_PBS/(FF_PBS_targ+FF_PBS);
                tmp=(FF_PBS+FF_AFP)/(FF_PBS_targ + FF_AFP_targ + FF_PBS + FF_AFP);
                Index4=Learn_Div_Learn-tmp;

                
                % ====== OUTPUT ALL INDICES
                Index1_AllSyls=[Index1_AllSyls Index1];
                Index1_MP_AllSyls=[Index1_MP_AllSyls Index1_MP];
                Index1_quotient_AllSyls =[ Index1_quotient_AllSyls Index1_quotient];
                Index3_AllSyls=[Index3_AllSyls Index3];
                Index4_AllSyls=[Index4_AllSyls Index4];
                
                % ====== GET STATS ABOUT SYLS
                similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
                SameType_AllSyls=[SameType_AllSyls similar];
                
                
                % === COLLECT FOR MEAN IN THIS EXPERIMENT
                if similar
                    Index1_SameTypes=[Index1_SameTypes Index1];
                    Index1_SameTypes_Quotient=[Index1_SameTypes_Quotient Index1_quotient];
                    Index4_SameTypes=[Index4_SameTypes Index4];
                end
                
                
                % ==== display output for user:
                if similar==1
                    disp(' ');
                    disp([birdname '-' exptname '-' syl]);
                    disp(['target: ' num2str(FF_AFP_targ) '/' num2str(FF_PBS_targ)]);
                    disp(['nontarg: ' num2str(FF_AFP) '/' num2str(FF_PBS)]);
                    disp(['Index1: ' num2str(Index1)]);
                    disp(['Index3: ' num2str(Index3)]);
                    disp(['Index4: ' num2str(Index4)]);
                    
                end
            end
            
            % === GET MEAN FOR SAME-TYPE SYLS FOR THIS EXPERIMENT
            if ~isempty(Index1_SameTypes)
            Index1_MeanOverSameType_AllExpt=[Index1_MeanOverSameType_AllExpt mean(Index1_SameTypes)];
            Index1_Quotient_MeanOverSameType_AllExpt=[Index1_Quotient_MeanOverSameType_AllExpt mean(Index1_SameTypes_Quotient)];
            
            Index4_MeanOverSameType_AllExpt=[Index4_MeanOverSameType_AllExpt mean(Index4_SameTypes)];
            end
            
        end
    end
end
    
    
    %% ================= PLOT HISTOGRAMS OF INDICES
    
%     lt_figure; hold on;
%     inds=SameType_AllSyls==1;
%     title('ind 1, difference');
%     lt_plot_histogram(Index1_AllSyls(inds), '', 1, 0)
% %     lt_plot_cdf(Index1_AllSyls(inds))
    
%     lt_figure; hold on;
%     inds=SameType_AllSyls==1;
%     title('ind 1, difference');
% %     lt_plot_histogram(Index1_AllSyls(inds), '', 1, 0)
%     lt_plot_cdf(Index1_AllSyls(inds))
%     line([0 0], ylim)

    lt_figure; hold on;
    title('index 1, same type');
    inds=SameType_AllSyls==1;
    plot(1, Index1_AllSyls(inds), 'o');
    p=signrank(Index1_AllSyls(inds));
    if p<0.15;
        lt_plot_annotation(1, ['p=' num2str(p)], 'b');
    end

    lt_figure; hold on;
    title('index 3, same type');
    inds=SameType_AllSyls==1;
    plot(1, Index3_AllSyls(inds), 'ok');
    p=signrank(Index3_AllSyls(inds));
        if p<0.15;
        lt_plot_annotation(1, ['p=' num2str(p)], 'b');
    end


    lt_figure; hold on;
    title('index 4, same type');
    inds=SameType_AllSyls==1;
    plot(1, Index4_AllSyls(inds), 'ok');
    p=signrank(Index4_AllSyls(inds));
        if p<0.15;
        lt_plot_annotation(1, ['p=' num2str(p)], 'b');
    end

    
%     % ----
%     lt_figure; hold on;
%     inds=SameType_AllSyls==1;
%     
%     plot(1, Index1_MP_AllSyls(inds), 'o');
%     lt_plot_bar(1, mean(Index1_MP_AllSyls))
    
        
    % ----
%     lt_figure; hold on;
%     inds=SameType_AllSyls==1;
%     title('quotient, all same types');
%     lt_plot_cdf(Index1_quotient_AllSyls(inds), 'k');
%     p=signrank(Index1_quotient_AllSyls(inds));
%     if p<0.1
%        lt_plot_annotation(1, ['signrank p=' num2str(p)], 'k');
%     end
    
    
    % ----
    lt_figure; hold on;
    inds=SameType_AllSyls==1;
    title('ind 1, quotient');
    plot(1, Index1_quotient_AllSyls(inds), 'o');
%     lt_plot_bar(1, mean(Index1_quotient_AllSyls))
    
%     % ----
%     lt_figure; hold on;
%     inds=SameType_AllSyls==1;
%     title('quotient');
%         lt_plot_histogram(Index1_quotient_AllSyls(inds), '', 1);
% %     lt_plot_bar(1, mean(Index1_quotient_AllSyls))
%     line([1 1], ylim);
%     p=signrank(Index1_quotient_AllSyls(inds))
    
    % ----
%     lt_figure; hold on;
%     inds=SameType_AllSyls==1;
%     xlabel('quotient');
%     ylabel('difference');
%   
%     plot(Index1_quotient_AllSyls(inds), Index1_AllSyls(inds), 'o')
%     line(xlim, [0 0]);
%     line([1 1], ylim);
    
    
% ======================== FIGURES, ONE DOT FOR EACH TARGET
lt_figure; hold on;
title('one dot per expt; index 1, diff');
    plot(1, Index1_MeanOverSameType_AllExpt, 'ok');
    line(xlim, [0 0]);
    
    p=signrank(Index1_MeanOverSameType_AllExpt);
%     lt_plot_bar(1, mean(Index1_MP_AllSyls))
        if p<0.15;
        lt_plot_annotation(1, ['p=' num2str(p)], 'b');
    end


lt_figure; hold on;
title('one dot per expt; index 1, quotient');
    plot(1, Index1_Quotient_MeanOverSameType_AllExpt, 'ok');
    line(xlim, [1 1]);
    
    p=signrank(Index1_Quotient_MeanOverSameType_AllExpt);
%     lt_plot_bar(1, mean(Index1_MP_AllSyls))
    if p<0.15;
        lt_plot_annotation(1, ['p=' num2str(p)], 'b');
    end

lt_figure; hold on;
title('one dot per expt; index 4');
    plot(1, Index4_MeanOverSameType_AllExpt, 'ok');
    line(xlim, [0 0]);
    
    p=signrank(Index4_MeanOverSameType_AllExpt);
    
%     lt_plot_bar(1, mean(Index1_MP_AllSyls))

    if p<0.15;
        lt_plot_annotation(1, ['p=' num2str(p)], 'b');
    end

    
    
    %                 % --- other stuff
    %                 Y_PosRelTarg=[Y_PosRelTarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).distance_from_targ];
    %                 Y_FFmean_pbs=[Y_FFmean_pbs FF_PBS];
    %                 Y_FFmean_musc=[Y_FFmean_musc FF_MUSC];
    %                 Y_FFsem_pbs=[Y_FFsem_pbs FFsem_PBS];
    %                 Y_FFsem_musc=[Y_FFsem_musc FFsem_MUSC];
    %                 Y_AFP_bias=[Y_AFP_bias FF_PBS-FF_MUSC];
    %
    %                 Y_syls=[Y_syls, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
    %                 Y_similar_diff=[Y_similar_diff SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
    %                 Y_istarg=[Y_istarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
    %                 Y_presimilar=[Y_presimilar SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl];
    % %                 Y_AcousticDist=[Y_AcousticDist SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];
    % %                 targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
    % %                 Y_Corr=[Y_Corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl)];
    %
    %                 Yexpt=[Yexpt exptname(end-3:end)];
    %                 Ybird=[Ybird birdname(1:4)];
    %
    %                 Expt_count_all=[Expt_count_all expt_count];
    %
    %
    %                 Y_Generalization_MP=[Y_Generalization_MP Generalization_MP];
    %                 Y_Generalization_AFP=[Y_Generalization_AFP Generalization_AFP];
    %                 Y_Generalization_Learn=[Y_Generalization_Learn Generalization_Learn];
    %
    % %
    %             end
    %
    %
    %             % ================= Flip sign if learning at targsyl is negative
    %             if Y_FFmean_pbs(Y_istarg==1)<0; % negative learning
    %                 Y_FFmean_pbs=-1.*Y_FFmean_pbs;
    %                 Y_FFmean_musc=-1.*Y_FFmean_musc;
    %                 Y_AFP_bias=-1.*Y_AFP_bias;
    %             end
    %
    %             % ========= Normalize by targsyl if desired (PBS learning
    %             % by taergsyl)
    %             if norm_by_targsyl==1;
    %                 learning_by_targ=Y_FFmean_pbs(Y_istarg==1);
    %
    %                 Y_FFmean_pbs=Y_FFmean_pbs./learning_by_targ;
    %                 Y_FFmean_musc=Y_FFmean_musc./learning_by_targ;
    %                 Y_AFP_bias=Y_AFP_bias./learning_by_targ;
    %             end
    %
    %
    %
    %             % ============================ COLLECT DATA TO PLOT FOR ALL
    %             % EXPERIMENTS
    %             if any(~isnan(Y_FFmean_pbs)); % if any are not nan.
    %
    %                 % -- cv stuff
    %                 cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS=[cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS cvRatio_MUSCoverPBS_usingAllVals_ALLSYLS];
    %                 cvRatio_pvalue_UsingAllVals_ALLEXPTS=[cvRatio_pvalue_UsingAllVals_ALLEXPTS cvRatio_pvalue_UsingAllVals_ALLSYLS];
    %
    %                 CVofAllDays_UsingValsDividedByDayMean_PBS_ALLEXPTS=[CVofAllDays_UsingValsDividedByDayMean_PBS_ALLEXPTS CVofAllDays_UsingValsDividedByDayMean_PBS_ALLSYLS];
    %                 CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLEXPTS=[CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLEXPTS CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLSYLS];
    %
    %                 cvPBS_alldays_ALLEXPTS=[cvPBS_alldays_ALLEXPTS cvPBS_alldays_ALLSYLS];
    %                 cvMUSC_alldays_ALLEXPTS=[cvMUSC_alldays_ALLEXPTS  cvMUSC_alldays_ALLSYLS];
    %
    %
    %                 % -- other stuff
    %                 LearningPBS_all=[LearningPBS_all Y_FFmean_pbs];
    %                 MPbias_all=[MPbias_all Y_FFmean_musc];
    %                 AFPbias_all=[AFPbias_all Y_AFP_bias];
    %                 SimDiff_all=[SimDiff_all Y_similar_diff];
    %                 TargStatus_all=[TargStatus_all Y_istarg];
    %                 PreSylSimilar_all=[PreSylSimilar_all Y_presimilar];
    %
    %                 Generalization_AFP_all=[Generalization_AFP_all Y_Generalization_AFP];
    %                  Generalization_MP_all=[Generalization_MP_all Y_Generalization_MP];
    %                  Generalization_Learn_all=[Generalization_Learn_all Y_Generalization_Learn];
    %
    %
    %
    %                 Yexpt_all=[Yexpt_all Yexpt];
    %                 Ybird_all=[Ybird_all Ybird];
    %                 Y_PosRelTarg_All=[Y_PosRelTarg_All Y_PosRelTarg];
    %
    %
    %
    %                 expt_count=expt_count+1;
    %
    %             end
    %         end
    %     end
    % end
    
    
