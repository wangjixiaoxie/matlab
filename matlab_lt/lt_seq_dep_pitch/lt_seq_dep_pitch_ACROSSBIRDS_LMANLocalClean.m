function [PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANlearning(SeqDepPitch_AcrossBirds, PARAMS, norm_by_targsyl, epochfield_input, UseBaselineForCV)
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


%% === list names of all types of syls
if (0)
disp(' --- ');

for i=1:length(SeqDepPitch_AcrossBirds.birds);
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                continue
            end
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ==1 ...
                    & SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl==0;
                
                
                disp([birdname '-' exptname '-' syl]);
            end
        end
        
        
        
    end
end
end

%% ==== list birds and expts
disp(' =============== birds/expts');
for i=1:length(SeqDepPitch_AcrossBirds.birds);
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;

        disp([birdname '-' exptname]);

    end
end


%% [ COLLECT DATA] - PLOT AFP BIAS VS. MP LEARNING, DISTRIBUTIONS, AND OTHER THINGS
% lt_figure; 
% hold on;
% title('
% line([0 0.05], [0 0.05]);

epochfield=epochfield_input;

LearningPBS_all=[];

MPbias_all=[];
AFPbias_all=[];
SimDiff_all=[];
TargStatus_all=[];
PreSylSimilar_all=[];
Expt_count_all=[];
Yexpt_all={};
Ybird_all={};
Y_PosRelTarg_All=[];

Generalization_MP_all=[];
Generalization_AFP_all=[];
Generalization_Learn_all=[];

cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS=[];
cvRatio_pvalue_UsingAllVals_ALLEXPTS=[];

CVofAllDays_UsingValsDividedByDayMean_PBS_ALLEXPTS=[];
CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLEXPTS=[];

cvPBS_alldays_ALLEXPTS={};
cvMUSC_alldays_ALLEXPTS={};

TargLearnDirAll=[];


expt_count=1;

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        
        % ==== ONE FIGURE PER EXPT
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        % ========== COLLECT DATA
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, epochfield);
            
        disp(['COLLLECTED DATA FOR : ' birdname '-' exptname]);
            
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
            Yexpt={};
            Ybird={};
            Y_PosRelTarg=[];
            
            Y_Generalization_MP=[];
            Y_Generalization_AFP=[];
             Y_Generalization_Learn=[];
           
            
            % -- for CV stuff
            cvRatio_MUSCoverPBS_usingAllVals_ALLSYLS=[];
            cvRatio_pvalue_UsingAllVals_ALLSYLS=[];
            
            CVofAllDays_UsingValsDividedByDayMean_PBS_ALLSYLS=[];
            CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLSYLS=[];
            
            cvPBS_alldays_ALLSYLS={};
            cvMUSC_alldays_ALLSYLS={};
            
            % ===== STATS AT TARGET
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            FF_PBS_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl).meanFF_pbs;
            FF_MUSC_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl).meanFF_musc;
            FF_AFP_targ=FF_PBS_targ-FF_MUSC_targ;
            
            targlearndir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
            
            for j=1:length(SylsUnique);
                syl=SylsUnique{j};
                
                % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                % MUSC)
                FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                
                FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                
                
                % ====== CALCULATE AFP AND MUSC RELATIVE TO
                % THOSE VALUES AT TARGET (I.E. GENERALIZATIONS)
                FF_AFP=FF_PBS-FF_MUSC;
                
                Generalization_MP=FF_MUSC/FF_MUSC_targ;
                Generalization_AFP=FF_AFP/FF_AFP_targ;
                Generalization_Learn=FF_PBS/FF_PBS_targ;
                
                
                
                % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                % ======================== calculate CV PBS and MUSC (for each day of
                if UseBaselineForCV==1;
                TvalsPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).Tvals_WithinTimeWindow;
                TvalsMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).Tvals_WithinTimeWindow;
                    
                FFvalsPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_WithinTimeWindow;
                FFvalsMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).rawFF_WithinTimeWindow;
                    
                else
                % inactivation, and mean across days)
                TvalsPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).TvalsWithinWindow_PBS;
                TvalsMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).TvalsWithinWindow_MUSC;
                
                FFvalsPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).FFminusBaseWithinWindow_PBS_vals;
                FFvalsMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).FFminusBaseWithinWindow_MUSC_vals;
                % convert FFvals to actual values, not diff from base
                FFvalsPBS=FFvalsPBS+SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
                FFvalsMUSC=FFvalsMUSC+SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF_WithinTimeWindow;
                end
                
                % make sure tvals correspond to ffvals
                assert(length(TvalsPBS)==length(FFvalsPBS) & length(TvalsMUSC)==length(FFvalsMUSC), 'PRoblem - tvals and ffvals dont match');
                
                
                % -- for each day, calculate CV
                TvalsPBS_round=floor(TvalsPBS);
                TvalsMUSC_round=floor(TvalsMUSC);
                
                days=unique(TvalsPBS_round);
                
                cvPBS_alldays=[];
                cvMUSC_alldays=[];
                ffvals_DivideDayMean_AllDays_PBS=[];
                ffvals_DivideDayMean_AllDays_MUSC=[];
                % ffvals_MinusDayMean_AllDays_PBS=[];
                % ffvals_MinusDayMean_AllDays_MUSC=[];
                
                if isempty(days)
                    continue
                end
                
                for k=days
                    % for each day, get PBS and MUSC cv
                    
                    % --- PBS
                    inds=TvalsPBS_round==k; % only this day
                    ffvals=FFvalsPBS(inds);
                    cvPBS=std(ffvals)/abs(mean(ffvals));
                    
                    %                    ffvals_MinusMean=ffvals-mean(ffvals);
                    %                    ffvals_MinusDayMean_AllDays_PBS=[ffvals_MinusDayMean_AllDays_PBS ffvals_MinusMean];
                    
                    ffvals_DivideDayMean_AllDays_PBS=[ffvals_DivideDayMean_AllDays_PBS ffvals/mean(ffvals)];
                    
                    
                    
                    % --- MUSC
                    inds=TvalsMUSC_round==k;
                    ffvals=FFvalsMUSC(inds);
                    cvMUSC=std(ffvals)/abs(mean(ffvals));
                    
                    %                    ffvals_MinusMean=ffvals-mean(ffvals);
                    %                    ffvals_MinusDayMean_AllDays_MUSC=[ffvals_MinusDayMean_AllDays_MUSC ffvals_MinusMean];
                    
                    ffvals_DivideDayMean_AllDays_MUSC=[ffvals_DivideDayMean_AllDays_MUSC ffvals/mean(ffvals)];
                    
                    
                    % ====== save cv vals
                    cvPBS_alldays=[cvPBS_alldays cvPBS];
                    cvMUSC_alldays=[cvMUSC_alldays cvMUSC];
                    
                end
                
                % === get 1) mean of day CVs, and 2) CV over all days
                % (detrended)
%                 MeanOfDayCVs_PBS=mean(cvPBS_alldays);
%                 MeanOfDayCVs_MUSC=mean(cvMUSC_alldays);
                
                CVofAllDays_UsingValsDividedByDayMean_PBS=std(ffvals_DivideDayMean_AllDays_PBS);
                CVofAllDays_UsingValsDividedByDayMean_MUSC=std(ffvals_DivideDayMean_AllDays_MUSC);
                
                % ratio of MUSC CV to PBS CV, and whether that is
                % significant
                cvRatio_MUSCoverPBS_usingAllVals=CVofAllDays_UsingValsDividedByDayMean_MUSC/CVofAllDays_UsingValsDividedByDayMean_PBS;
                [~, p]=vartest2(ffvals_DivideDayMean_AllDays_PBS, ffvals_DivideDayMean_AllDays_MUSC, 'tail', 'right');
                if DispEachSylCVpval==1
                disp([birdname '-' exptname '-' syl ': ' num2str(cvRatio_MUSCoverPBS_usingAllVals), '; p=' num2str(p)]);
                end
%                 plot(CVofAllDays_UsingValsDividedByDayMean_MUSC, CVofAllDays_UsingValsDividedByDayMean_PBS, 'o');
                


                % ===== OUTPUT DATA
                % --- cv related stuff
                cvRatio_MUSCoverPBS_usingAllVals_ALLSYLS=[cvRatio_MUSCoverPBS_usingAllVals_ALLSYLS cvRatio_MUSCoverPBS_usingAllVals];
                cvRatio_pvalue_UsingAllVals_ALLSYLS=[cvRatio_pvalue_UsingAllVals_ALLSYLS p];
                
                CVofAllDays_UsingValsDividedByDayMean_PBS_ALLSYLS=[CVofAllDays_UsingValsDividedByDayMean_PBS_ALLSYLS CVofAllDays_UsingValsDividedByDayMean_PBS];
                CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLSYLS=[CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLSYLS CVofAllDays_UsingValsDividedByDayMean_MUSC];
                
                cvPBS_alldays_ALLSYLS=[cvPBS_alldays_ALLSYLS cvPBS_alldays];
                cvMUSC_alldays_ALLSYLS=[cvMUSC_alldays_ALLSYLS  cvMUSC_alldays];
                % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                
                
                
                % --- other stuff              
                Y_PosRelTarg=[Y_PosRelTarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).distance_from_targ];
                Y_FFmean_pbs=[Y_FFmean_pbs FF_PBS];
                Y_FFmean_musc=[Y_FFmean_musc FF_MUSC];
                Y_FFsem_pbs=[Y_FFsem_pbs FFsem_PBS];
                Y_FFsem_musc=[Y_FFsem_musc FFsem_MUSC];
                Y_syls=[Y_syls, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
                Y_similar_diff=[Y_similar_diff SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
                Y_istarg=[Y_istarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
                Y_AFP_bias=[Y_AFP_bias FF_PBS-FF_MUSC];
                Y_presimilar=[Y_presimilar SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl];
%                 Y_AcousticDist=[Y_AcousticDist SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];                
%                 targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
%                 Y_Corr=[Y_Corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl)];
                
                Yexpt=[Yexpt exptname(end-3:end)];
                Ybird=[Ybird birdname(1:4)];

                Expt_count_all=[Expt_count_all expt_count];
                
                
                Y_Generalization_MP=[Y_Generalization_MP Generalization_MP];
                Y_Generalization_AFP=[Y_Generalization_AFP Generalization_AFP];
                Y_Generalization_Learn=[Y_Generalization_Learn Generalization_Learn];
                
                TargLearnDirAll=[TargLearnDirAll targlearndir];
                
            end
            
            
            % ================= Flip sign if learning at targsyl is negative
            if Y_FFmean_pbs(Y_istarg==1)<0; % negative learning
                Y_FFmean_pbs=-1.*Y_FFmean_pbs;
                Y_FFmean_musc=-1.*Y_FFmean_musc;
                Y_AFP_bias=-1.*Y_AFP_bias;
            end
            
            % ========= Normalize by targsyl if desired (PBS learning
            % by taergsyl)
            if norm_by_targsyl==1;
                learning_by_targ=Y_FFmean_pbs(Y_istarg==1);
                
                Y_FFmean_pbs=Y_FFmean_pbs./learning_by_targ;
                Y_FFmean_musc=Y_FFmean_musc./learning_by_targ;
                Y_AFP_bias=Y_AFP_bias./learning_by_targ;
            end
            
            
            
            % ============================ COLLECT DATA TO PLOT FOR ALL
            % EXPERIMENTS
            if any(~isnan(Y_FFmean_pbs)); % if any are not nan.
                
                % -- cv stuff
                cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS=[cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS cvRatio_MUSCoverPBS_usingAllVals_ALLSYLS];
                cvRatio_pvalue_UsingAllVals_ALLEXPTS=[cvRatio_pvalue_UsingAllVals_ALLEXPTS cvRatio_pvalue_UsingAllVals_ALLSYLS];
                
                CVofAllDays_UsingValsDividedByDayMean_PBS_ALLEXPTS=[CVofAllDays_UsingValsDividedByDayMean_PBS_ALLEXPTS CVofAllDays_UsingValsDividedByDayMean_PBS_ALLSYLS];
                CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLEXPTS=[CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLEXPTS CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLSYLS];
                
                cvPBS_alldays_ALLEXPTS=[cvPBS_alldays_ALLEXPTS cvPBS_alldays_ALLSYLS];
                cvMUSC_alldays_ALLEXPTS=[cvMUSC_alldays_ALLEXPTS  cvMUSC_alldays_ALLSYLS];

                
                % -- other stuff
                LearningPBS_all=[LearningPBS_all Y_FFmean_pbs];
                MPbias_all=[MPbias_all Y_FFmean_musc];
                AFPbias_all=[AFPbias_all Y_AFP_bias];
                SimDiff_all=[SimDiff_all Y_similar_diff];
                TargStatus_all=[TargStatus_all Y_istarg];
                PreSylSimilar_all=[PreSylSimilar_all Y_presimilar];
                
                Generalization_AFP_all=[Generalization_AFP_all Y_Generalization_AFP];
                 Generalization_MP_all=[Generalization_MP_all Y_Generalization_MP];
                 Generalization_Learn_all=[Generalization_Learn_all Y_Generalization_Learn];
               
                
                
                Yexpt_all=[Yexpt_all Yexpt];
                Ybird_all=[Ybird_all Ybird];
                Y_PosRelTarg_All=[Y_PosRelTarg_All Y_PosRelTarg];
                
                
                
                expt_count=expt_count+1;
                
            end
        end
    end
end


%% ==== display num expts/birds for both targ and nontarg

% TARG
disp(' ==== TARG SYLS');
disp(['numexpts: ' num2str(max(Expt_count_all))]);
disp(['numbirds: ' num2str(length(unique(Ybird_all)))]);

% NONTARG (SAME)
disp(' ==== SAME NONTARG')
inds=find(TargStatus_all==0 & SimDiff_all==1);

exptnumtmp=[];
birdnametmp={};
for i=1:length(inds)
    ind=inds(i);
    
    expt=Expt_count_all(ind);
    bname=Ybird_all{ind};
    
    exptnumtmp=[exptnumtmp expt];
    birdnametmp=[birdnametmp bname];
end

disp(['n: ' num2str(length(inds))])
disp(['numexpts: ' num2str(length(unique(exptnumtmp)))]);
disp(['numbirds: ' num2str(length(unique(birdnametmp)))]);

%% ==== BAR PLOTS OF LEARNING VS. MUSC [2]

Learning_thresh=-10000;

lt_figure;
hold on;
title('PBS vs. MUSC');

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PRESIM
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PREDIFF
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% -------- DIFF/PRESIM
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];

% -------- DIFF/PREDIFF
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];

% ---------- ALL NONTARGS
inds=TargStatus_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];

% -------- ALL (EXCEPT DIFF TYPE, DIFF SEQ)
inds=TargStatus_all==0 & ~(SimDiff_all==0 & PreSylSimilar_all==0) & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% ---------- SAME TYPE
inds=TargStatus_all==0 & SimDiff_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];

% ---------- DIFF TYPE
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];

% ---------- PRESIM
inds=TargStatus_all==0 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ---------- PREDIFF
inds=TargStatus_all==0 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ================= PLOT (RAW)
X=1:length(Ylearn_raw);

for i=1:length(X);
    xtmp=[i-0.2, i+0.2];
    ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
    plot(xtmp, ytmp, '-k')
    
    plot(i-0.2, Ylearn_raw{i}, 'ok');
    plot(i+0.2, YMP_raw{i}, 'or');
end


% ================== PLOT MEAN
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type, same-seq','same-type, diff-seq', 'diff-type same seq', 'diff type, diff seq', 'ALL NONTARGS', 'ALL[EXCEPT DTDS]', 'SAMETYPE', 'DIFFTYPE', 'PRESIM[ALL]', 'PREDIFF[ALL]'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)




%% ==== BAR PLOTS OF LEARNING VS. MUSC

Learning_thresh=-10000;

count=1;
SubplotsPerFig=8;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PRESIM
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PREDIFF
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% -------- DIFF
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% ================= PLOT (RAW)
X=1:length(Ylearn_raw);

for i=1:length(X);
    xtmp=[i-0.2, i+0.2];
    ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
    plot(xtmp, ytmp, '-k')
    
    plot(i-0.2, Ylearn_raw{i}, 'ok');
    plot(i+0.2, YMP_raw{i}, 'or');
end


% ================== PLOT MEAN
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type, same-seq','same-type, diff-seq', 'diff-type'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)



%% ==== BAR PLOTS OF LEARNING VS. MUSC [JUST TARG, STSQ, AND OTHERS]

Learning_thresh=-10000;

count=1;
SubplotsPerFig=8;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PRESIM
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PREDIFF + DIFF
inds=TargStatus_all==0 & ((SimDiff_all==1 & PreSylSimilar_all==0)...
    | SimDiff_all==0) ...
    & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ================= PLOT (RAW)
X=1:length(Ylearn_raw);

for i=1:length(X);
    xtmp=[i-0.2, i+0.2];
    ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
    plot(xtmp, ytmp, '-k')
    
    plot(i-0.2, Ylearn_raw{i}, 'ok');
    plot(i+0.2, YMP_raw{i}, 'or');
end


% ================== PLOT MEAN
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type, same-seq','others'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)




%% ==== BAR PLOTS OF LEARNING VS. MUSC [

Learning_thresh=-10000;

count=1;
SubplotsPerFig=8;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PRESIM + DTSS
inds=TargStatus_all==0 & (SimDiff_all==1 ...
    | (SimDiff_all==0 & PreSylSimilar_all==1)) ...
    & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ---- DIFF TYPE, DIFF SEQ
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% ================= PLOT (RAW)
X=1:length(Ylearn_raw);

for i=1:length(X);
    xtmp=[i-0.2, i+0.2];
    ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
    plot(xtmp, ytmp, '-k')
    
    plot(i-0.2, Ylearn_raw{i}, 'ok');
    plot(i+0.2, YMP_raw{i}, 'or');
end


% ================== PLOT MEAN
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','SameType & DTSS','Diff Type Diff Trans'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)




%% ==== BAR PLOTS OF LEARNING VS. MUSC [NO INDIVIDUAL DATAPOINTS]

Learning_thresh=-10000;


[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];




% ------- SIMILAR/ PRESIM
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PREDIFF
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% -------- DIFF
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% ================= PLOT (RAW)
X=1:length(Ylearn_raw);

% for i=1:length(X);
%     xtmp=[i-0.2, i+0.2];
%     ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
%     plot(xtmp, ytmp, '-k')
%     
%     plot(i-0.2, Ylearn_raw{i}, 'ok');
%     plot(i+0.2, YMP_raw{i}, 'or');
% end


% ================== PLOT MEAN
for i=1:length(Ylearn_raw);
    Ylearn_mean=mean(Ylearn_raw{i});
Ymp_mean=mean(YMP_raw{i});

Ylearn_sem=lt_sem(Ylearn_raw{i});
Ymp_sem=lt_sem(YMP_raw{i});

% -- highlight if learning is significant
p = signrank(Ylearn_raw{i});

if p<0.05;
    lt_plot_bar(i-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
    hold on;
    lt_plot_bar(i+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});
else
    hbar=lt_plot_bar(i-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
    set(hbar, 'LineStyle','--');
    hold on;
    hbar=lt_plot_bar(i+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});
    set(hbar, 'LineStyle','--');
end
end

% 
% if p<0.05;
%     lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
%     hold on;
%     lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});
% else
%     hbar=lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
%     set(hbar, 'LineStyle','--');
%     hold on;
%     lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});
%     set(hbar, 'LineStyle','--');
% end
% 

% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i-0.1, 1.1*mean(Ylearn_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, 1.1*mean(Ylearn_raw{i}), '***', 'r', 15);
    elseif p<0.005
        lt_plot_text(i, 1.1*mean(Ylearn_raw{i}), '**', 'r', 15);
    elseif p<0.05
        lt_plot_text(i, 1.1*mean(Ylearn_raw{i}), '*', 'r', 15);
    end

end

% ---------- LEARNING VS 0
% for i=1:length(Ylearn_raw);
%    
%     p = signrank(Ylearn_raw{i});
%     
%     if p<0.0005
%         lt_plot_text(i-0.2, 0, '***', 'c', 15);
%     elseif p<0.005
%         lt_plot_text(i-0.2, 0, '**', 'c', 15);
%     elseif p<0.05
%         lt_plot_text(i-0.2, 0, '*', 'c', 15);
%     end
% end


% % ---------- MP VS 0
% for i=1:length(YMP_raw);
%     
%     p = signrank(YMP_raw{i});
%     
%     if p<0.0005
%         lt_plot_text(i+0.2, 0, '***', 'c', 15);
%     elseif p<0.005
%         lt_plot_text(i+0.2, 0, '**', 'c', 15);
%     elseif p<0.05
%         lt_plot_text(i+0.2, 0, '*', 'c', 15);
%     end
% 
% end
% 


Xlabels={'Targets','same-type, same-seq','same-type, diff-seq', 'diff-type'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)


%% PLOT RATIOS

%% == PLOT DISTRIBUTIONS OF AFP BIAS
plotPDF=0;

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('distributions of AFP bias');

hbar_all=[];

% ------ TARGETS
inds=TargStatus_all==1;

afp=AFPbias_all(inds);

[~,~,hbar]=lt_plot_histogram(afp, '', 1, plotPDF, '', 1, 'k');
hbar_all=[hbar_all hbar];


% ----- SAME TYPE , SAME SEQ
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1;

afp=AFPbias_all(inds);

[~,~,hbar]=lt_plot_histogram(afp, '', 1, plotPDF, '', 1, 'b');
hbar_all=[hbar_all hbar];


% ----- SAME TYPE , DIFF SEQ
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0;

afp=AFPbias_all(inds);

[~,~,hbar]=lt_plot_histogram(afp, '', 1, plotPDF, '', 1, 'c');
hbar_all=[hbar_all hbar];


% ----- DIFF TYPE
inds=TargStatus_all==0 & SimDiff_all==0;

afp=AFPbias_all(inds);

[~,~,hbar]=lt_plot_histogram(afp, '', 1, plotPDF, '', 1,'r');
hbar_all=[hbar_all hbar];


% =====
Legendlabel={'Targets','Same type, same seq', 'Same type, diff seq', 'Diff type'};
legend(hbar_all, Legendlabel);

line([0 0], ylim, 'Color','k', 'LineStyle','--')


%% +++++++++++++++++++ AFP DISTRUBTIONS (SEPARATE PLOTS)
count=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


hsplots=[];
plotPDF=1;

% get Xcenters
tmp=[AFPbias_all MPbias_all];
[~,Xcenters]=lt_plot_histogram(tmp, '', 0, plotPDF, '', 1, 'k');
Xcenters=linspace(Xcenters(1), Xcenters(end),  21); 

%% == PLOT DISTRIBUTIONS OF AFP BIAS [targ]
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[TARG] distributions of AFP bias');
hsplots=[hsplots hsplot];

% ------ TARGETS
inds=TargStatus_all==1;
afp=AFPbias_all(inds);
[~,~,hbar]=lt_plot_histogram(afp, Xcenters, 1, plotPDF, '', 1, 'k');
line([0 0], ylim, 'Color', 'k','LineStyle','--');


%% == PLOT DISTRIBUTIONS OF AFP BIAS [SAME TYPE]

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[SAME TYPE] distributions of AFP bias');
hsplots=[hsplots hsplot];

% ------ 
inds=TargStatus_all==0 & SimDiff_all==1;
color='b';
afp=AFPbias_all(inds);
[~,~,hbar]=lt_plot_histogram(afp, Xcenters, 1, plotPDF, '', 1, color);
line([0 0], ylim, 'Color', 'k','LineStyle','--');


%% == PLOT DISTRIBUTIONS OF AFP BIAS [DIFF TYPE]

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[DIFF TYPE] distributions of AFP bias');
hsplots=[hsplots hsplot];

% ------ 
inds=TargStatus_all==0 & SimDiff_all==0;
color='r';
afp=AFPbias_all(inds);
[~,~,hbar]=lt_plot_histogram(afp, Xcenters, 1, plotPDF, '', 1, color);
line([0 0], ylim, 'Color', 'k','LineStyle','--');


%% +++++++++++++++++++ MP DISTRUBTIONS (SEPARATE PLOTS)

%% == PLOT DISTRIBUTIONS OF MP BIAS [targ]
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[TARG] distributions of MP bias');
hsplots=[hsplots hsplot];

inds=TargStatus_all==1;
color='k';

mp=MPbias_all(inds);
[~,~,hbar]=lt_plot_histogram(mp, Xcenters, 1, plotPDF, '', 1, color);
line([0 0], ylim, 'Color', 'k','LineStyle','--');

%% == PLOT DISTRIBUTIONS OF MP BIAS [SAME]
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[SAME] distributions of MP bias');
hsplots=[hsplots hsplot];

inds=TargStatus_all==0 & SimDiff_all==1;
color='b';

mp=MPbias_all(inds);
[~,~,hbar]=lt_plot_histogram(mp, Xcenters, 1, plotPDF, '', 1, color);
line([0 0], ylim, 'Color', 'k','LineStyle','--');

%% == PLOT DISTRIBUTIONS OF MP BIAS [SAME]
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[DIFF] distributions of MP bias');
hsplots=[hsplots hsplot];

inds=TargStatus_all==0 & SimDiff_all==0;
color='r';

mp=MPbias_all(inds);
[~,~,hbar]=lt_plot_histogram(mp, Xcenters, 1, plotPDF, '', 1, color);
line([0 0], ylim, 'Color', 'k','LineStyle','--');


linkaxes(hsplots, 'xy');


%% ==== Plot distributions of AFP bias, bar plot form [ALL (ADJACENCY]

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('AFP bias');

hbar_all=[];

% ------ TARGETS
x=1;
inds=TargStatus_all==1;
color='k';

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});



% ----- SAME TYPE , SAME SEQ
x=2; color='b';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1;

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});



% ----- SAME TYPE , DIFF SEQ
x=3; color='c';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0;

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});

% ----- DIFF TYPE
x=4;
color='r';
inds=TargStatus_all==0 & SimDiff_all==0;

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});



% =====
Legendlabel={'Targets','Same type, same seq', 'Same type, diff seq', 'Diff type'};
set(gca, 'XTick', 1:length(Legendlabel));
set(gca, 'XTickLabel', Legendlabel);
rotateXLabels(gca, 90)


%% ==== Plot distributions of AFP bias, bar plot form [ADJACENT ONLY]

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('AFP bias [ADJACENT]');

hbar_all=[];

% ------ TARGETS
x=1;
inds=TargStatus_all==1;
color='k';

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});



% ----- SAME TYPE , SAME SEQ
x=2; color='b';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1  & (Y_PosRelTarg_All==-1 | Y_PosRelTarg_All==1);
if any(inds)
afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
end


% ----- SAME TYPE , DIFF SEQ
x=3; color='c';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0 & (Y_PosRelTarg_All==-1 | Y_PosRelTarg_All==1);

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});

% ----- DIFF TYPE
x=4;
color='r';
inds=TargStatus_all==0 & SimDiff_all==0 & (Y_PosRelTarg_All==-1 | Y_PosRelTarg_All==1);

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});



% =====
Legendlabel={'Targets','Same type, same seq', 'Same type, diff seq', 'Diff type'};
set(gca, 'XTick', 1:length(Legendlabel));
set(gca, 'XTickLabel', Legendlabel);
rotateXLabels(gca, 90)


%% ==== Plot distributions of AFP bias, bar plot form [NOT ADJACENT ONLY]

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('AFP bias [NOT ADJACENT]');

hbar_all=[];

% ------ TARGETS
x=1;
inds=TargStatus_all==1;
color='k';

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});



% ----- SAME TYPE , SAME SEQ
x=2; color='b';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1  & ~(Y_PosRelTarg_All==-1 | Y_PosRelTarg_All==1);

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});



% ----- SAME TYPE , DIFF SEQ
x=3; color='c';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0 & ~(Y_PosRelTarg_All==-1 | Y_PosRelTarg_All==1);

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});

% ----- DIFF TYPE
x=4;
color='r';
inds=TargStatus_all==0 & SimDiff_all==0 & ~(Y_PosRelTarg_All==-1 | Y_PosRelTarg_All==1);

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});



% =====
Legendlabel={'Targets','Same type, same seq', 'Same type, diff seq', 'Diff type'};
set(gca, 'XTick', 1:length(Legendlabel));
set(gca, 'XTickLabel', Legendlabel);
rotateXLabels(gca, 90)



%% ====== PLOT DISTRUBTIONS OF MP BIAS
plotPDF=0;

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('distributions of MP bias');

hbar_all=[];

% ------ TARGETS
inds=TargStatus_all==1;

afp=MPbias_all(inds);

[~,~,hbar]=lt_plot_histogram(afp, '', 1, plotPDF, '', 1, 'k');
hbar_all=[hbar_all hbar];



% ----- SAME TYPE , SAME SEQ
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1;

afp=MPbias_all(inds);

[~,~,hbar]=lt_plot_histogram(afp, '', 1, plotPDF, '', 1, 'b');
hbar_all=[hbar_all hbar];


% ----- SAME TYPE , DIFF SEQ
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0;

afp=MPbias_all(inds);

[~,~,hbar]=lt_plot_histogram(afp, '', 1, plotPDF, '', 1, 'c');
hbar_all=[hbar_all hbar];


% ----- DIFF TYPE
inds=TargStatus_all==0 & SimDiff_all==0;

afp=MPbias_all(inds);

[~,~,hbar]=lt_plot_histogram(afp, '', 1, plotPDF, '', 1,'r');
hbar_all=[hbar_all hbar];

% =====
Legendlabel={'Targets','Same type, same seq', 'Same type, diff seq', 'Diff type'};
legend(hbar_all, Legendlabel);
line([0 0], ylim, 'Color','k', 'LineStyle','--')



%% ==== Plot distributions of MP bias, bar plot form

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('MP bias');

hbar_all=[];

% ------ TARGETS
x=1;
inds=TargStatus_all==1;
color='k';

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});



% ----- SAME TYPE , SAME SEQ
x=2; color='b';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1;

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});



% ----- SAME TYPE , DIFF SEQ
x=3; color='c';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0;

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});

% ----- DIFF TYPE
x=4;
color='r';
inds=TargStatus_all==0 & SimDiff_all==0;

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});



% =====
Legendlabel={'Targets','Same type, same seq', 'Same type, diff seq', 'Diff type'};
set(gca, 'XTick', 1:length(Legendlabel));
set(gca, 'XTickLabel', Legendlabel);
rotateXLabels(gca, 90)



%% ==== Plot distributions of MP bias, bar plot form [adjacent

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('MP bias [ADJACENT]');

hbar_all=[];

% ------ TARGETS
x=1;
inds=TargStatus_all==1;
color='k';

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});



% ----- SAME TYPE , SAME SEQ
x=2; color='b';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1  & (Y_PosRelTarg_All==-1 | Y_PosRelTarg_All==1);

afp=MPbias_all(inds);
if any(inds)
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
end


% ----- SAME TYPE , DIFF SEQ
x=3; color='c';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0  & (Y_PosRelTarg_All==-1 | Y_PosRelTarg_All==1);

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});

% ----- DIFF TYPE
x=4;
color='r';
inds=TargStatus_all==0 & SimDiff_all==0  & (Y_PosRelTarg_All==-1 | Y_PosRelTarg_All==1);

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});



% =====
Legendlabel={'Targets','Same type, same seq', 'Same type, diff seq', 'Diff type'};
set(gca, 'XTick', 1:length(Legendlabel));
set(gca, 'XTickLabel', Legendlabel);
rotateXLabels(gca, 90)


%% ==== Plot distributions of MP bias, bar plot form [adjacent

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('MP bias [NOT ADJACENT]');

hbar_all=[];

% ------ TARGETS
x=1;
inds=TargStatus_all==1;
color='k';

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});



% ----- SAME TYPE , SAME SEQ
x=2; color='b';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1  & ~(Y_PosRelTarg_All==-1 | Y_PosRelTarg_All==1);

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});



% ----- SAME TYPE , DIFF SEQ
x=3; color='c';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0  & ~(Y_PosRelTarg_All==-1 | Y_PosRelTarg_All==1);

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});

% ----- DIFF TYPE
x=4;
color='r';
inds=TargStatus_all==0 & SimDiff_all==0  & ~(Y_PosRelTarg_All==-1 | Y_PosRelTarg_All==1);

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});



% =====
Legendlabel={'Targets','Same type, same seq', 'Same type, diff seq', 'Diff type'};
set(gca, 'XTick', 1:length(Legendlabel));
set(gca, 'XTickLabel', Legendlabel);
rotateXLabels(gca, 90)


%% ==== Plot distributions of CV REDUCTION

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('CV reduction (filled = not sign.; (signi/total))');

hbar_all=[];

% ------ TARGETS
x=1;
inds=TargStatus_all==1;
color='k';

cv_ratio=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);
cv_ratio_pval=cvRatio_pvalue_UsingAllVals_ALLEXPTS(inds);

% --- if pval <0.05, plot open
inds2=cv_ratio_pval<0.05;
plot(x+0.1, cv_ratio(inds2), 'o', 'Color', color)
% --- if pval>=0.05, plot solid
inds2=cv_ratio_pval>=0.05;
if any(inds2)
plot(x+0.1, cv_ratio(inds2), 'o', 'Color', color, 'MarkerFaceColor',color)
end

lt_plot_bar(x, mean(cv_ratio), {'Errors', lt_sem(cv_ratio), 'Color',color});
% --- overlay number that are themselves significant
NumberSigni=sum(cv_ratio_pval<0.05);
NumberTot=length(cv_ratio_pval);
lt_plot_text(x-0.2,  max(cv_ratio)+0.06, ['(' num2str(NumberSigni) '/' num2str(NumberTot) ')'], color);




% ----- SAME TYPE , SAME SEQ
x=2; color='b';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1;

cv_ratio=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);
cv_ratio_pval=cvRatio_pvalue_UsingAllVals_ALLEXPTS(inds);

% --- if pval <0.05, plot open
inds2=cv_ratio_pval<0.05;
plot(x+0.1, cv_ratio(inds2), 'o', 'Color', color)
% --- if pval>=0.05, plot solid
inds2=cv_ratio_pval>=0.05;
if any(inds2)
plot(x+0.1, cv_ratio(inds2), 'o', 'Color', color, 'MarkerFaceColor',color)
end
lt_plot_bar(x, mean(cv_ratio), {'Errors', lt_sem(cv_ratio), 'Color',color});

% --- overlay number that are themselves significant
NumberSigni=sum(cv_ratio_pval<0.05);
NumberTot=length(cv_ratio_pval);
lt_plot_text(x-0.2,  max(cv_ratio)+0.06, ['(' num2str(NumberSigni) '/' num2str(NumberTot) ')'], color);



% ----- SAME TYPE , DIFF SEQ
x=3; color='c';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0;

cv_ratio=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);
cv_ratio_pval=cvRatio_pvalue_UsingAllVals_ALLEXPTS(inds);

% --- if pval <0.05, plot open
inds2=cv_ratio_pval<0.05;
if any(inds2)>0;
plot(x+0.1, cv_ratio(inds2), 'o', 'Color', color)
end
% --- if pval>=0.05, plot solid
inds2=cv_ratio_pval>=0.05;
if any(inds2)
plot(x+0.1, cv_ratio(inds2), 'o', 'Color', color, 'MarkerFaceColor',color)
end

lt_plot_bar(x, mean(cv_ratio), {'Errors', lt_sem(cv_ratio), 'Color',color});
% --- overlay number that are themselves significant
NumberSigni=sum(cv_ratio_pval<0.05);
NumberTot=length(cv_ratio_pval);
lt_plot_text(x-0.2,  max(cv_ratio)+0.06, ['(' num2str(NumberSigni) '/' num2str(NumberTot) ')'], color);



% ----- DIFF TYPE
x=4;
color='r';
inds=TargStatus_all==0 & SimDiff_all==0;

cv_ratio=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);
cv_ratio_pval=cvRatio_pvalue_UsingAllVals_ALLEXPTS(inds);

% --- if pval <0.05, plot open
inds2=cv_ratio_pval<0.05;
plot(x+0.1, cv_ratio(inds2), 'o', 'Color', color)
% --- if pval>=0.05, plot solid
inds2=cv_ratio_pval>=0.05;
if any(inds2)
plot(x+0.1, cv_ratio(inds2), 'o', 'Color', color, 'MarkerFaceColor',color)
end

lt_plot_bar(x, mean(cv_ratio), {'Errors', lt_sem(cv_ratio), 'Color',color});
% --- overlay number that are themselves significant
NumberSigni=sum(cv_ratio_pval<0.05);
NumberTot=length(cv_ratio_pval);
lt_plot_text(x-0.2,  max(cv_ratio)+0.06, ['(' num2str(NumberSigni) '/' num2str(NumberTot) ')'], color);


% =====
Legendlabel={'Targets','Same type, same seq', 'Same type, diff seq', 'Diff type'};
set(gca, 'XTick', 1:length(Legendlabel));
set(gca, 'XTickLabel', Legendlabel);
rotateXLabels(gca, 90)

line(xlim, [1 1], 'LineStyle' ,'--','Color','k')

line(xlim, [0.72 0.72], 'LineStyle','--','Color','b')

% ==== make sure no group's cv is diff from targ's
% 1) targ vs. same-type, same-seq
inds=TargStatus_all==1;
cv_ratio1=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);

inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1;
cv_ratio2=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);

p=ranksum(cv_ratio1, cv_ratio2);
disp(['Targ vs. SametypeSameseq; p=' num2str(p)]);


% 1) targ vs. same-type, diff-seq
inds=TargStatus_all==1;
cv_ratio1=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);

inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0;
cv_ratio2=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);

p=ranksum(cv_ratio1, cv_ratio2);
disp(['Targ vs. SametypeDiffseq; p=' num2str(p)]);

% 1) targ vs. difftype
inds=TargStatus_all==1;
cv_ratio1=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);

inds=TargStatus_all==0 & SimDiff_all==0;
cv_ratio2=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);

p=ranksum(cv_ratio1, cv_ratio2);
disp(['Targ vs. Difftype; p=' num2str(p)]);



%% ==== Plot distributions of CV REDUCTION

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('CV reduction (filled = not sign.; (signi/total))');

hbar_all=[];

% ------ TARGETS
x=1;
inds=TargStatus_all==1;
color='k';

cv_ratio=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);
cv_ratio_pval=cvRatio_pvalue_UsingAllVals_ALLEXPTS(inds);

% --- if pval <0.05, plot open
inds2=cv_ratio_pval<0.05;
plot(x+0.1, cv_ratio(inds2), 'o', 'Color', color)
% --- if pval>=0.05, plot solid
inds2=cv_ratio_pval>=0.05;
if any(inds2)
plot(x+0.1, cv_ratio(inds2), 'o', 'Color', color, 'MarkerFaceColor',color)
end

lt_plot_bar(x, mean(cv_ratio), {'Errors', lt_sem(cv_ratio), 'Color',color});
% --- overlay number that are themselves significant
NumberSigni=sum(cv_ratio_pval<0.05);
NumberTot=length(cv_ratio_pval);
lt_plot_text(x-0.2,  max(cv_ratio)+0.06, ['(' num2str(NumberSigni) '/' num2str(NumberTot) ')'], color);
% -- sign rank test
p=signrank(cv_ratio);
lt_plot_text(x,  max(cv_ratio)+1.1, ['p=' num2str(p)])



% ----- SAME TYPE 
x=2; color='b';
inds=TargStatus_all==0 & SimDiff_all==1;

cv_ratio=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);
cv_ratio_pval=cvRatio_pvalue_UsingAllVals_ALLEXPTS(inds);

% --- if pval <0.05, plot open
inds2=cv_ratio_pval<0.05;
plot(x+0.1, cv_ratio(inds2), 'o', 'Color', color)
% --- if pval>=0.05, plot solid
inds2=cv_ratio_pval>=0.05;
if any(inds2)
plot(x+0.1, cv_ratio(inds2), 'o', 'Color', color, 'MarkerFaceColor',color)
end
lt_plot_bar(x, mean(cv_ratio), {'Errors', lt_sem(cv_ratio), 'Color',color});

% --- overlay number that are themselves significant
NumberSigni=sum(cv_ratio_pval<0.05);
NumberTot=length(cv_ratio_pval);
lt_plot_text(x-0.2,  max(cv_ratio)+0.06, ['(' num2str(NumberSigni) '/' num2str(NumberTot) ')'], color);
% -- sign rank test
p=signrank(cv_ratio);
lt_plot_text(x,  max(cv_ratio)+1.1, ['p=' num2str(p)])


% ----- DIFF TYPE
x=3;
color='r';
inds=TargStatus_all==0 & SimDiff_all==0;

cv_ratio=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);
cv_ratio_pval=cvRatio_pvalue_UsingAllVals_ALLEXPTS(inds);

% --- if pval <0.05, plot open
inds2=cv_ratio_pval<0.05;
plot(x+0.1, cv_ratio(inds2), 'o', 'Color', color)
% --- if pval>=0.05, plot solid
inds2=cv_ratio_pval>=0.05;
if any(inds2)
plot(x+0.1, cv_ratio(inds2), 'o', 'Color', color, 'MarkerFaceColor',color)
end

lt_plot_bar(x, mean(cv_ratio), {'Errors', lt_sem(cv_ratio), 'Color',color});
% --- overlay number that are themselves significant
NumberSigni=sum(cv_ratio_pval<0.05);
NumberTot=length(cv_ratio_pval);
lt_plot_text(x-0.2,  max(cv_ratio)+0.06, ['(' num2str(NumberSigni) '/' num2str(NumberTot) ')'], color);
% -- sign rank test
p=signrank(cv_ratio);
lt_plot_text(x,  max(cv_ratio)+1.1, ['p=' num2str(p)])


% =====
Legendlabel={'Targets','Same type', 'Diff type'};
set(gca, 'XTick', 1:length(Legendlabel));
set(gca, 'XTickLabel', Legendlabel);
rotateXLabels(gca, 90)

line(xlim, [1 1], 'LineStyle' ,'--','Color','k')

line(xlim, [0.72 0.72], 'LineStyle','--','Color','b')

% ==== make sure no group's cv is diff from targ's
% 1) targ vs. same-type
inds=TargStatus_all==1;
cv_ratio1=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);

inds=TargStatus_all==0 & SimDiff_all==1;
cv_ratio2=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);

p=ranksum(cv_ratio1, cv_ratio2);
disp(['Targ vs. Sametype; p=' num2str(p)]);


% 1) targ vs. difftype
inds=TargStatus_all==1;
cv_ratio1=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);

inds=TargStatus_all==0 & SimDiff_all==0;
cv_ratio2=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);

p=ranksum(cv_ratio1, cv_ratio2);
disp(['Targ vs. Difftype; p=' num2str(p)]);



%% ========== SEPARATE cv reduction for each experiment
lt_figure; hold on;

NumExpts=max(Expt_count_all);
XlabelsAll={};

for i=1:NumExpts
    
    inds=find(Expt_count_all==i);
    
    % --- collect CV ratio for all syls
    cvratios=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);
    
    plot(i, cvratios, 'ok');
    lt_plot_bar(i, mean(cvratios), {'Errors',lt_sem(cvratios)});
   
    % --- lower than 1?
    p=signrank(cvratios, 1, 'tail', 'left');
    lt_plot_text(i, max(cvratios)+0.01, ['p=' num2str(p, '%3.2g')], 'r');
    
    % --- collect name
    XlabelsAll=[XlabelsAll [Ybird_all{inds(1)} '-' Yexpt_all{inds(1)}]];
    
end
    
   line(xlim, [0.72 0.72], 'LineStyle','--','Color','b')
   line(xlim, [1 1], 'LineStyle','--','Color','b')

set(gca, 'XTick', 1:length(XlabelsAll), 'XTickLabel', XlabelsAll)
rotateXLabels(gca, 90);


%% ========= greater reversion if greater cv reduction?
lt_figure; hold on;

% ------ TARGETS (consolidiation)
lt_subplot(2,2,1); hold on;
title('targets');
xlabel('cv(MUSC)/cv(PBS)');
ylabel('consolidation');

inds=TargStatus_all==1;
color='k';

cv_ratio=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);
learning=LearningPBS_all(inds);
mpbias=MPbias_all(inds);
consolidation = mpbias./learning;

x=cv_ratio;
y=consolidation;

lt_regress(y, x, 1, 0, 1, 1, color);
xlim([-0.1 1.1]); ylim([-0.1 1.1]); line([0 0], ylim); line([1 1], ylim); line(xlim, [0 0]); line(xlim, [1 1]);


% ------ TARGETS (reversion hz, )
lt_subplot(2,2,2); hold on;
title('targets');
xlabel('cv(MUSC)/cv(PBS)');
ylabel('reversion (hz)');

inds=TargStatus_all==1;
color='k';

cv_ratio=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);
learning=LearningPBS_all(inds);
mpbias=MPbias_all(inds);
reversion = sign(learning).*(learning-mpbias);

x=cv_ratio;
y=reversion;

lt_regress(y, x, 1, 0, 1, 1, color);
xlim([-0.1 1.1]); ylim([-40 150]); line([0 0], ylim); line([1 1], ylim); line(xlim, [0 0])

% ------ SAME (reversion)
lt_subplot(2,2,3); hold on;
title('same type');
xlabel('cv(MUSC)/cv(PBS)');
ylabel('reversion (hz) (in dir opposite targ shift');

inds=TargStatus_all==0 & SimDiff_all==1;
color='b';

cv_ratio=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);
learning=LearningPBS_all(inds);
mpbias=MPbias_all(inds);
reversion = (learning-mpbias);

x=cv_ratio;
y=reversion;

lt_regress(y, x, 1, 0, 1, 1, color);
xlim([-0.1 1.1]); ylim([-40 150]); line([0 0], ylim); line([1 1], ylim); line(xlim, [0 0])

% ------ SAME (reversion)
lt_subplot(2,2,4); hold on;
title('same type');
xlabel('cv(MUSC)/cv(PBS)');
ylabel('reversion (hz) (in dir opposite own learning dir');

inds=TargStatus_all==0 & SimDiff_all==1;
color='b';

cv_ratio=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);
learning=LearningPBS_all(inds);
mpbias=MPbias_all(inds);
reversion = sign(learning).*(learning-mpbias);

x=cv_ratio;
y=reversion;

lt_regress(y, x, 1, 0, 1, 1, color);
xlim([-0.1 1.1]); ylim([-40 150]); line([0 0], ylim); line([1 1], ylim); line(xlim, [0 0])

if (0)
% ------ DIFF
lt_subplot(4,2,5); hold on;
title('diff type');
xlabel('cv(MUSC)/cv(PBS)');
ylabel('reversion (hz) (in dir opposite targ shift');

inds=TargStatus_all==0 & SimDiff_all==0;
color='r';

cv_ratio=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);
learning=LearningPBS_all(inds);
mpbias=MPbias_all(inds);
reversion = (learning-mpbias);

x=cv_ratio;
y=reversion;

lt_regress(y, x, 1, 0, 1, 1, color);


% ------ 
lt_subplot(4,2,6); hold on;
title('diff');
xlabel('cv(MUSC)/cv(PBS)');
ylabel('reversion (hz) (in dir opposite own learning dir');

inds=TargStatus_all==0 & SimDiff_all==0;
color='r';

cv_ratio=cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS(inds);
learning=LearningPBS_all(inds);
mpbias=MPbias_all(inds);
reversion = sign(learning).*(learning-mpbias);

x=cv_ratio;
y=reversion;

lt_regress(y, x, 1, 0, 1, 1, color);
end

%% ===== AFP BIAS [BAR PLOTS]
count=1;
subplotrows=2;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('AFP bias');

hbar_all=[];


% ------ TARGETS
x=1;
inds=TargStatus_all==1;
color='k';

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});

[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end


% ----- SAME TYPE , SAME SEQ
x=2; color='b';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1;

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end



% ----- SAME TYPE , DIFF SEQ
x=3; color='c';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0;

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end


% ----- DIFF TYPE, same seq
x=4;
color='r';
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==1;

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end

% ----- DIFF TYPE, diff seq
x=5;
color='m';
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==0;

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end


% ---------- ALL NONTARGS
inds=TargStatus_all==0;
x=6;
color='k';


afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end



% -------- ALL (EXCEPT DIFF TYPE, DIFF SEQ)
inds=TargStatus_all==0 & ~(SimDiff_all==0 & PreSylSimilar_all==0) & LearningPBS_all>Learning_thresh;
x=7;
color='k';


afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end



% ---------- SAME TYPE
inds=TargStatus_all==0 & SimDiff_all==1 & LearningPBS_all>Learning_thresh;
x=8;
color='k';

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end




% ---------- DIFF TYPE
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>Learning_thresh;
x=9;
color='k';

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end



% ---------- PRESIM
inds=TargStatus_all==0 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh;
x=10;
color='k';

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end




% ---------- PREDIFF
inds=TargStatus_all==0 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh;
x=11;
color='k';

afp=AFPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end





% =====
Legendlabel={'Targets', 'Same type, same seq', 'Same type, diff seq', 'Diff type, same seq', 'Diff type, diff seq', 'ALL NONTARGS', 'ALL (EXCEPT DTDS)', 'SAME TYPE', 'DIFF TYPE', 'PRESIM', 'PREDIFF'};
set(gca, 'XTick', 1:length(Legendlabel));
set(gca, 'XTickLabel', Legendlabel);
rotateXLabels(gca, 90)


%% ==== Plot distributions of MP bias, bar plot form [ALL THINGS]

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('MP bias');

hbar_all=[];

% ------ TARGETS
x=1;
inds=TargStatus_all==1;
color='k';

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end



% ----- SAME TYPE , SAME SEQ
x=2; color='b';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1;

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end



% ----- SAME TYPE , DIFF SEQ
x=3; color='c';
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0;

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});




% ----- DIFF TYPE , same seq
x=4;
color='r';
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==1;

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end


% ----- DIFF TYPE , diff seq
x=5;
color='m';
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==0;

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end




% ---------- ALL NONTARGS
inds=TargStatus_all==0;
x=6;
color='k';


afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end



% -------- ALL (EXCEPT DIFF TYPE, DIFF SEQ)
inds=TargStatus_all==0 & ~(SimDiff_all==0 & PreSylSimilar_all==0) & LearningPBS_all>Learning_thresh;
x=7;
color='k';


afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end



% ---------- SAME TYPE
inds=TargStatus_all==0 & SimDiff_all==1 & LearningPBS_all>Learning_thresh;
x=8;
color='k';

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end




% ---------- DIFF TYPE
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>Learning_thresh;
x=9;
color='k';

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end





% ---------- PRESIM
inds=TargStatus_all==0 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh;
x=10;
color='k';

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end




% ---------- PREDIFF
inds=TargStatus_all==0 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh;
x=11;
color='k';

afp=MPbias_all(inds);
plot(x, afp, 'o', 'Color', color)
lt_plot_bar(x, mean(afp), {'Errors', lt_sem(afp), 'Color',color});
[~, p]=ttest(afp);
if p<0.005;
    lt_plot_text(x+0.2, 1.1*max(afp), '**', 'r');
elseif p<0.05
    lt_plot_text(x+0.2, 1.1*max(afp), '*', 'r');
end





% =====
Legendlabel={'Targets', 'Same type, same seq', 'Same type, diff seq', 'Diff type, same seq', 'Diff type, diff seq', 'ALL NONTARGS', 'ALL (EXCEPT DTDS)', 'SAME TYPE', 'DIFF TYPE', 'PRESIM', 'PREDIFF'};
set(gca, 'XTick', 1:length(Legendlabel));
set(gca, 'XTickLabel', Legendlabel);
rotateXLabels(gca, 90)
    

%% +++++++++++++++++++= BAR PLOTS OF LEARNING VS MUSC (TAKE INTO ACCOUNT ADJACENCY TO TARGET)
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


%% ==== ALL COMBINED

Learning_thresh=-10000;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hold on;
title('All Nontargs');

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PRESIM
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PREDIFF
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% -------- DIFF/PRESIM
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];

% -------- DIFF/PREDIFF
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ================= PLOT (RAW)
X=1:length(Ylearn_raw);

for i=1:length(X);
    xtmp=[i-0.2, i+0.2];
    ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
    plot(xtmp, ytmp, '-k')
    
    plot(i-0.2, Ylearn_raw{i}, 'ok');
    plot(i+0.2, YMP_raw{i}, 'or');
end


% ================== PLOT MEAN
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type, same-trans','same-type, diff-trans', 'diff-type same-trans', 'diff-type, diff-trans'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)


%% ==== [NOT ADJACENT TO TARG]

Learning_thresh=-10000;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hold on;
title('Not adjacent to targ');

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PRESIM
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh ...
    & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PREDIFF
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh ...
    & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% -------- DIFF/PRESIM
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh ...
    & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];

% -------- DIFF/PREDIFF
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh ...
    & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ================= PLOT (RAW)
X=1:length(Ylearn_raw);

for i=1:length(X);
    xtmp=[i-0.2, i+0.2];
    ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
    plot(xtmp, ytmp, '-k')
    
    plot(i-0.2, Ylearn_raw{i}, 'ok');
    plot(i+0.2, YMP_raw{i}, 'or');
end


% ================== PLOT MEAN
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type, same-trans','same-type, diff-trans', 'diff-type same-trans', 'diff-type, diff-trans'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)


%% ==== [ADJACENT TO TARG]

Learning_thresh=-10000;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hold on;
title('Adjacent to targ');

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PRESIM
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh ...
    & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

if ~any(inds);
   Ylearn_raw=[Ylearn_raw {nan}];
    YMP_raw=[YMP_raw {nan}];
else

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];
end


% ------- SIMILAR/ PREDIFF
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh ...
    & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% -------- DIFF/PRESIM
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh ...
    & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];

% -------- DIFF/PREDIFF
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh ...
    & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ================= PLOT (RAW)
X=1:length(Ylearn_raw);

for i=1:length(X);
    xtmp=[i-0.2, i+0.2];
    ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
    plot(xtmp, ytmp', '-k')
    
    plot(i-0.2, Ylearn_raw{i}, 'ok');
    plot(i+0.2, YMP_raw{i}, 'or');
end


% ================== PLOT MEAN
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end



% musc vs pbs
for i=1:length(Ylearn_raw);
    
    if isnan(Ylearn_raw{i})
        continue
    end
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
        if isnan(Ylearn_raw{i})
        continue
        end

    
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end



% ---------- MP VS 0
for i=1:length(YMP_raw);
            if isnan(Ylearn_raw{i})
        continue
        end

    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type, same-trans','same-type, diff-trans', 'diff-type same-trans', 'diff-type, diff-trans'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)



%% ==== [ADJACENT TO TARG] [same type, diff type]

Learning_thresh=-10000;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hold on;
title('adjacent to targ');

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR
inds=TargStatus_all==0 & SimDiff_all==1 & LearningPBS_all>Learning_thresh ...
    & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% -------- DIFF
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>Learning_thresh ...
    & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% ================= PLOT (RAW)
X=1:length(Ylearn_raw);

for i=1:length(X);
    xtmp=[i-0.2, i+0.2];
    ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
    plot(xtmp, ytmp', '-k')
    
    plot(i-0.2, Ylearn_raw{i}, 'ok');
    plot(i+0.2, YMP_raw{i}, 'or');
end


% ================== PLOT MEAN
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type', 'diff-type'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)



%% ==== [ALL NONTARG] [same type, diff type] [NO BALLS]

Learning_thresh=-10000;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hold on;
title('All nontarg');

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR
inds=TargStatus_all==0 & SimDiff_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% -------- DIFF
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% ================= PLOT (RAW)
X=1:length(Ylearn_raw);

for i=1:length(X);
    xtmp=[i-0.2, i+0.2];
    ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
    plot(xtmp, ytmp', '-k')
end


% ================== PLOT MEAN
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type', 'diff-type'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)



Learning_thresh=-10000;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hold on;
title('All nontarg');

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR
inds=TargStatus_all==0 & SimDiff_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% -------- DIFF
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% ================= PLOT (RAW)
X=1:length(Ylearn_raw);

for i=1:length(X);
    xtmp=[i-0.2, i+0.2];
    ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
    plot(xtmp, ytmp', '-k')
end


% ================== PLOT MEAN
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type', 'diff-type'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)

%% ==== [NOT ADJACENT TO TARG] [same type, diff type]

Learning_thresh=-10000;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hold on;
title('Not adjacent to targ');

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR
inds=TargStatus_all==0 & SimDiff_all==1 & LearningPBS_all>Learning_thresh ...
    & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% -------- DIFF
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>Learning_thresh ...
    & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% ================= PLOT (RAW)
X=1:length(Ylearn_raw);

for i=1:length(X);
    xtmp=[i-0.2, i+0.2];
    ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
    plot(xtmp, ytmp, '-k')
    
    plot(i-0.2, Ylearn_raw{i}, 'ok');
    plot(i+0.2, YMP_raw{i}, 'or');
end


% ================== PLOT MEAN
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type', 'diff-type'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)

%% ==== [ALL NONTARG] [same type, diff type] [NO BALLS] [GOOD]
lt_figure; hold on;

Learning_thresh=-10000;
title('All nontarg');

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR
inds=TargStatus_all==0 & SimDiff_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% -------- DIFF
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% ================= PLOT (RAW)
X=1:length(Ylearn_raw);

for i=1:length(X);
    xtmp=[i-0.2, i+0.2];
    ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
    plot(xtmp, ytmp', '-k')
end


% ================== PLOT MEAN
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});


% === overlay reversion magnitude
for i=1:length(Ylearn_raw);
    
    ReversionAll=YMP_raw{i}-Ylearn_raw{i};
    revMean=mean(ReversionAll);
    revSEM=lt_sem(ReversionAll);
    lt_plot_text(i-0.1, 1.3*max(YMP_raw{i}), ['ReverMean(sem)= ' num2str(revMean) '(' num2str(revSEM) ')'], 'k')
    lt_plot_text(i-0.1, 1.4*max(YMP_raw{i}), ['PBS= ' num2str(mean(Ylearn_raw{i})) '(' num2str(lt_sem(Ylearn_raw{i})) ')'], 'k')
    lt_plot_text(i-0.1, 1.5*max(YMP_raw{i}), ['MUSC= ' num2str(mean(YMP_raw{i})) '(' num2str(lt_sem(YMP_raw{i})) ')'], 'k')
    
end

% ==== compare reversion magnitudes (raw reversion)
OutputString=[];
for i=1:length(Ylearn_raw)
    Reversion1=YMP_raw{i}-Ylearn_raw{i};
    for ii=i+1:length(Ylearn_raw);
        Reversion2=YMP_raw{ii}-Ylearn_raw{ii};
        
        p=ranksum(Reversion1, Reversion2);
        
        OutputString=[OutputString ['| ' num2str(i) ' vs. ' num2str(ii) ': ' num2str(p)]];
        
    end
end

lt_plot_annotation(1, OutputString, 'm')


% ==== overlay sample size
for i=1:length(Ylearn_raw);
    
    N=numel(Ylearn_raw{i});
    lt_plot_text(i-0.1, 1.2*max(YMP_raw{i}), ['N=' num2str(N)], 'g')
end


% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), ['p=' num2str(p)], 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    else
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), ['p=' num2str(p)], 'b', 12);
        
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type', 'diff-type'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)




%% +++++++++++++++++++= BAR PLOTS OF LEARNING VS MUSC (TAKE INTO ACCOUNT ADJACENCY TO TARGET) [JUST MEANS]
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


%% ==== ALL COMBINED

Learning_thresh=-10000;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hold on;
title('All Nontargs');

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PRESIM
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PREDIFF
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% -------- DIFF/PRESIM
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];

% -------- DIFF/PREDIFF
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% % ================= PLOT (RAW)
% % X=1:length(Ylearn_raw);
% 
% for i=1:length(X);
%     xtmp=[i-0.2, i+0.2];
%     ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
%     plot(xtmp, ytmp, '-k')
%     
%     plot(i-0.2, Ylearn_raw{i}, 'ok');
%     plot(i+0.2, YMP_raw{i}, 'or');
% end
% 

% ================== PLOT MEAN
X=1:length(Ylearn_raw);

Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type, same-trans','same-type, diff-trans', 'diff-type same-trans', 'diff-type, diff-trans'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)


%% ==== [NOT ADJACENT TO TARG]

Learning_thresh=-10000;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hold on;
title('Not adjacent to targ');

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PRESIM
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh ...
    & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PREDIFF
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh ...
    & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% -------- DIFF/PRESIM
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh ...
    & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];

% -------- DIFF/PREDIFF
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh ...
    & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% % ================= PLOT (RAW)
% X=1:length(Ylearn_raw);
% 
% for i=1:length(X);
%     xtmp=[i-0.2, i+0.2];
%     ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
%     plot(xtmp, ytmp, '-k')
%     
%     plot(i-0.2, Ylearn_raw{i}, 'ok');
%     plot(i+0.2, YMP_raw{i}, 'or');
% end
% 
% 
% ================== PLOT MEAN
X=1:length(Ylearn_raw);
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type, same-trans','same-type, diff-trans', 'diff-type same-trans', 'diff-type, diff-trans'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)


%% ==== [ADJACENT TO TARG]

Learning_thresh=-10000;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hold on;
title('Adjacent to targ');

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR/ PRESIM
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh ...
    & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

if ~any(inds);
   Ylearn_raw=[Ylearn_raw {nan}];
    YMP_raw=[YMP_raw {nan}];
else

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];
end


% ------- SIMILAR/ PREDIFF
inds=TargStatus_all==0 & SimDiff_all==1 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh ...
    & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% -------- DIFF/PRESIM
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==1 & LearningPBS_all>Learning_thresh ...
    & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];

% -------- DIFF/PREDIFF
inds=TargStatus_all==0 & SimDiff_all==0 & PreSylSimilar_all==0 & LearningPBS_all>Learning_thresh ...
    & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% % ================= PLOT (RAW)
% X=1:length(Ylearn_raw);
% 
% for i=1:length(X);
%     xtmp=[i-0.2, i+0.2];
%     ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
%     plot(xtmp, ytmp', '-k')
%     
%     plot(i-0.2, Ylearn_raw{i}, 'ok');
%     plot(i+0.2, YMP_raw{i}, 'or');
% end


% ================== PLOT MEAN
X=1:length(Ylearn_raw);

Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end



% musc vs pbs
for i=1:length(Ylearn_raw);
    
    if isnan(Ylearn_raw{i})
        continue
    end
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
        if isnan(Ylearn_raw{i})
        continue
        end

    
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end



% ---------- MP VS 0
for i=1:length(YMP_raw);
            if isnan(Ylearn_raw{i})
        continue
        end

    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type, same-trans','same-type, diff-trans', 'diff-type same-trans', 'diff-type, diff-trans'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)



%% ==== [ADJACENT TO TARG] [same type, diff type]

Learning_thresh=-10000;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hold on;
title('adjacent to targ');

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR
inds=TargStatus_all==0 & SimDiff_all==1 & LearningPBS_all>Learning_thresh ...
    & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% -------- DIFF
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>Learning_thresh ...
    & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% 
% % ================= PLOT (RAW)
% X=1:length(Ylearn_raw);
% 
% for i=1:length(X);
%     xtmp=[i-0.2, i+0.2];
%     ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
%     plot(xtmp, ytmp', '-k')
%     
%     plot(i-0.2, Ylearn_raw{i}, 'ok');
%     plot(i+0.2, YMP_raw{i}, 'or');
% end


% ================== PLOT MEAN
X=1:length(Ylearn_raw);
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type', 'diff-type'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)




%% ==== [NOT ADJACENT TO TARG] [same type, diff type]

Learning_thresh=-10000;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hold on;
title('Not adjacent to targ');

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR
inds=TargStatus_all==0 & SimDiff_all==1 & LearningPBS_all>Learning_thresh ...
    & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% -------- DIFF
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>Learning_thresh ...
    & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% % ================= PLOT (RAW)
% X=1:length(Ylearn_raw);
% 
% for i=1:length(X);
%     xtmp=[i-0.2, i+0.2];
%     ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
%     plot(xtmp, ytmp, '-k')
%     
%     plot(i-0.2, Ylearn_raw{i}, 'ok');
%     plot(i+0.2, YMP_raw{i}, 'or');
% end


% ================== PLOT MEAN
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);

X=1:length(Ylearn_raw);

lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type', 'diff-type'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)




%% ==== [ALL NONTARG] [same type, diff type]

lt_figure; hold on;
title('All nontarg');

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR
inds=TargStatus_all==0 & SimDiff_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% -------- DIFF
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% % ================= PLOT (RAW)
% X=1:length(Ylearn_raw);
% 
% for i=1:length(X);
%     xtmp=[i-0.2, i+0.2];
%     ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
%     plot(xtmp, ytmp', '-k')
%     
%     plot(i-0.2, Ylearn_raw{i}, 'ok');
%     plot(i+0.2, YMP_raw{i}, 'or');
% end


% ================== PLOT MEAN
X=1:length(Ylearn_raw);
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});



% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
    
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '***', 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type', 'diff-type'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)

%% +++++++++++++++++++++ AFP AND MP GENERALIZATION
%% =======

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


%% ====== [NOT ADJACENT]

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
disp('GENERALIZATION (afp -- mp -- learning) [NOT ADJACENT]');
title('GENERALIZATION (afp -- mp -- learning) [NOT ADJACENT]');


% ====  SIMILAR
x=1;
color='b';
inds=TargStatus_all==0 & SimDiff_all==1 & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);
GenLearn=Generalization_Learn_all(inds);


X=[x-0.2 x x+0.2];
Y=[GenAFP' GenMP' GenLearn'];

plot(X, Y, 'o-','Color',color);
% -- plot mean
Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(X, Ymean, {'Errors', Ysem});
% -- stats
[~, p]=ttest(Y(:,1), Y(:,2)); % paired diff
disp(['1 vs 2; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,1), Y(:,3)); % paired diff
disp(['1 vs 3; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,2), Y(:,3)); % paired diff
disp(['2 vs 3; p=' num2str(p, '%3.2g')]);

[~, p]=ttest(Y(:,1)); % diff from 0
if p<0.1; 
    lt_plot_text(x-0.2, 1.1*max(Y(:,1)), ['p=' num2str(p, '%3.2g')]);
end
[~, p]=ttest(Y(:,2)); % diff from 0
if p<0.1; 
    lt_plot_text(x, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
 [~, p]=ttest(Y(:,3)); % diff from 0
if p<0.1; 
    lt_plot_text(x+0.2, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
   

% === DIFF
x=2;
color='r';
inds=TargStatus_all==0 & SimDiff_all==0 & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);
GenLearn=Generalization_Learn_all(inds);


X=[x-0.2 x x+0.2];
Y=[GenAFP' GenMP' GenLearn'];

plot(X, Y, 'o-','Color',color);
% -- plot mean
Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(X, Ymean, {'Errors', Ysem});
% -- stats
[~, p]=ttest(Y(:,1), Y(:,2)); % paired diff
disp(['1 vs 2; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,1), Y(:,3)); % paired diff
disp(['1 vs 3; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,2), Y(:,3)); % paired diff
disp(['2 vs 3; p=' num2str(p, '%3.2g')]);

[~, p]=ttest(Y(:,1)); % diff from 0
if p<0.1; 
    lt_plot_text(x-0.2, 1.1*max(Y(:,1)), ['p=' num2str(p, '%3.2g')]);
end
[~, p]=ttest(Y(:,2)); % diff from 0
if p<0.1; 
    lt_plot_text(x, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
 [~, p]=ttest(Y(:,3)); % diff from 0
if p<0.1; 
    lt_plot_text(x+0.2, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end



% == ALL
x=3;
color='k';
inds=TargStatus_all==0 & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);
GenLearn=Generalization_Learn_all(inds);


X=[x-0.2 x x+0.2];
Y=[GenAFP' GenMP' GenLearn'];

plot(X, Y, 'o-','Color',color);
% -- plot mean
Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(X, Ymean, {'Errors', Ysem});
% -- stats
[~, p]=ttest(Y(:,1), Y(:,2)); % paired diff
disp(['1 vs 2; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,1), Y(:,3)); % paired diff
disp(['1 vs 3; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,2), Y(:,3)); % paired diff
disp(['2 vs 3; p=' num2str(p, '%3.2g')]);

[~, p]=ttest(Y(:,1)); % diff from 0
if p<0.1; 
    lt_plot_text(x-0.2, 1.1*max(Y(:,1)), ['p=' num2str(p, '%3.2g')]);
end
[~, p]=ttest(Y(:,2)); % diff from 0
if p<0.1; 
    lt_plot_text(x, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
 [~, p]=ttest(Y(:,3)); % diff from 0
if p<0.1; 
    lt_plot_text(x+0.2, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end



%% ====== [NOT ADJACENT and ADJACENT]

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
disp('GENERALIZATION (afp -- mp -- learning) [ALL]');
title('GENERALIZATION (afp -- mp -- learning) [ALL]');


% ====  SIMILAR
x=1;
color='b';
inds=TargStatus_all==0 & SimDiff_all==1;

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);
GenLearn=Generalization_Learn_all(inds);


X=[x-0.2 x x+0.2];
Y=[GenAFP' GenMP' GenLearn'];

plot(X, Y, 'o-','Color',color);
% -- plot mean
Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(X, Ymean, {'Errors', Ysem});
% -- stats
[~, p]=ttest(Y(:,1), Y(:,2)); % paired diff
disp(['1 vs 2; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,1), Y(:,3)); % paired diff
disp(['1 vs 3; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,2), Y(:,3)); % paired diff
disp(['2 vs 3; p=' num2str(p, '%3.2g')]);

[~, p]=ttest(Y(:,1)); % diff from 0
if p<0.1; 
    lt_plot_text(x-0.2, 1.1*max(Y(:,1)), ['p=' num2str(p, '%3.2g')]);
end
[~, p]=ttest(Y(:,2)); % diff from 0
if p<0.1; 
    lt_plot_text(x, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
 [~, p]=ttest(Y(:,3)); % diff from 0
if p<0.1; 
    lt_plot_text(x+0.2, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
   

% === DIFF
x=2;
color='r';
inds=TargStatus_all==0 & SimDiff_all==0;

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);
GenLearn=Generalization_Learn_all(inds);


X=[x-0.2 x x+0.2];
Y=[GenAFP' GenMP' GenLearn'];

plot(X, Y, 'o-','Color',color);
% -- plot mean
Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(X, Ymean, {'Errors', Ysem});
% -- stats
[~, p]=ttest(Y(:,1), Y(:,2)); % paired diff
disp(['1 vs 2; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,1), Y(:,3)); % paired diff
disp(['1 vs 3; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,2), Y(:,3)); % paired diff
disp(['2 vs 3; p=' num2str(p, '%3.2g')]);

[~, p]=ttest(Y(:,1)); % diff from 0
if p<0.1; 
    lt_plot_text(x-0.2, 1.1*max(Y(:,1)), ['p=' num2str(p, '%3.2g')]);
end
[~, p]=ttest(Y(:,2)); % diff from 0
if p<0.1; 
    lt_plot_text(x, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
 [~, p]=ttest(Y(:,3)); % diff from 0
if p<0.1; 
    lt_plot_text(x+0.2, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end



% == ALL
x=3;
color='k';
inds=TargStatus_all==0;

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);
GenLearn=Generalization_Learn_all(inds);


X=[x-0.2 x x+0.2];
Y=[GenAFP' GenMP' GenLearn'];

plot(X, Y, 'o-','Color',color);
% -- plot mean
Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(X, Ymean, {'Errors', Ysem});
% -- stats
[~, p]=ttest(Y(:,1), Y(:,2)); % paired diff
disp(['1 vs 2; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,1), Y(:,3)); % paired diff
disp(['1 vs 3; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,2), Y(:,3)); % paired diff
disp(['2 vs 3; p=' num2str(p, '%3.2g')]);

[~, p]=ttest(Y(:,1)); % diff from 0
if p<0.1; 
    lt_plot_text(x-0.2, 1.1*max(Y(:,1)), ['p=' num2str(p, '%3.2g')]);
end
[~, p]=ttest(Y(:,2)); % diff from 0
if p<0.1; 
    lt_plot_text(x, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
 [~, p]=ttest(Y(:,3)); % diff from 0
if p<0.1; 
    lt_plot_text(x+0.2, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end




%% ====== [ADJACENT]

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');


% ====  SIMILAR
x=1;
color='b';
inds=TargStatus_all==0 & SimDiff_all==1 & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);
GenLearn=Generalization_Learn_all(inds);


X=[x-0.2 x x+0.2];
Y=[GenAFP' GenMP' GenLearn'];

plot(X, Y, 'o-','Color',color);
% -- plot mean
Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(X, Ymean, {'Errors', Ysem});
% -- stats
[~, p]=ttest(Y(:,1), Y(:,2)); % paired diff
disp(['1 vs 2; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,1), Y(:,3)); % paired diff
disp(['1 vs 3; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,2), Y(:,3)); % paired diff
disp(['2 vs 3; p=' num2str(p, '%3.2g')]);

[~, p]=ttest(Y(:,1)); % diff from 0
if p<0.1; 
    lt_plot_text(x-0.2, 1.1*max(Y(:,1)), ['p=' num2str(p, '%3.2g')]);
end
[~, p]=ttest(Y(:,2)); % diff from 0
if p<0.1; 
    lt_plot_text(x, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
 [~, p]=ttest(Y(:,3)); % diff from 0
if p<0.1; 
    lt_plot_text(x+0.2, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
   

% === DIFF
x=2;
color='r';
inds=TargStatus_all==0 & SimDiff_all==0 & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);
GenLearn=Generalization_Learn_all(inds);


X=[x-0.2 x x+0.2];
Y=[GenAFP' GenMP' GenLearn'];

plot(X, Y, 'o-','Color',color);
% -- plot mean
Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(X, Ymean, {'Errors', Ysem});
% -- stats
[~, p]=ttest(Y(:,1), Y(:,2)); % paired diff
disp(['1 vs 2; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,1), Y(:,3)); % paired diff
disp(['1 vs 3; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,2), Y(:,3)); % paired diff
disp(['2 vs 3; p=' num2str(p, '%3.2g')]);

[~, p]=ttest(Y(:,1)); % diff from 0
if p<0.1; 
    lt_plot_text(x-0.2, 1.1*max(Y(:,1)), ['p=' num2str(p, '%3.2g')]);
end
[~, p]=ttest(Y(:,2)); % diff from 0
if p<0.1; 
    lt_plot_text(x, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
 [~, p]=ttest(Y(:,3)); % diff from 0
if p<0.1; 
    lt_plot_text(x+0.2, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end



% == ALL
x=3;
color='k';
inds=TargStatus_all==0 & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);
GenLearn=Generalization_Learn_all(inds);


X=[x-0.2 x x+0.2];
Y=[GenAFP' GenMP' GenLearn'];

plot(X, Y, 'o-','Color',color);
% -- plot mean
Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(X, Ymean, {'Errors', Ysem});
% -- stats
[~, p]=ttest(Y(:,1), Y(:,2)); % paired diff
disp(['1 vs 2; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,1), Y(:,3)); % paired diff
disp(['1 vs 3; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,2), Y(:,3)); % paired diff
disp(['2 vs 3; p=' num2str(p, '%3.2g')]);

[~, p]=ttest(Y(:,1)); % diff from 0
if p<0.1; 
    lt_plot_text(x-0.2, 1.1*max(Y(:,1)), ['p=' num2str(p, '%3.2g')]);
end
[~, p]=ttest(Y(:,2)); % diff from 0
if p<0.1; 
    lt_plot_text(x, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
 [~, p]=ttest(Y(:,3)); % diff from 0
if p<0.1; 
    lt_plot_text(x+0.2, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end


%% ====== [NOT ADJACENT]

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
disp('GENERALIZATION (afp -- mp -- learning) [NOT ADJACENT]');
title('GENERALIZATION (afp -- mp -- learning) [NOT ADJACENT]');


% ====  SIMILAR
x=1;
color='b';
inds=TargStatus_all==0 & SimDiff_all==1 & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);
GenLearn=Generalization_Learn_all(inds);


X=[x-0.2 x x+0.2];
Y=[GenAFP' GenMP' GenLearn'];

% -- plot mean
Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(X, Ymean, {'Errors', Ysem});
% -- stats
[~, p]=ttest(Y(:,1), Y(:,2)); % paired diff
disp(['1 vs 2; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,1), Y(:,3)); % paired diff
disp(['1 vs 3; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,2), Y(:,3)); % paired diff
disp(['2 vs 3; p=' num2str(p, '%3.2g')]);

[~, p]=ttest(Y(:,1)); % diff from 0
if p<0.1; 
    lt_plot_text(x-0.2, 1.1*max(Y(:,1)), ['p=' num2str(p, '%3.2g')]);
end
[~, p]=ttest(Y(:,2)); % diff from 0
if p<0.1; 
    lt_plot_text(x, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
 [~, p]=ttest(Y(:,3)); % diff from 0
if p<0.1; 
    lt_plot_text(x+0.2, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
   

% === DIFF
x=2;
color='r';
inds=TargStatus_all==0 & SimDiff_all==0 & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);
GenLearn=Generalization_Learn_all(inds);


X=[x-0.2 x x+0.2];
Y=[GenAFP' GenMP' GenLearn'];

% -- plot mean
Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(X, Ymean, {'Errors', Ysem});
% -- stats
[~, p]=ttest(Y(:,1), Y(:,2)); % paired diff
disp(['1 vs 2; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,1), Y(:,3)); % paired diff
disp(['1 vs 3; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,2), Y(:,3)); % paired diff
disp(['2 vs 3; p=' num2str(p, '%3.2g')]);

[~, p]=ttest(Y(:,1)); % diff from 0
if p<0.1; 
    lt_plot_text(x-0.2, 1.1*max(Y(:,1)), ['p=' num2str(p, '%3.2g')]);
end
[~, p]=ttest(Y(:,2)); % diff from 0
if p<0.1; 
    lt_plot_text(x, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
 [~, p]=ttest(Y(:,3)); % diff from 0
if p<0.1; 
    lt_plot_text(x+0.2, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end



% == ALL
x=3;
color='k';
inds=TargStatus_all==0 & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);
GenLearn=Generalization_Learn_all(inds);


X=[x-0.2 x x+0.2];
Y=[GenAFP' GenMP' GenLearn'];

% -- plot mean
Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(X, Ymean, {'Errors', Ysem});
% -- stats
[~, p]=ttest(Y(:,1), Y(:,2)); % paired diff
disp(['1 vs 2; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,1), Y(:,3)); % paired diff
disp(['1 vs 3; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,2), Y(:,3)); % paired diff
disp(['2 vs 3; p=' num2str(p, '%3.2g')]);

[~, p]=ttest(Y(:,1)); % diff from 0
if p<0.1; 
    lt_plot_text(x-0.2, 1.1*max(Y(:,1)), ['p=' num2str(p, '%3.2g')]);
end
[~, p]=ttest(Y(:,2)); % diff from 0
if p<0.1; 
    lt_plot_text(x, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
 [~, p]=ttest(Y(:,3)); % diff from 0
if p<0.1; 
    lt_plot_text(x+0.2, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end




%% ====== [ADJACENT]

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');


% ====  SIMILAR
x=1;
color='b';
inds=TargStatus_all==0 & SimDiff_all==1 & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);
GenLearn=Generalization_Learn_all(inds);


X=[x-0.2 x x+0.2];
Y=[GenAFP' GenMP' GenLearn'];

% -- plot mean
Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(X, Ymean, {'Errors', Ysem});
% -- stats
[~, p]=ttest(Y(:,1), Y(:,2)); % paired diff
disp(['1 vs 2; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,1), Y(:,3)); % paired diff
disp(['1 vs 3; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,2), Y(:,3)); % paired diff
disp(['2 vs 3; p=' num2str(p, '%3.2g')]);

[~, p]=ttest(Y(:,1)); % diff from 0
if p<0.1; 
    lt_plot_text(x-0.2, 1.1*max(Y(:,1)), ['p=' num2str(p, '%3.2g')]);
end
[~, p]=ttest(Y(:,2)); % diff from 0
if p<0.1; 
    lt_plot_text(x, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
 [~, p]=ttest(Y(:,3)); % diff from 0
if p<0.1; 
    lt_plot_text(x+0.2, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
   

% === DIFF
x=2;
color='r';
inds=TargStatus_all==0 & SimDiff_all==0 & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);
GenLearn=Generalization_Learn_all(inds);


X=[x-0.2 x x+0.2];
Y=[GenAFP' GenMP' GenLearn'];

% -- plot mean
Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(X, Ymean, {'Errors', Ysem});
% -- stats
[~, p]=ttest(Y(:,1), Y(:,2)); % paired diff
disp(['1 vs 2; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,1), Y(:,3)); % paired diff
disp(['1 vs 3; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,2), Y(:,3)); % paired diff
disp(['2 vs 3; p=' num2str(p, '%3.2g')]);

[~, p]=ttest(Y(:,1)); % diff from 0
if p<0.1; 
    lt_plot_text(x-0.2, 1.1*max(Y(:,1)), ['p=' num2str(p, '%3.2g')]);
end
[~, p]=ttest(Y(:,2)); % diff from 0
if p<0.1; 
    lt_plot_text(x, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
 [~, p]=ttest(Y(:,3)); % diff from 0
if p<0.1; 
    lt_plot_text(x+0.2, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end



% == ALL
x=3;
color='k';
inds=TargStatus_all==0 & ~(Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);
GenLearn=Generalization_Learn_all(inds);


X=[x-0.2 x x+0.2];
Y=[GenAFP' GenMP' GenLearn'];

% -- plot mean
Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(X, Ymean, {'Errors', Ysem});
% -- stats
[~, p]=ttest(Y(:,1), Y(:,2)); % paired diff
disp(['1 vs 2; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,1), Y(:,3)); % paired diff
disp(['1 vs 3; p=' num2str(p, '%3.2g')]);
[~, p]=ttest(Y(:,2), Y(:,3)); % paired diff
disp(['2 vs 3; p=' num2str(p, '%3.2g')]);

[~, p]=ttest(Y(:,1)); % diff from 0
if p<0.1; 
    lt_plot_text(x-0.2, 1.1*max(Y(:,1)), ['p=' num2str(p, '%3.2g')]);
end
[~, p]=ttest(Y(:,2)); % diff from 0
if p<0.1; 
    lt_plot_text(x, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end
 [~, p]=ttest(Y(:,3)); % diff from 0
if p<0.1; 
    lt_plot_text(x+0.2, 1.1*max(Y(:,2)), ['p=' num2str(p, '%3.2g')]);
end

%% ===== SCATTER [AFP VS. MP]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');

title('AFP vs. MP generalization [not adjacent]')
ylabel('AFP generalization');
xlabel('MP generalization');
% ====  SIMILAR
color='b';
inds=TargStatus_all==0 & SimDiff_all==1 & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);

lt_plot_45degScatter(GenMP, GenAFP, color)


% ====  DIFF
color='r';
inds=TargStatus_all==0 & SimDiff_all==0 & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

GenAFP=Generalization_AFP_all(inds);
GenMP=Generalization_MP_all(inds);

lt_plot_45degScatter(GenMP, GenAFP, color)



xlim([-2 2]);
ylim([-2 2]);


%% ===== SCATTER [AFP VS. LEARNING] [NOT ADJACENT]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');

title('AFP vs. LEARNING generalization [not adjacent]')
ylabel('AFP generalization');
xlabel('LEARNING generalization');
% ====  SIMILAR
color='b';
inds=TargStatus_all==0 & SimDiff_all==1 & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

GenAFP=Generalization_AFP_all(inds);
GenLearn=Generalization_Learn_all(inds);

lt_plot_45degScatter(GenLearn, GenAFP, color)


% ====  DIFF
color='r';
inds=TargStatus_all==0 & SimDiff_all==0 & (Y_PosRelTarg_All~=-1 & Y_PosRelTarg_All~=1);

GenAFP=Generalization_AFP_all(inds);
GenLearn=Generalization_Learn_all(inds);

lt_plot_45degScatter(GenLearn, GenAFP, color)



xlim([-2 2]);
ylim([-2 2]);



%% ===== SCATTER [AFP VS. MP] [targ]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
hsplots=[];
hsplots=[hsplots hsplot];

title('targ [all]')
ylabel('AFP bias');
xlabel('MP bias');

color='k';
inds=TargStatus_all==1;

X=MPbias_all(inds);
Y=AFPbias_all(inds);

lt_plot_45degScatter(X, Y, color)


%% ===== SCATTER [AFP VS. MP] [SAME]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
hsplots=[hsplots hsplot];

title('same [all]')
ylabel('AFP bias');
xlabel('MP bias');

color='b';
inds=TargStatus_all==0 & SimDiff_all==1;

X=MPbias_all(inds);
Y=AFPbias_all(inds);

lt_plot_45degScatter(X, Y, color)

%% ===== SCATTER [AFP VS. MP] [DIFF]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
hsplots=[hsplots hsplot];

title('Diff [all]')
ylabel('AFP bias');
xlabel('MP bias');

color='r';
inds=TargStatus_all==0 & SimDiff_all==0;

X=MPbias_all(inds);
Y=AFPbias_all(inds);

lt_plot_45degScatter(X, Y, color)

linkaxes(hsplots, 'xy');


%% +++++++++++++++++++++++++++++++++++++++++ AFP BIAS VS. LEARNING
%% ===== SCATTER [AFP VS. MP] [targ]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
hsplots=[];
hsplots=[hsplots hsplot];

title('targ [all]')
ylabel('AFP bias');
xlabel('pitch shift');

color='k';
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

lt_plot_45degScatter(X, Y, color)


%% ===== SCATTER [AFP VS. MP] [SAME]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
hsplots=[hsplots hsplot];

title('same [all]')
ylabel('AFP bias');
xlabel('pitch shift');

color='b';
inds=TargStatus_all==0 & SimDiff_all==1;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

lt_plot_45degScatter(X, Y, color)

%% ===== SCATTER [AFP VS. MP] [DIFF]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
hsplots=[hsplots hsplot];

title('Diff [all]')
ylabel('AFP bias');
xlabel('pitch shift');

color='r';
inds=TargStatus_all==0 & SimDiff_all==0;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

lt_plot_45degScatter(X, Y, color)

linkaxes(hsplots, 'xy');


%% +++++++++++++++++++++++++++++++++++++++++ MP BIAS VS. LEARNING
%% ===== SCATTER [AFP VS. MP] [targ]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
hsplots=[];
hsplots=[hsplots hsplot];

title('targ [all]')
ylabel('MP bias');
xlabel('pitch shift');

color='k';
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot_45degScatter(X, Y, color)

%% ===== SCATTER [AFP VS. MP] [SAME]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
hsplots=[hsplots hsplot];

title('same [all]')
ylabel('MP bias');
xlabel('pitch shift');

color='b';
inds=TargStatus_all==0 & SimDiff_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot_45degScatter(X, Y, color)

%% ===== SCATTER [AFP VS. MP] [DIFF]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
hsplots=[hsplots hsplot];

title('Diff [all]')
ylabel('MP bias');
xlabel('pitch shift');

color='r';
inds=TargStatus_all==0 & SimDiff_all==0;
X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_plot_45degScatter(X, Y, color)

linkaxes(hsplots, 'xy');


%% +++++++++++++++++++++++++++++++++++++++++ MP BIAS VS. LEARNING [regression
%% ===== SCATTER [AFP VS. MP] [targ]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
hsplots=[];
hsplots=[hsplots hsplot];

title('targ [all]')
ylabel('MP bias');
xlabel('pitch shift');

color='k';
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, color)
% line(xlim, [0 0], 'Color','k','LineStyle','--');
% line([0 0], ylim, 'Color','k','LineStyle','--')

%% ===== SCATTER [AFP VS. MP] [SAME]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
hsplots=[hsplots hsplot];

title('same [all]')
ylabel('MP bias');
xlabel('pitch shift');

color='b';
inds=TargStatus_all==0 & SimDiff_all==1;

X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, color)

%% ===== SCATTER [AFP VS. MP] [DIFF]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
hsplots=[hsplots hsplot];

title('Diff [all]')
ylabel('MP bias');
xlabel('pitch shift');

color='r';
inds=TargStatus_all==0 & SimDiff_all==0;
X=LearningPBS_all(inds);
Y=MPbias_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, color)

linkaxes(hsplots, 'xy');


%% +++++++++++++++++++++++++++++++++++++++++ AFP BIAS VS. LEARNING [regression
%% ===== SCATTER [AFP VS. MP] [targ]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
hsplots=[hsplots hsplot];

title('targ [all]')
ylabel('AFP bias');
xlabel('pitch shift');

color='k';
inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, color)
% line(xlim, [0 0], 'Color','k','LineStyle','--');
% line([0 0], ylim, 'Color','k','LineStyle','--')

%% ===== SCATTER [AFP VS. MP] [SAME]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
hsplots=[hsplots hsplot];

title('same [all]')
ylabel('AFP bias');
xlabel('pitch shift');

color='b';
inds=TargStatus_all==0 & SimDiff_all==1;

X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, color)

%% ===== SCATTER [AFP VS. MP] [DIFF]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% disp('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
% title('GENERALIZATION (afp -- mp -- learning) [ADJACENT]');
hsplots=[hsplots hsplot];

title('Diff [all]')
ylabel('AFP bias');
xlabel('pitch shift');

color='r';
inds=TargStatus_all==0 & SimDiff_all==0;
X=LearningPBS_all(inds);
Y=AFPbias_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, color)

linkaxes(hsplots, 'xy');



%% ======================== USE ANCOVA TO ASK WHETHER AFP BIAS IS DIFFERENT (AFTER TAKING INTO ACCOUNT PREDICTED CORR BETWEEN AFP BIAS AND LEARNING)
% ALSO ASK WHETHER SLOPES ARE DIFFERENT, ETC
keyboard;

inds=TargStatus_all==1 | (TargStatus_all==0 & SimDiff_all==1); % both targets and same types


X=LearningPBS_all(inds);
Xname='learning';
Y=-AFPbias_all(inds);
Yname='reversion (hz)';

Zgroups=TargStatus_all(inds)==0 & SimDiff_all(inds)==1;
% Zgroups(Zgroups~=1)=100;
Zname='nontarget?';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname, 'on', 'separate lines');
[c, m, h, nms]=multcompare(stats, 'alpha', 0.05, 'estimate', 'slope');

% [h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname, 'on', 'parallel lines');
% [c, m, h, nms]=multcompare(stats, 'estimate', 'pmm');

% ==== PLOT
lt_figure; hold on;
xlabel('laerning');
ylabel('reversion');
title('targ + same-type');

inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Y=-AFPbias_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'k')

inds=TargStatus_all==0 & SimDiff_all==1;

X=LearningPBS_all(inds);
Y=-AFPbias_all(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b')

%% [ALL NONTARGS] ======================== USE ANCOVA TO ASK WHETHER AFP BIAS IS DIFFERENT (AFTER TAKING INTO ACCOUNT PREDICTED CORR BETWEEN AFP BIAS AND LEARNING)
% ALSO ASK WHETHER SLOPES ARE DIFFERENT, ETC
keyboard;

% ===== correlation between sequential position (collapse into +1 and others) and correlations?
inds=TargStatus_all==1 | (TargStatus_all==0); % both targets and same types


X=LearningPBS_all(inds);
Xname='learning';
Y=AFPbias_all(inds);
Yname='afp bias';

Zgroups=TargStatus_all(inds)==0;
% Zgroups(Zgroups~=1)=100;
Zname='nontarget?';

[h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname, 'on', 'separate lines');
[c, m, h, nms]=multcompare(stats, 'alpha', 0.05, 'estimate', 'slope');

% [h, a, c, stats]=aoctool(X, Y, Zgroups, '', Xname, Yname, Zname, 'on', 'parallel lines');
% [c, m, h, nms]=multcompare(stats, 'estimate', 'pmm');

% ==== PLOT
lt_figure; hold on;
xlabel('laerning');
ylabel('afp bias');
title('targ + allnontargs');

inds=TargStatus_all==1;

X=LearningPBS_all(inds);
Xname='learning';
Y=AFPbias_all(inds);
Yname='afp bias';

lt_regress(Y, X, 1, 0, 1, 1, 'k')

inds=TargStatus_all==0;

X=LearningPBS_all(inds);
Xname='learning';
Y=AFPbias_all(inds);
Yname='afp bias';

lt_regress(Y, X, 1, 0, 1, 1, 'b')

