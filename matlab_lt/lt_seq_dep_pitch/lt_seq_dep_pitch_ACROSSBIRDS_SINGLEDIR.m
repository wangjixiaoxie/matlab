function [SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_SINGLEDIR(SeqDepPitch_AcrossBirds, PARAMS)
%% LT 7/21/15 - Plots single dir learning experiments


%% START ANALYSIS ON SUBSET OF EXPERIMENTS (e.g. those with only one target)
%% SORT OUT ONLY THE EXPERIMENTS WITH ONLY ONE TARG, OR STARTED OUT WITH ONLY ONE TARG


NumBirds=length(SeqDepPitch_AcrossBirds.birds);

% copy strcuture, save backup.
if ~exist('SeqDepPitch_AcrossBirds_ORIGINAL', 'var'); % do not save a new backup if it already exists.
    SeqDepPitch_AcrossBirds_ORIGINAL=SeqDepPitch_AcrossBirds;
end

% remove experiments - first figure out what inds to remove
expts_to_remove=[];
for i=1:NumBirds;
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    expts_to_remove{i}=[];
    for ii=1:NumExperiments;
        numtargs=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs;
        
        if numtargs==2 | numtargs==9;
            % mult targets. note down index and remove later
            expts_to_remove{i}=[expts_to_remove{i} ii];
        end
    end
end

% actually remove them
for i=1:length(expts_to_remove);
    if ~isempty(expts_to_remove{i});
        SeqDepPitch_AcrossBirds.birds{i}.experiment(expts_to_remove{i})=[];
        disp(['TO GET ONLY ONE TARGETS: removed: bird: ' num2str(i) '; expt: ' num2str(expts_to_remove{i})]);
    end
end


PARAMS.global.did_take_only_onetarg_expts=1;
NumBirds=length(SeqDepPitch_AcrossBirds.birds);


%% MAKE A NEW STRUCTURE - each datapoint is one syl (across all birds and experiments)
try
SeqDepPitch_AcrossBirds=rmfield(SeqDepPitch_AcrossBirds, 'AllSyllables');
catch err
end

Learning_all=[];
Learning_TargDir_All=[];
LearningRelTarg_all=[];
Similar_all=[];
PreSimilar_all=[];
TwoBackSimilar_all=[];
Target_all=[];
AcousticDist_all=[];
Corr_motif_all=[];
Corr_motifsubsong_all=[];
Corr_song_all=[];
BirdNum_all=[];
ExptNum_all=[];
Syl_all={};
SylSingle_all={};
OverallExptNum_all=[];
DayVals_DurWN_zscore={};
DayVals_Baseline_zscore={};
DayVals_DurWN_zscore_targdirsign={};
DayVals_Baseline_zscore_targdirsign={};

baselineDriftTargDir_All=[];
baselineDrift_RelTarg_All=[];

AdjacentToTarg_All=[];


count=0;

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
                
        
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;

        count=count+1;
        for j=1:length(syls_unique);
            syl=syls_unique{j};

            % ========== EXTRACT DATA
            
            % --- Learning
            Learning_all=[Learning_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean];
            
            % learning rel target learning dir
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            target_learning_sign=sign(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean);
            learning_targsylDir=target_learning_sign*SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            Learning_TargDir_All=[Learning_TargDir_All learning_targsylDir];
            
            
            % --- Learning rel targ
            LearningRelTarg_all=[LearningRelTarg_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ];
            
            if LearningRelTarg_all<-2
                keyboard
            end
            
            
            % ==== baseline drift
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl), 'DRIFT_UsingTwoDays');
                baselineDrift_z=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.zscore_mean;
                baselineDriftTargDir=target_learning_sign*baselineDrift_z;
                baselineDrift_RelTarg=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.zscore_rel_targ;
            else
                baselineDrift_z=nan;
                baselineDriftTargDir=nan;
                baselineDrift_RelTarg=nan;
            end
          
            baselineDriftTargDir_All=[baselineDriftTargDir_All baselineDriftTargDir];
            baselineDrift_RelTarg_All=[baselineDrift_RelTarg_All baselineDrift_RelTarg];
    
            
            % --
            Similar_all = [Similar_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
            
            % --
            Target_all = [Target_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
            
            % ---
            AcousticDist_all = [AcousticDist_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];
            
            % --- adjacent to target (pre or post)
            PosRelTarg=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).distance_from_targ;
            if PosRelTarg==-1 | PosRelTarg==1;
                AdjacentToTarg=1;
            else
                AdjacentToTarg=0;
            end
            AdjacentToTarg_All=[AdjacentToTarg_All AdjacentToTarg];
            
            
            try
            if ~isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).distance_from_targ)
                Corr_motif_all=[Corr_motif_all ...
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(targsyl)];
                
                Corr_motifsubsong_all=[Corr_motifsubsong_all ...
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED.corrcoeff_vs.(targsyl)];
            else
                Corr_motif_all=[Corr_motif_all nan];
                Corr_motifsubsong_all=[Corr_motifsubsong_all nan];
            end
            catch err
                Corr_motif_all=[Corr_motif_all nan];
                Corr_motifsubsong_all=[Corr_motifsubsong_all nan];
            end
            
            
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl), 'CORRELATIONS');
                corrsong=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl);
            else
            corrsong=nan;
            end
            
            Corr_song_all=[Corr_song_all corrsong];
            
            % --
            BirdNum_all=[BirdNum_all i];
            
            % --
                ExptNum_all=[ExptNum_all ii];
                
                %
                Syl_all=[Syl_all syl];
                
                % single syl
                SylSingle_all=[SylSingle_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
               
                
                % preceding syllable
                presimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
                PreSimilar_all=[PreSimilar_all presimilar];
                
                
                % two back
                TwoBackSim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back_same_as_targ;
                TwoBackSimilar_all=[TwoBackSimilar_all TwoBackSim];

                
                % 
                OverallExptNum_all=[OverallExptNum_all count];
            
                % day by day pitch (zscore) vals - during WN
                WNonInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
                WNoffInd=min([length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_zscore), ...
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd]); % min([length(data) WNoffInd]);
                dayvals1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_zscore(WNonInd:WNoffInd);
                
                DayVals_DurWN_zscore=[DayVals_DurWN_zscore dayvals1];
                
                % day by day pitch (zscore) vals - during baseline
                WNonInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
                dayvals2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_zscore(1:WNonInd-1);
                
                DayVals_Baseline_zscore=[DayVals_Baseline_zscore dayvals2];
                
                % day by day pitch (zscore) vals - during WN [learning dir
                % sign]
                DayVals_DurWN_zscore_targdirsign=[DayVals_DurWN_zscore_targdirsign target_learning_sign*dayvals1];
                
                % day by day pitch (zscore) vals - during baseline
                % [learning dir sign]
                DayVals_Baseline_zscore_targdirsign=[DayVals_Baseline_zscore_targdirsign target_learning_sign*dayvals2];
                
                
                
                % ++++++++++++++++ REPORT FEATURES THAT INTERESTED IN
                % 1) syls that are 3+ back
                if Target_all(end)==0 & TwoBackSim==1 & presimilar==1
                    disp(['3 back syl: ' syl '('  birdname '-' exptname '); targ: ' targsyl]);
                end
              
                
        end
    end
end

% ----- OUTPUT
SeqDepPitch_AcrossBirds.AllSyllables.Learning_all=Learning_all;
SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All=Learning_TargDir_All;
SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all=LearningRelTarg_all;
SeqDepPitch_AcrossBirds.AllSyllables.Similar_all=Similar_all;
SeqDepPitch_AcrossBirds.AllSyllables.PreSimilar_all=PreSimilar_all;
SeqDepPitch_AcrossBirds.AllSyllables.TwoBackSimilar_all=TwoBackSimilar_all;
SeqDepPitch_AcrossBirds.AllSyllables.Target_all=Target_all;
SeqDepPitch_AcrossBirds.AllSyllables.AcousticDist_all=AcousticDist_all;
SeqDepPitch_AcrossBirds.AllSyllables.BirdNum_all=BirdNum_all;
SeqDepPitch_AcrossBirds.AllSyllables.ExptNum_all=ExptNum_all;
SeqDepPitch_AcrossBirds.AllSyllables.Syl_all=Syl_all;
SeqDepPitch_AcrossBirds.AllSyllables.SylSingle_all=SylSingle_all;
SeqDepPitch_AcrossBirds.AllSyllables.OverallExptNum_all=OverallExptNum_all;
SeqDepPitch_AcrossBirds.AllSyllables.Corr_motif_all=Corr_motif_all;
SeqDepPitch_AcrossBirds.AllSyllables.Corr_motifsubsong_all=Corr_motifsubsong_all;
SeqDepPitch_AcrossBirds.AllSyllables.Corr_song_all=Corr_song_all;

SeqDepPitch_AcrossBirds.AllSyllables.AdjacentToTarg_All=AdjacentToTarg_All;

SeqDepPitch_AcrossBirds.AllSyllables.baselineDriftTargDir_Z=baselineDriftTargDir_All;
SeqDepPitch_AcrossBirds.AllSyllables.baselineDriftRelTarg=baselineDrift_RelTarg_All;



SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore=DayVals_DurWN_zscore;
SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore=DayVals_Baseline_zscore;
SeqDepPitch_AcrossBirds.AllSyllables.DayVals_DurWN_zscore_targdirsign=DayVals_DurWN_zscore_targdirsign;
SeqDepPitch_AcrossBirds.AllSyllables.DayVals_Baseline_zscore_targdirsign=DayVals_Baseline_zscore_targdirsign;






%% MAKE A NEW STRUCTURE - each datapoint is one experiment
count_expt=1; % experiments
try
    SeqDepPitch_AcrossBirds=rmfield(SeqDepPitch_AcrossBirds, 'AllExpt');
catch err
end

SeqDepPitch_AcrossBirds.AllExpt=[];

for i=1:NumBirds
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        % === GET INDS FOR THIS EXPT
        inds_expt=SeqDepPitch_AcrossBirds.AllSyllables.BirdNum_all==i & ...
            SeqDepPitch_AcrossBirds.AllSyllables.ExptNum_all==ii;
        
        
        % ==== 1) EXTRACT DATA, SEPARATED INTO 1) TARGET, 2) NONTARGETS, 3)
        % SIMILAR, 4) DIFF;
        
        % --- TARGET
        inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==1 & inds_expt;
        DatStruct_tmp=lt_structure_subsample_all_fields(SeqDepPitch_AcrossBirds.AllSyllables, inds);
        
        SeqDepPitch_AcrossBirds.AllExpt(count_expt).Target=DatStruct_tmp;
        
        
        % --- NONTARGETS
        inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & inds_expt;
        DatStruct_tmp=lt_structure_subsample_all_fields(SeqDepPitch_AcrossBirds.AllSyllables, inds);
        
        SeqDepPitch_AcrossBirds.AllExpt(count_expt).Nontargets=DatStruct_tmp;
        
        % --- SIMILAR
        inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & inds_expt;
        DatStruct_tmp=lt_structure_subsample_all_fields(SeqDepPitch_AcrossBirds.AllSyllables, inds);
        
        SeqDepPitch_AcrossBirds.AllExpt(count_expt).Similar=DatStruct_tmp;
        
        % --- DIFFERENT
        inds=SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0 & SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==0 & inds_expt;
        DatStruct_tmp=lt_structure_subsample_all_fields(SeqDepPitch_AcrossBirds.AllSyllables, inds);
        
        SeqDepPitch_AcrossBirds.AllExpt(count_expt).Different=DatStruct_tmp;
        
        
        % ==== GET VECTOR OF ACOUSTIC DISTANCES FOR ALL SIMILAR AND DIFF
        % TYPE SYLS
        AcousticDist_vect_Similar=[];
        AcousticDist_vect_Diff=[];
        
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                continue;
            end
            ac_dist=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore;
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ==1;
                AcousticDist_vect_Similar=[AcousticDist_vect_Similar ac_dist];
            elseif SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ==0;
                AcousticDist_vect_Diff=[AcousticDist_vect_Diff ac_dist];
            end
        end
        
        % -- OUTPUT
        SeqDepPitch_AcrossBirds.AllExpt(count_expt).AcousticDist_vect_Similar=AcousticDist_vect_Similar;
        SeqDepPitch_AcrossBirds.AllExpt(count_expt).AcousticDist_vect_Diff=AcousticDist_vect_Diff;
        
        % ==== ASSIGN EACH SIMILAR TYPE SYL A NUMBER - how many diff types
        % are closer to target than you are?
        vect_of_acousticscores=[];
        for j=1:length(AcousticDist_vect_Similar);
            ac_dist_sim=AcousticDist_vect_Similar(j);
            
            vect_of_acousticscores=[vect_of_acousticscores sum(ac_dist_sim>AcousticDist_vect_Diff)];
        end
        % --- OUTPUT
        SeqDepPitch_AcrossBirds.AllExpt(count_expt).ForEachSimSyl_HowManyDiffTypeCloserToTarg=vect_of_acousticscores;
        
        
        % ==== ASSIGN EACH DIFF TYPE SYL A NUMBER - how many diff types
        % are closer to target than you are?
        vect_of_acousticscores=[];
        for j=1:length(AcousticDist_vect_Diff);
            ac_dist_diff=AcousticDist_vect_Diff(j);
            
            vect_of_acousticscores=[vect_of_acousticscores sum(ac_dist_diff>AcousticDist_vect_Similar)];
        end
        % --- OUTPUT
        SeqDepPitch_AcrossBirds.AllExpt(count_expt).ForEachDiffSyl_HowManySimTypeCloserToTarg=vect_of_acousticscores;
        
        
        
        
        % =================================== OTHER OUTPUTS (for this expt)
        SeqDepPitch_AcrossBirds.AllExpt(count_expt).birdnum=i;
        SeqDepPitch_AcrossBirds.AllExpt(count_expt).exptnum=ii;
        
        count_expt=count_expt+1;
        
        
        
    end
end
        
        

%% ==== manual control
keyboard



%% ==== linear model




%% ===== PLOT ForEachSimSyl_HowManyDiffTypeCloserToTarg;
lt_figure; hold on;
title('for each syl, count the number of opposite-type syls that are closer to the target than it is');
xlabel('number of opposite-type syls that are closer to the target');

Scores_all_sim=[SeqDepPitch_AcrossBirds.AllExpt.ForEachSimSyl_HowManyDiffTypeCloserToTarg];
Scores_all_diff=[SeqDepPitch_AcrossBirds.AllExpt.ForEachDiffSyl_HowManySimTypeCloserToTarg];
Xcenters=0:max([Scores_all_sim Scores_all_diff])+1;

% Each datapoint is similar type
lt_plot_histogram(Scores_all_sim, Xcenters, 1, 1, '',1);
    
% Each datapoint is diff type
[~, ~, hbar]=lt_plot_histogram(Scores_all_diff, Xcenters, 1, 1, '',1);
set(hbar, 'Color','r');


xlim([-0.5 max(Xcenters)+0.5]);


%% HIT RATE
close all;
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_HITRATE(SeqDepPitch_AcrossBirds, PARAMS);

pause;


%% [SEPARATION] CHANGE OVER TIME OF GENERALIZATION
close all;
OnlyAnalyzeIfHasAllDaysData=1;
ThrowOutIfDriveMore=1;
DaysToPlot=1:4;
ThrowOut_pu11_exptwithbadday7=0;
[SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_Separation(SeqDepPitch_AcrossBirds, PARAMS, OnlyAnalyzeIfHasAllDaysData, DaysToPlot,ThrowOutIfDriveMore,ThrowOut_pu11_exptwithbadday7);

%% LEARNING TAKING INTO ACCOUNT DRIFT/BASELINE

% ==== 1) PLOT LEARNING IN COMPARISON TO BASELINE IN VARIOUS WAYS, AND 2) BASELINE DRIFT ANALYSIS (SIGNIFICANCE TESTS, MORE COMPLICATED
% MODELS, ETC)
[SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_DRIFT(SeqDepPitch_AcrossBirds, PARAMS);

% ===== 2) DRIFT : TAKE EQUAL NUMBER OF BASELINE DAYS AS NUMBER OF DAYS
% DURING WN
close all;
plot_text=0;
[SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_DRIFT2(SeqDepPitch_AcrossBirds, PARAMS, plot_text);

% ===== 3) DRIFT : TAKE EQUAL NUMBER OF BASELINE DAYS AS NUMBER OF DAYS
% DURING WN [same as v2, but modified to take any input for which WN days
% to look at]
close all;
WNdaynum_start=4; % i.e. if WN starts on ind 8, then if this is 3, will starr count learning as day 10
WNdaynum_end=4; % if this is 4, then would stop counting learning on day 11. baseline would be last 2 baseline days

lt_seq_dep_pitch_ACROSSBIRDS_DRIFT3(SeqDepPitch_AcrossBirds, PARAMS, WNdaynum_start, WNdaynum_end);




%% Plot distributions of learning

close all;
lt_seq_dep_pitch_ACROSSBIRDS_LearningDistribution(SeqDepPitch_AcrossBirds, PARAMS)

%% Plot distributions v2 (take into account hand label vs. computer labeled syls)
close all
acoustThresh=PARAMS.SylClassify.SylDistCutoff; % for computer label.
lt_seq_dep_pitch_ACROSSBIRDS_LearDistribv2(SeqDepPitch_AcrossBirds, PARAMS, acoustThresh)


%% === plot distribution v3 [lines for expts, plot raw val]
% ---- separates by experiment. IDEALLY PLOTS ALL EXPTS, INCLUDING THOSE
% THAT DON'T PASS LEARNING THRESHOLD. 
close all
lt_seq_dep_pitch_ACROSSBIRDS_LearningDistribution3(SeqDepPitch_AcrossBirds, PARAMS)

%% Does effect of presyl hold even if account for acoustic similarity?
close all;
inds1=SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;

lt_seq_dep_pitch_ACROSSBIRDS_PresylAcoustic(SeqDepPitch_AcrossBirds, PARAMS, inds1)

%% Does effect of presyl hold even if account for Corr?
close all;
inds1=SeqDepPitch_AcrossBirds.AllSyllables.Similar_all==1 & SeqDepPitch_AcrossBirds.AllSyllables.Target_all==0;
lt_seq_dep_pitch_ACROSSBIRDS_PresylCorr(SeqDepPitch_AcrossBirds, PARAMS, inds1)



%% ===== ANTIGENERALiATION due to use-based learning?
close all;
lt_seq_dep_pitch_ACROSSBIRDS_Antigen(SeqDepPitch_AcrossBirds, PARAMS);


%% PLOTS OF LEARNING BASED ON MOTIF DISTANCE etc (USING EARLY CONSOL FF)

% [IGNORE!! does not use learning metric];
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_MOTIFS(SeqDepPitch_AcrossBirds, PARAMS);
pause;

% PLOTS BASED ON MOTIF DISTANCE etc (USING LEARNING METRIC)
close all;
use_learning_metric=1;
plot_hit_rate=0;
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_MOTIFS(SeqDepPitch_AcrossBirds, PARAMS, use_learning_metric, plot_hit_rate);
pause;

%% PLOT VARIOUST THINGS LIKE - corrleations, motifs, acoustic distance, learning.
% ======== Things like correlations, etc
close all;
% All Analyses in here [USING FFVALS CONSOL START]
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PLOT(SeqDepPitch_AcrossBirds, PARAMS);
pause;
close all;

% All Analyses in here[ USING LEARNING METRIC]
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PLOT_v2(SeqDepPitch_AcrossBirds, PARAMS);
pause;
close all;

% === DIRTY VERSION OF ABOVE THAT PLOTS CORR AND ACOUSTIC RELATIVE TO MOTIF
% (REL TARGET) [USING FFVALS CONSOL START]
% POSITION - must clean it up. (it hijacks the learning code)
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PLOT_DIRTY_reltarg(SeqDepPitch_AcrossBirds, PARAMS);
pause;

% === DIRTY VERSION OF ABOVE THAT PLOTS CORR AND ACOUSTIC RELATIVE TO MOTIF
% (ALL PAIRWISE) [USING FFVALS CONSOL START]
% POSITION - must clean it up. (it hijacks the learning code)
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PLOT_DIRTY_allpairs(SeqDepPitch_AcrossBirds, PARAMS);
pause;

%% PLOT LEARNING SORTING SYLLABLES (filter based on one category only)
% this works only for categorical things (including multicategory, like
% syllable distance)

% [ALL USING LEARNING METRIC, EXCEPT CONSOLIDATION STUFF];

close all;
% PARAMS:
% -- consolidation, what data to throw out?
PARAMS.FilterPlot.MinConsolidDur=5; % in days (inclusize);
PARAMS.FilterPlot.MaxDayFromWNStart=10; % i.e. number of days between WN starting and consolid start (e.g. WN on day 5, consolid day 7, is 3 days);


% ==== CHOOSE FILTER
% PARAMS.FilterPlot.sylID_filter{1}='similar_to_targ'; % i.e. the name of the category to filter by
% PARAMS.FilterPlot.sylID_filter{2}={0, 1}; % the possible classes of the category

% filter based on preceding syllable similarity to target
% PARAMS.FilterPlot.sylID_filter{1}='presyl_similar_to_targ_presyl';
% PARAMS.FilterPlot.sylID_filter{2}={0,1}; %

% filter based on post syllable similarity to target post
PARAMS.FilterPlot.sylID_filter{1}='postsyl_similar_to_targ_postsyl';
PARAMS.FilterPlot.sylID_filter{2}={0,1}; %

% ==== ADDITIONAL FILTER
PARAMS.FilterPlot.extra_filter=1; % if 1, then filter similar vs diff
PARAMS.FilterPlot.similar_only=1; % if 1, takes only similar, if 0, takes only different; only applies if extra_filter is on. 

% ===== What days to plot for across days analysis?
PARAMS.FilterPlot.FirstDayToPlot=1 ; % lock start to what WN day? (1st WN day would be 1);
PARAMS.FilterPlot.LastDayToPlot=[]; % lock end to what WN day? (5 would be 5th WN day, 'end' would be end of WN, [] would be openended);


[FILTERED_DATA, SeqDepPitch_AcrossBirds]=lt_seq_dep_pitch_ACROSSBIRDS_FilterPlot(PARAMS, SeqDepPitch_AcrossBirds);

    


