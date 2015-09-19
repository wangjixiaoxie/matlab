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
        
        if numtargs>1;
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
Target_all=[];
AcousticDist_all=[];
BirdNum_all=[];
ExptNum_all=[];
Syl_all={};
SylSingle_all={};
OverallExptNum_all=[];

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
            
            % --
            Similar_all = [Similar_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
            
% --
            Target_all = [Target_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
            
                % ---
                AcousticDist_all = [AcousticDist_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];
                    
                % -- 
                BirdNum_all=[BirdNum_all i];
                
                % -- 
                ExptNum_all=[ExptNum_all ii];
                
                %
                Syl_all=[Syl_all syl];
                
                % single syl
                SylSingle_all=[SylSingle_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
               
                % 
                OverallExptNum_all=[OverallExptNum_all count];
            
                
                
        end
    end
end

% ----- OUTPUT
SeqDepPitch_AcrossBirds.AllSyllables.Learning_all=Learning_all;
SeqDepPitch_AcrossBirds.AllSyllables.Learning_TargDir_All=Learning_TargDir_All;
SeqDepPitch_AcrossBirds.AllSyllables.LearningRelTarg_all=LearningRelTarg_all;
SeqDepPitch_AcrossBirds.AllSyllables.Similar_all=Similar_all;
SeqDepPitch_AcrossBirds.AllSyllables.Target_all=Target_all;
SeqDepPitch_AcrossBirds.AllSyllables.AcousticDist_all=AcousticDist_all;
SeqDepPitch_AcrossBirds.AllSyllables.BirdNum_all=BirdNum_all;
SeqDepPitch_AcrossBirds.AllSyllables.ExptNum_all=ExptNum_all;
SeqDepPitch_AcrossBirds.AllSyllables.Syl_all=Syl_all;
SeqDepPitch_AcrossBirds.AllSyllables.SylSingle_all=SylSingle_all;
SeqDepPitch_AcrossBirds.AllSyllables.OverallExptNum_all=OverallExptNum_all;

%% HIT RATE

[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_HITRATE(SeqDepPitch_AcrossBirds, PARAMS);

pause;


%% LEARNING TAKING INTO ACCOUNT DRIFT/BASELINE

% ==== 1) PLOT LEARNING IN COMPARISON TO BASELINE IN VARIOUS WAYS, AND 2) BASELINE DRIFT ANALYSIS (SIGNIFICANCE TESTS, MORE COMPLICATED
% MODELS, ETC)
[SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_DRIFT(SeqDepPitch_AcrossBirds, PARAMS);

% ===== 2) DRIFT : TAKE EQUAL NUMBER OF BASELINE DAYS AS NUMBER OF DAYS
% DURING WN
[SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_DRIFT2(SeqDepPitch_AcrossBirds, PARAMS);




%% PLOTS OF LEARNING BASED ON MOTIF DISTANCE etc (USING EARLY CONSOL FF)

% [IGNORE!! does not use learning metric];
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_MOTIFS(SeqDepPitch_AcrossBirds, PARAMS);
pause;

% PLOTS BASED ON MOTIF DISTANCE etc (USING LEARNING METRIC)
use_learning_metric=1;
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_MOTIFS(SeqDepPitch_AcrossBirds, PARAMS, use_learning_metric);
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
PARAMS.FilterPlot.sylID_filter{1}='presyl_similar_to_targ_presyl';
PARAMS.FilterPlot.sylID_filter{2}={0,1}; %

% filter based on post syllable similarity to target post
% PARAMS.FilterPlot.sylID_filter{1}='postsyl_similar_to_targ_postsyl';
% PARAMS.FilterPlot.sylID_filter{2}={0,1}; %

% ==== ADDITIONAL FILTER
PARAMS.FilterPlot.extra_filter=1; % if 1, then filter similar vs diff
PARAMS.FilterPlot.similar_only=0; % if 1, takes only similar, if 0, takes only different; only applies if extra_filter is on. 

% ===== What days to plot for across days analysis?
PARAMS.FilterPlot.FirstDayToPlot=1 ; % lock start to what WN day? (1st WN day would be 1);
PARAMS.FilterPlot.LastDayToPlot=[]; % lock end to what WN day? (5 would be 5th WN day, 'end' would be end of WN, [] would be openended);


[FILTERED_DATA, SeqDepPitch_AcrossBirds]=lt_seq_dep_pitch_ACROSSBIRDS_FilterPlot(PARAMS, SeqDepPitch_AcrossBirds);

    


