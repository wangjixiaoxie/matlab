%% LT 4/13/15 - Compiles data and plots across experiments/birds.  Need to specific dir of saved structures
%% TO DO:
% 1) remove repeat experiments - i.e. average over them.


%% INPUT PARAMS
clear all; close all;

PARAMS.global.use_zscore_learning=1; % not important. keep as 1.
PARAMS.global.ConsolBinDaySize=2;  % days at start and end of consolid
PARAMS.global.learning_metric='zscore';
% options: 'zscore' (zscore in epoch first pass threshold, defined
% below)
%               'ff_consolstart' (mean ff, at consol start, defined in cell)
PARAMS.global.remove_bad_syls=1; % e.g. those hard to quantify
PARAMS.global.SylsToRemove=...
    {'gr41gr90','SeqDepPitchShift',{'d'}, ...
    'gr41gr90','SeqDepPitchShift2',{'d'}, ...
    'gr41gr90','SeqDepPitchLMAN',{'d'}, ...
    'gr41gr90','SeqDepPitchLMAN2',{'d'}, ... % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
    'rd12pu6','SeqDepPitch',{'h'}}; % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});


%% INPUT DATA STRUCTURES
% from PlotLearning: AllDays_PlotLearning and Params

% Cell holding experiments {birdname, experient ID, NumTargets(e.g. 1, 2, 0(means 1 and more in diff epochs) 9(means no learning, just baseline), directory of plotlearning, ...
% ConsolidationStart, ConsolidationEnd, RegularExpressionsCell, {One day before start bidir, day of last bidir, {bidir syl 1, bidir syl 2, multidir syl 3...}}, LMANinactivation};

% Note:
% "bidir" is actually multidir, just enter additional syls if you want
% multidir

% for example:
% RegularExpressionsCell={'abbccbb', 'dccbb'}; OPTIONAL, only if will have
% to perform reg expr analysis.
% LMANinactivation= 1 or 0; (if does not exist, then is 0)

ExperimentList{1}={'pu53wh88','SyntaxDepPitchShift_abDowncbUp',2,'/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SyntaxDepPitchShift/SeqFilterCompile',...
    '16Jul2014','21Jul2014', {'abbb', 'accbb'}, {'07Jul2014','16Jul2014', 'aB', 'cB'}, 0};

% ExperimentList{22}={'pu53wh88','SyntaxDepPitchShift_cbUp',1,'/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SyntaxDepPitchShift_cbUp/SeqFilterCompile',...
%     '04Aug2014','06Aug2014', {'abbb', 'accbb'}, {}, 0}; % THROW OUT AS MAX LEARNING ONLY 50HZ

ExperimentList{2}={'pu53wh88','SeqDepPitchShift',1,'/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SeqDepPitchShift/SeqFilterCompile',...
    '06Nov2014','11Nov2014', {'abbb', 'accbb'}, {}, 0};

ExperimentList{3}={'pu53wh88','SeqDepPitchShift2',1','/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SeqDepPitchShift2/SeqFilterCompile',...
    '04Dec2014','11Dec2014', {'abbb', 'accbb'}, {}, 0};

ExperimentList{4}={'pu53wh88','SeqDepPitchShift3',0,'/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SeqDepPitchShift3/SeqFilterCompile',...
    '21Feb2015','26Feb2015', {'abbb', 'accbb'}, {'01Mar2015', '07Mar2015', 'abB', 'abbB'}, 0};


ExperimentList{5}={'pu11wh87','SyntaxDepPitchShift_abUP',1,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SyntaxDepPitchShift_abUP/SeqFilterCompile',...
    '12Jul2014','15Jul2014', {'abbccbb', 'dccbb'}, {}, 0};

ExperimentList{6}={'pu11wh87','SeqDepPitchShift2',0,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchShift2/SeqFilterCompile',...
    '03Dec2014','11Dec2014', {'abbccbb', 'dccbb'}, {'11Dec2014', '20Dec2014','dccB', 'bccB'}, 0};

ExperimentList{7}={'pu11wh87','SeqDepPitchShift3',0,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchShift3/SeqFilterCompile',...
    '26Feb2015','06Mar2015', {'abbccbb', 'dccbb'}, {}, 0}; % note, no bidir, actually switched target after end of unidir learning.

ExperimentList{8}={'pu11wh87','SyntaxDepPitchShift_cbDOWN',1,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SyntaxDepPitchShift_cbDOWN/SeqFilterCompile',...
    '28Jul2014','03Aug2014', {'abbccbb', 'dccbb'}, {}, 0}; % note, only analyzing abbccbb motif (ignore dccbb data entirely)

ExperimentList{9}={'pu11wh87','SeqDepPitchShift',2,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchShift/SeqFilterCompile',...
    '02Nov2014','06Nov2014', {'abbccbb', 'dccbb'},{'29Oct2014', '02Nov2014', 'aB', 'bccB'}, 0}; % note, only analyzing abbccbb motif, ignore dccbb motif entirely (they look the same)



ExperimentList{10}={'pu37wh20','SeqDepPitchShift',1,'/bluejay3/lucas/birds/pu37wh20/seq_dep_pitch_SeqDepPitchShift/SeqFilterCompile',...
    '05Nov2014','06Nov2014', {'abbb', 'abccbb'}, {}, 0};


ExperimentList{11}={'pu64bk13','SeqDepPitchShift',0,'/bluejay4/lucas/birds/pu64bk13/seq_dep_pitch_SeqDepPitchShift/SeqFilterCompile',...
    '17Dec2014','21Dec2014', {'dbbb'}, {'21Dec2014', '28Dec2014', 'dbB', 'dbbB'}, 0};

ExperimentList{12}={'pu64bk13','SeqDepPitchShift2',0,'/bluejay4/lucas/birds/pu64bk13/seq_dep_pitch_SeqDepPitchShift2/SeqFilterCompile',...
    '17Feb2015','22Feb2015', {'dbbb'}, {'28Feb2015', '07Mar2015', 'dB', 'dbB'}, 0};


ExperimentList{13}={'gr41gr90','SeqDepPitchShift',0,'/bluejay3/lucas/birds/gr41gr90/seq_dep_pitch_SeqDepPitchShift/SeqFilterCompile',...
    '15Feb2015','24Feb2015', {'bbacbb'}, {'28Feb2015', '07Mar2015', 'cbB', 'bBa'}, 0}; % NOTE - need to separate all gr41 expts to j motif only

ExperimentList{14}={'gr41gr90','SeqDepPitchShift2',0,'/bluejay3/lucas/birds/gr41gr90/seq_dep_pitch_SeqDepPitchShift2/SeqFilterCompile',...
    '18Apr2015','23Apr2015', {'bbacbb'}, {'10May2015', '19May2015', 'cB', 'cbB', 'Bba', 'bBa'}, 0};


ExperimentList{15}={'rd23gr89','SeqDepPitch',1,'/bluejay4/lucas/birds/rd23gr89/seq_dep_pitch_SeqDepPitch/SeqFilterCompile',...
    '09Jun2015','14Jun2015', {'ah', 'dbbgcb'}, {}, 0};


ExperimentList{16}={'pu35wh17','SeqDepPitch',1,'/bluejay4/lucas/birds/pu35wh17/seq_dep_pitch_SeqDepPitch/SeqFilterCompile',...
    '10Jun2015','14Jun2015', {'abbb', 'jbbb'}, {}, 0};

ExperimentList{17}={'pu35wh17','SeqDepPitch2',1,'/bluejay4/lucas/birds/pu35wh17/seq_dep_pitch_SeqDepPitch2/SeqFilterCompile',...
    '29Jun2015','30Jun2015', {'abbb', 'jbbb'}, {}, 0};

ExperimentList{18}={'pu35wh17','SeqDepPitch3',1,'/bluejay4/lucas/birds/pu35wh17/seq_dep_pitch_SeqDepPitch3/SeqFilterCompile',...
    '12Jul2015','13Jul2015', {'abbb', 'jbbb'}, {}, 0};


% ---------------------- LMAN EXPERIMENTS
ExperimentList{19}={'pu11wh87','SeqDepPitchLMAN', 0,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchLMAN/SeqFilterCompile',...
    '27May2015','08Jun2015', {'abbccbb', 'dccbb'}, {'08Jun2015', '16Jun2015', 'bccB', 'dccB'}, 1};
    % NEED TO USE EARLIER DAY THAN what I called consolidation start (based
    % on stability, about 8 days post 1st musc day (which is about 8 days post WN start), although they look
    % similar
    % - need to clean up labeling, as, likely because of WN, I get some
    % outliers on each day
    % - NEED TO LABEL CATCH SONG for bidir (that's when started catch song
    % method.) ONce do that, turn off skipping of post syls.
    % NOT CATCH

ExperimentList{20}={'pu11wh87','SeqDepPitchLMAN2', 0,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchLMAN2/SeqFilterCompile',...
    '07Jul2015','12Jul2015', {'abbccbb', 'dccbb'}, {'15Jul2015', '23Jul2015', 'aB', 'abB'}, 1}; % LABEL/MUSC GOOD - but need to account for WN during single dir.  
    % remove post syl (how many? abB is target. I will start just removing bC). SOME DAYS could use more labels, but important days
    % are relatively well labeled. Consol start is earliest MUSC day
    % rise: 8 days
    % hold: 8
    % bidir rise: 3
    % bidir hold: 8
    % NOT CATCH (SINGLE DIR)
    
ExperimentList{21}={'pu11wh87','SeqDepPitchLMAN3', 1,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchLMAN3/SeqFilterCompile',...
    '','', {'abbccbb', 'dccbb'}, {}, 1}; 

ExperimentList{22}={'pu53wh88','SeqDepPitchLMAN', 1,'/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SeqDepPitchLMAN/SeqFilterCompile',...
    '23Jun2015','24Jun2015', {'abbb', 'accbb'}, {}, 1}; % CONFIRM GOOD LABELING/INACTIVATION [to do: day 19 and 21 look potentially good.  need more labels late (after 3.5 hr) for those days.
%     the late days look good earlier (e.g. 2hrs post) (with higher musc).
%     Either way bird learned slowly and inactivation days are ~wk and more
%     during learning, so hard to compare to other experiments.

ExperimentList{23}={'gr41gr90','SeqDepPitchLMAN', 1,'/bluejay3/lucas/birds/gr41gr90/seq_dep_pitch_SeqDepPitchLMAN/SeqFilterCompile',...
    '','', {'jbbacbb'}, {}, 1}; % NOTE - LABEL/MUSC checked good. PROBLEM: FP on couple other B was high.

ExperimentList{24}={'gr41gr90','SeqDepPitchLMAN2', 1,'/bluejay3/lucas/birds/gr41gr90/seq_dep_pitch_SeqDepPitchLMAN2/SeqFilterCompile',...
    '','', {'jbbacbb', 'dbbacbb'}, {}, 1}; 


% ==== REPEATS - only select these if you want to analyze repeats
Params.global.repeats_analysis=0; 

% ExperimentList{1}={'pu53wh88','SyntaxDepPitchShift_abDowncbUp',2,'/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SyntaxDepPitchShift/SeqFilterCompile',...
%     '16Jul2014','21Jul2014', {'abbb', 'accbb'}, {'07Jul2014','16Jul2014', 'aB', 'cB'}, 0};




% ==================== PARAMS
PARAMS.global.timestamp_analy=lt_get_timestamp(0);
PARAMS.global.ExperimentList=ExperimentList;


%% PREPROCESS
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PREPROCESS(ExperimentList, PARAMS);

%% NOT ENOUGH SAMPLE FOR ACOUSTIC? 
% Display sample sizes for all expts and syls for structure stats (passing
% duration threshold), and all plot learning (all)
NumBirds=length(SeqDepPitch_AcrossBirds.birds);
disp(' ');
disp(' ======================================= ');
disp('N_acoustic/N_all');
disp('------------');

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments =length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        disp(' ');
        disp([birdname '-' exptname]);
        
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            N_acoustic=size(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_all_rends,1);
            N_all=numel(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF);
            
            if N_acoustic<N_all;
                
                disp([syl ':  ' num2str(N_acoustic) '/' num2str(N_all)]);
            end
        end
    end
end
        






%% ZSCORE DATA 
PARAMS.zscore.zthresh_learning=1; % take 1st N days after target passes this value of shift (will take max of default day (see below) or this day)
PARAMS.zscore.default_min_day=4; 
PARAMS.zscore.max_day=4; % inclusive, limit days to within min and max days (wn start =1) then throw out this experiment
PARAMS.zscore.N_days_post_thresh=1; % how many days to use to quantify learning
PARAMS.zscore.account_for_NoData=1; % If early WN days don't have data, then the window will be moved forward
% the number of no-data days (assumes no singing early days)
[SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_ZSCORE(SeqDepPitch_AcrossBirds, PARAMS);


%% DEFINE LEARNING - all future analysis will be done with this metric
% use z-score or FF.  which epoch?
NumBirds=length(SeqDepPitch_AcrossBirds.birds);


for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            
            % ============ DEFINE THE LEARNING METRIC
            if PARAMS.global.learning_metric=='zscore';
                % then use z-score, with epoch defined in ZSCORE function above
                if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh,'DATA'); % some don't have z-score defined becasue not enough elarning.
                    learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DATA.(syl).mean_Zscore;
                    learning_sem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DATA.(syl).sem_Zscore;
                    learning_DayInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DayInds;
                else
                    learning=nan;
                    learning_sem=nan;
                    learning_DayInds=nan;
                end
                
                
            elseif PARAMS.global.learning_metric=='ff_consolstart';
                % then use mean FFval at time of consolidation start
                
                learning = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.consolid_start_meanFF_minbase;
                learning_sem=lt_sem(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).ConsolStart_FFvals_minbase);
                consol_start=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolStartInd;
                consol_start_numdays=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.Consol_DayBins;
                learning_DayInds=consol_start:consol_start+consol_start_numdays-1;

            else
                disp('PROBLEM - LEARNING METRIC NOT DEFINED PROPERLY, ANALYSES WILL NOT WORK');
                pause;
            end
            
            % ---- PUT INTO OUTPUT STRUCTURE
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean=learning;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.sem=learning_sem;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.dayIndsUsed=learning_DayInds;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.LearningMetric.dayIndsUsed=learning_DayInds;
            
        end
        
        
        % ======== GET REALTUIVE TO TARGSYL
        % --- targsyl info
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        learning_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            learning_rel_targ=learning/learning_targ;
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ=learning_rel_targ;
            
        end
    end
end


%% 

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




%% SINGLE DIR LEARNING ANALYSIS


% ===== OPTIONS - removing experiments before analysis
% Option 1) REMOVE ANY EXPERIMENTS THAT DONT HAVE GOOD LEARNING
filter='good_learning';
[SeqDepPitch_AcrossBirds_GoodLearning, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);
% Remove LMAN experiments
filter='notLMAN';
[SeqDepPitch_AcrossBirds_GoodLearning, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds_GoodLearning, filter);
% Remove experiments where learning metric is "nan";
filter='learning_metric';
[SeqDepPitch_AcrossBirds_filtered, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);


% ====== LEARNING ANALYSIS [NOTED INSIDE WHAT IS USING CONSOLD START VS.
% LEARNING METRIC]
[SeqDepPitch_AcrossBirds_SINGLEDIR, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_SINGLEDIR(SeqDepPitch_AcrossBirds_filtered, PARAMS);




% Option 2) ==== RUN THIS IF DON'T WANT TO FILTER
% ====== LEARNING ANALYSIS
[SeqDepPitch_AcrossBirds_SINGLEDIR, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_SINGLEDIR(SeqDepPitch_AcrossBirds, PARAMS);


%% MULTIDIR LEARNING ANALYSES [HAVE NOT IMPLEMENTED LEARNING METRIC]
close all;
PARAMS.global.MULTIDIR.DayBinSize=2; % num days to take at start and end.

[SeqDepPitch_AcrossBirds_MULTIDIR, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_MULTIDIR(SeqDepPitch_AcrossBirds, PARAMS);


%% LMAN LEARNING ANALYSIS

% 1) SORT OUT ONLY THE THE EXPEIRMENTS THAT HAVE LMAN INACTIVATION DATA
% copy strcuture, save backup.
SeqDepPitch_AcrossBirds_ORIGINAL=SeqDepPitch_AcrossBirds;
filter = 'LMAN';
[SeqDepPitch_AcrossBirds_LMAN, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);


% 2) RUN VARIOUS ANALYSIS
% ==== BASELINE CORRELATIONS
close all;
[SeqDepPitch_AcrossBirds_LMAN, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMAN(SeqDepPitch_AcrossBirds_LMAN, PARAMS);


% ===== PLOT EFFECT OF MUSCIMOL ON LEARNING
% THROW OUT BAD SYLS? i.e. those immediately post WN since labeled
% not-catch songs. (will keep for plots summarizing experiments, but will
% throw out for actual analyses);
close all;
remove_bad_syls=1;
SylsToRemove_SingleDir={'pu11wh87','SeqDepPitchLMAN',{'bccbB'}, ...
    'pu11wh87','SeqDepPitchLMAN2', {'bC'}, ...
    'pu53wh88','SeqDepPitchLMAN', {'abbB'}}; % triplets {'birdname','exptname',syls(cell)}, eg. {'pu11wh87','SeqDepPitchShift',{'dccB','aB'});
SylsToRemove_BiDir={'pu11wh87','SeqDepPitchLMAN',{'bccbB', 'dccB'}};

PARAMS.global.LMAN.SylsToRemove_SingleDir=SylsToRemove_SingleDir;
PARAMS.global.LMAN.SylsToRemove_BiDir=SylsToRemove_BiDir;

% ===== SINGLE DIR LEARNING
norm_by_targsyl=1;
[PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANlearning(SeqDepPitch_AcrossBirds_LMAN, PARAMS, norm_by_targsyl);

% ==== PLOT FOR BIDIR LEARNING
close all;
[PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANbidir(SeqDepPitch_AcrossBirds_LMAN, PARAMS);

% ==== SEQUENCE LEARNING - EFFECT OF MUSC
close all; 
[PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANsquence(SeqDepPitch_AcrossBirds_LMAN, PARAMS);





% ==== [DIRTY!!!!] PLOT STUFF RELATIVE TO MOTIF DISTANCE
close all;
plotLMAN_musc=1; % if 0, then plots LMAN experiments, but PBS data
% FIRST, change the name of the field so the MUSC data will be where normal
% data are
SeqDepPitch_AcrossBirds_TEMP=SeqDepPitch_AcrossBirds_LMAN;
if plotLMAN_musc==1;
    for i=1:length(SeqDepPitch_AcrossBirds_TEMP.birds);
        numexperiments=length(SeqDepPitch_AcrossBirds_TEMP.birds{i}.experiment);
        
        for ii=1:numexperiments;
            SeqDepPitch_AcrossBirds_TEMP.birds{i}.experiment{ii}.Syl_ID_Dimensions=SeqDepPitch_AcrossBirds_TEMP.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN;
        end
    end
end
% Things like correlations, learning etc.
[SeqDepPitch_AcrossBirds_TEMP, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PLOT(SeqDepPitch_AcrossBirds_TEMP, PARAMS);

% Rel Target correaltions/acoustic vs. motif pos
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PLOT_DIRTY_reltarg(SeqDepPitch_AcrossBirds_TEMP, PARAMS);

% All pairs corre/acoustic vs. motif post
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PLOT_DIRTY_allpairs(SeqDepPitch_AcrossBirds_TEMP, PARAMS);

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% SAVING
savedir='/bluejay4/lucas/across_birds/seq_dep_pitch';
currdir=pwd;

try
    cd(savedir);
catch err
    mkdir(savedir);
    cd(savedir);
end

save('SeqDepPitch_AcrossBirds', 'SeqDepPitch_AcrossBirds')
try
    save('SeqDepPitch_AcrossBirds_ONETARG', 'SeqDepPitch_AcrossBirds_ONETARG')
catch err
end

cd(currdir)

