%% LT 4/13/15 - Compiles data and plots across experiments/birds.  Need to specific dir of saved structures
function lt_seq_dep_pitch_ACROSSBIRDS
%% TO DO:

%% PARAMS TO ENTER
%% NEED DATA STRUCTURES
clear all; close all;

PARAMS.global.use_zscore_learning=1; % then uses zscore for all start of consolid learning stuff. [IN PROGRESS]
PARAMS.global.ConsolBinDaySize=2;  % days at start and end of consolid


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


ExperimentList{11}={'pu64bk13','SeqDepPitchShift',0,'/bluejay3/lucas/birds/pu64bk13/seq_dep_pitch_SeqDepPitchShift/SeqFilterCompile',...
    '17Dec2014','21Dec2014', {'dbbb'}, {'21Dec2014', '28Dec2014', 'dbB', 'dbbB'}, 0};

ExperimentList{12}={'pu64bk13','SeqDepPitchShift2',0,'/bluejay3/lucas/birds/pu64bk13/seq_dep_pitch_SeqDepPitchShift2/SeqFilterCompile',...
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

ExperimentList{20}={'pu11wh87','SeqDepPitchLMAN2', 0,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchLMAN2/SeqFilterCompile',...
    '07Jul2015','12Jul2015', {'abbccbb', 'dccbb'}, {'15Jul2015', '23Jul2015', 'aB', 'abB'}, 1}; % NEED TO UPDATE DATES ONCE LABEL MORE, ALSO CLEAN UP, REMOVED DAYS WITH LITTLE LMAN DATA.

ExperimentList{21}={'pu53wh88','SeqDepPitchLMAN', 1,'/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SeqDepPitchLMAN/SeqFilterCompile',...
    '23Jun2015','24Jun2015', {'abbb', 'accbb'}, {}, 1};

ExperimentList{22}={'gr41gr90','SeqDepPitchLMAN', 9,'/bluejay3/lucas/birds/gr41gr90/seq_dep_pitch_SeqDepPitchLMAN/SeqFilterCompile',...
    '','', {'jbbacbb'}, {}, 1}; % NOTE - many FP, do not use for single dir. Also, remove "d"

% ==================== PARAMS
PARAMS.global.timestamp_analy=lt_get_timestamp(0);
PARAMS.global.ExperimentList=ExperimentList;


%% EXTRACT EXPERIMENT DATA STRUCTURES
NumExpt=length(ExperimentList);
SeqDepPitch_AcrossBirds=struct('birds',[]); % holding structure

for i=1:length(ExperimentList);
    birdname=ExperimentList{i}{1};
    ExptID=ExperimentList{i}{2};
    NumTargs=ExperimentList{i}{3};
    PlotLearningDir=ExperimentList{i}{4};
    ConsolStartDate=ExperimentList{i}{5};
    ConsolEndDate=ExperimentList{i}{6};
    RegularExpressionsList=ExperimentList{i}{7};
    LMANinactivated=ExperimentList{i}{9};
    
    % Multidir stuff
    if ~isempty(ExperimentList{i}{8});
        MultiDir_OneDayBeforeStart=ExperimentList{i}{8}{1};
        MultiDir_LastDay=ExperimentList{i}{8}{2};
        
        num_multidir_syls=length(ExperimentList{i}{8})-2;
        MultiDirSyls={};
        for j=1:num_multidir_syls;
            MultiDirSyls{j}=ExperimentList{i}{8}{j+2};
        end
        
        multidir_exists=1;
    else
        multidir_exists=0;
    end
    
    
    % Slot into structure
    if isempty(SeqDepPitch_AcrossBirds.birds);
        % then no birds slotted yet
        % start bird 1
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.PlotLearningDir=PlotLearningDir;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.NumTargs=NumTargs;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.ExptID=ExptID;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.DATES.ConsolStartDate=ConsolStartDate;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.DATES.ConsolEndDate=ConsolEndDate;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.INFORMATION.RegularExpressionsList=RegularExpressionsList;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.INFORMATION.LMANinactivated=LMANinactivated;
        SeqDepPitch_AcrossBirds.birds{1}.birdname=birdname;
        
        
        if multidir_exists==1;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.DATES.MultiDir_OneDayBeforeStart=MultiDir_OneDayBeforeStart;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.DATES.MultiDir_LastDay=MultiDir_LastDay;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.INFORMATION.num_multidir_syls=num_multidir_syls;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.INFORMATION.MultiDirSyls=MultiDirSyls;
        end
    else
        % find the right bird to slot data into
        c=0;
        for j=1:length(SeqDepPitch_AcrossBirds.birds);
            if strcmp(SeqDepPitch_AcrossBirds.birds{j}.birdname,birdname);
                % then this is right bird
                % enter new data, new expt
                ind=length(SeqDepPitch_AcrossBirds.birds{j}.experiment)+1;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.PlotLearningDir=PlotLearningDir;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.NumTargs=NumTargs;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.ExptID=ExptID;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.DATES.ConsolStartDate=ConsolStartDate;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.DATES.ConsolEndDate=ConsolEndDate;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.INFORMATION.RegularExpressionsList=RegularExpressionsList;              
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.INFORMATION.LMANinactivated=LMANinactivated;
                SeqDepPitch_AcrossBirds.birds{j}.birdname=birdname;
                
                        if multidir_exists==1;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.DATES.MultiDir_OneDayBeforeStart=MultiDir_OneDayBeforeStart;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.DATES.MultiDir_LastDay=MultiDir_LastDay;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.INFORMATION.num_multidir_syls=num_multidir_syls;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.INFORMATION.MultiDirSyls=MultiDirSyls;
                        end
                
                % change c to show that bird found.
                c=1;
                
                continue % quit loop
            end
        end
        if c==0
            % then no bird was found, start new bird
            ind=length(SeqDepPitch_AcrossBirds.birds)+1;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.PlotLearningDir=PlotLearningDir;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.NumTargs=NumTargs;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.ExptID=ExptID;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.DATES.ConsolStartDate=ConsolStartDate;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.DATES.ConsolEndDate=ConsolEndDate;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.INFORMATION.RegularExpressionsList=RegularExpressionsList;              
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.INFORMATION.LMANinactivated=LMANinactivated;
            SeqDepPitch_AcrossBirds.birds{ind}.birdname=birdname;
            
                    if multidir_exists==1;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.DATES.MultiDir_OneDayBeforeStart=MultiDir_OneDayBeforeStart;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.DATES.MultiDir_LastDay=MultiDir_LastDay;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.INFORMATION.num_multidir_syls=num_multidir_syls;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.INFORMATION.MultiDirSyls=MultiDirSyls;
                    end
                    
        end
    end
end



%-- LOAD AND COMPILE STRUCTURES
for i=1:length(SeqDepPitch_AcrossBirds.birds);
    for ii=1:length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        % load
        tmp=load([SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir '/Params']);
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params=tmp.Params;
        tmp=load([SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir '/AllDays_PlotLearning']);
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning=tmp.AllDays_PlotLearning;
        
    end
end


% -- HOW many birds? 
NumBirds=length(SeqDepPitch_AcrossBirds.birds);


% -- MAKE RANDOM EDITS (ad hoc)
 
% pu37wh20, % remove bccB as this was actually also a target (along with dccB);
% Don't need to, not in data nymore...
% for i=1:NumBirds;
%     if strcmp(SeqDepPitch_AcrossBirds.birds{i}.birdname,'pu37wh20');
%     
%         for ii=1:length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%             if strcmp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID,'SeqDepPitchShift');
%                 
% %                 % remove stuff
% %                 ii=find(strcmp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.SylsSame,'bccB'))
%                 
%             end
%     
%     end
% end


% CHECK DIRECTION OF PITCH SHIFT


%% == MAKE SURE Unique syls are defined for all expts (and if not, then do so)

SeqDepPitch_AcrossBirds=lt_seq_dep_pitch_ACROSSBIRDS_UniqueSyls(SeqDepPitch_AcrossBirds);




% ===== PLOT SUMMARY OF METADATA OF EXPERIMENTS (E.G. SAMPLE SIZE)
% IN PROGRESS
% e.g. bar plot, each col is a bird, with height as num experiments of
% different ytpes)


%% Correlations and Acoustic distance analyses.
close all;
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_CORR_and_PCA(SeqDepPitch_AcrossBirds, PARAMS);




%% PLOT AND EXTRACT LEARNING PARAMS, LEARNING, SUMMARY PLOTS FOR EACH EXPERIMENT
close all;
DayBinSize=PARAMS.global.ConsolBinDaySize;
[SeqDepPitch_AcrossBirds]=lt_seq_dep_pitch_ACROSSBIRDS_SummaryPlot(SeqDepPitch_AcrossBirds, DayBinSize);

%% MAKE SURE ALL LMAN EXPERIMENTS HAVE LMAN RAW DATA (e.g. regexpr) - IF NOT, RUN ANALYSES

currdir=pwd;
close all;

for i=1:NumBirds;
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        % SKIP IF THIS IS NOT LMAN EXPERIMENT
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==0;
            continue;
        end
        
        
        % ==== FIRST, make sure that regexpr has baselkine data for both
        % musc (must already have pbs, because checked in previous code)
        if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr, 'baseline_MUSC');
            % Then need to go and do those analysis
            
            close all;
            disp(' ')
            disp('NEED TO DO ANALYSIS REGEXPR for MUSC');
            disp('Running...')
            
            % 1) go to dir
            cd(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir);
            
            % 2) run RegExpr and RegExprPlot analyses
            clear Params;
            clear AllDays_RegExpr;
            clear AllDays_RawDatStruct;
            clear AllDays_PlotLearning
            
            load('Params');
            load('AllDays_RawDatStruct');
            load('AllDays_PlotLearning');
            
            % ------ RUN REG EXPR
            Params.RegExpr.expressions=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.RegularExpressionsList;
            saveON=1;
            DoLMAN=1;
            [Params, AllDays_RegExpr] = lt_seq_dep_pitch_RegExpr(Params, AllDays_RawDatStruct, saveON, DoLMAN, AllDays_PlotLearning);
            
            Params.PlotRegExpr.plotWNdays=1;
            saveON=1;
            LMANon=1;
            SuppressQueries=1;
            [Params, AllDays_RegExpr]=lt_seq_dep_pitch_PlotRegExpr(Params, AllDays_RegExpr ,saveON, LMANon, SuppressQueries);
            
            % ---------- SAVE OUTPUTS
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params=Params;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr=AllDays_RegExpr;
            
        end
        
        % ==== SECOND, make sure has baseline data for PBS
        if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr, 'baseline');
            % Then need to go and do those analysis (does not remove any
            % anlsyis of LMAN data done above).
            
            close all;
            disp(' ')
            disp('NEED TO DO ANALYSIS REGEXPR for PBS');
            disp('Running...')
            
            % 1) go to dir
            cd(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir);
            
            % 2) run RegExpr and RegExprPlot analyses
            clear Params;
            clear AllDays_RegExpr;
            
            load('Params');
            load('AllDays_RegExpr');
            
            % ------ RUN REG EXPR PLOR
            Params.PlotRegExpr.plotWNdays=1;
            saveON=1;
            LMANon=0;
            SuppressQueries=1;
            [Params, AllDays_RegExpr]=lt_seq_dep_pitch_PlotRegExpr(Params, AllDays_RegExpr ,saveON, LMANon, SuppressQueries);
            
            % ---------- SAVE OUTPUTS
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params=Params;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr=AllDays_RegExpr;
            
        end
        
    end
end
disp('DONE!...');
cd(currdir);


%% FILTER - To give each syl a sylID (its features)
close all; 
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_FILTER(SeqDepPitch_AcrossBirds, PARAMS);


%% [LMAN] FILTER - To give each syl a sylID (its features) (LMAN EXPERIMENTS)
% Modified to allow to extract for LMAN data too.  Only does one thing -
% takes the output Syl_ID_Dimensions field and adds that as a new field
% Syl_ID_Dimensions_LMAN to the structure

for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    for ii=1:numexperiments;
        
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1;
            % Then this has LMAN data
            
            
            % === FIRST, extract only that experiment (don't overwrite
            % struct) (replace bird 1 expt 1 with the current expt)
            SeqDepPitch_AcrossBirds_ORIGINAL=SeqDepPitch_AcrossBirds;
            
            SeqDepPitch_AcrossBirds.birds{1}.experiment{1}=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii};
            
            filter = 'LMAN';
            [SeqDepPitch_AcrossBirds_LMAN, ~]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);
            
            
            % ==== SECOND, change the name of the field so that the next
            % program analyzes LMAN data, instead of PBS. (also remove the
            % SylID field that was previosly done).
            SeqDepPitch_AcrossBirds_LMAN.birds{1}.experiment{1}=rmfield(SeqDepPitch_AcrossBirds_LMAN.birds{1}.experiment{1}, ...
                'Syl_ID_Dimensions');
            
            disp(' ');
            disp('NOTE: LMAN Data filter will work for: ');
            % ---- To get correlations
            SeqDepPitch_AcrossBirds_LMAN.birds{1}.experiment{1}.Data_RegExpr.AllDays_RegExpr.baseline=...
                SeqDepPitch_AcrossBirds_LMAN.birds{1}.experiment{1}.Data_RegExpr.AllDays_RegExpr.baseline_MUSC;
            
            disp('Correlations...');
            
            disp('DOES NOT YET WORK FOR: ');
            disp('Acoustic distance');
            disp('Learning_by_target (not important)')
            disp(' ');
            
            % ==== THIRD, RUN FILTER on LMAN DATA
            [SeqDepPitch_AcrossBirds_LMAN, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_FILTER(SeqDepPitch_AcrossBirds_LMAN, PARAMS);
            
            
            
            % ===== FOURTH, Take the SylID field from filtered LMAN data and put it
            % intot he global all experimtns structure
            SeqDepPitch_AcrossBirds_ORIGINAL.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN=SeqDepPitch_AcrossBirds_LMAN.birds{1}.experiment{1}.Syl_ID_Dimensions;
            
            % ====== RESTORE THE ORIGINAL SeqDepPitch_AcrossBirds, with one
            % change (Syl_ID_Dimensions_LMAN is added for one expt/bird);
            SeqDepPitch_AcrossBirds=SeqDepPitch_AcrossBirds_ORIGINAL;
        end
    end
end

%% BASELINE DATA - MAKE SURE ALL HAVE FFVALS AND TVALS IN BASELINE DATA (I DID NOT CODE THIS INTO THE ORIGINAL PLOT LEARNING FUNCTION)

for i=1:NumBirds;
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        baseline_days=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.BaselineDays;
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            % ==== GET VECTOR OF BASELINE FF AND TIMES
            ffvals_all=[];
            tvals_all=[];
            for k=baseline_days;
                ffvals_all=[ffvals_all cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})];
                tvals_all=[tvals_all cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals{k})];
            end
               
            % ---- make sure is same as what was done before
            n_old=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).n;
            meanFF_old=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF;
            
            n_new=length(ffvals_all);
            meanFF_new=mean(ffvals_all);
            
            if n_old~=n_new || meanFF_old~=meanFF_new;
                disp('PROBLEM - NEW BASELINE ANALYSIS DOES NOT MATCH OLD BASELINE');
            end
            
            % ----- UPDATE DATA 
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).Tvals=tvals_all; % Tvals
            disp([num2str(i) ' ' num2str(ii) ' ' syl]);
        end
    end
end


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



%% MULTIDIR LEARNING ANALYSES
close all;
PARAMS.global.MULTIDIR.DayBinSize=2; % num days to take at start and end.

[SeqDepPitch_AcrossBirds_MULTIDIR, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_MULTIDIR(SeqDepPitch_AcrossBirds, PARAMS);


%% SINGLE DIR LEARNING ANALYSIS


% ===== OPTIONS - removing experiments before analysis
% Option 1) REMOVE ANY EXPERIMENTS THAT DONT HAVE GOOD LEARNING
filter='good_learning';
[SeqDepPitch_AcrossBirds_GoodLearning, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);
% ===== BASELINE DRIFT ANALYSIS
[SeqDepPitch_AcrossBirds_GoodLearning, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_DRIFT(SeqDepPitch_AcrossBirds_GoodLearning, PARAMS);
% ====== LEARNING ANALYSIS
[SeqDepPitch_AcrossBirds_SINGLEDIR, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_SINGLEDIR(SeqDepPitch_AcrossBirds_GoodLearning, PARAMS);


% Option 2) ==== RUN THIS IF DON'T WANT TO FILTER
% ===== BASELINE DRIFT ANALYSIS
[SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_DRIFT(SeqDepPitch_AcrossBirds, PARAMS);
% ====== LEARNING ANALYSIS
[SeqDepPitch_AcrossBirds_SINGLEDIR, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_SINGLEDIR(SeqDepPitch_AcrossBirds, PARAMS);

%% LMAN LEARNING ANALYSIS
close all;
[SeqDepPitch_AcrossBirds_LMAN, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMAN(SeqDepPitch_AcrossBirds, PARAMS);


% ===== PLOT EFFECT OF MUSCIMOL ON LEARNING
close all;
norm_by_targsyl=1;
[SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANlearning(SeqDepPitch_AcrossBirds, PARAMS, norm_by_targsyl);

% ==== PLOT FOR BIDIR LEARNING
close all;
[SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANbidir(SeqDepPitch_AcrossBirds, PARAMS);



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

