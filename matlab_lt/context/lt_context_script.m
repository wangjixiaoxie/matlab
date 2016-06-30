%% ++++++++++++++++++++++++++++++++++++INTRO
% Method 1 and Method 2 below



%% +++++++++++++++++++++++++++++++++ METHOD 1 - autolabels using evtafv4 some - uses tameplate matching FF for context analysis
% === useful when have ginormous data, e.g. many days with rap[id context
% switching.

%% LT 4/26/15 - general script with directions for how to perform context analyses
% This requires using EvTAF_v4_LT_v1 or higher - i.e. has note group
% functionality.  It allows automatic switching of note groups (a note
% group is a set of notes that are imposed during a given epoch).


%% Single day data extraction
% extracts data using evtafsim, also extracts notegroup information. puts
% all into one structure and saves into a common dir.
lt_make_batch(1);

% -- 1) what batch and config file to use to run evtafsim.  if needed, make a batch.
Params.batch='batch.keep';
% Params.config= '/bluejay4/lucas/birds/rd23gr89/config_042115_fixedfreqbounds.evconfig2'; % for CtxtDepPitch (before 4/26)
% Params.config= '/bluejay4/lucas/birds/rd23gr89/042715_CtxtDepPitch2_TargB2_day2/config.evconfig2'; % for CtxtDepPitch2 (4/26 +)
Params.config= '/bluejay4/lucas/birds/rd23gr89/042815_CtxtDepPitch2_TargB2__day3/config.evconfig2'; % for CtxtDepPitch2 (4/28 +) (better hits)

Params.expt_phrase='Stimulation'; % folder where will save. leave blankt o get automatically from dir name
Params.dataID='All'; % e.g. id of data (e.g. 'All' for all data in data). if blank, uses batch name.

% -- Run evtafsim and save information
AllData=lt_context_ExtractDayData(Params);



%% DATA EXTRACTION - ACROSS DAYS
clear all; close all;
% metadata params
Params.metadata.experiment = 'LightDepPitch';
Params.metadata.condition='';
Params.metadata.notes='';
Params.metadata.date_range={'05May2015','05May2015'};
Params.metadata.only_labeled_dirs=0;

% analysis params
Params.analysis.batch='batch.keep';
% last version config (template used at very end of SeqDepPitch)
Params.analysiALL AL [EDGES FP REMOVED], edges checked s.config = '/bluejay4/lucas/birds/pu35wh17/061215_SeqDepPitch_durWN_day9/config_061115.evconfig2';
Params.analysis.dataID='All'; % e.g. id of data (e.g. 'All' for all data in data). if blank, uses batch name.
Params.analysis.Make_BatchKeep=1; % will make a new batch.keep for each day


Dirs_Missing_notegroups=lt_context_ExtractData_AcrossDays(Params);

%% BELOW, STUFF TO COLLECT DATA AND PLOT - NEED TO CONSOLIDATE INTO VARIOUS FUNCTIONS





%% COLLECT AllData Structures from multiple days
clear all; close all;

% Params
Params_alldays.CollectAllData=1; % if 1, then collects all data starting from first day
Params_alldays.firstday='05May2015';
Params_alldays.lastday='14May2015';

Params_alldays.NoteToPlot=2; % this is the note whose detects we will analyze (i.e. this note should get all renditions of the syl)
Params_alldays.RunBin=10;

Params_alldays.BoundaryTimes={'05May2014-1423', '08May2014-1423'}; % in format of e.g. 05May2014-1423, these are times of switching in experiment (e.g. turning WN off and on, changing pitch contingency, etc)

Params_alldays.Edge_Num_Rends = 20; % num rends to call "edges" (defualt: queries)
Params_alldays.Probe_CSplus=[1 2]; % [from to] (actual NG nums) (e.g. from no light --> light on(probe))
Params_alldays.Probe_CSminus=[1 3]; % [from to] (actual NG nums) (e.g. from no light --> no light (probe))

Params_alldays.PhaseToCompare1=4; % e.g. [light + WN up] phase
Params_alldays.PhaseToCompare2=5; % e.g. [light + WN dn] phase 

Params_alldays.throw_out_if_epoch_diff_days=1; % throws out any transitions that overlap with O/N (potentially 2 per O/N)

lt_context_CompileAndPlot(Params_alldays);




%% ++++++++++++++++++++++++++++++++ METHOD 2 - hand labeled - exctracts data, uses pitch contour

%% ======================= Extracting data across days
% --- day directories must be in format
% [date]_[experimentname]_[context_name].
% -- Run this in bird folder to extract all days data to subfolder

clear all; close all;

Params.syl_list={'b'}; % single syls
Params.freq_range_list={[3000 4000]};
Params.pc_dur_list=[0.11];
Params.pc_harms_list=[1];

Params.batch='batch.labeled.all.edges';
Params.experiment = 'CtxtDepPitch';

date_range={'15Nov2015','18Nov2015'}; % e.g. {'20Apr2015','20May2015'}. leave blank ('') for all days

lt_extract_AllDaysPC(Params, date_range)


%% %% Compiling data - go to "extract" folder first
clear all; close all;
Params_global.CompilePC.PC_window_list={'b', [40 70]}; % syl, value pairs [single syls]
Params_global.CompilePC.FirstDay='';
Params_global.CompilePC.LastDay='';
plotON=1; % pitch contours, all days, all syls
saveON=1;

% Regular expressions - first calculates FF etc, then performs regular
% expressions
Params_global.CompilePC.regexp_list={'c(b)'}; % e.g. {'dcc(b)', 'ab+(g)'} : dcc(b) means match dccb, and give me ind of b in parantheses.  ab+g means match ab... (anly length repeat), then g. give me ind of g
Params_global.CompilePC.regexp_fieldnames={'dccB','bccB'}; % whatever
% want to call structure field (if this cell array not defined, then will
% attempt to use the regexp names.
    
[ALLDATSTRUCT, Params_global]= lt_extract_CompilePC(plotON, Params_global, saveON);


%% Convert to context1 format
% --- TO BE ABLE TO RUN USING CONTEXT PLOT SAME AS FOR METHOD 1
Params_global.ConvertToContext1.NoteGroupNum_codes={'away', 1, 'home', 2, 'awayProbe', 3}; % {NoteGroup_name, NoteGroupNum} pairs - name must match what is in "condition" field. NoteGroupNum can be anything (keep it from 1, 2, ...);
% Params_global.ConvertToContext1.NoteNum_codes={'dcc_b_', 1, 'bcc_b_', 2}; % {notestring, notenum} pairs - notestring either single syl (e.g. 'a') or regexp, using underscores (e.g. 'dcc_b_')
Params_global.ConvertToContext1.NoteNum_codes={'c_b_', 1}; % {notestring, notenum} pairs - notestring either single syl (e.g. 'a') or regexp, using underscores (e.g. 'dcc_b_')
% syl='b';


[AllSongsDataMatrix, Params_alldays]= lt_context2_ConvertToContext1(ALLDATSTRUCT, Params_global);


%% PLOT
% USING SAME CONTEXT PLOT CODE FROM METHOD 1
close all;
Params_alldays.NoteToPlot=1;
Params_alldays.RunBin=10;

Params_alldays.BoundaryTimes={'15Nov2015-0000'}; % in format of e.g. 05May2014-1423, these are times of switching in experiment (e.g. turning WN off and on, changing pitch contingency, etc)
Params_alldays.Edge_Num_Rends = 40; % num rends to call "edges"

Params_alldays.throw_out_if_epoch_diff_days=0; % throws out any transitions that overlap with O/N (potentially 2 per O/N)
one_switch_a_day=1; % manual switching experiemnts.

lt_context_PLOT(AllSongsDataMatrix, Params_alldays, one_switch_a_day);



