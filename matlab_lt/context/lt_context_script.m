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
Params.analysis.config = '/bluejay4/lucas/birds/pu35wh17/061215_SeqDepPitch_durWN_day9/config_061115.evconfig2';
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

Params_input.Edge_Num_Rends = 20; % num rends to call "edges" (defualt: queries)
Params_alldays.Probe_CSplus=[1 2]; % [from to] (actual NG nums) (e.g. from no light --> light on(probe))
Params_alldays.Probe_CSminus=[1 3]; % [from to] (actual NG nums) (e.g. from no light --> no light (probe))

Params_alldays.PhaseToCompare1=4; % e.g. [light + WN up] phase
Params_alldays.PhaseToCompare2=5; % e.g. [light + WN dn] phase 

Params_alldays.throw_out_if_epoch_diff_days=1; % throws out any transitions that overlap with O/N (potentially 2 per O/N)

lt_context_CompileAndPlot(Params_alldays);









