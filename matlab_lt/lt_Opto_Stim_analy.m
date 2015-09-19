

%% LT 12/2014 - to analyze opto stim - both acute effects, and effects on learning - i.e. a syllable some delay away
% Can work with catch trials and multiple notes (i.e. evtafv4)

% PROCEDURE (in progress):
% 1) get data for all syl
% 2) get info on timing of nearest triggers (for all triggers)
% 3) get info on catch (trial)
% 3) process data as desired - e.g. pitch, entropy

% ASSUMPTIONS:
% 1) Stim (not WN, but e.g. laser stim) have no frequency criteria, so all trigs will be recorded online

% close all; clear all;

%% PARAMS

clear Params
clear StatsStruct
    Params.batch='batch.labeled.all';
    Params.Fs=32000;
%         WHERE TO ALIGN SONGS
    Params.SylTarg='b'; % aligns to onset of this
    Params.PreNote='c';
    
    Params.PreDur=0.25; % sec, how much data to take before targ
    Params.PostSylOnsetDur=0.35; % sec, data post to take
    Params.TargNoteNum=1; % template note directed to target
    Params.TargTrigCatchRate=0; % catch trial fraction at target

%         STIM
%     Params.NoteAfterStim=Params.SylTarg; % pick a note that happens after stim, will define trials as stim catch or not based on the stim occuring before this.
%     Params.NoteAfterStim_pre=Params.PreNote; % note before (i.e. btaf)
    Params.notenum_stim=0; % of Stim - CONFIRMED THIS WORKS
    Params.StimDelay=0; % in ms, delay from trigger to stim on
    Params.StimDur=230; % duration of Stim, in ms
    Params.StimType='pulse'; % either 'constant' or 'pulse'
    
%     	for pitch contour
    Params.PC_freqrange=[2450 4050]; % for both pitch contour and find note
    Params.pc_harms=1;

%     	To sort into stim vs. notstim trials
    Params.StimLatencyWindow=[-250, 0]; % [a b] a and b are times (in negative ms) preceding syl onset. If stim is within this window, then will call a "stim trial."

    
    % To note information about data
    Params.DataInfo='UnDir'; % 'Dir' means this is directed song.

    
    % FOR WN SONG
%     Params.batch='batch.labeled.all';
%     Params.Fs=32000;
% %         WHERE TO ALIGN SONGS
%     Params.SylTarg='c'; % aligns to onset of this
%     Params.PreNote='';
%     
%     Params.PreDur=0.2; % sec, how much data to take before targ
%     Params.PostSylOnsetDur=0.45; % sec, data post to take
%     Params.TargNoteNum=1; % template note directed to target
%     Params.TargTrigCatchRate=0.5; % catch trial fraction at target
% 
% %         STIM
% %     Params.NoteAfterStim=Params.SylTarg; % pick a note that happens after stim, will define trials as stim catch or not based on the stim occuring before this.
% %     Params.NoteAfterStim_pre=Params.PreNote; % note before (i.e. btaf)
%     Params.notenum_stim=0; % of Stim - CONFIRMED THIS WORKS
%     Params.StimDelay=0; % in ms, delay from trigger to stim on
%     Params.StimDur=230; % duration of Stim, in ms
%     Params.StimType='pulse'; % either 'constant' or 'pulse'
%     
% %     	for pitch contour
%     Params.PC_freqrange=[2450 4050]; % for both pitch contour and find note
%     Params.pc_harms=1;
% 
% %     	To sort into stim vs. notstim trials
%     Params.StimLatencyWindow=[-150, 0]; % [a b] a and b are times (in negative ms) preceding syl onset. If stim is within this window, then will call a "stim trial."
% 
%     
%     % To note information about data
%     Params.DataInfo='UnDir'; % 'Dir' means this is directed song.


% RUN DATA EXTRACT --------------------------------------------

[DatStructCompiled, Params]=lt_Opto_Stim_analy_EXTRACTDATA(Params);

% clear Params
% clear StatsStruct
%     Params.batch='batch.labeled.notcatch';
%     Params.Fs=32000;
% %         WHERE TO ALIGN SONGS
%     Params.SylTarg='c'; % aligns to onset of this
%     Params.PreNote='-';
%     
%     Params.PreDur=0.1; % sec, how much data to take before targ
%     Params.PostSylOnsetDur=0.5; % sec, data post to take
%     Params.TargNoteNum=0; % template note directed to target
%     Params.TargTrigCatchRate=0; % catch trial fraction at target
% 
% %         STIM
% %     Params.NoteAfterStim=Params.SylTarg; % pick a note that happens after stim, will define trials as stim catch or not based on the stim occuring before this.
% %     Params.NoteAfterStim_pre=Params.PreNote; % note before (i.e. btaf)
%     Params.notenum_stim=1; % of Stim - CONFIRMED THIS WORKS
%     Params.StimDelay=0; % in ms, delay from trigger to stim on
%     Params.StimDur=230; % duration of Stim, in ms
%     Params.StimType='pulse'; % either 'constant' or 'pulse'
%     
% %     	for pitch contour
%     Params.PC_freqrange=[2450 4050]; % for both pitch contour and find note
%     Params.pc_harms=1;
% 
% %     	To sort into stim vs. notstim trials
%     Params.StimLatencyWindow=[0, 200]; % [a b] a and b are times (in negative ms) preceding syl onset. If stim is within this window, then will call a "stim trial."
% 
% 
% % RUN DATA EXTRACT --------------------------------------------
% 
% [DatStructCompiled, Params]=lt_Opto_Stim_analy_EXTRACTDATA(Params);

%% (OPTIONAL) - SIMULATE CATCH BY TAKING RANDOM SAMPLES
% Is useful if using catch songs - there all trials are called catch (in .rec file), even
% if you have catch trials as well. Therefore, if want to get stim catch
% and notcatch (trial-wise) take subsample.

Params.SimCatchFrac=0.5;

[DatStructCompiled, Params]=lt_Opto_Stim_analy_SimulateCatch(DatStructCompiled,Params);





%% PROCESS AND PLOT -
% For only a single type of trial then only enter one field
% For 2 or more, will plot comaprisons, so enter all fields names
% Params.FieldsToCheck{1}='COMBO_StimCatch_TargAllCatch';
% Params.FieldsToCheck{2}='COMBO_StimNotCatch_TargAllCatch';
Params.FieldsToCheck{1}='StimCatch';
Params.FieldsToCheck{2}='StimNotCatch';

% Params.FieldsToCheck{1}='COMBO_StimCatch_TargAllCatch';
% Params.FieldsToCheck{2}='COMBO_StimNotCatch_TargAllCatch';

% Params.FieldsToCheck{1}='SIMULATE_Catch';
% Params.FieldsToCheck{2}='SIMULATE_NotCatch';

% Params.FieldsToCheck{1}='All';
% Params.FieldsToCheck{1}='COMBO_TargAllCatch';

% Params.FieldsToCheck{1}='COMBO_TargAllCatch';
 

% RUN
[StatsStruct, Params]=lt_Opto_Stim_analy_PLOT_Compare2(DatStructCompiled,Params);


%% TO PLOT TRIGGERED BY STIM TIME


lt_Opto_Stim_analy_PLOT_StimTriggered


%% TO PLOT HEAT MAP OF ALL TRIALS 



%% TO LOOK AT RAW DATA AND PULL OUT SPECIFIC SONGS TO LOOK AT



%% LOOK AT CERTAIN TIME WINDOW - EXTRACT AND PLOT AVERAGES
% Must first run lt_Opto_Stim_analy_PLOT_Compare2 to extract PC data. enter
% outputs from that into function here.
close all

% targetting b1, with WN (Experiment 1)
Params.TimeWindowList{1}=[70 115];
Params.TimeWindowList{2}=[178 210];
Params.TimeWindowList{3}=[255 266];


% just target b2 with WN (OLD) (1bin)
% Params.TimeWindowList{1}=[70 115];
% Params.TimeWindowList{2}=[173 205];
% Params.TimeWindowList{3}=[245 267];
% Params.TimeWindowList{4}=[325.5 326.5];

% % NEW - targetting b2 (2bins)
% Params.TimeWindowList{1}=[70 115];
% Params.TimeWindowList{2}=[173 205];
% Params.TimeWindowList{3}=[245 267];
% Params.TimeWindowList{4}=[325 331];

% NEW - targetting b2 (2bins) (changed to this0
Params.TimeWindowList{1}=[70 115];
Params.TimeWindowList{2}=[173 205];
Params.TimeWindowList{3}=[245 267];
Params.TimeWindowList{4}=[328 331];


% UNDIRECTED - (when I was stimming in WN off)
Params.TimeWindowList{1}=[70 115];
Params.TimeWindowList{2}=[173 205];
Params.TimeWindowList{3}=[245 267];
Params.TimeWindowList{4}=[328 348];

% DIRECTED
% Params.TimeWindowList{1}=[74 118];
% Params.TimeWindowList{2}=[175 206];
% Params.TimeWindowList{3}=[244 266];
% Params.TimeWindowList{4}=[328 347];


% during WN
% Params.TimeWindowList{1}=[70 115];
% Params.TimeWindowList{2}=[173 205];
% Params.TimeWindowList{3}=[253 262];

    for i=1:length(Params.TimeWindowList);
        TimeWind = Params.TimeWindowList{i};
        if ceil(Params.TimeWindowList{i})~=Params.TimeWindowList{i}; % means that there is fraction
            Params.TimeField{i}=['Wind' num2str(floor(TimeWind(1))) 'fr_' num2str(floor(TimeWind(2))) 'frms']; % put "fr" at end to signify fraction
        else
        Params.TimeField{i}=['Wind' num2str(TimeWind(1)) '_' num2str(TimeWind(2)) 'ms'];
        end
    end

    % RUN
[StatsStruct, Params]=lt_Opto_Stim_analy_PLOT_TimeWindow(StatsStruct,Params);

% DISPLAY SUMMARY STATS
lt_Opto_Stim_analy_PLOT_TimeWindow_Summary(StatsStruct,Params);


%% Heat map (z-scored pitch, relative to nonstim), ordered by trigger timing 
load Params
load StatsStruct
clear DATSTRUCT

close all
% PARAMS:
Params.HEATMAP.stim_delay=80; % in ms.(delay from trig to light onset
Params.HEATMAP.runningbin=10; % num trials to run over (sorted by stim onset)
% Params.HEATMAP.risecross=9; % 1st cross rising amplitude
% Params.HEATMAP.tmin=0.11e4; % index at start of window of interest
% Params.HEATMAP.tmax=0.32e4; % 


fieldname=['delay' num2str(Params.HEATMAP.stim_delay)]; % used for SUMMARY below

[Params, DATSTRUCT]=lt_Opto_Stim_analy_PLOT_HeatMap(Params,StatsStruct);


% SUMMARY - HEAT MAP OVER DIFFERENT BATCH FILES
lt_Opto_Stim_analy_SUMMARY_waterfall;


%% PLOT OVER TIME - to look at learning
% Takes output from lt_Opto_Stim_analy_PLOT_TimeWindow and plots that over
% time - scatter and running average;
close all;
Params.SmthBin=20; % smooth # rends

[StatsStruct,Params]=lt_Opto_Stim_analy_PLOT_TimeWindow_OverTime(StatsStruct,Params,1);



%% PERFORM STATISTICS on time window data 

[StatsStruct,Params]=lt_Opto_Stim_analy_PLOT_TimeWindow_Statistics(StatsStruct,Params);



