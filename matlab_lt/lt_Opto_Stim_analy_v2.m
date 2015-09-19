%% LT 4/7/15 - v2: changed to overwrite StatsStruct each time run a subprogram, instead of making new subdirectory,
% each struct was about 100MB, takes up space way too fast.
% all subprograms will be called v2


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

%%

    Params.batch='batch.labeled.all';
    Params.ExptID='Stim'; % an identifier for the data in this batch file. leave as '' to use 'All' as default. ('Stim' is for stim epochs).
    Params.DataInfo='UnDir'; % To note information about data, 'Dir' and 'UnDir' means this is directed or undirected song.
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
    Params.StimDur=240; % duration of Stim, in ms
    Params.StimType='const'; % either 'constant' or 'pulse'
    
%     	for pitch contour
    Params.PC_freqrange=[2450 4050]; % for both pitch contour and find note
    Params.pc_harms=1;

%     	To sort into stim vs. notstim trials
    Params.StimLatencyWindow=[-250, 0]; % [a b] a and b are times (in negative ms) preceding syl onset. If stim is within this window, then will call a "stim trial."


    
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

[DatStructCompiled, Params]=lt_Opto_Stim_analy_EXTRACTDATA_v2(Params);

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
% Params.FieldsToCheck{1}='NotStim';
% Params.FieldsToCheck{2}='Stim';



% Params.FieldsToCheck{1}='COMBO_StimCatch_TargAllCatch';
% Params.FieldsToCheck{2}='COMBO_StimNotCatch_TargAllCatch';

% Params.FieldsToCheck{1}='SIMULATE_Catch';
% Params.FieldsToCheck{2}='SIMULATE_NotCatch';

% Params.FieldsToCheck{1}='All';
% Params.FieldsToCheck{1}='COMBO_TargAllCatch';

% Params.FieldsToCheck{1}='COMBO_TargAllCatch';
 

% RUN
[StatsStruct, Params]=lt_Opto_Stim_analy_PLOT_Compare2_v2(DatStructCompiled,Params);


%% TO PLOT TRIGGERED BY STIM TIME


lt_Opto_Stim_analy_PLOT_StimTriggered


%% TO PLOT HEAT MAP OF ALL TRIALS 



%% TO LOOK AT RAW DATA AND PULL OUT SPECIFIC SONGS TO LOOK AT



%% LOOK AT CERTAIN TIME WINDOW - EXTRACT AND PLOT AVERAGES
% Must first run lt_Opto_Stim_analy_PLOT_Compare2 to extract PC data. enter
% outputs from that into function here.
close all

% pre 1/20, all syls
% Params.TimeWindowList{1}=[70 115];
% Params.TimeWindowList{2}=[178 210];
% Params.TimeWindowList{3}=[255 275];
% Params.TimeWindowList{4}=[338 358];
% 


% just target b2 with WN (OLD) (1bin) (i.e. for 3/27)
% Params.TimeWindowList{1}=[70 115];
% Params.TimeWindowList{2}=[173 205];
% Params.TimeWindowList{3}=[245 267];
% Params.TimeWindowList{4}=[325.5 326.5];

% % NEW - targetting b2 with WN (2bins)
% Params.TimeWindowList{1}=[70 115];
% Params.TimeWindowList{2}=[173 205];
% Params.TimeWindowList{3}=[245 267];
% Params.TimeWindowList{4}=[325 331];

% NEW - targetting b2 with WN(2bins) (changed to this0
Params.TimeWindowList{1}=[70 115];
Params.TimeWindowList{2}=[173 205];
Params.TimeWindowList{3}=[245 267];
Params.TimeWindowList{4}=[328 331];


% UNDIRECTED - noWN (when I was stimming in WN off)
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
    RunStats=1;
    KeepOutliers=0; % for running stats and plotting.
[StatsStruct, Params]=lt_Opto_Stim_analy_PLOT_TimeWindow_v2(StatsStruct,Params,RunStats,KeepOutliers);

% % DISPLAY SUMMARY STATS
% lt_Opto_Stim_analy_PLOT_TimeWindow_Summary(StatsStruct,Params);


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
KeepOutliers=0;

[StatsStruct,Params]=lt_Opto_Stim_analy_PLOT_TimeWindow_OverTime_v2(StatsStruct,Params,1,KeepOutliers);



%% PERFORM STATISTICS on time window data 

KeepOutliers=0;
[StatsStruct,Params]=lt_Opto_Stim_analy_PLOT_TimeWindow_Statistics_v2(StatsStruct,Params,KeepOutliers);





%% ================================== 
%% ============== SUMMARY ANALYSES

%% EXCTRACTING DATA ACROSS DAYS
% RUN THIS IN BIRD FOLDER.
clear all; close all;

% =======, to find directories
Params_metadata.experiment='Reversion1'; % 1st underscore ...
Params_metadata.condition='';
Params_metadata.notes='';
Params_metadata.date_range={'08Jul2015', '08Jul2015'};
Params_metadata.only_labeled_dirs=1;


% ===== For opto analysis
Params_glob.DataInfo='UnDir'; % To note information about data, 'Dir' and 'UnDir' means this is directed or undirected song.
Params_glob.Fs=32000;

%         WHERE TO ALIGN SONGS
Params_glob.SylTarg='c'; % aligns to onset of this
Params_glob.PreNote='';

Params_glob.PreDur=0.10; % sec, how much data to take before targ
Params_glob.PostSylOnsetDur=0.4; % sec, data post to take
Params_glob.TargNoteNum=1; % template note directed to target
Params_glob.TargTrigCatchRate=0; % catch trial fraction at target

%         STIM
Params_glob.notenum_stim=0; % of Stim - CONFIRMED THIS WORKS
Params_glob.StimDelay=55; % in ms, delay from trigger to stim on
Params_glob.StimDur=100; % duration of Stim, in ms
Params_glob.StimType='const'; % either 'constant' or 'pulse'

%     	for pitch contour
Params_glob.PC_freqrange=[2400 3650]; % for both pitch contour and find note
Params_glob.pc_harms=1;

%     	To sort into stim vs. notstim trials
Params_glob.StimLatencyWindow=[0, 100]; % [a b] a and b are times (in negative ms) preceding syl onset. If stim is within this window, then will call a "stim trial."

% Time Window List
% DURING WN (on window 3)
Params_glob.TimeWindowList{1}=[115 155];
Params_glob.TimeWindowList{2}=[197 240];
Params_glob.TimeWindowList{3}=[264 272];

% Plotting over time
Params_glob.SmthBin=20; % smooth # rends

% Do you want to delete old opto analysis folder if it exists?
Params_glob.Delete_old_analysis_folder=1;


% =============================== RUN
lt_Opto_Stim_analy_SUMMARY_MultDayAnaly_v3(Params_metadata, Params_glob);


%% PLOTTING DATA ACROSS DAYS - One analysis for each experiment

% ==== Params for analysis
BirdDir='/bluejay3/lucas/birds/wh73pk61/';
TimeFieldsOfInterest = 1:4; % i.e. time windows in pitch contour
statfield='ffvals';
BaselineDays=1:3;
plotStimEpochs=1; % if 1, then separates all data to stim epochs (even if multiple in one day)

% ====== ENTER DIRECTORIES OF INTEREST 2 WAYS - manually and with auto
% metadata collection (can do both)
% ==== METHOD 1) MANUALLY
ListOfDirs1=...
    {'030815_HVCChR2_XStim_StimOFF/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2009/PLOT_All_12Mar2015_2010/TimeWindow_12Mar2015_2010/OverTime_12Mar2015_2010',...
    '030815_HVCChR2_XStim_StimON_pulse/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2010/PLOT_StimCatch_StimNotCatch_12Mar2015_2012/TimeWindow_12Mar2015_2012/OverTime_12Mar2015_2012',...
    '030715_HVCChR2_XStim_StimOFF/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2015/PLOT_All_12Mar2015_2016/TimeWindow_12Mar2015_2016/OverTime_12Mar2015_2016',...
    '030715_HVCChR2_XStim_StimON_pulse/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2016/PLOT_StimCatch_StimNotCatch_12Mar2015_2018/TimeWindow_12Mar2015_2018/OverTime_12Mar2015_2018',...
    '030615_HVCChR2_XStim_StimON_pulse/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2022/PLOT_StimCatch_StimNotCatch_12Mar2015_2025/TimeWindow_12Mar2015_2025/OverTime_12Mar2015_2025',...
    '032615_HVCChR2_NoteGroups_StimONWnON_day6/lt_Opto_Stim_analy_batch.labeled.all_28Mar2015_1809X/PLOT_StimCatch_StimNotCatch_28Mar2015_1811/TimeWindow_28Mar2015_1813/OverTime_28Mar2015_1814/',...
    '032615_HVCChR2_NoteGroups_StimONWnOFF_day6/lt_Opto_Stim_analy_batch.labeled.all_28Mar2015_1828X/PLOT_StimCatch_StimNotCatch_28Mar2015_1830/TimeWindow_28Mar2015_1833/OverTime_28Mar2015_1833/',...
    '032715_HVCChR2_NoteGroups_StimOFFWnON_Targb2only_day1/lt_Opto_Stim_analy_v2/All',...
    '032715_HVCChR2_NoteGroups_StimONWnON_Targb2only_day1/lt_Opto_Stim_analy_v2/Stim',...
    '032815_HVCChR2_NoteGroups_StimONWnON_Targb2only_day2/lt_Opto_Stim_analy_batch.labeled.all_28Mar2015_1718X/PLOT_StimCatch_StimNotCatch_28Mar2015_1719/TimeWindow_30Mar2015_1232/OverTime_30Mar2015_1232',...
    '032815_HVCChR2_NoteGroups_StimOFFWnON_Targb2only_day2/lt_Opto_Stim_analy_v2/All',...    
    '032815_HVCChR2_NoteGroups_StimONWnON_Targb2only_day_stim2/lt_Opto_Stim_analy_batch.labeled.all_29Mar2015_1832X/PLOT_StimCatch_StimNotCatch_29Mar2015_1833/TimeWindow_30Mar2015_1234/OverTime_30Mar2015_1234',...
    '032915_HVCChR2_NoteGroups_StimONWnON_Targb2only_day3/lt_Opto_Stim_analy_batch.labeled.all_29Mar2015_1910X/PLOT_StimCatch_StimNotCatch_29Mar2015_1913/TimeWindow_30Mar2015_1238/OverTime_30Mar2015_1238',...
    '032915_HVCChR2_NoteGroups_StimOFFWnON_Targb2only_day3/lt_Opto_Stim_analy_v2/All',...
    '033015_HVCChR2_NoteGroups_StimOFFWnON_Targb2only_day4/lt_Opto_Stim_analy_v2/All',...
    '033015_HVCChR2_NoteGroups_StimONWnON_Targb2only_day4/lt_Opto_Stim_analy_batch.labeled.all_30Mar2015_1853X/PLOT_StimCatch_StimNotCatch_30Mar2015_1856/TimeWindow_30Mar2015_1857/OverTime_30Mar2015_1858/',...
    '033115_HVCChR2_NoteGroups_StimOFFWnON_Targb2only_day5/lt_Opto_Stim_analy_batch.labeled.all_01Apr2015_1306X/PLOT_All_01Apr2015_1311/TimeWindow_01Apr2015_1313/OverTime_01Apr2015_1314',...
    '033115_HVCChR2_NoteGroups_StimONWnON_Targb2only_day5/lt_Opto_Stim_analy_batch.labeled.all_01Apr2015_1315X/PLOT_StimCatch_StimNotCatch_01Apr2015_1316/TimeWindow_01Apr2015_1335/OverTime_01Apr2015_1335',...
    '040115_HVCChR2_NoteGroups_StimONWnON_Targb2only_day6/lt_Opto_Stim_analy_v2/Stim',...
    '040115_HVCChR2_NoteGroups_StimOffWnON_Targb2only_day6/lt_Opto_Stim_analy_v2/All'};
    
% ==== METHOD 2) METADATA automatically
MetadataStruct=lt_metadata_collect;

experiment = '';
condition='';
notes='';
date_range={'20Apr2015','20May2015'};
only_labeled_dirs=1;

ListOfDirs2=lt_metadata_find_dirs(MetadataStruct, experiment, condition, notes, date_range, only_labeled_dirs);

% get the correct subdir that contains opto stats
c=1;
for i=1:length(ListOfDirs2);
    cd(ListOfDirs2{i});
    
    tmp=[];
    try
        cd lt_Opto_Stim_analy_v2
        try cd 'All';
            tmp='All';
        catch err
            try cd 'Stim';
                tmp='Stim';
            end
        end
    catch err
        disp(['error - no lt_Opto_Stim_analy_v2 data on: ' ListOfDirs2{i} ' - throwing day out']);
    cd(BirdDir);        
        continue
    end
    
    % go back up to main dir
    cd(BirdDir);
    
    % modify name to add opto analysis name
    ListOfDirs2_modified{c}=[ListOfDirs2{i} '/lt_Opto_Stim_analy_v2/' tmp];
    c=c+1;
end


% ===== COMBINE DIRS
ListOfDirs_all=[ListOfDirs1 ListOfDirs2_modified];


lt_Opto_Stim_analy_SUMMARY_PlotOverTime(BirdDir, ListOfDirs_all,TimeFieldsOfInterest,statfield,BaselineDays,plotStimEpochs)
