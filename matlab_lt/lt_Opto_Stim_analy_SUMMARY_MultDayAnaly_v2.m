%% LT 7/9/15 - MODIFIED TO AUTOMATICALLY FIND DIR NAMES (uses lt_metadata analysis)
clear all; close all;

% PARAMS, to find directories
experiment='Reversion1'; % 1st underscore ...
condition='';
notes='';
date_range={'01Jul2015', '09Jul2015'};
only_labeled_dirs=1;
prog_vers=2; % version of program to run - DONT CHANGE

% FIRST, RUN IN FOLDER CONTAINING ALL DAY DIRS\
MetadataStruct=lt_metadata_collect;

% SECOND, extracts dirs that meet criteria given above.
[MetadataStruct_filtered]=lt_metadata_find_dirs(MetadataStruct, experiment, condition, notes, date_range, only_labeled_dirs, prog_vers);


%% goes thru days and runs lt_Opto_Stim_analy_v2 if required; (5/16/15)

clear all; close all
BaseDir='/bluejay3/lucas/birds/wh73pk61/';

ListOfDirs=...
    {'032715_HVCChR2_NoteGroups_StimOFFWnON_Targb2only_day1', ...
    '032815_HVCChR2_NoteGroups_StimOFFWnON_Targb2only_day2',...
    
    '032915_HVCChR2_NoteGroups_StimOFFWnON_Targb2only_day3',...
    '051315_ReversionWNup_StimOFFWnON', ...
    '051415_ReversionWNup_StimOFFWnON',...
    '051415_ReversionWNup_StimONWnON',...
    '051515_ReversionWNup_StimOFFWnON', ...
    '051615_ReversionWNup_StimOFFWnON',...
    '051615_ReversionWNup_StimONWnON',...
    '051715_ReversionWNup_StimOFFWnOFF',...
    '051815_ReversionWNup_StimOFFWnOFF',...
    '051915_ReversionWNup_StimOFFWnOFF',...
    '052015_ReversionWNup_StimONWnOFF',...
    '052015_ReversionWNup_StimOFFWnOFF'};

% ListOfBatchNames={'batch.labeled.all',...
%     'batch.labeled.all',...
%     'batch.labeled.catch',...
%     'batch.labeled.catch',...
%     'batch.labeled.all',...
%     'batch.labeled.catch',...
%     'batch.labeled.all',...
%     'batch.labeled.all',...
%     'batch.labeled.all'};

ListOfTrialTypes=[0 0 0 0 0 1 0 0 1 0 0 0 1 0];
% if 0, then 'All'; if 1 then 'StimCatch' and 'StimNotCatch' (i.e is stim
% epoch)


%% RUN
for i=1:length(ListOfDirs);
    
    close all
    dirname=[BaseDir ListOfDirs{i}];
    cd(dirname);
    
    
    clear Params
    clear StatsStruct
    
    %% Extract
    % -- get batch
    lt_v2_db_transfer_calls(2);
    Params.batch='batch.labeled.all';
    %     Params.batch=ListOfBatchNames{i};
    
    if ListOfTrialTypes(i)==0;
        % then is no stim epoch
        Params.ExptID='All'; % an identifier for the data in this batch file. leave as '' to use 'All' as default. ('Stim' is for stim epochs).
    elseif ListOfTrialTypes(i)==1;
        % then is stim epoch
        Params.ExptID='Stim';
    end
    
    Params.DataInfo='UnDir'; % To note information about data, 'Dir' and 'UnDir' means this is directed or undirected song.
    Params.Fs=32000;
    %         WHERE TO ALIGN SONGS
    Params.SylTarg='b'; % aligns to onset of this
    Params.PreNote='';
    
    Params.PreDur=0.25; % sec, how much data to take before targ
    Params.PostSylOnsetDur=0.35; % sec, data post to take
    Params.TargNoteNum=1; % template note directed to target
    Params.TargTrigCatchRate=0; % catch trial fraction at target
    
    %         STIM
    Params.notenum_stim=0; % of Stim - CONFIRMED THIS WORKS
    Params.StimDelay=0; % in ms, delay from trigger to stim on
    Params.StimDur=240; % duration of Stim, in ms
    Params.StimType='const'; % either 'constant' or 'pulse'
    
    %     	for pitch contour
    Params.PC_freqrange=[2450 4050]; % for both pitch contour and find note
    Params.pc_harms=1;
    
    %     	To sort into stim vs. notstim trials
    Params.StimLatencyWindow=[-250, 0]; % [a b] a and b are times (in negative ms) preceding syl onset. If stim is within this window, then will call a "stim trial."
    
    
    % RUN DATA EXTRACT --------------------------------------------
    
    [DatStructCompiled, Params]=lt_Opto_Stim_analy_EXTRACTDATA_v2(Params);
    
    
    
    
    %% PROCESS AND PLOT -
    if ListOfTrialTypes(i)==0;
        Params.FieldsToCheck{1}='All';
    elseif ListOfTrialTypes(i)==1;
        Params.FieldsToCheck{1}='StimCatch';
        Params.FieldsToCheck{2}='StimNotCatch';
    end
    
    % RUN
    [StatsStruct, Params]=lt_Opto_Stim_analy_PLOT_Compare2_v2(DatStructCompiled,Params);
    
    
    
    
    %% LOOK AT CERTAIN TIME WINDOW - EXTRACT AND PLOT AVERAGES
    % Must first run lt_Opto_Stim_analy_PLOT_Compare2 to extract PC data. enter
    % outputs from that into function here.
    close all
    Params.TimeWindowList{1}=[70 115];
    Params.TimeWindowList{2}=[173 205];
    Params.TimeWindowList{3}=[245 267];
    Params.TimeWindowList{4}=[328 331];
    
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
    
    
    
    
    %% PLOT OVER TIME - to look at learning
    % Takes output from lt_Opto_Stim_analy_PLOT_TimeWindow and plots that over
    % time - scatter and running average;
    close all;
    Params.SmthBin=20; % smooth # rends
    KeepOutliers=0;
    
    [StatsStruct,Params]=lt_Opto_Stim_analy_PLOT_TimeWindow_OverTime_v2(StatsStruct,Params,1,KeepOutliers);
    
    
end