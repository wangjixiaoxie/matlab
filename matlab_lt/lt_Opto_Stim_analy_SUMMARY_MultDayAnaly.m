clear all; close all
BaseDir='/bluejay3/lucas/birds/wh73pk61/';

ListOfDirs=...
    {'050615_ReversionWNup_StimOFFWnOFF/',...
    '050615_ReversionWNup_StimOFFWnON',...
    '050715_ReversionWNup_StimOFFWnON',...
    '050915_ReversionWNup_StimOFFWnON',...
    '051015_ReversionWNup_StimOFFWnON',...
    '051015_ReversionWNup_StimONWnON',...
    '051115_ReversionWNup_StimOFFWnON',...
    '051215_ReversionWNup_StimOFFWnON',...
    '051215_ReversionWNup_StimONWnON'};

% ListOfBatchNames={'batch.labeled.all',...
%     'batch.labeled.all',...
%     'batch.labeled.catch',...
%     'batch.labeled.catch',...
%     'batch.labeled.all',...
%     'batch.labeled.catch',...
%     'batch.labeled.all',...
%     'batch.labeled.all',...
%     'batch.labeled.all'};

ListOfTrialTypes=[0 0 0 0 0 1 0 0 1];
% if 0, then 'All'; if 1 then 'StimCatch' and 'StimNotCatch' (i.e is stim
% epoch)


%% RUN
for i=1:length(ListOfDirs);
    
    close all
    dirname=[BaseDir ListOfDirs{i}];
    cd(dirname);
    
    
    clear Params
    clear StatsStruct
    
    % -- get batch
    lt_v2_db_transfer_calls(2);
    Params.batch='batch.labeled.all';
    %     Params.batch=ListOfBatchNames{i};
    
    % -- other params
    Params.Fs=32000;
    %         WHERE TO ALIGN SONGS
    Params.SylTarg='b'; % aligns to onset of this
    Params.PreNote='';
    
    Params.PreDur=0.25; % sec, how much data to take before targ
    Params.PostSylOnsetDur=0.35; % sec, data post to take
    Params.TargNoteNum=1; % template note directed to target
    Params.TargTrigCatchRate=0; % catch trial fraction at target
    
    %         STIM
    %     Params.NoteAfterStim=Params.SylTarg; % pick a note that happens after stim, will define trials as stim catch or not based on the stim occuring before this.
    %     Params.NoteAfterStim_pre=Params.PreNote; % note before (i.e. btaf)
    Params.notenum_stim=0; % of Stim - CONFIRMED THIS WORKS
    Params.StimDelay=[]; % in ms, delay from trigger to stim on
    Params.StimDur=240; % duration of Stim, in ms
    Params.StimType=[]; % either 'constant' or 'pulse'
    
    %     	for pitch contour
    Params.PC_freqrange=[2450 4050]; % for both pitch contour and find note
    Params.pc_harms=1;
    
    %     	To sort into stim vs. notstim trials
    Params.StimLatencyWindow=[-250, 0]; % [a b] a and b are times (in negative ms) preceding syl onset. If stim is within this window, then will call a "stim trial."
    
    
    % RUN DATA EXTRACT --------------------------------------------
    
    [DatStructCompiled, Params]=lt_Opto_Stim_analy_EXTRACTDATA(Params);
    
    
    
    
    %% PROCESS AND PLOT -
    if ListOfTrialTypes(i)==0;
        Params.FieldsToCheck{1}='All';
    elseif ListOfTrialTypes(i)==1;
        Params.FieldsToCheck{1}='StimCatch';
        Params.FieldsToCheck{2}='StimNotCatch';
    end
    
    % RUN
    [StatsStruct, Params]=lt_Opto_Stim_analy_PLOT_Compare2(DatStructCompiled,Params);
    
    
    
    
    %% LOOK AT CERTAIN TIME WINDOW - EXTRACT AND PLOT AVERAGES
    % Must first run lt_Opto_Stim_analy_PLOT_Compare2 to extract PC data. enter
    % outputs from that into function here.
    close all
    Params.TimeWindowList{1}=[70 115];
    Params.TimeWindowList{2}=[173 205];
    Params.TimeWindowList{3}=[245 267];
    Params.TimeWindowList{4}=[328 348];
    for i=1:length(Params.TimeWindowList);
        TimeWind = Params.TimeWindowList{i};
        Params.TimeField{i}=['Wind' num2str(TimeWind(1)) '_' num2str(TimeWind(2)) 'ms'];
    end
    
    [StatsStruct, Params]=lt_Opto_Stim_analy_PLOT_TimeWindow(StatsStruct,Params);
    
    
    
    
    %% PLOT OVER TIME - to look at learning
    % Takes output from lt_Opto_Stim_analy_PLOT_TimeWindow and plots that over
    % time - scatter and running average;
    close all;
    Params.SmthBin=5; % smooth # rends
    
    [StatsStruct,Params]=lt_Opto_Stim_analy_PLOT_TimeWindow_OverTime(StatsStruct,Params,1);
    
end