function [DatStructCompiled, Params]=lt_Opto_Stim_analy_EXTRACTDATA_v2(Params)


%% LT 4/7/15 - v2 - modified saving so overwrites instead of making new subdir


%% LT 1/15/15 - called from lt_Opto_Stim_analy script.
% run from day folder.
% extracts sound data and sorts for stim/WN - for opto combined with WN experiments

% e.g. Params entry:
% Params.batch='batch.labeled.all';
% Params.Fs=32000;
% WN TARG
% Params.SylTarg='b'; WN target
% Params.PreNote='c';
% Params.PreDur=0.25; % sec, how much data to take before targ
% Params.PostSylOnsetDur=0.3; % sec, data post to take
% Params.TargNoteNum=1; % template note directed to target
% Params.TargTrigCatchRate=0.5; % catch trial fraction at target
% STIM
% Params.notenum_stim=0; % of Stim - CONFIRMED THIS WORKS
% Params.StimDur=250; % duration of Stim, in ms
%
% for pitch contour
% Params.PC_freqrange=[2450 4050]; % for both pitch contour and find note
% Params.pc_harms=1;

% To sort into stim vs. notstim trials
% Params.StimLatencyWindow=[-250, 0]; % [a b] a and b are times (in negative ms) preceding syl onset. If stim is within this window, then will call a "stim trial."



%% AUTO PARAMS
PreTime=-Params.PreDur;
DatDur=(Params.PreDur+Params.PostSylOnsetDur)*Params.Fs; % how much data (time wise) to get?

% experiment ID default is "All"
if strcmp(Params.ExptID,'');
    Params.ExptID='All';
end


%% EXTRACT SOUND DATA and trigger data.

SylDat=findwnote2tw_v4_LT(Params.batch, Params.SylTarg,Params.PreNote,PreTime,...
    Params.PC_freqrange,DatDur,1,1,'obs0',0);


NumSyls=size(SylDat,2);

%% GET TIMING OF NEAREST TRIGGERS
% for desired note, get the nearest trigger.

% OLD VERSION - uses catch of 1st desired trig in song... bad
% if (1)
% for i=1:NumSyls;
%     SylOns=SylDat(i).ons(SylDat(i).NotePos); % onset time of syl
%
%     % see if preceded by trigger
%     if ~isempty(SylDat(i).ttimes); % i.e. if triggers exist in this song
%
%         % filter trigs and catch statuses to only look at desired note
%         TrigsOfDesiredNote=SylDat(i).ttimes(SylDat(i).NoteNums==Params.notenum_stim); % which triggers are of the desired note?
%         CatchStatusesOfDesiredNotes=SylDat(i).CatchList(SylDat(i).NoteNums==Params.notenum_stim);
%
%         if isempty(TrigsOfDesiredNote); % no trigger of wanted note detected
%
%             TrigToOnsetDur=[];
%             TrigCatch=[];
%
%         else
%
%             dummy=SylOns-TrigsOfDesiredNote; % times between trigs to syl onset (positive = trig occured first)
%             dummy=dummy(dummy>0); % only want trigs before syl
%             TrigToOnsetDur=min(dummy); % most recent trig to syl.
%
% %             dummy=TrigsOfDesiredNote-SylOns; % times between trigs to syl onset (negative = trig occured first)
% %             dummy=dummy(dummy>Params.StimLatencyWindow(1) & dummy<Params.StimLatencyWindow(2)); % only want trigs occuring within the window of interest
%
%             [~,ind]=min(abs(dummy)); % find ind of closest trig
%             TrigToOnsetDur=dummy(ind); % most recent trig to syl.
%
%             % was it catch?
%             Ind=find(dummy==TrigToOnsetDur);
%             TrigCatch=CatchStatusesOfDesiredNotes(Ind);
%         end
%     end
%
%     SylDat(i).MostRecentTrig.(['Note' num2str(Params.notenum_stim)]).TimeSince=TrigToOnsetDur;
%     SylDat(i).MostRecentTrig.(['Note' num2str(Params.notenum_stim)]).TrigCatch=TrigCatch;
%
% end
% end

% NEW VERSION - fixed problem with catch status
if (0)
    for i=1:NumSyls;
        SylOns=SylDat(i).ons(SylDat(i).NotePos); % onset time of syl
        
        % see if preceded by trigger
        if ~isempty(SylDat(i).ttimes); % i.e. if triggers exist in this song
            
            % filter trigs and catch statuses to only look at desired note
            TrigsOfDesiredNote=SylDat(i).ttimes(SylDat(i).NoteNums==Params.notenum_stim); % which triggers are of the desired note?
            CatchStatusesOfDesiredNotes=SylDat(i).CatchList(SylDat(i).NoteNums==Params.notenum_stim);
            
            if isempty(TrigsOfDesiredNote); % no trigger of wanted note detected
                
                TrigToOnsetDur=[];
                TrigCatch=[];
                
            else
                
                dummy=SylOns-TrigsOfDesiredNote; % times between trigs to syl onset (positive = trig occured first)
                dummy2=dummy(dummy>0); % only want trigs before syl
                TrigToOnsetDur=min(dummy2); % most recent trig to syl.
                
                %             dummy=TrigsOfDesiredNote-SylOns; % times between trigs to syl onset (negative = trig occured first)
                %             dummy=dummy(dummy>Params.StimLatencyWindow(1) & dummy<Params.StimLatencyWindow(2)); % only want trigs occuring within the window of interest
                
                TrigToOnsetDur=min(dummy2); % most recent trig to syl.
                
                % was it catch?
                Ind=find(dummy==TrigToOnsetDur);
                TrigCatch=CatchStatusesOfDesiredNotes(Ind);
            end
        end
        
        SylDat(i).MostRecentTrig.(['Note' num2str(Params.notenum_stim)]).TimeSince=TrigToOnsetDur;
        SylDat(i).MostRecentTrig.(['Note' num2str(Params.notenum_stim)]).TrigCatch=TrigCatch;
        
    end
end

%% GET TIMING OF NEAREST TRIGGERS - and decide if is stim or not.
% for desired note, get the nearest trigger.

if (1)
    StimInds=[];
    NoStimInds=[];
    
    
    for i=1:NumSyls;
        SylOns=SylDat(i).ons(SylDat(i).NotePos); % onset time of syl
        
        TrigToOnsetDur=[];
        TrigCatch=[];
        
        % see if preceded by trigger
        if isempty(SylDat(i).ttimes); % i.e. if triggers exist in this song
            
            % THEN IS NOT STIM - no triggers in song
            NoStimInds=[NoStimInds i];
        else
            % filter trigs and catch statuses to only look at desired note
            TrigsOfDesiredNote=SylDat(i).ttimes(SylDat(i).NoteNums==Params.notenum_stim); % which triggers are of the desired note?
            
            CatchStatuses=SylDat(i).CatchList;
            CatchStatusesOfDesiredNotes=CatchStatuses(SylDat(i).NoteNums==Params.notenum_stim);
            
            if isempty(TrigsOfDesiredNote); % no trigger of wanted note detected
                
                % THEN IS NOT STIM - no trigger of desired note
                NoStimInds=[NoStimInds i];
                
                
                
                
                
            else
                
                %             dummy=SylOns-TrigsOfDesiredNote; % times between trigs to syl onset (positive = trig occured first)
                %             dummy=dummy(dummy>0); % only want trigs before syl
                %             TrigToOnsetDur=min(dummy); % most recent trig to syl.
                
                dummy=SylOns-TrigsOfDesiredNote; % times between trigs to syl onset (positive = trig occured first)
                dummy2=dummy(-dummy>Params.StimLatencyWindow(1) & -dummy<Params.StimLatencyWindow(2)); % only want trigs occuring within the window of interest
                
                
                if length(dummy2)==0;
                    % THEN IS NOT STIM - trigger does not fall within desired
                    % window.
                    
                    NoStimInds=[NoStimInds i];
                    
                else
                    % THEN IS STIM TRIAL
                    StimInds=[StimInds i];
                    
                    
                    
                    
                    [~,ind]=min(abs(dummy2)); % find ind of closest trig
                    TrigToOnsetDur=dummy2(ind); % most recent trig to syl.
                    
                    % was it catch?
                    Ind=find(dummy==TrigToOnsetDur);
                    TrigCatch=CatchStatusesOfDesiredNotes(Ind);
                end
            end
            
            
            
        end
        
        SylDat(i).MostRecentTrig.(['Note' num2str(Params.notenum_stim)]).TimeSince=TrigToOnsetDur;
        SylDat(i).MostRecentTrig.(['Note' num2str(Params.notenum_stim)]).TrigCatch=TrigCatch;
        
    end
end

%% DECIDE IF A TRIAL IS CLOSE ENOUGH TO STIM TO COUNT AS "STIM TRIAL"
% i.e. using the time range specified in inputs.
% (all stim, even catch trial is called stim)

if (0)
    StimInds=[];
    NoStimInds=[];
    
    for i=1:NumSyls;
        TrigLaten=SylDat(i).MostRecentTrig.(['Note' num2str(Params.notenum_stim)]).TimeSince;
        
        if -TrigLaten>Params.StimLatencyWindow(1) & -TrigLaten<Params.StimLatencyWindow(2); % then is a stim trial
            StimInds=[StimInds i];
        else
            NoStimInds=[NoStimInds i];
        end
    end
    
    StimInds;
    NoStimInds;
end


% SORT STIM TRIALS BASED ON CATCH
% FIRST
% Want stimulation trials categorized into:
% 1) catch song and catch trial
% 2) catch song and notcatch trial
% Choosing catch Params.batch takes care of song catch criteria.

StimCatchInds=[];
StimNotCatchInds=[];

for i=1:NumSyls;
    WasItCatch=SylDat(i).MostRecentTrig.(['Note' num2str(Params.notenum_stim)]).TrigCatch;
    if WasItCatch==1;
        StimCatchInds=[StimCatchInds i];
    elseif WasItCatch==0;
        StimNotCatchInds=[StimNotCatchInds i];
    end
end

StimCatchInds;
StimNotCatchInds;


%% ASSIGN STIM IDs:
% StimID=0; not stim
% StimID=1; Stim, catch
% StimID=2; Stim, notcatch

StimID=nan(1,NumSyls);

for i=1:NumSyls;
    if sum(StimInds==i)==1; % then is stim trial
        if sum(StimCatchInds==i)==1; % then is a catch stim
            StimID(i)=1;
        elseif sum(StimCatchInds==i)==0; % then is not catch
            StimID(i)=2;
        end
    elseif sum(StimInds==i)==0; % then is not stim trial
        StimID(i)=0;
    end
end


% Assign Structures based on stim ID:
SylDatNoStim=SylDat(StimID==0);
SylDatStim=SylDat(StimID==1 | StimID==2);
SylDatStimCatch=SylDat(StimID==1);
SylDatStimNotCatch=SylDat(StimID==2);





%% EXTRACT TRIALS BASED ON WHETHER TARGET WAS TRIGGERED (ACTUAL TRIGGER)
% TargTrigID=-1; not trig, simulated catch
% TargTrigID=0; not trig, simlated not catch
% TargTrigID=1; trig, catch
% TargTrigID=2; trig, notcatch

% NOTE: Not trig includes both those purposely missed (i.e. pitch
% threshold) and those template failed to match.
% - That is fine, if my goal is to compare stim to not stim.
% - catch trial: for those triggered, catch trial is encoded into .rec.
%     for other syls, catch trial is not encoded, so I simulate using rand.


TargTrigID=nan(1,NumSyls);

for i=1:NumSyls;
    if (SylDat(i).TRIG==1 && SylDat(i).NoteNum==Params.TargNoteNum); % yes triggered, and is correct trigger note
        
        if SylDat(i).CatchTrial==1; % then is catch
            TargTrigID(i)=1;
        else % is not catch
            TargTrigID(i)=2;
        end
        
        
    else % not trig
        if rand<Params.TargTrigCatchRate; % Simulate catch for not trigs
            TargTrigID(i)=-1; % then is not trig, simluated catch
        else % then is not trig, simulated not catch
            TargTrigID(i)=0;
        end
    end
end

% Assign Structures based on TargTrigID:
SylDatTargNoTrig=SylDat(TargTrigID==0 | TargTrigID==-1);
SylDatTargNoTrigCatch=SylDat(TargTrigID==-1);
SylDatTargNoTrigNotCatch=SylDat(TargTrigID==0);

SylDatTargTrig=SylDat(TargTrigID==1 | TargTrigID==2);
SylDatTargTrigCatch=SylDat(TargTrigID==1);
SylDatTargTrigNotCatch=SylDat(TargTrigID==2);



%% EXTRACT VARIOUS COMBINATIONS OF TRIG AND STIM
% Do this using the stim and trig IDs


% i.e. TrigCatch and StimCatch

% i.e. TrigCatch and StimNotCatch

% i.e. NoTrig and StimCatch

% i.e. NoTrig and StimNotCatch

% simply trig catch (don't care what stim was) (takes a subset of syls,
% unbiased by whether or not template matched it)
SylDat_TargAllCatch=SylDat(TargTrigID==-1 | TargTrigID==1); % Stim catch, and target catch (whether or not triggered)


% i.e. StimNotCatch and AllTargCatch(TrigCatch + noTrig subset)
SylDat_StimCatch_TargAllCatch=SylDat(StimID==1 & (TargTrigID==-1 | TargTrigID==1)); % Stim catch, and target catch (whether or not triggered)

SylDat_StimNotCatch_TargAllCatch=SylDat(StimID==2 & (TargTrigID==-1 | TargTrigID==1)); % Stim catch, and target catch (whether or not triggered)

%


%% COMPILE ALL DATA STRUCTURES from above.

DatStructCompiled.All=SylDat;

DatStructCompiled.Stim=SylDatStim;
DatStructCompiled.NotStim=SylDatNoStim;
DatStructCompiled.StimCatch=SylDatStimCatch;
DatStructCompiled.StimNotCatch=SylDatStimNotCatch;

DatStructCompiled.TargNotTrig=SylDatTargNoTrig;
DatStructCompiled.TargNotTrig_Catch=SylDatTargNoTrigCatch;
DatStructCompiled.TargNotTrig_NotCatch=SylDatTargNoTrigNotCatch;

DatStructCompiled.TargTrig=SylDatTargTrig;
DatStructCompiled.TargTrigCatch=SylDatTargTrigCatch;
DatStructCompiled.TargTrigNotCatch=SylDatTargTrigNotCatch;

DatStructCompiled.COMBO_StimCatch_TargAllCatch=SylDat_StimCatch_TargAllCatch;
DatStructCompiled.COMBO_StimNotCatch_TargAllCatch=SylDat_StimNotCatch_TargAllCatch;
DatStructCompiled.COMBO_TargAllCatch=SylDat_TargAllCatch;

% DISPLAY sample sizes
DatFields=fieldnames(DatStructCompiled);
for i=1:length(DatFields);
    disp([DatFields{i} ' has ' num2str(length(DatStructCompiled.(DatFields{i}))) ' trials.']);
end


%% ADD ALL DETECTS -

if (0) % skip, as this problem (of not recording all detects online) is not a problem if
    %     I assume my template matches the target 100% of the time.  Then, if that trial is not designated as TRIG,
    %     then I know it had to have escaped
    
    % As it stands, SylDat has incomplete info on trigs - i.e. evtafv4 saves to
    % .rec file only trigs that are actually triggered (e.g. pass freq
    % trheshodl), so to get detects(without trigger), need to run offline.
    
    % However, a nuisance is that offline detection will therefore not have
    % catch statuses.  Therefore If I want to look at catch trials, I need to
    % take i) catch trial of online trigs + ii) % subset of offline detects
    % reflecting catch %.
    
    % So: here I will run offline to get all trigs, then add a SylDat field
    % with those ttimes
    
    % 1) First, run evtafsim offline to get X.rec files (since those not passing freq
    % threshold will not be on .rec)
    EvTAFv4Sim_LT(Params.batch,config_v4,'obs0');
    
    % 2) Using those X.rec files, get all ttimes
    ADDX=1;
    SylDatOffline=findwnote2tw_v4_LT(Params.batch, Params.SylTarg,Params.PreNote,PreTime,...
        frequency_range,DatDur,1,1,'obs0',0,ADDX);
    
    % 3) Update data structure
    for i=1:NumSyls;
        SylDat(i).ttimesAllDetect=SylDatOffline(i).ttimes; % all detects
    end
    
end


%% ============ CONVERT ANALOG LASER PULSES INTO DIGITAL (OPTIONAL)

if isfield(SylDat, 'dat_otherchans') && strcmp(Params.ExptID, 'Stim')
    
    % ========================= CONVERSION
    AllDifferences_LaserOnset = []; % sanity check - rec file and pulse show same times?
    AllDifferences_LaserOffset = [];
    
    LaserOnsets_CatchTrials = []; % sanity check - this should be empty
    for i=1:length(SylDat)
        
        %     figure; hold on; lt_plot_histogram(SylDat(i).dat_otherchans{1}, '', 1, 1, '', 0, 'k');
        
        tmp = diff(SylDat(i).dat_otherchans{1});
        
        %     figure; hold on;
        %     lt_plot_histogram(tmp, '', 1, 0, '', '', 'k')
        
        if max(tmp) > 100*median(abs(tmp)) & min(tmp) < -100*median(abs(tmp));
            % then this is stim trial
            ThrCrosses = midcross(SylDat(i).dat_otherchans{1}, Params.Fs);
%             assert(mod(length(ThrCrosses), 2)==0, 'THR CROSSES FOR LASER NOT EVEN!');
            if mod(length(ThrCrosses), 2)~=0
                % then 'THR CROSSES FOR LASER NOT EVEN!. assume that first
                % cross is onset of first stim...
            LaserOnsets_secFromSegmentOn = ThrCrosses(1:2:end-2);
            LaserOffsets_secFromSegmentOn = ThrCrosses(2:2:end-1);
            else
            LaserOnsets_secFromSegmentOn = ThrCrosses(1:2:end-1);
            LaserOffsets_secFromSegmentOn = ThrCrosses(2:2:end);
            end
        else
            LaserOnsets_secFromSegmentOn = [];
            LaserOffsets_secFromSegmentOn = [];
        end
        
        SylDat(i).LaserOnsets_secFromSegmentOn = LaserOnsets_secFromSegmentOn;
        SylDat(i).LaserOffsets_secFromSegmentOn = LaserOffsets_secFromSegmentOn;
        
        
        % ================= CONFIRM THAT MATCHES TRIG ONSET EXTRACTED FROM REC
        % FILE
        lasernote = ['Note' num2str(Params.notenum_stim)];
        laser_catch = SylDat(i).MostRecentTrig.(lasernote).TrigCatch;
        laser_time_onset = 1000*Params.PreDur - SylDat(i).MostRecentTrig.(lasernote).TimeSince; % milseconds
        
        if ~isempty(laser_catch) & ~isempty(LaserOnsets_secFromSegmentOn);
            if laser_catch == 1
                % then there should be no pulse recorded
                LaserOnsets_CatchTrials =[LaserOnsets_secFromSegmentOn'];
                
            else
                pulseonset_recfile = laser_time_onset+Params.StimDelay;
                AllDifferences_LaserOnset = [AllDifferences_LaserOnset ...
                    pulseonset_recfile-1000*LaserOnsets_secFromSegmentOn(1)];
                
                pulseoffset_recfile = [laser_time_onset+Params.StimDelay + Params.StimDur];
                AllDifferences_LaserOffset = [AllDifferences_LaserOffset ...
                    pulseoffset_recfile-1000*LaserOffsets_secFromSegmentOn(end)];
            end
        else
            % then no laser at all;
            LaserOnsets_CatchTrials =[LaserOnsets_secFromSegmentOn'];
        end
        
        
        
        % ==== TROUBLESHOOTING - overlay raw data, trigger time, and actual
        % pulses (and digitally extracted pulses)
        if (0)
            figure; hold on;
            
            lt_subplot(3,1,1); hold on;
            plot(SylDat(i).datt, 'k');
            plot(SylDat(i).dat_otherchans{1}, 'r');
            
            
            
            lt_subplot(3,1,2); hold on;
            lt_plot_spectrogram(SylDat(i).datt, Params.Fs, 1, 1);
            plot(1000*(1:length(SylDat(i).dat_otherchans{1}))/Params.Fs, SylDat(i).dat_otherchans{1}, 'g');
            
            %     LaserOnsetTime = SylDat(i). % in sec, relative to window on
            %    SylDat.ttimes
            lasernote = ['Note' num2str(Params.notenum_stim)];
            laser_catch = SylDat(i).MostRecentTrig.(lasernote).TrigCatch;
            laser_time_onset = 1000*Params.PreDur - SylDat(i).MostRecentTrig.(lasernote).TimeSince; % milseconds
            
            
            if laser_catch==0
                line([laser_time_onset laser_time_onset], ylim, 'Color', 'g');
                lt_plot_text(laser_time_onset+20, 1000, ['STIM ON'], 'g');
            else
                line([laser_time_onset laser_time_onset], ylim, 'Color', 'b');
                lt_plot_text(laser_time_onset+20, 1000, ['stim off (catch)'], 'b');
                
            end
            
            % ---- overlay extracted digital
            if isempty(SylDat(i).LaserOnsets_secFromSegmentOn)
                lt_plot_text(laser_time_onset+20, 2000, ['NO DIG stim extracted'], 'b');
            else
                % -- plot extracted digitals
                for j=1:length(SylDat(i).LaserOnsets_secFromSegmentOn)
                    line(1000*[SylDat(i).LaserOnsets_secFromSegmentOn(j) ...
                        SylDat(i).LaserOffsets_secFromSegmentOn(j)], [5000 5000], 'Color', 'b',...
                        'LineWidth', 2);
                end
            end
            
            pause; close all;
        end
        
        
    end
    
    lt_figure; hold on;
    lt_subplot(3,1,1); hold on;
    xlabel('onsets (rec file, accounting for delay minus dig input) (ms)');
    if ~isempty(AllDifferences_LaserOnset)
    lt_plot_histogram(AllDifferences_LaserOnset, '', 1, 0, '', '' ,'k');
    end
    
    lt_subplot(3,1,2); hold on;
    xlabel('offsets (rec file, accounting for delay minus dig input) (ms)');
        if ~isempty(AllDifferences_LaserOffset)
lt_plot_histogram(AllDifferences_LaserOffset, '', 1, 0, '', '' ,'k');
        end
        
        lt_subplot(3,1,3); hold on;
        xlabel('THIS SHOULD BE EMPTY - laser onsets on catch trials(rec file) (ms)');
    if ~isempty(LaserOnsets_CatchTrials)
        lt_plot_histogram(LaserOnsets_CatchTrials, '', 1, 0, '', '' ,'k');
    end
    
    
    % ==== REMOVE RAW SYLDAT
    SylDat = rmfield(SylDat, 'dat_otherchans');
    
end

%% SAVE
sdir=pwd;
savetime=lt_get_timestamp(0);
% savefolder=[sdir '/lt_Opto_Stim_analy_' Params.batch '_' savetime 'X'];
savefolder=[sdir '/lt_Opto_Stim_analy_v2'];

% PUT THINGS INTO PARAMS
Params.StimID=StimID;
Params.TargTrigID=TargTrigID;


% 1st dir
try
    cd(savefolder);
catch err
    mkdir(savefolder);
    cd(savefolder);
end

% 2nd dir
try 
    cd(Params.ExptID);
catch err
    mkdir(Params.ExptID); 
    cd(Params.ExptID);
end


% PUT folder name into params
Params.savefolder=[savefolder '/' Params.ExptID] ;



% if there is already a save file here, move it to old folder.
FilesInDir=dir('*.mat');

if length(FilesInDir)>0; 
    
    try 
        cd('OldAnalysis');
        
    catch err
        mkdir('OldAnalysis');
        cd('OldAnalysis');

    end
    
    mkdir(['Moved_' savetime]);
    cd ../
    
    eval(['!mv * OldAnalysis/Moved_' savetime]) % moves old save files
end

% save
save('DatStructCompiled.mat','DatStructCompiled');
save('Params.mat','Params');

% write a text file that tells you when files were made
fid1=fopen(['DONE_ExtractData_' savetime '.txt'],'w');
fclose(fid1);

mkdir('FIGURES')
cd('FIGURES')
mkdir('ExtractData');
lt_save_all_figs
cd ..









