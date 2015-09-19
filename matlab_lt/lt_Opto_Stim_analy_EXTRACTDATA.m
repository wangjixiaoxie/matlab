function [DatStructCompiled, Params]=lt_Opto_Stim_analy_EXTRACTDATA(Params)
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


%% SAVE
sdir=pwd;
savetime=lt_get_timestamp(0);
savefolder=[sdir '/lt_Opto_Stim_analy_' Params.batch '_' savetime 'X'];

mkdir(savefolder);
cd(savefolder);

% params
Params.StimID=StimID;
Params.TargTrigID=TargTrigID;
Params.savefolder=savefolder;

% save
save('DatStructCompiled.mat','DatStructCompiled');
save('Params.mat','Params');