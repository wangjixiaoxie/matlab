function [AllData, error_NoNoteGroup]=lt_context_ExtractDayData(Params, saveON, OnlyContinueIfNoteGroups)
%% LT 4/26/15 - Extract single day data for context shift experiments

% INSTRUCTIONS:
% Run this in the day folder. Must have created a batch with songs you want
% to analyze. Must have note_group_log.ltrec2 files that have note group
% information for all songs in that batch.  Could be multiple note group
% logs, this code will automatically go through all log files.  If missing
% any song from note group logs, this code will still run, but those songs
% will miss NG data.  Can't have one song in multiple note group logs, this
% code will also skip those songs. But that is unlikely to occur.

% Inputs e.g.:
% Params.batch = 'batch.catch.keep';
% Params.config=
% '/bluejay4/lucas/birds/rd23gr89/config_042115_fixedfreqbounds.evconfig2';
% Params.exptID = 'light_context'; ID that identifies the data used in this experiment. if
% left blank (''), will use batch name
% saveON  = 1, saves in birdfolder, all days in one folder. if 0, then no
% save
% OnlyContinueIfNoteGroups = 1, then will ask user whether continue if finds that song is
% missing from note group. default = 0;

% Outputs;
% AllData
% error_NoNoteGroup - 1 if failed to find notegroup. 0 if found (outputs
% whether or not OnlyContinueIfNoteGroups is 1)

% DESCRIPTION: takes an evtafv4 config file. it will be used to run
% evtafsim on all songs in the batch file.  It will output a data
% structure containing offline detects and online hits and metadata for
% each song in batch. That strucutre is also saved in a subfolder. It will
% perform evtafsim for all notes in the config file (even if all you care
% about is a single note). Takes me about a second per song file.


%% PARAMS
if ~exist('OnlyContinueIfNoteGroups','var');
    OnlyContinueIfNoteGroups=0;
end

if ~exist('saveON','var');
    saveON=1;
end


% Extract info about this day;
if isfield(Params, 'expt_phrase');
    [Params.birdname, Params.bluejaynum, date, ~]=lt_get_birdname_date_from_dir(1);
else
    [Params.birdname, Params.bluejaynum, date, Params.expt_phrase]=lt_get_birdname_date_from_dir(1);
end

Params.date=date{2}; % date in string form;
Params.day_directory=pwd;
Params.timestamp_analysis=lt_get_timestamp(0);


%% First check to make sure this day folder contains note group log.

% -- COllect note group data
disp('Collecting note group data...')

% get list of notegroup files in current day folder
NoteGroupFiles=dir('*.ltrec2');

% Collect note group log data from multiple files into one structure.
AllNoteGroupData=struct('SongFname_cbin',[],'NotesInGroup',[], 'CurrNoteGroup',[]);
for i=1:length(NoteGroupFiles);
    disp(['found ' NoteGroupFiles(i).name]);
    
    tmp=lt_read_ltrec2(NoteGroupFiles(i).name,0);
    
    % slot into global structure
    AllNoteGroupData.SongFname_cbin=[AllNoteGroupData.SongFname_cbin tmp.SongFname_cbin];
    AllNoteGroupData.NotesInGroup=[AllNoteGroupData.NotesInGroup tmp.NotesInGroup];
    AllNoteGroupData.CurrNoteGroup=[AllNoteGroupData.CurrNoteGroup tmp.CurrNoteGroup];
end


% -- Check that all songs in back are present in the log
disp(' ')
fid=fopen(Params.batch);
songname=fgetl(fid);
c=1;
check=0;
while ischar(songname);
    
    if ~any(strcmp(AllNoteGroupData.SongFname_cbin,songname))
        % then this song does not have note group log
        check=1;
        disp(['PROBLEM: song num ' num2str(c) ' ('  songname ') is missing from note group logs (skipping)']);
    end
    
    % update song
    c=c+1;
    songname=fgetl(fid);
end
disp(' ');
disp(' ');

% --- output error if no missing notegroup
error_NoNoteGroup=check;

% --- detected any missing song from note group log?
if check==1 && OnlyContinueIfNoteGroups==1;
    continue_input=input('found missing song from note group log. continue? (y or n) ', 's');
    
    if continue_input=='n';
        dafdsf; % this causes error
    end
end


%% Run simulation of evtaf to get hit, detect, song metadata, info
disp('Running evtaf sim...')
AllData=EvTAFv4Sim_LT_temp(Params.batch,Params.config,'obs0');

disp(['Found ' num2str(length(AllData)) ' songs today']);

%% Collect note group information from note group log/s and put into AllData structure
disp('Collecting note group data...')

% get list of notegroup files in current day folder
NoteGroupFiles=dir('*.ltrec2');

% Collect note group log data from multiple files into one structure.
AllNoteGroupData=struct('SongFname_cbin',[],'NotesInGroup',[], 'CurrNoteGroup',[]);
for i=1:length(NoteGroupFiles);
    disp(['found ' NoteGroupFiles(i).name]);
    
    tmp=lt_read_ltrec2(NoteGroupFiles(i).name,0);
    
    % slot into global structure
    AllNoteGroupData.SongFname_cbin=[AllNoteGroupData.SongFname_cbin tmp.SongFname_cbin];
    AllNoteGroupData.NotesInGroup=[AllNoteGroupData.NotesInGroup tmp.NotesInGroup];
    AllNoteGroupData.CurrNoteGroup=[AllNoteGroupData.CurrNoteGroup tmp.CurrNoteGroup];
end

% For each song, look through note group data to ask what note group it was
for i=1:length(AllData);
    
    ind=find(strcmp(AllNoteGroupData.SongFname_cbin,AllData(i).filename));
    
    if isempty(ind);
        % then this song has data, but is not annotated in note groups.
        % tell user
        disp(['Problem: song num ' num2str(i) ' ('  AllData(i).filename ') is missing from note group logs (skipping)']);
        continue
        
    elseif length(ind)>1;
        % then for some reason this song has been logged in multiple note
        % group logs. tell user:
        disp(['Problem: song num ' num2str(i) ' ('  AllData(i).filename ') is annotated in >1 note group logs (skipping)']);
        continue
        
    else
        % add data to AllData
        AllData(i).NoteGroup=AllNoteGroupData.CurrNoteGroup(ind);
        AllData(i).NotesInNoteGroup=AllNoteGroupData.NotesInGroup(ind);
    end    
end


%% SAVE 
if saveON==1;
disp('Saving...');

% find or make dir;
cd ..
try 
cd('context_data');
catch err
    mkdir('context_data');
cd('context_data');
end

% make subdir for this experiment
try 
    cd(Params.expt_phrase);
catch err
    mkdir(Params.expt_phrase);
    cd(Params.expt_phrase);
end

% save day's data
save(['Params_' Params.date],'Params');
save(['AllData_' Params.date],'AllData');

end

disp('DONE!');
















