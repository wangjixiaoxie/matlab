function RecDatStruct=lt_read_ltrec2(fname_ltrec,makebatch,getrecFB)
%% LT 6/5/15 - added default option - run with no arguments
% looks for all .ltrec2 files files in dir and compiles data for all.

%% LT 3/16/15 - Loads .ltrec2 files (they are logs of note groups in EvTAFv4_LT).
% .ltrec2 has one line for each song - it logs what note group was active
% for that song. This is one log file that conitnously appends.

% Instructions: run in same directory containing the .ltrec2 file.  This will
% save the output structure in same directory.

% inputs example:
% fname_ltrec='note_group_log.ltrec2';
% makebatch=1; % will make a batch separating songs based on what note
% group they were; if 0, then will not do so (batch will only contain songs
% that are actual .cbin files in this dir)
% getrecFB=1; % processes those batches to get songs wtih FB in rec file.



%%
RunDefaults=0;
if ~exist('fname_ltrec','var');
    RunDefaults=1;
end

if ~exist('makebatch','var');
    makebatch=1;
end

if ~exist('getrecFB','var');
    getrecFB=1;
end


%% Extract rec data into cell array

disp('Extracting .ltrec2 data into data structure...');

RecDatCell={[],[],[]};
if RunDefaults==1;
    % === FIND ALL .ltrec2 files in dir
    NoteGroupFiles=dir('*.ltrec2');
    
    % == extract data for each one
    for i=1:length(NoteGroupFiles);
        % for each ltrec2 file, extract info to cell array
        disp(['found ' NoteGroupFiles(i).name]);
        
        fname_ltrec=NoteGroupFiles(i).name;
        fid1=fopen(fname_ltrec);
        
        tmp=textscan(fid1,'%s %s %s');
        RecDatCell{1}=[RecDatCell{1}; tmp{1}];
        RecDatCell{2}=[RecDatCell{2}; tmp{2}];
        RecDatCell{3}=[RecDatCell{3}; tmp{3}];
    end
    
elseif RunDefaults==0; % then open file entered by user
    fid1=fopen(fname_ltrec);
    RecDatCell=textscan(fid1,'%s %s %s');
end

% NOTE: Output is RecDatCell, contains info for all songs.

%% Put different rec data into separate fields of output struct
RecDatStruct=struct;

for i=1:length(RecDatCell{1});
    % File names
    tmp=strfind(RecDatCell{1}{i},'\'); % want just cbin name, without path
    RecDatStruct.SongFname_cbin{i}=RecDatCell{1}{i}(tmp(end)+1:end);
    
    % Note group
    currnotegroup=str2double(RecDatCell{2}{i}(14));
    
    RecDatStruct.CurrNoteGroup(i)=currnotegroup;
    
    % Members of note group
    tmp=strfind(RecDatCell{3}{i},'_');
    
    for ii=1:length(tmp);
        RecDatStruct.NotesInGroup{i}(ii)=str2double(RecDatCell{3}{i}(tmp(ii)+1));
    end
end


%% SAVE

if RunDefaults==1;
    savename='NoteGroupLog_all_structure'; % take name string excluding '.ltrec' part
    save(savename, 'RecDatStruct');
    
elseif RunDefaults==0; % then file name was specified
    savename=[fname_ltrec(1:end-7) '_structure']; % take name string excluding '.ltrec' part
    save(savename, 'RecDatStruct');

    disp(['savename: ' savename]);
end

disp('Done, saved!')

%% OPTIONAL - make batch

    % then make one batch for each note group log
    if makebatch==1;
        
        % get list of songs in current directory - will only put in batch the
        % songs in the dir
        SongsInDir=dir('*.cbin');
        
        % MAKE BATCH FILE
        AllNoteGroups=unique(RecDatStruct.CurrNoteGroup); % note groups
        NumNoteGroups=length(AllNoteGroups);
        
        % open batch files
        for i=1:NumNoteGroups;
            fid(i)=fopen(['batch.NoteGroup_' num2str(AllNoteGroups(i))],'w');
        end
        
        % go through all songs and see which note group it fits into
        for i=1:length(RecDatStruct.SongFname_cbin);
            
            if ~isempty(findstr([SongsInDir.name],RecDatStruct.SongFname_cbin{i})); % true if this song exists in todays folder
                
                ind=find(AllNoteGroups==RecDatStruct.CurrNoteGroup(i)); % Index of fid
                
                fprintf(fid(ind),'%s\n',RecDatStruct.SongFname_cbin{i});
            end
        end
        
        % close all batch
        for i=1:NumNoteGroups;
            fclose(fid(i));
        end
        
        
        
        
        
        %% OF THOSE BATCHES, GET SONGS THAT HAVE REC WITH FB (I.E. REAL SONGS)
        if getrecFB==1;
            for i=1:NumNoteGroups;
                try
                    batchf=['batch.NoteGroup_' num2str(AllNoteGroups(i))];
                    lt_rec_files_find_FB_v3(batchf);
                    close all;
                catch err
                    continue
                end
            end
        end
        
    end



