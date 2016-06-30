function [output_cell_or_struct]=lt_metadata_find_dirs(MetadataStruct, experiment, condition, notes, date_range, only_labeled_dirs, prog_vers)
%% LT 7/9/15 - v2 option
% activate by putting 2 in prog_vers argument. default is 1 (v1)
% version 2: outputs not just dir names, but also the parameters that
% Metadatsstruct puts in (e.g. experiment, condition, notes, date, ...) -
% outputs a structure, not a cell.

%% LT 5/19/15 - v1- after using lt_metadata_collect run this to sort out dirs that match your desired "experiment", "conditions", and/or "notes"

% MetadataStruct - output of lt_metadata_collect
% experiment - name of experiment (after 1st underscore)
% condition - phrase between 2nd and 3rd uscore
% notes - phrase btw 3rd and 4th uscore
% date_range={'05Nov2014', '10Nov2014'};
% only_labeled_dirs=1 (only takes dirs with either label or autolabeled songs) or 0 (takes all dirs) (default - takes all)

% leave any blank(''), or don't define, if you don't care what it is.

% output: cell containing dirnames matching your filter

% FOR EXAMPLE:
% MetadataStruct=lt_metadata_collect;
% 
% experiment = '';
% condition='';
% notes='';
% date_range={'20Apr2015','20May2015'};
% only_labeled_dirs=1;
% 
% ListOfDirs2=lt_metadata_find_dirs(MetadataStruct, experiment, condition, notes, date_range, only_labeled_dirs, 2);

%% DEFAULTS
if ~exist('prog_vers', 'var');
    prog_vers=1;
end

if ~exist('experiment', 'var');
    experiment='';
end

if ~exist('condition', 'var');
    condition='';
end

if ~exist('notes', 'var');
    notes='';
end

if ~exist('date_range', 'var');
    date_range={}; % collect all dates
end

if ~exist('only_labeled_dirs', 'var');
    only_labeled_dirs=0;
end

%% RUN
%% Version 1
if prog_vers==1;
    output_cell_or_struct={};
    for i=1:length(MetadataStruct);
        
        % skip if empty
        if isempty(MetadataStruct(i).dirname);
            continue
        end
        
        % == check if matches conditions desired, if fails anything, then continue to
        % next directory
        if ~isempty(experiment)
            if strcmp(MetadataStruct(i).experiment, experiment)==0;
                continue;
            end
        end
        
        if ~isempty(condition)
            if strcmp(MetadataStruct(i).condition, condition)==0;
                continue;
            end
        end
        
        if ~isempty(notes)
            if strcmp(MetadataStruct(i).notes, notes)==0;
                continue;
            end
        end
        
        if ~isempty(date_range);
            if datenum(MetadataStruct(i).date,'ddmmmyyyy')<datenum(date_range{1}) || ...
                    datenum(MetadataStruct(i).date,'ddmmmyyyy')>datenum(date_range{2}); % then is out of range
                continue
            end
        end
        
        % check if cotnains labeled song
        if only_labeled_dirs==1;
            if MetadataStruct(i).labeled_songs_exist==0 && MetadataStruct(i).autolabeled_songs_exist==0; % i.e. lacking both labeled and autolabeled songs
                continue
            end
        end
        
        output_cell_or_struct=[output_cell_or_struct MetadataStruct(i).dirname];
    end
    
elseif prog_vers==2;
    
    Inds_to_keep=[];
    for i=1:length(MetadataStruct);
        
        % skip if empty
        if isempty(MetadataStruct(i).dirname);
            continue
        end
        
        % == check if matches conditions desired, if fails anything, then continue to
        % next directory
        if ~isempty(experiment)
            if strcmp(MetadataStruct(i).experiment, experiment)==0;
                continue;
            end
        end
        
        if ~isempty(condition)
            if strcmp(MetadataStruct(i).condition, condition)==0;
                continue;
            end
        end
        
        if ~isempty(notes)
            if strcmp(MetadataStruct(i).notes, notes)==0;
                continue;
            end
        end
        
        if ~isempty(date_range);
            if datenum(MetadataStruct(i).date,'ddmmmyyyy')<datenum(date_range{1}) || ...
                    datenum(MetadataStruct(i).date,'ddmmmyyyy')>datenum(date_range{2}); % then is out of range
                continue
            end
        end
        
        % check if cotnains labeled song
        if only_labeled_dirs==1;
            if MetadataStruct(i).labeled_songs_exist==0 && MetadataStruct(i).autolabeled_songs_exist==0; % i.e. lacking both labeled and autolabeled songs
                continue
            end
        end
        
        Inds_to_keep=[Inds_to_keep i];
    end
    
    if ~isempty(Inds_to_keep);
        output_cell_or_struct=MetadataStruct(Inds_to_keep);
    else
        output_cell_or_struct=struct;
    end
end




