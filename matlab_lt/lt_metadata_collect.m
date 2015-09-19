%% LT 3/13/15 - Creates data structure that contains metadata of all directories of a bird
% Data divided analogous to directories - first by experiment, then by
% index (corresponding to day) then by conditions performed on that day.
% Requires consistent naming of directories.

% Run this in folder containing all day directories

function MetadataStruct=lt_metadata_collect

%% Initiate struct
MetadataStruct=struct;

%% Collect all dirs in current directory that are in correct format
[~, cellofnames, ~, ~]=lt_write_all_folder_contents_to_batch_v2('*','','', 'dir');

%% Go through all filenames and extract metadata
for i=1:length(cellofnames);
    Expt=[];
    Condition=[];
    Notes=[];
    
    dirname=cellofnames{i};
    
    if isempty(dirname); % empty is likley because cell is forced to be in order of dates.
        continue
    end
    
    % === Check format
    if ~all(ismember(dirname(1:6),'0123456789'));
        disp(['PROBLEM, a dir is not formated correctly; name: ' dirname]);
        continue
    end
    
    
    % ==== Extract meta data - up to four levels (separated by up to 3
    % underscores)
    uscores=strfind(dirname,'_'); % underscores separate metadata types
    
    if length(uscores)==1;
        Expt=dirname(uscores(1)+1:end);
    elseif length(uscores)>1;
        Expt=dirname(uscores(1)+1:uscores(2)-1);
    end
    
    if length(uscores)==2;
        Condition=dirname(uscores(2)+1:end);
    elseif length(uscores)>2;
        Condition=dirname(uscores(2)+1:uscores(3)-1);
    end
    
    % should have minimum 3 uscores. has more if have notes
    if length(uscores)==3;
        Notes=dirname(uscores(3)+1:end);
    elseif length(uscores)>3;
        Notes=dirname(uscores(3)+1:uscores(4)-1);
    end
    
    
    % ==== Put those metadata into one structure
    MetadataStruct(i).experiment=Expt;
    MetadataStruct(i).condition=Condition;
    MetadataStruct(i).notes=Notes;
    
    datedate=dirname(1:6);
    datedate=datestr(datenum(datedate, 'mmddyy'),'ddmmmyyyy');
    MetadataStruct(i).date=datedate;
    MetadataStruct(i).dirname=dirname;
    
    %
    %
    %
    %
    %     if ~isfield(MetadataStruct,Expt);
    %         ind=1;
    %     else
    %         ind=length(MetadataStruct.(Expt))
    %     end
    %
    
    
    % ===== DO SOME STUFF THAT REQUIRES GOING INTO THE DIR
    
    cd(MetadataStruct(i).dirname);
    % == Does this day have labeled songs?
    
    if exist('batch.labeled.all', 'file');
        MetadataStruct(i).labeled_songs_exist=1;
    else
        MetadataStruct(i).labeled_songs_exist=0;
    end
    
    
    
    % == Does this day have autolabeled songs?
    tmp=[];
    try
        tmp=ls('*.AutoLabeled*'); % find batch files with this suffix
    catch err
    end
    
    if isempty(tmp); % then no autolabeled
        MetadataStruct(i).autolabeled_songs_exist=0;
    else
        MetadataStruct(i).autolabeled_songs_exist=1;
    end
    
    cd ..
    
end


