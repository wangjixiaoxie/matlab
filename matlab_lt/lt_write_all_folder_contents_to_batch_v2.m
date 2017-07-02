%% LT 11/14/15 - arrayofdateinds added (1 is first day)

%% LT 5/16/15 - modified to either extract files of format [phrase]_[date], or directions (of format date_experiment_[optional: conditions]);

% argument: 
% filetype = 'dir' (song folder) or 'mat'; (default is mat)


%% LT - 5/6/15 - additional argument, which is present forces only to collect filenames that are within that date range.
% Files you are looking for have to be in formate: [phrase]_[date] - e.g. AllData_05May2015
% e.g.
% TargStr='AllData_*'
% Dates_first='05May2015'
% Dates_last='10May2015'
% if leave out either first or last (i.e. leave empty as ''), will get all in that direction.
% would take files in format: AllData_05May2015 and only get if within date range.
% output cell will start from day 1, so might have empty cells if no filename there.
% output batch will not be in order.

% If don't specificy Dates, then will get all filenames.


%% LT - 1/31/15 - v2 - this takes any string and searches for files that match.
% e.g.
% TargStr='RawDatStruct*' - looks for anything with that string, asterix is
% usable

%% LT - 9/30/13 - takes folder contents (any with .mat in the name) and write on
% separate lines in a new batch file.  this is useful to use
% right before lt_db_plot_over_experiment
% outputs the batch file name and syllable.

function [batch, cellofnames, Dates_first, Dates_last, arrayofdateinds, arrayofdatenums]=lt_write_all_folder_contents_to_batch_v2(TargStr,Dates_first,Dates_last, FileType, make_batch)

%% Dafaults

if ~exist('FileType','var');
FileType='mat';
end

if ~exist('Dates_first','var');
% then make blank, will fill in later
Dates_first=[];
Dates_last=[];
end

if ~exist('make_batch','var');
    make_batch=1;
end


%%

folder_contents=dir(TargStr);

% == sort out only the types of files you want
FilesToThrowOut=[];

if strcmp(FileType,'dir')==1;
    % only want directories holding songs
    for i=1:length(folder_contents);
        if isdir(folder_contents(i).name) && length(folder_contents(i).name)>6;
            % then is potentially a song folder
            try
                cd(folder_contents(i).name);
                
                % look for presence of song files
                cbinfiles=dir('*.cbin');
                if isempty(cbinfiles);
                    % then this dir does not contain songs, throw it out
                    FilesToThrowOut=[FilesToThrowOut i];
                end
                cd('..');
            catch err
                % then is not dir
                FilesToThrowOut=[FilesToThrowOut i];
            end
        else
            % not a dir, throw out
            FilesToThrowOut=[FilesToThrowOut i];
        end
        disp(i);
    end
elseif strcmp(FileType,'mat');
    % only keep .mat files
    for i=1:length(folder_contents);
        fname=folder_contents(i).name;
        
        if ~any(findstr(fname,'.mat'));
            
            FilesToThrowOut=[FilesToThrowOut i];
        end
    end
else
    disp([FileType ' is not undestood, should be either "mat" or "dir"']);
end

% == Throw out unwanted files
folder_contents(FilesToThrowOut)=[];


%% sort by date
earliest_datenum=[];
latest_datenum=[];
date_of_file={};
for i = 1:length(folder_contents)
    % find earliest and latest file dates
    
    % how to find date based on file name?
    if strcmp(FileType, 'dir');
        
        tmp=datenum(folder_contents(i).name(1:6),'mmddyy');
        date_of_file{i}=datestr(tmp,'ddmmmyyyy'); % convert datestr format
        
        
    elseif strcmp(FileType,'mat');
        ind=length(TargStr);
        
        date_of_file{i}=folder_contents(i).name(ind:ind+8);
    end
    
    earliest_datenum=min([earliest_datenum datenum(date_of_file{i},'ddmmmyyyy')]);
    latest_datenum=max([latest_datenum datenum(date_of_file{i},'ddmmmyyyy')]);
end


% == Use those dates to sort the folder contents
FileDatenums=datenum(date_of_file, 'ddmmmyyyy');
[FileDatenums, inds]=sort(FileDatenums);

% sort date_of_file and folder_contents
date_of_file=date_of_file(inds);
folder_contents=folder_contents(inds);


%% == Extract only the folder contents within dates you want
cellofnames={};
arrayofdateinds=[];
arrayofdatenums = [];

for i = 1:length(folder_contents)
    % === Do you care about dates?
%     if exist('Dates_first','var');
        
        % then fill in missing dates with extreme valules if they are not specified whether the other date is specified
        
        if isempty(Dates_first);
            % fill with extreme value to get everything at that end
            Dates_first=datestr(earliest_datenum,'ddmmmyyyy');
        end
        
        if isempty(Dates_last);
            Dates_last=datestr(latest_datenum,'ddmmmyyyy');
        end
        
        
        % --- check if this file is within the date range
        % get date of this file
        if datenum(date_of_file{i},'ddmmmyyyy')>=datenum(Dates_first,'ddmmmyyyy') ...
                && datenum(date_of_file{i},'ddmmmyyyy')<=datenum(Dates_last,'ddmmmyyyy');
            
%             % what is the index relative to 1st day (1st day is 1);
%             DateInd=datenum(date_of_file{i},'ddmmmyyyy')-datenum(Dates_first,'ddmmmyyyy')+1;

            % == save to output cell
            cellofnames=[cellofnames folder_contents(i).name];
            arrayofdateinds=[arrayofdateinds datenum(date_of_file{i},'ddmmmyyyy')-datenum(Dates_first,'ddmmmyyyy')+1];
            arrayofdatenums = [arrayofdatenums datenum(date_of_file{i},'ddmmmyyyy')];
        end
        
%     else
%         % used for sorting later.
%         Dates_first=datestr(earliest_datenum,'ddmmmyyyy');
%         Dates_last=datestr(latest_datenum,'ddmmmyyyy');
%         
%         % otherwise just save all filenames
%         cellofnames{i}=folder_contents(i).name;
%     end
    
end

%% check below
if (0) % redundant - sorts above
% sort cell of names by date order
cellofnames_sorted={};
for i=1:length(cellofnames);
    if ~isempty(cellofnames{i});
        ind=length(TargStr)+1;
        date_of_file=cellofnames{i}(ind:ind+8);
        
    dayind=datenum(date_of_file,'ddmmmyyyy')-datenum(Dates_first,'ddmmmyyyy')+1;
    
    % make new cell array, sorted
    cellofnames_sorted{dayind}=cellofnames{i};
    end
end

cellofnames=cellofnames_sorted;
end


% ===== write to batch file
if make_batch==1;
current_date=now;
day=datestr(current_date,'ddmmmyyyy');
time=datestr(current_date,'HHMM');

batch=['batch.folder_contents.' day '_' time];

fid = fopen(['batch.folder_contents.' day '_' time] ,'w');

for i=1:length(cellofnames);
    if ~isempty(cellofnames{i});
        % -- write this name to batch
        
        fprintf(fid,'%s\n',cellofnames{i});
    end
end


fclose(fid);
else
    batch='';
end

%% backup
% function [batch, cellofnames]=lt_write_all_folder_contents_to_batch_v2(TargStr);
%
%     folder_contents=dir([TargStr]);
%
%     current_date=datenum(now);
%     day=datestr(current_date,'ddmmmyyyy');
%     time=datestr(current_date,'HHMM');
%
%     batch=['batch.folder_contents.' day '_' time];
%
%     fid = fopen(['batch.folder_contents.' day '_' time] ,'w');
%
%     for i = 1:length(folder_contents)
%         fprintf(fid,'%s\n',folder_contents(i).name);
%
%         cellofnames{i}=folder_contents(i).name;
%     end
%
%     fclose(fid);
% end