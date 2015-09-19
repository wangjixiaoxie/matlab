function [] = db_seq_function_song_folders(data_type,motif,phrase,varargin)
%db_seq_function_song Goes through a bird folder and performs the function
%db_seq_func_day_save to each folder of 'data_type' with 'motif'
%   Input data_type (ex: 'amp'), motif (ex: {'ab' 'ac'}, and phrase 
%   ('' if you want all song folders, or specify a subset 'drug'). 
%   Can input up to three more variables: 
%       varargin{1} = varargin for db_seq_func_day_save function, if 'dur',
%       then specify {'begin_onset_offset' 'end_onset_offset'}. If 'amp'
%       then sepcify base syllable for ratio. NOTE: if you do this, run
%       'dur' and/or 'amp' first/second.
%       varargin{2} = directory (default is the current directory)  
%       varagin{3} = batchfile (default is 'batch.catch.keep')
%


%varargin{1} is the directory of interest. by default, it is the current
%directory (no input)
if length(varargin) > 1
    directory = varargin{2};
else
    directory = pwd;
end

%if data_type = 1, runs all seq_func_day
if data_type == 1
    data_type = {'dur' 'amp' 'frac' 'ent'};
else
end

%makes a list of song folders and a .mat with dates of folders
omit = 0;
db_list_song_folders(phrase,omit,directory);

%loads date file
load([directory '/song_folders_' phrase '.mat'])

%checks to see if motif is a cell, if not, converts it to a cell
if ischar(motif)
    motif = {motif};
end

%checks to see if data_type is a cell
if ischar(data_type)
    data_type = {data_type};
end

%gets batchfile (varargin{3} or 'batch.catch.keep')
if length(varargin) > 2
    batchfile = varargin{3};
else
    batchfile = 'batch.catch.keep';
end

%gets computer and birdname
slashes = strfind(directory,'/');
birdname = directory(slashes(end)+1:end);

%opens song_folders text file for name of song folders
fid = fopen([directory '/song_folders_' phrase '.txt'],'r');

%will go through each folder list in song_folders.txt and perform
%db_seq_func_day_save
for i = 1:length(song_folders)
    
    %gets the ith line of song_folders
    current_line = fgetl(fid);
    
    %changes the directory to that folder
    cd([directory '/' current_line])
    
    %loops through data_type
    for j = 1:max(size(data_type))
        %loops through motif
        for k = 1:max(size(motif))
            %displays what parameter and motif the program is analyzing
            display(['PARAMETER OF INTEREST: ' data_type{j} '    MOTIF: ' motif{k}])
            %runs this function to look at things like amplitude, entropy,
            %fraction of syllable sung, motif duration.
            if ~isempty(varargin) == 1
                try
                    %tries running with varargin{1}{j} if possible, else
                    %runs without varargin
                    db_seq_func_day_save(batchfile,directory,birdname,song_folders{i},data_type{j},motif{k},phrase,varargin{1}{j})
                catch err
                    db_seq_func_day_save(batchfile,directory,birdname,song_folders{i},data_type{j},motif{k},phrase)
                end
            else
                db_seq_func_day_save(batchfile,directory,birdname,song_folders{i},data_type{j},motif{k},phrase)
            end
        end
    end
    
    %changes back to original directory
    cd(directory)
end

%close song_folders.txt
fclose(fid);



end

