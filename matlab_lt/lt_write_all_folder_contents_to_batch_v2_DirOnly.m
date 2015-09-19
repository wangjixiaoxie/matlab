%% LT 3/13/15 - Modified from v2. only outputs directories (throws out files).
% savefile = 0; will not save batch file; 1 will save

%% LT - 1/31/15 - v2 - this takes any string and searches for files that match.
% e.g.
% TargStr='RawDatStruct*' - looks for anything with that string, asterix is
% usable
    
%% LT - 9/30/13 - takes folder contents (any with .mat in the name) and write on 
% separate lines in a new batch file.  this is useful to use
% right before lt_db_plot_over_experiment
% outputs the batch file name and syllable.

function [batch, OutNames]=lt_write_all_folder_contents_to_batch_v2(TargStr,savefile)

    folder_contents=dir([TargStr]);
        
    % Initiate output cell array
    OutNames={};
    c=1;
    
    % Initiate save file
    current_date=now;
    day=datestr(current_date,'ddmmmyyyy');
    time=datestr(current_date,'HHMM');
    
    batch=['batch.folder_contents.' day '_' time];
    fid = fopen(['batch.folder_contents.' day '_' time] ,'w');
    
    % collect names
    for i = 1:length(folder_contents)
        if folder_contents(i).isdir==1;
            if ~strcmp(folder_contents(i).name,'.') && ~strcmp(folder_contents(i).name,'..'); % don't want '.' and '..'
            fprintf(fid,'%s\n',folder_contents(i).name);
            
            OutNames{c}=folder_contents(i).name;
            c=c+1;
            end
        end
    end
    
    fclose(fid);
    
    if savefile==0; % then don't save
        eval(['!rm ' batch]);
    end
end