% LT - 9/30/13 - takes folder contents (any with .mat in the name) and write on 
% separate lines in a new batch file.  this is useful to use
% right before lt_db_plot_over_experiment
% outputs the batch file name and syllable.

function [batch]=lt_write_all_folder_contents_to_batch(use_syl,syllable)

if use_syl==1
    folder_contents=dir(['*_' syllable '.mat*']);
    
    current_date=datenum(now);
    day=datestr(current_date,'ddmmmyyyy');
    time=datestr(current_date,'HHMM');
    
    batch=['batch.folder_contents.' day '_' time];
    
    fid = fopen(['batch.folder_contents.' day '_' time] ,'w');
    
    for i = 1:length(folder_contents)
        fprintf(fid,'%s\n',folder_contents(i).name);
    end
    
    fclose(fid);
    
else
    folder_contents=dir('*.mat');
    
    current_date=datenum(now);
    day=datestr(current_date,'ddmmmyyyy');
    time=datestr(current_date,'HHMM');
    
    batch=['batch.folder_contents.' day '_' time];
    
    fid = fopen(['batch.folder_contents.' day '_' time] ,'w');
    
    for i = 1:length(folder_contents)
        fprintf(fid,'%s\n',folder_contents(i).name);
    end
    
    fclose(fid);
end