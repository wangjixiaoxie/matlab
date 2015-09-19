function [ motif_timing, absolute_time ] = lt_db_motif_timing(cbin_not_mat_file, motif, begin_on_or_offset, end_on_or_offset )
%db_motif_timing Will find the length of time for a motif. Input file,
%motif, and whether to begin with the offset or onset and end with the
%offset or onset.
%   Will go through each cbin or cbin.not.mat file find the motif of interest. Will
%   then find the start time of the motif (either the onset or offset of
%   the first syllable) and the end time (either the onset or offset of the
%   last syllable).

%% Preparing the variables

%Calculates the number of songs to make a cell for timing
[s,song_number] = system(['wc -l ' cbin_not_mat_file]);
clear s
song_number = str2double(song_number(1:strfind(song_number, cbin_not_mat_file)-2));

%Makes a cell that will have all the timing
motif_timing = cell(song_number,1);
%Makes a cell that will have when the motif was sung
absolute_time = cell(song_number,1);

%Calculates length of motif
motif_length = length(motif)-1;


%% Calculating the duration of the motif
% Opening the cbin_not_mat file
fid = fopen(cbin_not_mat_file, 'r');

%Gets the first line of the input file
next_line = fgetl(fid);
i = 1;

while ischar(next_line)
    %Checks to see if file is cbin or cbin.not.mat and loads the
    %cbin.not.mat file
    if isempty(strfind(next_line, '.not.mat'))
        load([next_line '.not.mat'])
    elseif ~isempty(strfind(next_line, '.not.mat'))
        load(next_line)
    end
    
    
    %finds the timing of the motif
    occurence_motif = strfind(labels, motif);
    
    if isempty(occurence_motif)
        %if there is no motif in the song, it will say NaN
        motif_timing{i} = NaN;
        absolute_time{i} = NaN;
    else
        %Calculates absolute time of motif
        absolute_time{i} = db_timing4_file(next_line) + onsets(occurence_motif).*(1/1000).*(1/86400);
        
        %otherwise it will calculate the length of the motif
        if strcmpi(begin_on_or_offset, 'on') && strcmpi(end_on_or_offset, 'on')
            motif_timing{i} =  onsets(occurence_motif + motif_length) - onsets(occurence_motif);
        elseif strcmpi(begin_on_or_offset, 'on') && strcmpi(end_on_or_offset, 'off')
            motif_timing{i} =  offsets(occurence_motif + motif_length) - onsets(occurence_motif);
        elseif strcmpi(begin_on_or_offset, 'off') && strcmpi(end_on_or_offset, 'on')    
            motif_timing{i} =  onsets(occurence_motif + motif_length) - offsets(occurence_motif);
        elseif strcmpi(begin_on_or_offset, 'off') && strcmpi(end_on_or_offset, 'off') 
            motif_timing{i} =  offsets(occurence_motif + motif_length) - offsets(occurence_motif);
        end
    end
    
    i = i+1;
    next_line = fgetl(fid);
end

fclose(fid);

%Converts motif_timing to a matrix
motif_timing = cell2mat(motif_timing);
absolute_time = cell2mat(absolute_time);

%Gets rid of NaN
motif_timing = motif_timing(~isnan(motif_timing));
absolute_time = absolute_time(~isnan(absolute_time));

end

