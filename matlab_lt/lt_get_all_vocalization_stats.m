%% LT 12/9/13 - use this to look at all vocalizations (i.e. those in labeled songs that have ampl. cross detection threshold).
function [durations_all_songs, durations_per_song]=lt_get_all_vocalization_stats(cbin_not_mat_file)

% MODIFIED FROM db_motif_timing

%% Preparing the variables

%Calculates the number of songs to make a cell for timing
% 
% onsets
% offsets
% durations_per_song
% 
% cbin_not_mat_file= 'batch.catch.keep';
%% Calculating the durations
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
    
    onsets_per_song{i}=onsets;
    offsets_per_song{i}=offsets;
    
    
    i=i+1;
    next_line=fgetl(fid);
end

fclose(fid);

% get durations from onset and offsets, per song
for i=1:length(onsets_per_song);
    durations_per_song{i}=offsets_per_song{i}-onsets_per_song{i};
end

% get data over all songs
durations_all_songs=cell2mat(durations_per_song');
onsets_all_songs=cell2mat(onsets_per_song');
offsets_all_songs=cell2mat(offsets_per_song');