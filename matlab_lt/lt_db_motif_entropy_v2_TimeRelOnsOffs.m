function [ syl_we absolute_time ] = lt_db_motif_entropy_v2( cbin_not_mat_file, motif, pretime, posttime)
% LT 1/16/14 - modified so that the time mean of entropy is taken BEFORE
% the 

% LT 12/13 - MODIFIED to go back to db_spec_syllable, since the change made
% on 12/12 needs much more to work (i.e. now every spec gram I get might
% have different number of time bins, and so can't concatenate into a
% matrix. 

% LT-12/12 - MODIFIED line 100 to use lt_db_spec_syllable.  
% difference is that now syl is defined as a fixed distance before offset
% (as opossed to relative to onset). 
% USE positive number for posttime (which will be subtracted from offset). 

% LT-12/02 - to allow argument for pre- and post-time for gettign
% spectrogram. (in seconds)

%MODIFIED BY LT 082313 TO ADD CATCH TO LINE 37 AND 45 IN CASE NOT ALLL
%SONGS IN TEH BATCH FILE HAVE EXISTING NOTMAT FILES.

%db_motif_entropy Finds the wiener entropy of every syllable in motif
%   Input a batch file, and it will calculate the average WE for each
%   syllable in the motif for the entire duration of the syllable.

%% Preparing the variables

%Calculates the number of songs to make a cell for timing
[s,song_number] = system(['wc -l ' cbin_not_mat_file]);
clear s
song_number = str2double(song_number(1:strfind(song_number, cbin_not_mat_file)-2));

%Makes a cell that will have when the motif was sung
for j = 1:length(motif)
    absolute_time.(motif(j)) = cell(song_number,1);
    syl_duration.(motif(j)) = cell(song_number,1);
end


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
        try
        load([next_line '.not.mat'])
        catch
            next_line=fgetl(fid);
        continue
        end
        
    elseif ~isempty(strfind(next_line, '.not.mat'))
        try
            load(next_line)
        catch
            next_line=fgetl(fid);
            continue
        end
        
    end

    %finds the timing of the motif
    occurence_motif = strfind(labels, motif);
    
    for j = 1:length(motif)
        if isempty(occurence_motif)
            %if there is no motif in the song, it will say ''
            syl_duration.(motif(j)){i} = '';
            absolute_time.(motif(j)){i} = '';
        else
            %Calculates absolute time of motif
            absolute_time.(motif(j)){i} = db_timing4_file(next_line) + onsets(occurence_motif+j-1).*(1/1000).*(1/86400);
            
            %Calculates the syllable duration
            syl_duration.(motif(j)){i} = offsets(occurence_motif+j-1)-onsets(occurence_motif+j-1);
        end
    end

    
    i = i+1;
    next_line = fgetl(fid);
end
fclose(fid);


for j = 1:length(motif)
    %Gets rid of all empty cells and converts to a matrix
    absolute_time.(motif(j)) = absolute_time.(motif(j))(~cellfun('isempty',absolute_time.(motif(j))));
    absolute_time.(motif(j)) = cell2mat(absolute_time.(motif(j)));
    syl_duration.(motif(j)) = syl_duration.(motif(j))(~cellfun('isempty',syl_duration.(motif(j))));
    syl_duration.(motif(j)) = cell2mat(syl_duration.(motif(j)));
    
    if isempty(absolute_time.(motif(j))) == 1
        if j == 1
            display(['the motif "' motif '" did not occur on this day'])
        else
        end
        syl_we.(motif(j)) = [];
    else
        %gets the spectograms of syllable in motif for that day
        if isempty(strfind(cbin_not_mat_file, '.not.mat'))
            [spec time] = db_spec_syllable(cbin_not_mat_file,motif(j),pretime,posttime,motif(1:j-1),motif(j+1:end));
%             [spec time] = lt_db_spec_syllable(cbin_not_mat_file,motif(j),pretime,posttime,motif(1:j-1),motif(j+1:end)); % LT added 12/12/13
        elseif ~isempty(strfind(cbin_not_mat_file, '.not.mat'))
            [spec time] = db_spec_syllable(cbin_not_mat_file(1:end-8),motif(j),pretime,posttime,motif(1:j-1),motif(j+1:end));
%             [spec time] = lt_db_spec_syllable(cbin_not_mat_file(1:end-8),motif(j),pretime,posttime,motif(1:j-1),motif(j+1:end)); % LT added 12/12

        end
        
        syl_we.(motif(j)) = zeros(length(syl_duration.(motif(j))),1);
        for k = 1:size(spec,3)
            %calculates the wiener entropy for each rendition of the syllable
            [wiener_entropy] = db_calc_wiener_entropy( spec(:,:,k) );
            
            %keeps the wiener entropy values only when the syllable was sung
            wiener_entropy = wiener_entropy(time >= 0 & time <= syl_duration.(motif(j))(k));
            
            %averages the wiener_entropy over the syllable
            syl_we.(motif(j))(k) = mean(wiener_entropy);
        end
    end
end
    
    
    
end

