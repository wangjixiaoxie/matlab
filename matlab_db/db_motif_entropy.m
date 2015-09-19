function [ syl_we, absolute_time ] = db_motif_entropy( batchfile, motif)
%db_motif_entropy Finds the wiener entropy of every syllable in motif
%   Input a batch file, and it will calculate the average WE for each
%   syllable in the motif for the entire duration of the syllable.

%% Preparing the variables

%Calculates the number of songs to make a cell for timing
[s,song_number] = system(['wc -l ' batchfile]);
clear s
song_number = str2double(song_number(1:strfind(song_number, batchfile)-2));

%Makes a cell that will have when the motif was sung
for j = 1:length(motif)
    absolute_time.(motif(j)) = cell(song_number,1);
    syl_duration.(motif(j)) = cell(song_number,1);
end


%% Calculating the duration of the motif
% Opening the cbin_not_mat file
fid = fopen(batchfile, 'r');

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
        if isempty(strfind(batchfile, '.not.mat'))
            [spec, time] = db_spec_syllable(batchfile,motif(j),.2,.2,motif(1:j-1),motif(j+1:end));
        elseif ~isempty(strfind(batchfile, '.not.mat'))
            [spec, time] = db_spec_syllable(batchfile(1:end-8),motif(j),.2,.2,motif(1:j-1),motif(j+1:end));
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

