function [ boot_percentiles boot_frac fraction_per_song] = db_fraction_syl( batchfile, syl)
%db_fraction_syl Calculates the fraction of the syllable of interest over
%the total number of syllables sung (as determined by segmentation).
%   Will consider a file a song only if there is at least one label in the
%   file, so make sure to do so. Will pool all the syllables together then
%   find the fraction and bootstrap 10,000 times to get a distribution of
%   the fraction. Will also make a cell keeping the fraction per song.

%Calculates the number of songs to make a cell for timing
[s,song_number] = system(['wc -l ' batchfile]);
clear s
song_number = str2double(song_number(1:strfind(song_number, batchfile)-2));

%Makes a cell that will have when the motif was sung
fraction_per_song = cell(song_number,1);
pool_notes = [];

%% Calculating the fraction of the syllable
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
    
    if sum(labels ~= '-') >= 1
        fraction_per_song{i} = (length(strfind(labels,syl))+length(syl)-1)./length(labels);
        pool_notes = [pool_notes labels];
    end
    
    next_line = fgetl(fid);
    i = i+1;
end

%gets rids of empty cells for fraction_per_song and converts it to a matric
fraction_per_song = fraction_per_song(~cellfun('isempty',fraction_per_song));
fraction_per_song = cell2mat(fraction_per_song);

%bootstrap procedure to find distribution of fraction from pool
num_bootstrap = 10000;
boot_frac = cell(num_bootstrap,1);

if length(syl) == 1
    for j = 1:num_bootstrap
        sampling = randi(length(pool_notes),[1 length(pool_notes)]);
        boot_run = pool_notes(sampling);
        boot_frac{j} = length(strfind(boot_run,syl))./length(boot_run);
    end
else
    for j = 1:num_bootstrap
        sampling = randi(length(fraction_per_song),[1 length(fraction_per_song)]);
        boot_run = fraction_per_song(sampling);
        boot_frac{j} = mean(boot_run);
    end
end


boot_frac = cell2mat(boot_frac);

% 2.5%, 50%, and 97.5% percentiles for boot_frac
boot_percentiles = [prctile(boot_frac,2.5) prctile(boot_frac,50) prctile(boot_frac,97.5)];




    

end

