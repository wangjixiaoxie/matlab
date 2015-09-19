function [number_of_songs, syl_rendition_amount]=lt_calc_syl_rendition_amount(batchfile,syllables)
    
%% gets nunber of songs labeled per day and the number of renditions of each
% specified syllable in that day.

% Opening the cbin_not_mat file
fid = fopen(batchfile, 'r');

%Gets the first line of the input file
next_line = fgetl(fid);
i = 1;

% variables
number_of_songs =0;
syl_rendition_amount=zeros(length(syllables),1);
    

while ischar(next_line)
    %Checks to see if file is cbin or cbin.not.mat and loads the
    %cbin.not.mat file
        if isempty(strfind(next_line, '.not.mat'))
        load([next_line '.not.mat'])
        elseif ~isempty(strfind(next_line, '.not.mat'))
        load(next_line)
        end

        if sum(labels ~= '-') >= 1
            for k=1:length(syllables);
            temp_syl_amount=length(strfind(labels,syllables(k)));
            syl_rendition_amount(k) = syl_rendition_amount(k)+temp_syl_amount; % add syl renditions to tally
            clear temp_syl_amount % add 1 to song tally
            end
            number_of_songs=number_of_songs+1;
        end

next_line = fgetl(fid);
i = i+1;
end





