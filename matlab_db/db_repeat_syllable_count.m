function [ matrix_repeat_length, cell_repeat_length ] = db_repeat_syllable_count( batchfile, motif, not_repeat_syllable, rep_position )
%db_repeat_syllable_count Counts the number of repeats for a syllable in a
%specified motif
%   db_repeat_syllable_count(batchfile, motif, repeat_syllable, rep_position). Ex:
%   db_repeat_syllable_count('batch', '[acd]b', '[^b]') will count the
%   number of 'b' following a, c, or d in a song
%   db_repeat_syllable_count('batch', '[^b]b','[^b]') will count the number
%   of 'b' following any character but 'b' in a song
%   db_repeat_syllable_count('batch','[ab][cd]','[^cd]',1) will count the
%   number of c or d that follow a or b. That the repeating syllable is 1
%   position from the start of the motif.

if nargin < 4
    rep_position = length(regexp(motif,'\w'))-1;
end

fid = fopen(batchfile,'r');
tline = fgetl(fid);



j = 1;
while ischar(tline)
    try
        if isempty(strfind(tline, '.not.mat'))
            load([tline '.not.mat'])
        else
            load tline
        end
        
        motif_positions = regexp(labels, motif);
        
        not_repeat_syl = regexp(labels, not_repeat_syllable);
        
        if not_repeat_syl(end) ~= length(labels)
            not_repeat_syl(end+1) = length(labels)+1;
        end
        
        for i = 1:length(motif_positions)
            next_not_syl = not_repeat_syl(motif_positions(i) < not_repeat_syl);
            next_not_syl = next_not_syl(1);
            cell_repeat_length{j}(i,1) = next_not_syl - motif_positions(i) - rep_position;
        end
%         cell_repeat_length = cell_repeat_length';
        j = j+1;
        display(tline)
    catch err
        display([tline ' does not have a .not.mat file'])
    end
    
    tline = fgetl(fid);
end

matrix_repeat_length = cell2mat(cell_repeat_length');
        
fclose(fid);
end

