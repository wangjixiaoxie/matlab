function [matrix_gap_timing, cell_gap_timing ] = db_gap_timing( batchfile, motif, span_of_timing, position_in_motif )
%db_gap_timing Finds the time between two syllables
%   Input a batch file (with contents either cbin or cbin.not.mat), a motif
%   (can be a regexp, i.e. 'ab', 'a.', 'a[^c]', '-[abc]'). Optional inputs
%   is the span of timing (for adjacent syllable, it is 1 (default)),
%   between 'a' and 'c' in 'abc' is 2, etc. Second optional input is
%   position_in_motif. Default it 0, for first position 'a' in 'abc', but
%   can be 1 for 'b' in 'abc', etc.

if nargin < 3
    position_in_motif = 0;
    span_of_timing = 1;
elseif nargin < 4
    position_in_motif = 0;
end

fid = fopen(batchfile,'r');
tline = fgetl(fid);

j = 1;
while ischar(tline)
    try
        if isempty(strfind('.not.mat',tline))
            load([tline '.not.mat'])
        else
            load(tline)
        end
        display(tline)
        
        offset_positions = regexp(labels, motif) + position_in_motif;
        onset_positions = offset_positions + span_of_timing;
        
        cell_gap_timing{j} =  onsets(onset_positions) - offsets(offset_positions);
        
        j = j+1;
    catch err
        display([tline ' does not have a .not.mat file'])
    end
    
    tline = fgetl(fid);
end

matrix_gap_timing = cell2mat(cell_gap_timing');

fclose(fid);

end

