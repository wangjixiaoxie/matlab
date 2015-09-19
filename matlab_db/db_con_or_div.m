function [ con_or_div modified_motifs] = db_con_or_div( motifs )
%db_con_or_div Given a cell of motifs, says whether it is a convergent or
%divergent sequence and gives you the modified motifs (for 'div' just the
%last syllable in the motifs, for 'con' just the first syllable in the
%motifs).
%   If the last syllable of every motif is the same, it is a convergent
%   ('con') sequence. If the first syllable of every motif is the same, it
%   is a divergent ('div') sequence.


for i = 1:length(motifs)
    %creates an array with the first syllable of the motifs
    first_syl(i) = motifs{i}(1);
    
    %creates an array with the last syllable of motifs
    last_syl(i) = motifs{i}(end);
    
end

%checks to see if divergent or convergent by seeing if all first_syl or
%last_syl are the same
if length(strfind(first_syl,first_syl(1))) == length(first_syl)
    con_or_div = 'div';
    
    %creates a modified motif eliminating initial syllable to be used
    %when bootstrapping to calculate trans probabilities
    for jj = 1:length(motifs)
        modified_motifs{jj} = motifs{jj}(end);
    end
    
elseif length(strfind(last_syl,last_syl(1))) == length(last_syl)
    con_or_div = 'con';
    
    %creates a modified motif eliminating all but first syllable to be used
    %when bootstrapping to calculate trans probabilities
    for jj = 1:length(motifs)
        modified_motifs{jj} = motifs{jj}(1);
    end
    
else
    %breaks the program if the motifs do not match up
    display('Input motifs do not reflect a sequence transition')
    return
end

end

