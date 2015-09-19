function [ num_song ] = db_get_num_syl_batch( batchfile, syllable )
%db_get_num_syl_batch Input batch and syllable, gives you the number of
%that syllable in each song
%   Merely counts the number of labeled syllable in each song of a batch
%   file.

% Gets all syllables per song
syl_order = db_getlabels(batchfile,'n');

num_song = zeros(length(syl_order),1);
for ii  = 1:length(syl_order)
    num_song(ii) = sum(syl_order{ii} == syllable);
    
end


end

