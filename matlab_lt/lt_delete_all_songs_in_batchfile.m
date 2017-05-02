function lt_delete_all_songs_in_batchfile(batch_to_delete)
% batch_to_delete= 'batch.08Apr2017_1901.dcrd';

fid=fopen(batch_to_delete);
song=fgetl(fid);

try cd('DELETED_SONGS');
    cd ..
catch err
    mkdir('DELETED_SONGS');
end
    

while ischar(song);
    disp(['moving ' song ' to DELETED_SONGS']);
% delete(song);
eval(['!mv ' song ' DELETED_SONGS']) 
song=fgetl(fid);
end


