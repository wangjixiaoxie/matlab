batch_to_delete= 'batch.catch.keep.notrand';

fid=fopen(batch_to_delete);
song=fgetl(fid);

while ischar(song);
delete(song);
song=fgetl(fid);
end


