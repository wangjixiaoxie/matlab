function lt_make_batch(input)

batchname = 'batch';
% batchname = input('What is the name of the batch file?  ', 's');

%makes a batch file
db_write_batch(batchname)

if input==1;
cleandirAuto('batch',1000,4,4)
end
end
