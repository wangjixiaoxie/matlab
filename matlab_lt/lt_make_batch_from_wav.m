function lt_make_batch(input)

batchname = 'batch';
% batchname = input('What is the name of the batch file?  ', 's');

%makes a batch file
fid = fopen(batchname,'w');

filenames = dir('*.wav');

for i = 1:length(filenames)
    fprintf(fid,'%s\n',filenames(i).name);
end

fclose(fid);

% ===
if input==1;
    cleandirAuto('batch',1000,4,4)
end
end
