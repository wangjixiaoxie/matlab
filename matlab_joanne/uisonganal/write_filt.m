function  write_filt(filt_file, filtsong, Fs)

%write song data to filt_file
% file format:
% all data is written "big-endian"
% Fs: sample rate in kHz 1*(short)
% songdata: length(filtsong)*(short)


%check if file already exists, if so exit, otherwise open new file 
if exist(filt_file,'file');
     disp('write_filt: file already exists');
     return
else
     [fid, message] = fopen(filt_file, 'w', 'b');
end 

% if couldn't open output, exit
if fid == -1
    disp('write_filt: couldn''t open output file')
    disp(message)
    return
end

%write data
fwrite(fid, Fs, 'double');
fwrite(fid, filtsong ,'double'); 

fclose(fid);

