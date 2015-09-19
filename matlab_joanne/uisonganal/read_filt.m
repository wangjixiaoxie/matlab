function  [filtsong, Fs] = read_filt(filt_file)

%read songdata from filt_file
% file format:
% files are assumed to be written by write_filt
% all data is written "big-endian"
% Fs: sample rate (in kHz) 1*(short)
% filtsong: length(filtsong)*(short) 

%try to open file 
[fid, message] = fopen(filt_file, 'r', 'b');

% if couldn't open output, exit
if fid == -1
    disp('read_filt: couldn''t open song file')
    disp(message)
    return
end

%read data
Fs = fread(fid, 1, 'double');
filtsong = fread(fid, inf, 'double');

fclose(fid);


