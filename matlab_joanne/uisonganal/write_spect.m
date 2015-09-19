function  write_spect(spect_file, spect, nfft, spect_win, noverlap, t_min, t_max, f_min, f_max)

%write spectrogram data to spect_file
% file format:
% all data is written "big-endian"
% n_rows: number of columns 1*(short)
% n_cols: number of rows 1*(short)
% nfft: window length (in number of points) 1*(short)
% spect_win: window used for spectrogram length(nfft)*double
% noverlap: number of overlapping points 1*(short)
% t_min: starting time of spectrogram 1*(double)
% t_max: ending time of spectrogram 1*(double)
% f_min: min frequency of spectrogram 1*(short)
% f_max: max frequency of spectrogram 1*(short)
% spectrogram data n_cols * n_rows *(short) 

%check if file already exists, if so exit, otherwise open new file 
if exist(spect_file,'file');
     disp('write_spect: file already exists')
     return
else 
    fid = fopen(spect_file, 'w', 'b');
end 

% if couldn't open output, exit
if fid == -1
    disp('write_spect: couldn''t open output file')
    disp(message)
    return
end

% get spect size
[n_rows, n_cols] = size(spect);

% check that spect_win is indeed same as length of spect window
if length(spect_win) ~= nfft
    disp('write_spect: length of spectrogram window ~= nfft')
    fclose(fid);
    return
end

%write data
fwrite(fid, n_rows, 'short');
fwrite(fid, n_cols, 'short'); 
fwrite(fid, nfft,'short'); 
fwrite(fid, spect_win, 'double');
fwrite(fid, noverlap,'short'); 
fwrite(fid, t_min,'double'); 
fwrite(fid, t_max,'double'); 
fwrite(fid, f_min,'short'); 
fwrite(fid, f_max,'short'); 
fwrite(fid, spect,'short'); 

fclose(fid);

