function  [spect, nfft, spect_win, noverlap, t_min, t_max, f_min, f_max] = read_spect(spect_file)

%read spectrograsm data from spect_file
% file format
% data is assumed to have been written by write_spect
% all data is written "big-endian"
% n_cols: number of columns 1*(short)
% n_rows: number of rows 1*(short)
% nfft: window length (in number of points)1*(short)
% spect_win: window used for spectrogram length(nfft)*double
% noverlap: number of overlapping points 1*(short)
% t_min: starting time of spectrogram 1*(double)
% t_max: ending time of spectrogram 1*(double)
% f_min: min frequency of spectrogram 1*(short)
% f_max: max frequency of spectrogram 1*(short)
% spectrogram data n_cols * n_rows *(short) 

%open file if it exists 

[fid, message] = fopen(spect_file, 'r', 'b');

% if couldn't open output, exit
if fid == -1
    disp('read_spect: couldn''t open song file')
    disp(message)
    return
end

%read data
n_rows = fread(fid, 1, 'short');
n_cols = fread(fid, 1, 'short');
nfft = fread(fid, 1, 'short');
spect_win = fread(fid, nfft, 'double');
noverlap = fread(fid, 1, 'short');
t_min = fread(fid, 1, 'double');
t_max = fread(fid, 1, 'double');
f_min = fread(fid, 1, 'short');
f_max = fread(fid, 1, 'short');
spect = fread(fid, [n_rows, n_cols], 'short');

fclose(fid);


