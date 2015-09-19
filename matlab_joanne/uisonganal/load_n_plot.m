function load_n_plot(smwin,filt)

% function load_n_plot(smwin,filt)
% 
% This function calles loadraw and plotsongs.  It also plots the
% local loudness in db (where the local loudness is calculated with
% a smoothing window of SMWIN milliseconds. FILT is a file filter
% string  used
% by loadraw.

[song,fs,name] = loadraw(filt);
lvarsq = short_time_lvar(song,fs,smwin);

loglvarsq = log_lvar(lvarsq);
song = (song-mean(song))/std(song);
loglvarsq = (loglvarsq-mean(loglvarsq))/std(loglvarsq);
plotsongs(song,fs,loglvarsq,name)