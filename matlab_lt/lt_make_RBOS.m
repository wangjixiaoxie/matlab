%% LT 2/5/14 - takes a cbin as input and outputs .wav reverse song.
function [song_data]=lt_make_RBOS(fname)
% fname - file name string as .cbin
% output_name - a string
%%


[song_data, fs]=ReadCbinFile(fname);
song_data=song_data(:,1);
song_data=flipud(song_data);
song_data=song_data/max(song_data);

wavwrite(song_data,fs,16,[fname '_REV.wav']);

end

% 
% %%
% [y,Fs]=audioread('test_RBOS.wav');
% b=audioplayer(y,fs);
% play(b);
