%% LT 2/5/14 - takes a cbin as input and outputs .wav song.
function [song_data]=lt_make_WAV(fname)
% fname - file name string as .cbin
% output_name - a string
%%


[song_data, fs]=ReadCbinFile(fname);
song_data=song_data(:,1);
song_data=song_data/max(song_data);
wavwrite(song_data,fs,16,[fname '.wav']);

end

% 
% %%
% [y,Fs]=audioread('test_RBOS.wav');
% b=audioplayer(y,fs);
% play(b);
