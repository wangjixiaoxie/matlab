function [ bandpass_song, time, fs, normal_raw, raw, bandpass_raw ] = db_normalize_song_waveform( varargin )
%db_normalize_song_waveform Input a cbin file and it will return a
%bandpassed normalized waveform (from -1 to 1) and time in seconds as
%Will also make a 32-bit wav file if you want
%(second argument = 1). Will also plot the waveform if you want (third
%argument = 1). Will play sound file if you want (fourth argument = 1).
%
% 

[raw, fs] = ReadCbinFile(varargin{1});
normal_raw = raw./max(raw);

time = linspace(0, length(normal_raw)./fs, length(normal_raw))';

bandpass_song = bandpass(normal_raw,fs,500,10000);

bandpass_raw = bandpass(raw,fs,500,10000);

if max(size(varargin)) > 1
    if varargin{2} == 1
        wavwrite(normal_raw,fs,32,[varargin{1}(1:end-4) 'wav']);
    end
end

if max(size(varargin)) > 2
    if varargin{3} == 1
        figure, plot(time,bandpass_song);
        ylim([-1 1]);
        xlim([0 max(time)]);
        title(varargin{1}, 'Interpreter', 'none')
        ylabel('amplitude')
        xlabel('time (sec)')
    end
end

if max(size(varargin)) > 3
    if varargin{4} == 1
        sound(bandpass_song,fs)
    end
end




end

