%% see or100or7 analysis, continued writing there.

%%

% %% make a tone_silent_WN file.
% clear all
% tone_pitch=1000; %in hz
% dur_pretone_silent=1;
% dur_tone=0.05;
% dur_interval=0.1;
% dur_WN=0.065;
% % rel_ampl_tone_vs_WN
% sample_rate=44100;
% 
% %make WN
% white=rand(ceil(sample_rate*dur_WN),1)*2-1;
% 
% % make the tone
% tone=sin(linspace(0, dur_tone*tone_pitch*2*pi, round(dur_tone*sample_rate)));
% 
% % final sound file
% final_sound=[zeros(sample_rate*dur_pretone_silent,1); tone'; zeros(dur_interval*sample_rate,1); white];
% 
% 
% figure; plot(final_sound)
% file_name='test.wav';
% audiowrite(file_name,final_sound,sample_rate);

%% Make tone-WN, but using SYS WN file for my WN
clear all
file='/bluejay3/lucas/birds/sounds/WN/WN_SYS/wn60length200.wav';
[white Fs]=wavread(file);

tone_pitch=500; %in hz
dur_pretone_silent=1;
dur_tone=0.1;
dur_interval=0.1;
% dur_WN=0.065;
rel_ampl_tone_vs_WN=0.1;
sample_rate=44100;
upsweep=[5000 8000]; % if do upsweep [startFreq endFreq]

%make WN
% white=rand(ceil(sample_rate*dur_WN),1)*2-1;

if (0)
    % make the pure tone
    tone=rel_ampl_tone_vs_WN*sin(linspace(0, dur_tone*tone_pitch*2*pi, round(dur_tone*sample_rate)));
end

% make the tone dirty
if (0)
    length_tone=size(tone,2);
    scaling_vector = 0.8 + (1.2-0.8).*rand(length_tone,1);
    tone_dirty=rel_ampl_tone_vs_WN*tone.*scaling_vector';
end

% make an upsweep
if (1)
    t=0:1/sample_rate:dur_tone;
    tone_UpSweep=rel_ampl_tone_vs_WN*chirp(t,upsweep(1),1,upsweep(2));
end

% final sound file
try
    final_sound=[zeros(sample_rate*dur_pretone_silent,1); tone'; zeros(dur_interval*sample_rate,1); white];
catch
    try
        final_sound=[zeros(sample_rate*dur_pretone_silent,1); tone_dirty'; zeros(dur_interval*sample_rate,1); white]; % for dirty tone
    catch
        final_sound=[zeros(sample_rate*dur_pretone_silent,1); tone_UpSweep'; zeros(dur_interval*sample_rate,1); white]; % for dirty tone
    end
end


figure; plot(final_sound)
file_name='/bluejay3/lucas/birds/sounds/test3_sysWN_Up.wav';
audiowrite(file_name,final_sound,sample_rate);

% to listen
% clear all
% file='test2_sysWN.wav';
[y Fs]=wavread(file_name);
player=audioplayer(y,Fs);
%%
play(player)




%% To listen to something else

clear all
file='wn60length200.wav';
[y Fs]=wavread(file);
player=audioplayer(y,Fs);
play(player)


