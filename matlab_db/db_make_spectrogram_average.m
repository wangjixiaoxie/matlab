function [] = db_make_spectrogram_average( spec, time, freq, syllables )
%db_make_spectrogram_average Takes output from db_spec_syllable and creates
%an average spectrogram


figure()
hold on
title(['Average spectrogram of "' syllables '"'], 'Interpreter', 'None')
xlabel('Time (msec)')
ylabel('Frequency (Hz)')
ylim([min(freq) max(freq)])
xlim([min(time) max(time)])

imagesc(time, freq, mean(log(spec),3))
syn



end

