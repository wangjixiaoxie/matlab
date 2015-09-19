

% need to label first.  this is used by db_check_template_timing..., but it
% is stripped down for simplicity here.  

close all; clear all

syllable_freq.batchfile='batch.notcatch';
syllable_freq.freq_range=[2000 4000] % input('what is the freq range? (e.g. [2000 4000]) ');
syllable_freq.target_note= 'b' %input('what is the target syllable? (e.g. b) ','s');



% syllable_freq.values = evtaf_freq(syllable_freq.batchfile, syllable_freq.freq_range, syllable_freq.target_note, 128, 'obs0', 0);
% syllable_freq.values = lt_evtaf_freqall(syllable_freq.batchfile, syllable_freq.target_note, '','',0,syllable_freq.freq_range, 128, 0,'obs0', 0);
syllable_freq.values = evtaf_freqJC(syllable_freq.batchfile, syllable_freq.freq_range, syllable_freq.target_note, 128, 'obs0', 0,0);
syllable_freq.values = evtaf_freqtw(syllable_freq.batchfile, syllable_freq.freq_range, syllable_freq.target_note, 128, 'obs0', 0);



%% FIGURES

% histogram of frequency values over renditions
    figure(1)
    [f,x] = hist(syllable_freq.values(:,2),30);
    bar(x,f/trapz(x,f))
    title(['FF Probability density function for ' syllable_freq.target_note])
    xlabel('Frequency (Hz)')
    ylabel('Density')
%     saveas(figure(1), ['check_template_' check_pitch.parameters.today_date '/FF_PDF_' check_pitch.parameters.sofinterest{i} '_' check_pitch.parameters.today_date], 'fig')
 
% histogram of frequency values over renditions
    figure(1)
    [f,x] = hist(syllable_freq.values.fvirt,30);
    bar(x,f/trapz(x,f))
    title(['FF Probability density function for ' syllable_freq.target_note])
    xlabel('Frequency (Hz)')
    ylabel('Density')
%     saveas(figure(1), ['check_template_' check_pitch.parameters.today_date '/FF_PDF_' check_pitch.parameters.sofinterest{i} '_' check_pitch.parameters.today_date], 'fig')


% scatter plot of values throughout day
times_in_hours=lt_convert_datenum_to_hour(syllable_freq.values(:,1));
figure(2); scatter(times_in_hours, syllable_freq.values(:,2))
