
first_day=input('first day of analysis (e.g. 21Jul2013, no string)? ','s');
last_day=input('last day of analysis (e.g. 22Jul2013, no string)? ','s');
num_of_pre_days=input('how many pre-WN days are included? ', 's');
syllable = input('what syllable to analyze (e.g. b, no string)? ', 's');

days{1}=datenum(first_day);
days{2}=datenum(last_day);

directory = '/home/lucas/data/song/pu13bk43/all_days_fraction_randomWN/';

counter=1;
for j=days{1}:days{2};
try 
    load([directory 'pu13bk43' '_' datestr(j, 'ddmmmyyyy') '_fraction_' syllable '.mat']);
catch
    fraction_bootstrap_medians(counter)=nan;
    fraction_bootstrap_std(counter)=nan;
    counter=counter+1;
    continue
end

fraction_bootstrap_medians(counter)=percentiles(2);
fraction_bootstrap_std(counter)=std(bootstrap);
counter=counter+1;
end

%% plot
    figure; hold on;
%     plot((days{1}:days{2})-days{1}+1,fraction_bootstrap_medians,'o')
    errorbar((days{1}:days{2})-days{1}+1,fraction_bootstrap_medians,fraction_bootstrap_std,'o');
    title('mean (+/- STD) syllable fraction in song vs. day'); xlabel(['day; last day of baseline is day #' num2str(num_of_pre_days)]); ylabel('fraction of syllables per song that are your syllable')
    saveas(figure(gcf), [directory 'median_fraction_vs_day' first_day 'to' last_day],'fig');
    
    fprintf('done!')


% % 
% % save([multiple_pitch.parameters.curr_directory '/' multiple_pitch.parameters.nameofbird '_' ...
% %     multiple_pitch.parameters.phrase '_FF_data.mat'], '-struct', 'multiple_pitch', 'FF')

