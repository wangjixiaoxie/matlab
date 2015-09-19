
first_day=input('first day of analysis (e.g. 21Jul2013, no string)? ','s');
last_day=input('last day of analysis (e.g. 22Jul2013, no string)? ','s');
num_of_pre_days=input('how many pre-WN days are included? ', 's');
syllable = input('what syllable to analyze (e.g. b, no string)? ', 's');

days{1}=datenum(first_day);
days{2}=datenum(last_day);

directory = '/home/lucas/data/song/pu13bk43/all_days_timing_randomWN/'

for j=days{1}:days{2};
load([directory 'pu13bk43' '_' datestr(j, 'ddmmmyyyy') '_motifduration_' syllable '.mat']);

duration_values{j}(:,1)=data;
duration_values{j}(:,2)=time;


end
%%
for j=days{1}:days{2};
    day_relative_to_start=j+1-days{1};
    duration_medians(day_relative_to_start)=median(duration_values{j}(:,1));
    duration_std(day_relative_to_start)=std(duration_values{j}(:,1));
    duration_means(day_relative_to_start)=mean(duration_values{j}(:,1));

end

    figure; hold on;
    errorbar((days{1}:days{2})-days{1}+1,duration_medians,duration_std,'o');
    title('mean (+/- STD) syllable duration vs. day'); xlabel(['day; last day of baseline is day #' num2str(num_of_pre_days)]); ylabel('syllable duration (ms)')
    saveas(figure(gcf), [directory 'median_duration_vs_day','fig']);
    
    fprintf('done!')


% % 
% % save([multiple_pitch.parameters.curr_directory '/' multiple_pitch.parameters.nameofbird '_' ...
% %     multiple_pitch.parameters.phrase '_FF_data.mat'], '-struct', 'multiple_pitch', 'FF')

