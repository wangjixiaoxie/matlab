vars.first_day=input('first day of analysis (e.g. 21Jul2013, no string)? ','s');
vars.last_day=input('last day of analysis (e.g. 22Jul2013, no string)? ','s');
vars.num_of_pre_days=input('how many pre-WN days are included? ', 's');
vars.syllable = input('what syllable to analyze (e.g. b, no string)? ', 's');

days{1}=datenum(vars.first_day);
days{2}=datenum(vars.last_day);

directory = '/home/lucas/data/song/pu13bk43/all_days_entropy_randomWN/'

for j=days{1}:days{2};
load([directory 'pu13bk43' '_' datestr(j, 'ddmmmyyyy') '_entropy_' vars.syllable '.mat']);
day_relative_to_start=j+1-days{1};
vars.entropy_values_and_times{day_relative_to_start}(:,1)=eval(['data.' vars.syllable]);
vars.entropy_values_and_times{day_relative_to_start}(:,2)=eval(['time.' vars.syllable]);



end


for j=days{1}:days{2};
    day_relative_to_start=j+1-days{1};
    vars.entropy_medians(day_relative_to_start)=median(vars.entropy_values_and_times{day_relative_to_start}(:,1));
    vars.entropy_std(day_relative_to_start)=std(vars.entropy_values_and_times{day_relative_to_start}(:,1));
    vars.entropy_means(day_relative_to_start)=mean(vars.entropy_values_and_times{day_relative_to_start}(:,1));
    vars.entropy_COV(day_relative_to_start)=sqrt(exp(vars.entropy_std(day_relative_to_start)^2)-1); %from wiki - COV for a log-normal distribution, which it is by confirmation in the following cell (i.e. take all points from a day and plot histogram with or without log transform)
        
end

    figure(1); hold on;
    errorbar((days{1}:days{2})-days{1}+1,vars.entropy_means,vars.entropy_std,'o');
    title('median (+/- STD) Wiener entropy (onset to offset) vs. day'); xlabel(['day; last day of baseline is day #' num2str(vars.num_of_pre_days)]); ylabel('Wiener entropy (-inf,0])')

    figure(2); hold on;
    plot((days{1}:days{2})-days{1}+1, vars.entropy_COV);
    title('Coeff. of Variance of weiner entropy(onset to offset) vs day'); xlabel(['day; last day of baseline is day #' num2str(vars.num_of_pre_days)])

% % 
% % save([multiple_pitch.parameters.curr_directory '/' multiple_pitch.parameters.nameofbird '_' ...
% %     multiple_pitch.parameters.phrase '_FF_data.mat'], '-struct', 'multiple_pitch', 'FF')

% save structure

save_dir=['lt_plot_syllable_entropy_over_days_results_' datestr(now,'ddmmmyyyy_HHMM')];
mkdir(save_dir)
cd(save_dir);
save('entropy_data_all_days','-struct', 'vars');

    saveas(figure(1), ['mean_entropy_vs_day_' vars.first_day 'to' vars.last_day],'fig');
    saveas(figure(2), ['COV_entropy_vs_day_' vars.first_day 'to' vars.last_day],'fig');
    fprintf('done!')

    
    
% %% logarithimically or linealrly distributed values of entropy?
% figure; hist(vars.entropy_values_and_times{5}(:,1))
% figure; hist(10.^vars.entropy_values_and_times{5}(:,1))
% 
% close all

