disp('move manually to folder equivalent to /all_days_harmonics_summary/lt_plot_harmonics_structure_201389_1514 and load the data structure. then type return.');
keyboard

first_day=1;
last_day=15;

close all;


for i=1:15;
    
    all_harmonics_second_over_first_log2{i}=log2(all_harmonics_second_over_first{i});
    second_over_first_mean_log2(i)=mean(all_harmonics_second_over_first_log2{i});
    second_over_first_STD_log2(i)=std(all_harmonics_second_over_first_log2{i});

end

%     COV_SecondOverFirst_log2=second_over_first_STD_log2./second_over_first_mean_log2;
%     COV_SecondOverFirst_log2=second_over_first_STD./second_over_first_mean;
    COV_SecondOverFirst_log2=sqrt(exp(second_over_first_STD_log2.^2)-1);% according to wikipedia, for log-normal distributions.

figure(1); hold on;

for index=first_day:last_day;
    maxtime_temp=max(times{index});
    mintime_temp=min(times{index});
    day_duration(index)=maxtime_temp-mintime_temp;
    plot(times{index}+sum(day_duration(1:index-1))+(index-1)*12, all_harmonics_second_over_first_log2{index},'o');
    
    errorbar(10+(index-1)*12+sum(day_duration(1:index))-day_duration(index)/2,second_over_first_mean_log2(index),second_over_first_STD_log2(index),'xr');
    title('ratio of 2nd to 1st harmonic (peak value) for each rendition.  1est 3 days are baseline. Rest of days are 70% WN');
end

figure(2); hold on;
plot(first_day:last_day, COV_SecondOverFirst_log2,'g'); title('COV of ratio of 2nd to 1st harmonic, for each day');

saveas(figure(1),'ratio_second_over_first_all_renditions_log2','fig');
saveas(figure(2),'COV_2ndOverFirst_EachDay','fig');

% %% log normal distribution of ratios for each day?
% close all
% j=2 % day number
% figure;
% hist(all_harmonics_second_over_first_log2{j})
% 
% figure 
% hist(all_harmonics_second_over_first{j});





