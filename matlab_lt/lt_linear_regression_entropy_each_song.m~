% perform linear regression of the entropy of each song versus
% the order # of the rendition within each song.  is there a pattern where
% syllables within each song incerases or decreases in that measure?

%8/19 LT modified from lt_linear_regression_harmonics...
%% 8/19 - try doing it  with the all syllable entropy calculated by db_seq_dunc_day_save...

% after loading data from all_days_entropy... folder in bird folder.


disp('WARNING: only run this after running lt_plot_harmonic_structure_for_all_rendition.m first, as it makes the needed data structures')
disp('manually move to /home/lucas/data/song/pu13bk43/all_days_entropy_randomWN/ that contains your data');
disp('then type return')

keyboard
current_directory=pwd;
first_day=input('first day of analysis (e.g. 21Jul2013, no string)? ','s');
last_day=input('last day of analysis (e.g. 22Jul2013, no string)? ','s');

days{1}=datenum(first_day);
days{2}=datenum(last_day);

cumulated_r=[];
cumulated_slope=[];

for j=days{1}:days{2};
    index=j-days{1}+1;
    clear time0; clear time1; clear time2; clear transition_times

% load the vector with entropy and time values;
disp(['load the correct entropy .mat file for' datestr(j,'DDMMYYYY') 'then type return']);
keyboard
cd(pwd);      

time0=time.b;

%make vector with time differentials (within song should have smaller
%differential.  will use this to determine which syllables were sung in the same song).
time1=time0;
time1(length(time0)+1)=0;
time2(2:length(time0)+1)=time0;
time2(1)=0;
time_diff=time1-time2';
time_diff=time_diff(1:length(time0));
time_diff(1)=time_diff(2);

% figure; hold on; plot(time_diff);plot(time0,'r')

transition_times=find(time_diff>0.001); % 10 seconds in units of hours. transition times are times for 1st note.
transition_times(2:length(transition_times)+1)=transition_times;
transition_times(1)=1;




% % this is to check the above calculation; uncomment  this block (and add an extra end way out down)  to check. (i.e. remove one %)
% figure; hold on; plot(times{index}(transition_times),1.8,'o'); plot(times{index},2,'o'); title('top: all times; bottom: calculated transition times');
% ylim([0 3])
% end

    clear r_regression
    clear slope_regression
    clear offset_regression
    for i=1:length(transition_times)-1;
    

    song_indices{i}=transition_times(i):transition_times(i+1);

[r_regression(i), slope_regression(i), offset_regression(i)]=regression(1:length(song_indices{i}),data.b(song_indices{i})')

    

    end
    
    % compile the results into a vector covering all days.
    cumulated_r=[cumulated_r r_regression];
    cumulated_slope=[cumulated_slope slope_regression];
% plot for just this day, histogram
    figure(index); hold on; hist_bin_contents_15bins(index,:)=hist(slope_regression,15); bar(hist_bin_contents_15bins(index,:),'hist'); title(['histogram of regression slopes for day ' datestr(j,'ddmmmyyyy')]);
    % get histograms with various sized bins
    hist_bin_contents(index,:)=hist(slope_regression); %10 bins
    hist_bin_contents_20bins(index,:)=hist(slope_regression,20);
    
    hist_bin_contents_r_squared_values(index,:)=hist(r_regression.^2,20) % r-squared values in histogram
    
    
end

% PLOT histogram over all days
all_days_hist_bin_contents=sum(hist_bin_contents,1); %10 bins
figure(index+1); hold on; bar(all_days_hist_bin_contents); title('histogram of regression slopes for all days (10bins)');

all_days_hist_bin_contents_15bins=sum(hist_bin_contents_15bins,1); % 15 bins
figure(index+2); hold on; bar(all_days_hist_bin_contents_15bins); title('histogram of regression slopes for all days (15bins)');
 
% 20bins
all_days_hist_bin_contents_20bins=sum(hist_bin_contents_20bins,1); % 15 bins
figure(index+3); hold on; bar(all_days_hist_bin_contents_20bins); title('histogram of regression slopes for all days (20bins)');

% hist of r values for all days
all_days_hist_bin_contents_r_squared_value=sum(hist_bin_contents_r_squared_values,1); % 15 bins
figure(index+5); hold on; bar(all_days_hist_bin_contents_r_squared_value); title('histogram of r-squared values for all days (20bins)');


% % plot histogram of pre vs. post during WN.
hist_bin_contents_pre=sum(hist_bin_contents_15bins(1:3,:),1); % pre days
hist_bin_contents_duringWN=sum(hist_bin_contents_15bins(4:15,:),1); % during WN

figure(index+6); hold on; subplot(2,1,1); bar(hist_bin_contents_pre); title('pre-days sum histogram');
subplot(2,1,2); bar(hist_bin_contents_duringWN); title('sum of days during WN');

% is there a relationship between r-value and slope?
figure(index+4); scatter(cumulated_r.^2,cumulated_slope); title('all slopes vs. r-squared values');
xlabel('r-squared value'); ylabel('slope');

% %% Do all of the above but with permutation (permuting the notes within each variable, thus keeping all other correlations the same.
% 
% % permute r and slope
% 
% 
% 
% 

% saving stuff
suffix2=fix(clock);
cd(current_directory)
mkdir(['lt_linear_regression_entropy_' num2str(suffix2(1)) num2str(suffix2(2)) num2str(suffix2(3)) '_' num2str(suffix2(4)) num2str(suffix2(5))]);
cd(['lt_linear_regression_entropy_' num2str(suffix2(1)) num2str(suffix2(2)) num2str(suffix2(3)) '_' num2str(suffix2(4)) num2str(suffix2(5))])

% save(['harmonics_compiled_data_ratios' datestr(days{1},'ddmmmyyyy') '_to_' datestr(days{2},'ddmmmyyyy')], '-struct', 'harmonics_compiled_data_all_days');

for fig=1:index;
    handle=figure(fig);
    saveas(handle, ['figure_' num2str(handle)], 'fig')
end

saveas(figure(index+1),'histogram_LinRegr_slopes_all_days_10bins', 'fig') 
saveas(figure(index+2),'histogram_LinRegr_slopes_all_days_15bins', 'fig') 
saveas(figure(index+3),'histogram_LinRegr_slopes_all_days_20bins', 'fig') 
saveas(figure(index+4),'histogram_LinRegr_Rsquared_all_days', 'fig') 
saveas(figure(index+5),'LinRegr_slope_vs_Rsquared_all_days', 'fig') 
saveas(figure(index+6),'LinRegr_slope_histogram_Pre_vs_duringWN', 'fig') 


