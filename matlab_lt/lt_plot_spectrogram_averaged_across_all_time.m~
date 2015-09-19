% is there filter bias for different frequency bands?  go to make_temp
% within each day folder and another subfolder.  

% do lt_day process just with the solo b syllable, then get the make temp
% which has spec info and plot entire specs for temporal areas of the
% syllable and outside the syllable.  outside syllable should not be diff.
% between days.

%% use 
%
 MAYBE USE mplotcbin to do this. if you want the whole song. below i just use around specific sylalbles.

close all








%%

clear all;
load make_temp.mat
clear temp; 
%first see which time window to take
% disp('run lt_db_make_template') % usually take 90:120, where not much happens.
% keyboard

temp=squeeze(mean(make_temp.all_spectrograms{1}(1:2:256,90:120,:),2)); %1:2:256 becuase that's waht I used to downsample to get 2nd/1st measurements
temp=mean(temp,2);

spec_all_time_averaged{1}=temp;

%see whether the day-to-day differences in spectral power can explain the 2nd/1st results I got
% these are the x-coords that encompass the 3 harmonics which I used for 2nd/1st analysis in lt_count_syllable... : [22 32; 47 59; 76 84]
test_second_over_first_integral_OutsideSyllable=(trapz(temp(47:59)))/(trapz(temp(22:32))) % manually move from day to day and write down results.



%plotting
figure(9); hold on;
plot(temp,'k'); title(['spec averaged over time points 90:120, outside syllable (points, not real time)']);
xlabel('red=7/17, green=7/18. blue =7/19; black=7/20')


%SECOND
clear temp
temp=squeeze(mean(make_temp.all_spectrograms{1}(1:2:256,60:90,:),2));
temp=mean(temp,2);

spec_all_time_averaged{1}=temp;


%plotting
figure(10); hold on;
plot(temp,'k'); title(['spec averaged over time points 60:90, syllable (points, not real time)']);
xlabel('red=7/16, green=7/17, blue=7/18; black=7/19')

test_second_over_first_integral_SyllableOnly=(trapz(temp(47:59)))/(trapz(temp(22:32)))

% THIRD, plot over entire spec time window:
clear temp
temp=squeeze(mean(make_temp.all_spectrograms{1}(1:2:256,:,:),2));
temp=mean(temp,2);

spec_all_time_averaged{1}=temp;

test_second_over_first_integral_AllTimePoints=(trapz(temp(47:59)))/(trapz(temp(22:32)))

%plotting



figure(11); hold on;
plot(temp,'k'); title(['spec averaged over all time points in window (240 pts) (points, not real time); red=7/16, green=7/17, blue=7/18, black=7/19']);


% RATIO OF SYLLABLE 2ND/1ST TO CONTROL 2ND/1ST VALUES. IF THIS DOES NOT
% CHANGE OVER DAYS THEN FILTER ISSUE EXPLAINS WHAT i THOUGHT WAS A REAL
% RESULT OF 2ND HARM OVER 1ST INCRESAING.
syllable_over_outside=test_second_over_first_integral_SyllableOnly/test_second_over_first_integral_OutsideSyllable
syllable_over_all=test_second_over_first_integral_SyllableOnly/test_second_over_first_integral_AllTimePoints