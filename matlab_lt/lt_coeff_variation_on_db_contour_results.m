% do lt_db_contour_and_FF_analysis_over_time_v3 first
% close all
% cd /home/lucas/data/song/all_calls_analysis_pitch/pitch/

disp('go to results folder and load pu13bk43.mat, then type return')
keyboard

% multiple_pitch.FF.a.mean_FF{6}=[]; Use this if you want to ignore a
% specific day
% multiple_pitch.FF.a.sd_FF{6}=[];

% change the syllable to what you have.
 mean=cell2mat(multiple_pitch.FF.b.mean_FF);
 std=cell2mat(multiple_pitch.FF.b.sd_FF);
 
 COV=std./mean;
 
 day1= input('first day? (e.g. 15Jul2013) ','s');
 day2=input('last day? (e.g. 29Jul2013) ','s');
 
 figure(1); plot(COV);
 title(['coeff of variation of FF over days, from' day1 'to' day2])
 
 
 saveas(figure(1),['coeff_of_variation_' day1 '_' day2] ,'fig');
 
 
 
%  %%
%  saveas(figure(gcf), ['FF_over_time_' num2str(multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}+1) 'days_b_skipping04Aug'], 'fig')