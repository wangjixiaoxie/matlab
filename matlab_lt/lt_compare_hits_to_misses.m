

% 1) get hits vs. misses, and calculate 1) spectrogram, 2) harmonics.

clear all;
close(figure(1));

keep_combined=input(' if you want to keep hits and misses combined type "y" (i.e. you do not want to compare). Otherwise type "n": ', 's');

if keep_combined=='y';
   add_tag='C';
else 
    add_tag='';
end
    
disp('move to day folder manually, then type return')
keyboard

suffix=fix(clock);
mkdir(['lt_compare_hits_to_misses_results_' num2str(suffix(1)) num2str(suffix(2)) num2str(suffix(3)) '_' num2str(suffix(4)) num2str(suffix(5)) add_tag]);
save_dir=[pwd '/lt_compare_hits_to_misses_results_' num2str(suffix(1)) num2str(suffix(2)) num2str(suffix(3)) '_' num2str(suffix(4)) num2str(suffix(5)) add_tag];



% change missed labels if desired

if keep_combined == 'n';
    
compare_HitsToMisses.batch_file='batch.catch.keep';
compare_HitsToMisses.old_syllable='b';
compare_HitsToMisses.new_syllable='N';
   
db_change_miss_label(compare_HitsToMisses.batch_file, compare_HitsToMisses.old_syllable, ...
    compare_HitsToMisses.new_syllable);

disp('missed syllable labels changed from b to N. perform following analyses on both syllables at once');

end

disp('missed syllable not converted to new label.')


% compare hit vs. missed spectrograms and slices

lt_count_syllable_harmonics_v2_HitsVsMiss


% calculate hit % and overall average syllable. 


% create summary plot over days of hit vs. missed spectrograms

% save things
save([save_dir '/peak_statistics'], 'peak_statistics')
save([save_dir '/make_temp'],'make_temp');
saveas(figure(1), [save_dir '/day_figs'], 'fig')
saveas(figure(2), [save_dir '/slices_AreaAndNotNormalized'], 'fig')
saveas(figure(3), [save_dir '/harmonic_relative_ampl_AreaAndNotNormalized'], 'fig')
saveas(figure(4), [save_dir '/all_slices_AreaAndNotNormalized'], 'fig')

if keep_combined=='n';
% change syllable labels back to all b
db_change_miss_label(compare_HitsToMisses.batch_file, compare_HitsToMisses.new_syllable, ...
    compare_HitsToMisses.old_syllable);
end

cd(save_dir)
