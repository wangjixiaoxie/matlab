clear make_temp.all_spectrogram_spectral_means



for i = 1:length(make_temp.parameters.syllables);
       
% for each syllable, calculate spectral mean using the three harmonics (i.e. freq that is
        % midpoint of the three harmonics, weighted by their amplitudes.)
        % (e.g. gives range from 1-3 depending on which harmonic dominates.
        make_temp.all_spectrogram_spectral_means{i}=nan(1,size(make_temp.all_spectrograms{i},3));
        for ind=1:size(make_temp.all_spectrograms{i},3);
             make_temp.all_spectrogram_spectral_means{i}(ind)=(sum((make_temp.all_spectrogram_ampl_first_three_harmonics{i}(:,ind).*[1;2;3]),1))./(sum(make_temp.all_spectrogram_ampl_first_three_harmonics{i}(:,ind),1));
        end
        
        % plot spectral means, and get mean and std over trials.  
        make_temp.mean_over_renditions_spectral_mean{i}=mean(make_temp.all_spectrogram_spectral_means{i}(:));
        make_temp.STD_over_renditions_spectral_mean{i}=std(make_temp.all_spectrogram_spectral_means{i}(:))
        disp('press any key to continue')
        pause 
        
end
