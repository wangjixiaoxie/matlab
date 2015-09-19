% steps for saucy
% label files
% run SAUCY_beta(good file,channel (in single quotes), number of models)
%     number of models = 2, always
%     e.g., SAUCY_beta('o85pu54_vol7_2-070804.003.cbin','0',2)
% 
% run SAUCY_beta in batch mode to propagate the model
%     e.g., SAUCY_beta('o85pu54_vol7_2-070804.003.cbin','0','batchmode')
%     
% run SAUCY_beta in track mode to see which files pass muster -> only analyze files where <1.5% error
% 
% calculate spectral features
% figure; one_syl_spectrogram('b')
% use this just to see how good the spectrogram looks
% tweek syllable_params_by_bird.m for specific syllable
% 
% quantify pitch
% quantify_pitch(syllable,sound channel) %this also quantifies amplitude
%     e.g., quantify_pitch('b','obs0r'); % this is correct --- errors may arise from .rec files present w/o .cbin files
% 
% quantify spectral entropy
%     e.g., quantify_spectral_entropy('b','obs0r');
%     
%     
% calculate spike counts for syllables using neural_by_syl_seq_SAUCY_beta(sequence,neural channel)
%     e.g., figure; neural_by_syl_seq_SAUCY_beta('b',0)
%     aligns by first syllable
%     might need to add birdname information (see code)
%     this spits out a figure of the spectrogram, one trace, and a raster plot of spikes around that sequence
%     this saves a 'combined_data_PCA_CH__TH.mat' file that subsequent code access
% 
% run the regression using sound_neural_unique_context(files,syllable)
%     e.g., sound_neural_unique_context('all','b')
%     define premotor window in the code
%     if working with new bird, there are three places in which you need to add birdname