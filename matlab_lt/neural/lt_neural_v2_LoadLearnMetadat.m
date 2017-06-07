function LearningMetaDat = lt_neural_v2_LoadLearnMetadat()

cd('/bluejay5/lucas/analyses/neural/');
load('LearningMetaDat.mat');
eval('!cp LearningMetaDat.mat BACKUP/');  % -- save a backup

