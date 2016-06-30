function [PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANplotraw(SeqDepPitch_AcrossBirds, PARAMS, norm_by_targsyl)
%% LT 7/21/15 - Plots LMAN learning
% norm_by_targsyl = 1; for across birds plots will normalize all metrics
% (learning, AFP, MP) by target syl learning. Default: 1

close all;


%% PARAMS

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

if ~exist('norm_by_targsyl','var');
    norm_by_targsyl=1;
end

%% WARNING TO USER
disp('DID NOT REMOVE SYLLABLES THAT SHOULD BE REMOVED (FOR SINGLE DIR OR BIDIR)');


%% 1) PLOT specific day, all syllables, color coded by whether is musc or not



