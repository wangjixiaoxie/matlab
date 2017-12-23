function lt_seq_dep_pitch_ACROSSBIRDS_LMANbasePlotLrn(DATSTRUCT, ...
    SeqDepPitch_AcrossBirds_LMAN, PARAMS)
%% lt 12/17/17 -



%% pARAMS



%% more similar AFP bias at baseline, more generalization?
lt_figure; hold on;

% ------------- go thru all syls and collect
numsyls = length(DATSTRUCT.singlesyls.AcousVec_PBS);

% comapre to AFP bias at targ
% is this targ?
% initial learning (rel targ dir)
% learnig at targ (positive)
% 4 day learning (rel targ dir)
% 4 day learning targ (positive)
% corr with targ (PBS)
% corr with targ (MUSC)

for i=1:numsyls
    
   % ----------  
    
end


